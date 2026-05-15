# encoding: UTF-8
# =============================================================================
# Clausulas F0 -- Shiny App
# =============================================================================
library(shiny)
library(tidyverse)
library(zoo)
library(ggplot2)
library(plotly)
library(scales)
library(DT)
library(readr)

# =============================================================================
# FUNCIONES ANALITICAS
# =============================================================================

hz_a_st <- function(f1, f2) 12 * log2(f2 / f1)

slope_st_fn <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  pos <- seq_along(x)
  coef(lm(x ~ pos))[["pos"]] * 12 / (log(2) * mean(x, na.rm = TRUE))
}

slope_pct_fn <- function(x) {
  if (sum(!is.na(x)) < 2) return(NA_real_)
  pos <- seq_along(x)
  coef(lm(x ~ pos))[["pos"]] / mean(x, na.rm = TRUE) * 100
}

# Per-clause slope on the full interleaved trajectory [ini1,fin1,ini2,fin2,…].
# Multiplied by 2 so the result is in ST per group (comparable to umbral_st).
slope_st_interleaved <- function(f0_ini, f0_fin) {
  vals <- as.vector(rbind(f0_ini, f0_fin))
  if (sum(!is.na(vals)) < 2) return(NA_real_)
  pos <- seq_along(vals)
  coef(lm(vals ~ pos))[["pos"]] * 2 * 12 / (log(2) * mean(vals, na.rm = TRUE))
}

# Direction from a numeric value vs threshold
etiquetar <- function(v, u) {
  case_when(is.na(v) ~ NA_character_, v > u ~ "ascendente",
            v < -u ~ "descendente", TRUE ~ "plano")
}

# Direction from within-group movement (f0_ini -> f0_fin)
dir_propio <- function(f0_ini, f0_fin, u) {
  st <- 12 * log2(f0_fin / f0_ini)
  case_when(st > u ~ "ascendente", st < -u ~ "descendente", TRUE ~ "plano")
}

# Assign clausulaN labels to consecutive runs of same direction.
# force_split: logical vector; TRUE at position i forces a new clause even
# when direction matches the previous group (valley-pattern split).
asignar_clausulas <- function(dir, force_split = rep(FALSE, length(dir))) {
  out <- character(length(dir)); cnt <- 1L; prev <- NA_character_
  for (i in seq_along(dir)) {
    d <- dir[i]
    f <- !is.na(force_split[i]) && force_split[i]
    if (is.na(d)) {
      out[i] <- if (i == 1) "clausula1" else out[i - 1]
    } else if (is.na(prev) || d != prev || f) {
      out[i] <- paste0("clausula", cnt); cnt <- cnt + 1L; prev <- d
    } else {
      out[i] <- out[i - 1]
    }
  }
  out
}

# Running-sum progression detector (Option A)
# Returns acum (cumulative ST/pct) and etiq (progresion label) per group.
# A step contributes to the current run if it points the same direction as the
# existing run (or if no direction has been established yet).
# When a non-flat step reverses the current direction the accumulator resets.
progresion_fn <- function(reajuste, umbral_paso, umbral_total) {
  n    <- length(reajuste)
  acum <- rep(NA_real_, n)
  etiq <- rep("neutro",  n)
  s <- 0; d <- 0L
  for (i in seq_len(n)) {
    r <- reajuste[i]
    if (is.na(r)) { acum[i] <- s; etiq[i] <- "neutro"; next }
    step <- if (r >  umbral_paso) 1L else if (r < -umbral_paso) -1L else 0L
    if (step != 0L && d != 0L && step != d) {
      s <- r; d <- step          # direction reversal: reset accumulator
    } else {
      if (step != 0L) d <- step  # adopt direction of first non-flat step
      s <- s + r
    }
    acum[i] <- round(s, 3)
    etiq[i] <- if (s >=  umbral_total) "prog_asc"
               else if (s <= -umbral_total) "prog_desc"
               else "neutro"
  }
  list(acum = acum, etiq = etiq)
}

# =============================================================================
# PARSEO: formato ancho -> largo
# =============================================================================
parse_wide <- function(df) {
  result    <- list()
  data_cols <- names(df)[-1]
  n_grupos  <- floor(length(data_cols) / 3)
  for (i in seq_len(nrow(df))) {
    row       <- df[i, ]
    enunciado <- as.character(row[[1]])
    for (g in seq_len(n_grupos)) {
      ini_col   <- data_cols[(g - 1) * 3 + 1]
      fin_col   <- data_cols[(g - 1) * 3 + 2]
      pausa_col <- data_cols[(g - 1) * 3 + 3]
      f0_ini <- suppressWarnings(as.numeric(row[[ini_col]]))
      f0_fin <- suppressWarnings(as.numeric(row[[fin_col]]))
      pausa  <- suppressWarnings(as.numeric(row[[pausa_col]]))
      if (!is.na(f0_ini) && !is.na(f0_fin)) {
        result[[length(result) + 1]] <- tibble(
          enunciado = enunciado, grupo = g,
          f0_ini = f0_ini, f0_fin = f0_fin, pausa = pausa)
      }
    }
  }
  if (length(result) == 0) return(NULL)
  bind_rows(result)
}

# =============================================================================
# CALCULO PRINCIPAL
# =============================================================================
calcular <- function(dl, ventana, umbral_pct, umbral_st, umbral_local,
                     umbral_prog_st, umbral_prog_pct, consecutivo = TRUE) {

  # Global sequential index (always across all groups regardless of mode)
  dl <- dl |> mutate(x = row_number())

  # Helper: renumber local clause names (per-utterance) to globally sequential ones
  renumber_clauses <- function(enun, local_cl) {
    keys  <- paste0(enun, "|||", local_cl)
    ukeys <- unique(keys)
    mapping <- setNames(paste0("clausula", seq_along(ukeys)), ukeys)
    mapping[keys]
  }

  # Steps 1-3: per-group metrics, rolling slope, clause assignment.
  # When consecutivo = TRUE  → computed across all groups (existing behaviour).
  # When consecutivo = FALSE → computed within each utterance independently.
  grp <- if (consecutivo) dl else group_by(dl, enunciado)

  dl2 <- grp |>
    mutate(
      dir_grupo     = dir_propio(f0_ini, f0_fin, umbral_st),
      delta_hz      = f0_ini - lag(f0_ini),
      delta_pct     = delta_hz / lag(f0_ini) * 100,
      delta_st      = hz_a_st(lag(f0_ini), f0_ini),
      dir_local     = etiquetar(delta_st, umbral_local),
      reajuste_hz   = f0_ini - lag(f0_fin),
      reajuste_pct  = (f0_ini - lag(f0_fin)) / lag(f0_fin) * 100,
      reajuste_st   = hz_a_st(lag(f0_fin), f0_ini),
      inflexion_pct = (f0_fin - f0_ini) / f0_ini * 100,
      inflexion_st  = hz_a_st(f0_ini, f0_fin),
      # Valley-pattern split: current rises AND previous fell (asymmetric).
      force_split   = inflexion_st > umbral_local &
                      !is.na(lag(inflexion_st)) &
                      lag(inflexion_st) < -umbral_local,
      force_split   = coalesce(force_split, FALSE),
      slope_st_v    = rollapply(f0_ini, width = ventana, FUN = slope_st_fn,
                                fill = NA, align = "right", partial = TRUE),
      slope_pct_v   = rollapply(f0_ini, width = ventana, FUN = slope_pct_fn,
                                fill = NA, align = "right", partial = TRUE),
      dir_slope_st  = etiquetar(slope_st_v,  umbral_st),
      dir_slope_pct = etiquetar(slope_pct_v, umbral_pct),
      # In non-consecutive mode, fall back to within-group direction when the
      # rolling window is too short (typically the first group of each utterance).
      # This makes clause assignment sensitive to each group's own movement.
      .dir_cl_st  = if (consecutivo) dir_slope_st  else coalesce(dir_slope_st,  dir_grupo),
      .dir_cl_pct = if (consecutivo) dir_slope_pct else coalesce(dir_slope_pct, dir_grupo),
      .cl_st      = asignar_clausulas(.dir_cl_st, force_split),
      .cl_pct     = asignar_clausulas(.dir_cl_pct),
      .cl_loc     = asignar_clausulas(dir_local)
    ) |>
    ungroup() |>
    mutate(
      clausula_slope_st  = if (consecutivo) .cl_st  else renumber_clauses(enunciado, .cl_st),
      clausula_slope_pct = if (consecutivo) .cl_pct else renumber_clauses(enunciado, .cl_pct),
      clausula_local     = if (consecutivo) .cl_loc else renumber_clauses(enunciado, .cl_loc),
      clausula           = clausula_slope_st
    ) |>
    select(-.cl_st, -.cl_pct, -.cl_loc, -.dir_cl_st, -.dir_cl_pct)

  # Step 4: per-clause slope (fits regression only on groups within each clause)
  # This corrects direction for multi-group clauses where the rolling window
  # bleeds across clause boundaries.
  # Single-group clauses use the within-group direction (dir_grupo) instead.
  cl_stats <- dl2 |>
    group_by(clausula) |>
    summarise(
      n_cl         = n(),
      slope_cl_st  = if (n() >= 2) slope_st_interleaved(f0_ini, f0_fin)            else NA_real_,
      slope_cl_pct = if (n() >= 2) slope_pct_fn((f0_ini + f0_fin) / 2)             else NA_real_,
      .groups = "drop"
    ) |>
    mutate(
      dir_cl_st  = etiquetar(slope_cl_st,  umbral_st),
      dir_cl_pct = etiquetar(slope_cl_pct, umbral_pct)
    )

  dl3 <- dl2 |>
    left_join(cl_stats, by = "clausula") |>
    mutate(
      dir_slope_st = case_when(
        n_cl >= 2 & !is.na(dir_cl_st) ~ dir_cl_st,   # multi-group: per-clause slope
        n_cl == 1                      ~ dir_grupo,    # single-group: own movement
        TRUE                           ~ dir_slope_st  # fallback
      ),
      dir_slope_pct = case_when(
        n_cl >= 2 & !is.na(dir_cl_pct) ~ dir_cl_pct,
        n_cl == 1                       ~ dir_grupo,
        TRUE                            ~ dir_slope_pct
      ),
      slope_st_v  = coalesce(slope_cl_st,  slope_st_v),
      slope_pct_v = coalesce(slope_cl_pct, slope_pct_v)
    )

  # Step 5: enunciado-level slope (F0 media per utterance)
  dl_enun <- dl3 |>
    group_by(enunciado) |>
    summarise(f0_media_enun = mean(f0_ini, na.rm = TRUE), .groups = "drop") |>
    mutate(
      slope_st_enun  = rollapply(f0_media_enun, width = ventana, FUN = slope_st_fn,
                                  fill = NA, align = "right", partial = TRUE),
      slope_pct_enun = rollapply(f0_media_enun, width = ventana, FUN = slope_pct_fn,
                                  fill = NA, align = "right", partial = TRUE),
      dir_enun_st    = etiquetar(slope_st_enun,  umbral_st),
      dir_enun_pct   = etiquetar(slope_pct_enun, umbral_pct),
      clausula_enun_st  = asignar_clausulas(dir_enun_st),
      clausula_enun_pct = asignar_clausulas(dir_enun_pct)
    )

  final <- dl3 |>
    left_join(dl_enun |> select(enunciado, f0_media_enun,
                                 slope_st_enun, slope_pct_enun,
                                 dir_enun_st, dir_enun_pct,
                                 clausula_enun_st, clausula_enun_pct),
              by = "enunciado")

  prog_st  <- progresion_fn(final$reajuste_st,  umbral_local,  umbral_prog_st)
  prog_pct <- progresion_fn(final$reajuste_pct, umbral_pct,    umbral_prog_pct)

  final |> mutate(
    acum_st        = prog_st$acum,
    progresion_st  = prog_st$etiq,
    acum_pct       = prog_pct$acum,
    progresion_pct = prog_pct$etiq
  )
}

# =============================================================================
# GRAFICO
# =============================================================================
hacer_grafico <- function(res, ventana, umbral_st, umbral_local) {
  cl_dirs  <- res |> filter(!is.na(dir_slope_st)) |> distinct(clausula) |> pull()
  cl_nodir <- setdiff(unique(res$clausula), cl_dirs)
  pal <- c(setNames(hue_pal()(max(length(cl_dirs), 1)), cl_dirs),
           setNames(rep("grey65", length(cl_nodir)), cl_nodir))

  ancho <- 0.35
  conexiones <- res |>
    mutate(x_start = x, x_end = lead(x), y_start = f0_fin, y_end = lead(f0_ini)) |>
    filter(!is.na(x_end))

  y_max <- max(c(res$f0_ini, res$f0_fin), na.rm = TRUE)
  flechas <- res |>
    filter(!is.na(dir_slope_st)) |>
    group_by(clausula, dir_slope_st) |>
    summarise(x_ini = min(x), x_fin = max(x), .groups = "drop") |>
    arrange(x_ini) |>
    mutate(fila = row_number() %% 2,
           y_arr  = y_max * 1.07 + fila * y_max * 0.05,
           y_lbl  = y_arr + y_max * 0.055,
           x_ini2 = ifelse(x_ini == x_fin, x_ini - 0.25, x_ini),
           x_fin2 = ifelse(x_ini == x_fin, x_fin + 0.25, x_fin))

  ggplot(res) +
    geom_segment(data = conexiones,
      aes(x = x_start + ancho, xend = x_end - ancho, y = y_start, yend = y_end),
      color = "grey72", linewidth = 0.55, linetype = "dashed") +
    geom_segment(aes(x = x - ancho, xend = x + ancho, y = f0_ini, yend = f0_fin,
                     color = clausula), linewidth = 2, lineend = "round") +
    geom_point(aes(x = x - ancho, y = f0_ini, fill = clausula),
               shape = 21, size = 2.8, color = "white", stroke = 0.7) +
    geom_point(aes(x = x + ancho, y = f0_fin, fill = clausula),
               shape = 24, size = 2.8, color = "white", stroke = 0.7) +
    geom_text(aes(x = x, y = pmin(f0_ini, f0_fin),
                  label = paste0(enunciado, "\ng", grupo)),
              vjust = 2.5, size = 2.4, color = "grey40", fontface = "italic") +
    geom_segment(data = flechas,
      aes(x = x_ini2 - ancho, xend = x_fin2 + ancho, y = y_arr, yend = y_arr,
          color = clausula),
      linewidth = 1.05,
      arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
    geom_text(data = flechas,
      aes(x = (x_ini2 + x_fin2) / 2, y = y_lbl,
          label = paste0(clausula, "\n(", dir_slope_st, ")"), color = clausula),
      size = 2.4, lineheight = 0.85, show.legend = FALSE) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    scale_x_continuous(breaks = res$x, labels = paste0(res$enunciado, "\ng", res$grupo),
                       expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.42))) +
    labs(title = "Progresiones de F0 por clausula",
         subtitle = sprintf("o=F0 ini  ^=F0 fin  |  Ventana: %d  |  Umbral ST: %.2f  |  Umbral local: %.2f ST",
                            ventana, umbral_st, umbral_local),
         x = NULL, y = "F0 (Hz)", color = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey50", size = 9),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 8, color = "grey40"),
          legend.position = "bottom")
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{background:#f4f3f0;font-family:'Helvetica Neue',Arial,sans-serif;color:#1a1a1a;margin:0;padding:0 24px 40px}
    .main-title{font-family:monospace;font-weight:700;font-size:22px;letter-spacing:-0.5px;padding:24px 0 2px}
    .main-subtitle{font-size:12px;color:#777;margin-bottom:18px;font-family:monospace}
    .params-bar{background:#1a1a1a;border-radius:10px;padding:14px 22px;margin-bottom:16px;display:flex;align-items:flex-end;gap:18px;flex-wrap:wrap}
    .params-bar label{color:#aaa!important;font-size:11px!important;font-family:monospace;text-transform:uppercase;letter-spacing:0.5px}
    .params-bar .form-control{background:#2d2d2d!important;border:1px solid #444!important;color:#f0f0f0!important;font-family:monospace;font-size:13px;border-radius:6px;width:90px!important;height:34px}
    .params-bar .form-control:focus{border-color:#8b9eff!important;box-shadow:0 0 0 2px rgba(139,158,255,0.2)!important}
    .params-bar select.form-control{width:130px!important}
    .param-label{color:#8b9eff;font-family:monospace;font-size:10px;text-transform:uppercase;letter-spacing:1px;margin-bottom:3px}
    .main-box{background:#fff;border-radius:12px;box-shadow:0 2px 14px rgba(0,0,0,0.07);overflow:hidden}
    .nav-tabs{background:#f8f7f4;border-bottom:2px solid #e8e7e4!important;padding:0 8px;margin:0}
    .nav-tabs>li>a{font-family:monospace!important;font-size:12px!important;font-weight:700!important;color:#888!important;border:none!important;border-bottom:3px solid transparent!important;padding:12px 20px!important;border-radius:0!important;background:transparent!important;text-transform:uppercase;letter-spacing:0.5px;transition:color 0.15s}
    .nav-tabs>li.active>a,.nav-tabs>li>a:hover{color:#1a1a1a!important;border-bottom:3px solid #1a1a1a!important;background:transparent!important}
    .tab-content{padding:24px}
    .nav-pills>li>a{font-family:monospace;font-size:11px;font-weight:600;text-transform:uppercase;letter-spacing:0.5px;color:#555;background:#f0efe9;border-radius:20px;padding:6px 16px;margin-right:6px}
    .nav-pills>li.active>a{background:#1a1a1a;color:#fff}
    textarea.data-textarea{font-family:monospace;font-size:12px;background:#fafaf8;border:1.5px solid #e0dfd9;border-radius:8px;padding:12px;width:100%;box-sizing:border-box;resize:vertical;color:#1a1a1a;line-height:1.7}
    textarea.data-textarea:focus{outline:none;border-color:#1a1a1a;box-shadow:0 0 0 3px rgba(26,26,26,0.05)}
    .btn-parse{background:#1a1a1a;color:#fff;font-family:monospace;font-size:12px;font-weight:700;border:none;border-radius:7px;padding:9px 22px;cursor:pointer;text-transform:uppercase;letter-spacing:0.5px}
    .btn-parse:hover{background:#333;color:#fff}
    .btn-clear{background:transparent;color:#888;font-family:monospace;font-size:12px;border:1.5px solid #ddd;border-radius:7px;padding:9px 18px;cursor:pointer}
    .btn-clear:hover{border-color:#888;color:#333}
    .sec-label{font-family:monospace;font-size:11px;text-transform:uppercase;letter-spacing:1px;color:#999;margin-bottom:8px;margin-top:16px}
    .sec-label:first-child{margin-top:0}
    .info-box{background:#f0f4ff;border-radius:8px;padding:11px 15px;font-size:11.5px;color:#555;font-family:monospace;margin-bottom:14px;border-left:3px solid #8b9eff;line-height:1.7}
    .stat-card{background:#f8f7f4;border-radius:10px;padding:16px 20px;margin-bottom:12px;border-left:4px solid #1a1a1a}
    .stat-card .stat-val{font-family:monospace;font-size:32px;font-weight:700;color:#1a1a1a;line-height:1}
    .stat-card .stat-lbl{font-size:11px;color:#888;margin-top:5px;font-family:monospace;text-transform:uppercase;letter-spacing:0.5px}
    .stat-asc{border-left-color:#27ae60}.stat-desc{border-left-color:#e74c3c}.stat-plan{border-left-color:#f39c12}
    .dataTables_wrapper{font-family:sans-serif;font-size:13px}
    table.dataTable thead th{font-family:monospace;font-size:11px;text-transform:uppercase;letter-spacing:0.5px;background:#f8f7f4;color:#555;border-bottom:2px solid #e0dfd9!important}
  "))),

  div(class = "main-title", "Clausulas F0"),
  div(class = "main-subtitle", "Analisis prosodico -- slope y direccion entonativa"),

  # PARAMS BAR
  div(class = "params-bar",
    div(div(class="param-label","Ventana slope"),
        numericInput("ventana",      NULL, value=3,   min=2, max=10,  step=1,    width="90px")),
    div(div(class="param-label","Umbral % / enun."),
        numericInput("umbral_pct",   NULL, value=5,   min=0.1,        step=0.5,  width="90px")),
    div(div(class="param-label","Umbral ST / enun."),
        numericInput("umbral_st",    NULL, value=0.8, min=0.05,       step=0.05, width="90px")),
    div(div(class="param-label","Umbral local ST"),
        numericInput("umbral_local", NULL, value=0.5, min=0.05,       step=0.05, width="90px")),
    div(div(class="param-label","Nivel slope"),
        selectInput("nivel_slope", NULL,
                    choices  = c("Grupo"="grupo","Enunciado"="enunciado"),
                    selected = "grupo", width = "130px")),
    div(div(class="param-label","Prog. umbral ST"),
        numericInput("umbral_prog_st",  NULL, value=2.0, min=0.5, step=0.5, width="90px")),
    div(div(class="param-label","Prog. umbral %"),
        numericInput("umbral_prog_pct", NULL, value=15,  min=1,   step=1,   width="90px"))
  ),

  div(class = "main-box",
    tabsetPanel(id = "tabs",

      # TAB 1: Datos
      tabPanel("01 | Datos",
        div(class="sec-label","Pegar, escribir o subir datos"),
        div(class="info-box", HTML(
          "Formato (separado por tabulaciones):<br>
           <b>enunciados | F01 | F02 | pausa1 | F03 | F04 | pausa2 | ...</b><br>
           Cada trio (F0_ini, F0_fin, pausa) = un grupo tonal. Celdas vacias = grupo inexistente.")),
        tags$textarea(id="raw_text", class="data-textarea", rows="10",
          placeholder="enunciados\tF01\tF02\tpausa1\tF03\tF04\tpausa2\nenun1\t120\t135\t0.3\t128\t145\t0.5\nenun2\t138\t155\t0.4\t\t\t"),
        br(), br(),
        div(class="sec-label","O subir archivo (.txt / .csv / .tsv)"),
        fileInput("file_upload", NULL, accept=c(".txt",".csv",".tsv"),
                  buttonLabel="Elegir archivo", placeholder="ningun archivo"),
        br(),
        div(style="display:flex;gap:10px;align-items:center;",
          actionButton("btn_parse","Procesar datos", class="btn-parse"),
          actionButton("btn_clear","Limpiar",        class="btn-clear"),
          div(style="margin-left:6px;",
            checkboxInput("consecutivo", "Grupos consecutivos entre enunciados",
                          value = TRUE)))
      ),

      # TAB 2: Tabla
      tabPanel("02 | Tabla",
        div(class="sec-label","Datos transformados -- un grupo tonal por fila"),
        DTOutput("tabla_larga"),
        br(),
        div(class="sec-label","Descripcion de variables"),
        div(class="info-box", style="line-height:1.9;", HTML("
          <b>Enunciado</b>: identificador del enunciado (frase, turno u oracion). Unidad de agrupacion superior, tomada directamente del archivo de entrada.<br>
          <b>Grupo</b>: numero de grupo tonal dentro del enunciado (1, 2, 3&hellip;). Cada trio de columnas en el archivo (F0 ini, F0 fin, pausa) define un grupo.<br>
          <b>F0 ini</b>: frecuencia fundamental en Hz al inicio del grupo. Punto de anclaje de todos los calculos de transicion y slope.<br>
          <b>F0 fin</b>: frecuencia fundamental en Hz al final del grupo. Junto con F0 ini define el contorno interno.<br>
          <b>Pausa</b>: duracion en segundos de la pausa al final del grupo (0 = sin pausa). Usada en la pestana 05 para correlacionar con el reajuste.<br>
          <br>
          <b>Inflexion %</b>: movimiento interno del grupo en porcentaje [(F0_fin &minus; F0_ini) / F0_ini &times; 100]. Positivo = sube; negativo = baja.<br>
          <b>Inflexion ST</b>: mismo movimiento en semitonos [12 &times; log2(F0_fin / F0_ini)]. Escala perceptivamente uniforme. Determina <i>Dir. propio</i>.<br>
          <br>
          <b>Reajuste %</b>: salto de F0 en la juntura entre grupos en porcentaje [(F0_ini(n) &minus; F0_fin(n&minus;1)) / F0_fin(n&minus;1) &times; 100]. Reset de altura desde donde termino el grupo anterior.<br>
          <b>Reajuste ST</b>: mismo salto en semitonos [12 &times; log2(F0_ini(n) / F0_fin(n&minus;1))]. Variable principal de transicion; usada en la pestana 05 y en la progresion acumulada.<br>
          <b>Acum. ST</b>: suma acumulada de Reajuste ST desde el ultimo cambio de direccion. Se reinicia cuando el reajuste invierte la tendencia (supera el umbral local en sentido contrario).<br>
          <b>Progresion</b>: etiqueta derivada de Acum. ST. <i>prog_asc</i> / <i>prog_desc</i> cuando la acumulacion supera el umbral configurado; <i>neutro</i> si no se ha alcanzado.<br>
          <br>
          <b>Dir. propio</b>: direccion del contorno interno del grupo (ascendente / descendente / plano) segun Inflexion ST frente al umbral ST. Fallback para clausulas de un unico grupo.<br>
          <b>dST</b>: diferencia en ST entre F0 ini del grupo actual y F0 ini del anterior [12 &times; log2(F0_ini(n) / F0_ini(n&minus;1))]. Solo se usa para calcular Dir. local.<br>
          <b>Dir. local</b>: direccion asignada a partir de dST frente al umbral local. Alimenta <i>clausula_local</i> (asignacion alternativa, no la principal).<br>
          <br>
          <b>Slope ST</b>: pendiente de regresion lineal de F0_ini sobre los grupos de la ventana deslizante, en ST. Para clausulas de 2+ grupos se recalcula sobre todos sus grupos (corrige efecto de borde). Magnitud central del analisis de tendencia.<br>
          <b>Dir. slope</b>: direccion definitiva del grupo (ascendente / descendente / plano) segun Slope ST frente al umbral ST. Usada en el grafico, en los boxplots y en el analisis chi cuadrado.<br>
          <b>Clausula</b>: clausula entonativa asignada al grupo (clausula1, clausula2&hellip;). Agrupa runs consecutivos de igual Dir. slope. Variable de segmentacion principal de toda la app.
        "))
      ),

      # TAB 3: Grafico
      tabPanel("03 | Grafico",
        div(class="sec-label","Progresiones de F0 por clausula"),
        fluidRow(column(12,
          div(style="display:flex;align-items:flex-end;gap:16px;flex-wrap:wrap;margin-bottom:14px;",
            div(div(class="param-label","Primer grupo"),
                numericInput("graf_desde",NULL,value=1, min=1,step=1,width="80px")),
            div(div(class="param-label","Ultimo grupo"),
                numericInput("graf_hasta",NULL,value=5, min=1,step=1,width="80px")),
            div(style="display:flex;gap:8px;align-items:flex-end;padding-bottom:2px;",
              actionButton("btn_todo",    "Todos",     class="btn-clear"),
              actionButton("btn_primero5","Primeros 5",class="btn-clear"),
              uiOutput("lbl_n_grupos"))))),
        plotlyOutput("grafico", height="520px"),
        br(),
        div(class="sec-label","Valores del grafico"),
        DTOutput("tabla_graf")
      ),

      # TAB 4: Clausulas
      tabPanel("04 | Clausulas",
        div(class="sec-label","Distribucion de clausulas entonativas"),
        fluidRow(
          column(3, uiOutput("stat_total")),
          column(3, uiOutput("stat_asc")),
          column(3, uiOutput("stat_desc")),
          column(3, uiOutput("stat_plan"))
        ),
        br(),
        div(class="sec-label","Detalle por clausula"),
        DTOutput("tabla_clausulas"),
        br(),
        div(class="sec-label","Comparacion slope ST vs slope %"),
        DTOutput("tabla_comp")
      ),

      # TAB 5: Transiciones
      tabPanel("05 | Transiciones",
        div(class="info-box", HTML(
          "Reajuste: diferencia F0_ini(n) - F0_fin(n-1) entre grupos consecutivos (inicio del grupo actual menos fin del grupo anterior).<br>
           Pausa: duracion registrada al final del grupo anterior.<br>
           Sub-pestana A (Por clausula): distribucion por clausula entonativa (ascendente/descendente/plano).<br>
           Sub-pestana B (Global): media general de todos los grupos + correlacion reajuste ~ pausa.")),
        tabsetPanel(id="tabs_trans", type="pills",

          tabPanel("Por clausula", br(),
            fluidRow(
              column(6, div(class="sec-label","Reajuste F0 (Hz)"),   plotlyOutput("box_hz_g",  height="280px")),
              column(6, div(class="sec-label","Reajuste F0 (%)"),    plotlyOutput("box_pct_g", height="280px"))
            ),
            fluidRow(
              column(6, div(class="sec-label","Reajuste F0 (ST)"),   plotlyOutput("box_st_g",  height="280px")),
              column(6, div(class="sec-label","Pausa anterior (s)"), plotlyOutput("box_pau_g", height="280px"))
            ),
            br(), div(class="sec-label","Estadisticos por clausula"),
            DTOutput("tabla_trans_g"),
            br(),
            uiOutput("glosario_variables")
          ),

          tabPanel("Global", br(),
            div(class="sec-label","Media general -- reajuste y pausa (todos los grupos)"),
            DTOutput("tabla_global"),
            br(),
            div(class="sec-label","Correlacion reajuste ~ pausa"),
            DTOutput("tabla_corr"),
            br(),
            fluidRow(
              column(8, div(class="sec-label","Dispersion reajuste ST ~ pausa anterior"),
                plotlyOutput("plot_corr", height="340px"))
            )
          )
        )
      ),

      # TAB 6: Chi cuadrado
      tabPanel("06 | Chi cuadrado",
        uiOutput("tab6_infobox"),
        fluidRow(
          column(3, uiOutput("chi_stat_chi2")),
          column(3, uiOutput("chi_stat_p")),
          column(3, uiOutput("chi_stat_v")),
          column(3, uiOutput("chi_stat_mix"))
        ),
        br(),
        uiOutput("tab6_label_bar"),
        plotlyOutput("barplot_chi", height="500px"),
        br(),
        div(class="sec-label","Tabla de contingencia (conteos)"),
        DTOutput("tabla_chi"),
        br(),
        uiOutput("tab6_label_table2"),
        DTOutput("tabla_mixtas")
      )
    )
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {
  `%||%` <- function(a, b) if (!is.null(a) && !is.na(a)) a else b

  datos_largos <- reactiveVal(NULL)

  # Load file into textarea
  observeEvent(input$file_upload, {
    txt <- readLines(input$file_upload$datapath, warn=FALSE) |> paste(collapse="\n")
    updateTextAreaInput(session, "raw_text", value=txt)
  })

  # Parse
  observeEvent(input$btn_parse, {
    txt <- input$raw_text
    if (is.null(txt) || nchar(trimws(txt)) == 0) {
      showNotification("No hay datos.", type="warning"); return()
    }
    sep <- if (grepl("\t", txt)) "\t" else ";"
    df  <- tryCatch(
      read_delim(I(txt), delim=sep, show_col_types=FALSE,
                 col_types=cols(.default=col_character())),
      error=function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
      showNotification("Error al leer datos.", type="error"); return()
    }
    largo <- parse_wide(df)
    if (is.null(largo)) {
      showNotification("No se encontraron grupos validos.", type="error"); return()
    }
    datos_largos(largo)
    updateTabsetPanel(session, "tabs", selected="02 | Tabla")
    showNotification(sprintf("OK: %d grupos tonales.", nrow(largo)), type="message")
  })

  observeEvent(input$btn_clear, {
    updateTextAreaInput(session, "raw_text", value="")
    datos_largos(NULL)
  })

  # Computed result
  resultado <- reactive({
    req(datos_largos())
    calcular(datos_largos(), input$ventana, input$umbral_pct,
             input$umbral_st, input$umbral_local,
             input$umbral_prog_st, input$umbral_prog_pct,
             consecutivo = isTRUE(input$consecutivo))
  })

  # ---------------------------------------------------------------------------
  # Tab 2: Tabla
  # ---------------------------------------------------------------------------
  output$tabla_larga <- renderDT({
    req(resultado())
    resultado() |>
      select(enunciado, grupo, f0_ini, f0_fin, pausa,
             inflexion_pct, inflexion_st,
             reajuste_pct, reajuste_st,
             acum_st, progresion_st,
             dir_grupo, delta_st, dir_local, slope_st_v, dir_slope_st,
             clausula_slope_st) |>
      mutate(across(where(is.double), \(x) round(x, 3))) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,"F0 ini"=f0_ini,"F0 fin"=f0_fin,
             "Pausa"=pausa,
             "Inflexion %"=inflexion_pct,"Inflexion ST"=inflexion_st,
             "Reajuste %"=reajuste_pct,"Reajuste ST"=reajuste_st,
             "Acum. ST"=acum_st,"Progresion"=progresion_st,
             "Dir. propio"=dir_grupo,"dST"=delta_st,
             "Dir. local"=dir_local,"Slope ST"=slope_st_v,
             "Dir. slope"=dir_slope_st,"Clausula"=clausula_slope_st) |>
      datatable(options=list(pageLength=20, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  # ---------------------------------------------------------------------------
  # Tab 3: Grafico
  # ---------------------------------------------------------------------------
  n_grupos_total <- reactive({ req(resultado()); nrow(resultado()) })

  observeEvent(n_grupos_total(), {
    n <- n_grupos_total()
    updateNumericInput(session, "graf_desde", max=n)
    updateNumericInput(session, "graf_hasta", value=min(5,n), max=n)
  })
  observeEvent(input$btn_todo, {
    n <- n_grupos_total()
    updateNumericInput(session, "graf_desde", value=1)
    updateNumericInput(session, "graf_hasta", value=n)
  })
  observeEvent(input$btn_primero5, {
    n <- n_grupos_total()
    updateNumericInput(session, "graf_desde", value=1)
    updateNumericInput(session, "graf_hasta", value=min(5,n))
  })
  output$lbl_n_grupos <- renderUI({
    req(n_grupos_total())
    div(style="font-family:monospace;font-size:11px;color:#8b9eff;padding-bottom:6px;",
        paste0("/ ", n_grupos_total(), " grupos"))
  })

  resultado_graf <- reactive({
    req(resultado())
    desde <- max(1,   input$graf_desde %||% 1)
    hasta <- min(nrow(resultado()), input$graf_hasta %||% 5)
    if (desde > hasta) hasta <- min(5, nrow(resultado()))
    resultado()[desde:hasta, ] |> mutate(x = row_number())
  })

  pl_config <- function(p, filename = "corolario") {
    p |> config(
      displayModeBar = TRUE,
      modeBarButtonsToRemove = list("select2d","lasso2d","autoScale2d"),
      toImageButtonOptions = list(format="png", filename=filename,
                                  width=1400, height=700, scale=2),
      locale = "es"
    )
  }

  output$grafico <- renderPlotly({
    req(resultado_graf())
    nivel <- input$nivel_slope %||% "grupo"
    df    <- resultado_graf()
    if (nivel == "enunciado") {
      df <- df |> mutate(clausula = clausula_enun_st, dir_slope_st = dir_enun_st)
    }
    p <- hacer_grafico(df, input$ventana, input$umbral_st, input$umbral_local)
    ggplotly(p, tooltip = c("x","y","colour")) |> pl_config("clausulas_f0")
  })

  output$tabla_graf <- renderDT({
    req(resultado_graf())
    resultado_graf() |>
      select(enunciado, grupo, f0_ini, f0_fin, pausa,
             reajuste_hz, reajuste_pct, reajuste_st,
             clausula_slope_st, dir_slope_st) |>
      mutate(across(where(is.double), \(x) round(x, 3))) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,"F0 ini"=f0_ini,"F0 fin"=f0_fin,
             "Pausa"=pausa,"Reajuste Hz"=reajuste_hz,"Reajuste %"=reajuste_pct,
             "Reajuste ST"=reajuste_st,"Clausula"=clausula_slope_st,"Dir."=dir_slope_st) |>
      datatable(options=list(pageLength=20, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  # ---------------------------------------------------------------------------
  # Tab 4: Clausulas
  # ---------------------------------------------------------------------------
  resumen <- reactive({
    req(resultado())
    resultado() |>
      group_by(clausula_slope_st, dir_slope_st) |>
      summarise(n_grupos=n(), enunciados=paste(unique(enunciado),collapse=", "),
                grupos=paste(paste0(enunciado,"-g",grupo),collapse=", "),
                f0_media=round(mean(f0_ini,na.rm=TRUE),1),
                slope_medio=round(mean(slope_st_v,na.rm=TRUE),3), .groups="drop")
  })

  mk_stat <- function(val, lbl, cls="stat-card") {
    div(class=cls, div(class="stat-val",val), div(class="stat-lbl",lbl))
  }
  output$stat_total <- renderUI({ req(resumen()); mk_stat(n_distinct(resumen()$clausula_slope_st),"clausulas totales") })
  output$stat_asc   <- renderUI({ req(resumen()); mk_stat(sum(resumen()$dir_slope_st=="ascendente",na.rm=TRUE),"ascendentes","stat-card stat-asc") })
  output$stat_desc  <- renderUI({ req(resumen()); mk_stat(sum(resumen()$dir_slope_st=="descendente",na.rm=TRUE),"descendentes","stat-card stat-desc") })
  output$stat_plan  <- renderUI({ req(resumen()); mk_stat(sum(resumen()$dir_slope_st=="plano",na.rm=TRUE),"planas","stat-card stat-plan") })

  output$tabla_clausulas <- renderDT({
    req(resumen())
    resumen() |>
      rename("Clausula"=clausula_slope_st,"Direccion"=dir_slope_st,"N grupos"=n_grupos,
             "Enunciados"=enunciados,"Grupos"=grupos,"F0 media"=f0_media,"Slope ST"=slope_medio) |>
      datatable(options=list(pageLength=15, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  output$tabla_comp <- renderDT({
    req(resultado())
    resultado() |>
      select(enunciado, grupo, clausula_slope_st, dir_slope_st,
             clausula_slope_pct, dir_slope_pct, slope_st_v, slope_pct_v) |>
      mutate(across(where(is.double),\(x) round(x,3)),
             coincide = dir_slope_st == dir_slope_pct) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,
             "Clausula ST"=clausula_slope_st,"Dir. ST"=dir_slope_st,
             "Clausula %"=clausula_slope_pct,"Dir. %"=dir_slope_pct,
             "Slope ST"=slope_st_v,"Slope %"=slope_pct_v,"Coincide"=coincide) |>
      datatable(options=list(pageLength=15, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  # ---------------------------------------------------------------------------
  # Tab 5: Transiciones
  # ---------------------------------------------------------------------------
  PAL_DIR <- c(ascendente="#27ae60", descendente="#e74c3c",
               plano="#f39c12", "sin dir."="grey65")

  # Base transitions dataset
  trans_base <- reactive({
    req(resultado())
    resultado() |>
      mutate(
        trans_hz  = f0_ini - lag(f0_fin),
        trans_pct = (f0_ini - lag(f0_fin)) / lag(f0_fin) * 100,
        trans_st  = hz_a_st(lag(f0_fin), f0_ini),
        pausa_ant = lag(pausa)
      ) |>
      filter(!is.na(trans_hz)) |>
      mutate(
        dir = coalesce(dir_slope_st, "sin dir."),
        cl  = clausula_slope_st
      )
  })

  # Shared boxplot builder
  box_plot <- function(df, x_var, y_var, ylab, show_hline = TRUE) {
    req(nrow(df) > 0)
    p <- ggplot(df, aes(x = .data[[x_var]], y = .data[[y_var]], fill = .data[[x_var]]))
    if (show_hline)
      p <- p + geom_hline(yintercept=0, linetype="dashed", color="grey55", linewidth=0.45)
    p +
      geom_boxplot(alpha=0.75, outlier.shape=21, outlier.size=1.8,
                   outlier.fill="white", width=0.55) +
      stat_summary(fun=mean, geom="point", shape=23, size=3,
                   fill="white", color="#1a1a1a") +
      scale_fill_manual(values=PAL_DIR, drop=FALSE) +
      labs(x=NULL, y=ylab, fill=NULL) +
      theme_minimal(base_size=12) +
      theme(panel.grid.minor=element_blank(),
            panel.grid.major.x=element_blank(),
            legend.position="none",
            axis.text.x=element_text(size=9, angle=30, hjust=1),
            plot.background=element_rect(fill="white", color=NA))
  }

  # Stats table helper
  stats_tbl <- function(df, grp_col) {
    df |>
      group_by(across(all_of(grp_col))) |>
      summarise(
        N             = n(),
        Hz_media      = round(mean(trans_hz,   na.rm=TRUE), 2),
        Hz_SD         = round(sd(trans_hz,     na.rm=TRUE), 2),
        Hz_mediana    = round(median(trans_hz, na.rm=TRUE), 2),
        Pct_media     = round(mean(trans_pct,  na.rm=TRUE), 2),
        Pct_SD        = round(sd(trans_pct,    na.rm=TRUE), 2),
        Pct_mediana   = round(median(trans_pct,na.rm=TRUE), 2),
        ST_media      = round(mean(trans_st,   na.rm=TRUE), 3),
        ST_SD         = round(sd(trans_st,     na.rm=TRUE), 3),
        ST_mediana    = round(median(trans_st, na.rm=TRUE), 3),
        Pausa_media   = round(mean(pausa_ant,  na.rm=TRUE), 3),
        Pausa_SD      = round(sd(pausa_ant,    na.rm=TRUE), 3),
        Pausa_mediana = round(median(pausa_ant,na.rm=TRUE), 3),
        .groups = "drop"
      ) |>
      datatable(options=list(pageLength=15, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  }

  # -- Sub-tab A: Por clausula (x = direction / clausula type) ----------------
  as_plotly_box <- function(p, fname) ggplotly(p, tooltip=c("y","fill")) |> pl_config(fname)

  output$box_hz_g  <- renderPlotly({ req(trans_base()); as_plotly_box(box_plot(trans_base(),"dir","trans_hz","Reajuste F0 (Hz)"),  "reajuste_hz") })
  output$box_pct_g <- renderPlotly({ req(trans_base()); as_plotly_box(box_plot(trans_base(),"dir","trans_pct","Reajuste F0 (%)"), "reajuste_pct") })
  output$box_st_g  <- renderPlotly({ req(trans_base()); as_plotly_box(box_plot(trans_base(),"dir","trans_st","Reajuste F0 (ST)"),  "reajuste_st") })
  output$box_pau_g <- renderPlotly({
    req(trans_base())
    df <- trans_base() |> filter(!is.na(pausa_ant))
    p  <- if (nrow(df)==0)
            ggplot()+annotate("text",x=1,y=1,label="Sin datos de pausa",size=5,color="grey60")+theme_void()
          else box_plot(df,"dir","pausa_ant","Pausa anterior (s)", show_hline=FALSE)
    as_plotly_box(p, "pausa")
  })
  output$tabla_trans_g <- renderDT({ req(trans_base()); stats_tbl(trans_base(),"dir") })

  output$glosario_variables <- renderUI({
    div(class="info-box", style="margin-top:18px;",
      HTML("<b>Glosario de variables</b><br><br>
      <b>Inflexion % / ST</b>: movimiento interno del grupo tonal [F0_fin(n) - F0_ini(n)], expresado en porcentaje y en semitonos. Mide si el contorno propio del grupo sube, baja o es plano.<br><br>
      <b>Reajuste F0 (Hz / % / ST)</b>: salto de F0 en la juntura entre grupos [F0_ini(n) - F0_fin(n-1)]. Valores positivos = subida en la transicion; negativos = bajada. El ST es perceptivamente uniforme (1 ST ≈ umbral de percepcion tonal).<br><br>
      <b>Acum. ST</b>: acumulacion de reajuste_st desde el ultimo cambio de direccion. Se reinicia cuando la direccion se invierte. Cuando el valor absoluto supera el umbral configurado (Prog. umbral ST), el grupo se etiqueta como progresion.<br><br>
      <b>Progresion</b>: etiqueta de progresion acumulada. <i>prog_asc</i> = la acumulacion de reajustes ascendentes supera el umbral; <i>prog_desc</i> = idem en sentido descendente; <i>neutro</i> = todavia no ha alcanzado el umbral o el movimiento es mixto.<br><br>
      <b>Pausa anterior (s)</b>: duracion en segundos de la pausa registrada al final del grupo anterior (pausa previa a la transicion analizada).<br><br>
      <b>Dir. slope (clausula)</b>: direccion entonativa asignada al grupo segun el slope de regresion sobre F0_ini en la ventana deslizante. Puede ser <i>ascendente</i>, <i>descendente</i> o <i>plano</i>.<br><br>
      <b>N</b>: numero de transiciones entre grupos tonales disponibles en la categoria.")
    )
  })

  # -- Sub-tab B: Global (medias generales + correlacion) ---------------------
  output$tabla_global <- renderDT({
    req(trans_base(), resultado())
    df   <- trans_base()
    res  <- resultado()
    # Pausa = 0 counts (original pausa column, all groups)
    n_pausa_cero    <- sum(!is.na(res$pausa) & res$pausa == 0)
    n_pausa_total   <- sum(!is.na(res$pausa))
    # Transitions where preceding pause = 0
    n_trans_p0      <- sum(!is.na(df$pausa_ant) & df$pausa_ant == 0)
    n_trans_total   <- nrow(df)

    df_pos <- df |> filter(!is.na(pausa_ant) & pausa_ant > 0)

    bind_rows(
      tibble(
        Concepto = "Grupos con pausa = 0",
        N = n_pausa_cero,
        `% del total` = round(n_pausa_cero / n_pausa_total * 100, 1),
        Media = NA_real_, SD = NA_real_, Mediana = NA_real_
      ),
      tibble(
        Concepto = "Transiciones con pausa_ant = 0",
        N = n_trans_p0,
        `% del total` = round(n_trans_p0 / n_trans_total * 100, 1),
        Media = NA_real_, SD = NA_real_, Mediana = NA_real_
      ),
      tibble(
        Concepto        = c("Reajuste Hz (pausa>0)","Reajuste % (pausa>0)","Reajuste ST (pausa>0)","Pausa anterior (s)"),
        N               = c(sum(!is.na(df_pos$trans_hz)), sum(!is.na(df_pos$trans_pct)),
                            sum(!is.na(df_pos$trans_st)), sum(!is.na(df_pos$pausa_ant))),
        `% del total`   = NA_real_,
        Media           = c(round(mean(df_pos$trans_hz,  na.rm=TRUE),2),
                            round(mean(df_pos$trans_pct, na.rm=TRUE),2),
                            round(mean(df_pos$trans_st,  na.rm=TRUE),3),
                            round(mean(df_pos$pausa_ant, na.rm=TRUE),3)),
        SD              = c(round(sd(df_pos$trans_hz,  na.rm=TRUE),2),
                            round(sd(df_pos$trans_pct, na.rm=TRUE),2),
                            round(sd(df_pos$trans_st,  na.rm=TRUE),3),
                            round(sd(df_pos$pausa_ant, na.rm=TRUE),3)),
        Mediana         = c(round(median(df_pos$trans_hz,  na.rm=TRUE),2),
                            round(median(df_pos$trans_pct, na.rm=TRUE),2),
                            round(median(df_pos$trans_st,  na.rm=TRUE),3),
                            round(median(df_pos$pausa_ant, na.rm=TRUE),3))
      )
    ) |>
    datatable(options=list(pageLength=10, scrollX=TRUE, dom="t"),
      rownames=FALSE, class="stripe hover")
  })

  output$tabla_corr <- renderDT({
    req(trans_base())
    # Correlacion solo sobre transiciones con pausa > 0
    df <- trans_base() |> filter(!is.na(trans_hz), !is.na(pausa_ant), pausa_ant > 0)
    n_excl <- nrow(trans_base() |> filter(!is.na(trans_hz))) - nrow(df)
    if (nrow(df) < 3) {
      return(datatable(
        tibble(Nota=paste0("Insuficientes datos con pausa > 0 (excluidos ", n_excl, " con pausa = 0).")),
        rownames=FALSE))
    }
    ct_hz  <- cor.test(df$trans_hz,  df$pausa_ant, method="pearson")
    ct_pct <- cor.test(df$trans_pct, df$pausa_ant, method="pearson")
    ct_st  <- cor.test(df$trans_st,  df$pausa_ant, method="pearson")
    cs_st  <- cor.test(df$trans_st,  df$pausa_ant, method="spearman", exact=FALSE)
    tibble(
      Reajuste  = c("Hz","Pct","ST (Pearson)","ST (Spearman)"),
      r         = round(c(ct_hz$estimate, ct_pct$estimate, ct_st$estimate, cs_st$estimate), 3),
      p_valor   = round(c(ct_hz$p.value,  ct_pct$p.value,  ct_st$p.value,  cs_st$p.value), 4),
      IC_95_inf = round(c(ct_hz$conf.int[1], ct_pct$conf.int[1], ct_st$conf.int[1], NA_real_), 3),
      IC_95_sup = round(c(ct_hz$conf.int[2], ct_pct$conf.int[2], ct_st$conf.int[2], NA_real_), 3),
      N         = nrow(df),
      Excluidos_pausa0 = n_excl
    ) |>
    datatable(options=list(pageLength=10, scrollX=TRUE, dom="t"),
      rownames=FALSE, class="stripe hover")
  })

  output$plot_corr <- renderPlotly({
    req(trans_base())
    df <- trans_base() |> filter(!is.na(trans_st), !is.na(pausa_ant), pausa_ant > 0)
    if (nrow(df) < 3) {
      p <- ggplot() +
        annotate("text",x=1,y=1,label="Sin datos suficientes con pausa > 0",size=5,color="grey60") +
        theme_void()
      return(ggplotly(p) |> pl_config("correlacion"))
    }
    r_val <- round(cor(df$trans_st, df$pausa_ant, use="complete.obs"), 3)
    p <- ggplot(df, aes(x=pausa_ant, y=trans_st, color=dir)) +
      geom_hline(yintercept=0, linetype="dashed", color="grey70", linewidth=0.4) +
      geom_point(size=2.5, alpha=0.8) +
      geom_smooth(method="lm", se=TRUE, color="#1a1a1a", fill="#1a1a1a", alpha=0.1,
                  linewidth=0.9, inherit.aes=FALSE, aes(x=pausa_ant, y=trans_st)) +
      scale_color_manual(values=PAL_DIR) +
      labs(x="Pausa anterior (s)", y="Reajuste F0 (ST)",
           title=paste0("Reajuste ST ~ Pausa (solo pausa > 0)  |  r = ", r_val),
           color="Clausula") +
      theme_minimal(base_size=12) +
      theme(panel.grid.minor=element_blank(),
            legend.position="bottom",
            plot.title=element_text(face="bold", size=11),
            plot.background=element_rect(fill="white", color=NA))
    ggplotly(p, tooltip=c("x","y","colour")) |> pl_config("correlacion_reajuste_pausa")
  })

  # ---------------------------------------------------------------------------
  # Tab 6: Chi cuadrado
  # ---------------------------------------------------------------------------

  # Dynamic section labels and info box
  output$tab6_infobox <- renderUI({
    if (isTRUE(input$consecutivo)) {
      div(class = "info-box", HTML(
        "Cada clausula se clasifica en tres categorias:<br>
         <b>Pura</b>: todos sus grupos pertenecen al mismo enunciado.<br>
         <b>Mixta — cruce en frontera</b>: contiene el grupo final de un enunciado y el grupo inicial del siguiente.<br>
         <b>Mixta — sin cruce directo</b>: agrupa grupos de enunciados distintos sin estar en la frontera inmediata.<br>
         El test chi cuadrado evalua la asociacion global entre enunciados y clausulas."))
    } else {
      div(class = "info-box", HTML(
        "Modo no consecutivo: las clausulas estan contenidas dentro de cada enunciado.<br>
         Los enunciados se clasifican segun cuantas clausulas tonales integran en su interior:<br>
         <b>1 clausula</b>: correspondencia plena entre enunciado y clausula entonativa.<br>
         <b>2 clausulas</b>: el enunciado se divide en dos clausulas de distinta direccion.<br>
         <b>3+ clausulas</b>: el enunciado presenta tres o mas clausulas internas.<br>
         El test chi cuadrado evalua si la distribucion de direcciones por enunciado es homogenea."))
    }
  })
  output$tab6_label_bar <- renderUI({
    lbl <- if (isTRUE(input$consecutivo))
      "Clausulas segun numero de enunciados que las componen"
    else
      "Enunciados segun numero de clausulas en su interior"
    div(class = "sec-label", lbl)
  })
  output$tab6_label_table2 <- renderUI({
    lbl <- if (isTRUE(input$consecutivo))
      "Clausulas mixtas (grupos de mas de un enunciado)"
    else
      "Enunciados con varias clausulas internas"
    div(class = "sec-label", lbl)
  })

  # NON-CONSECUTIVE: utterances classified by number of internal clauses
  cat_data_noc <- reactive({
    req(resultado())
    res <- resultado()

    enun_detail <- res |>
      group_by(enunciado) |>
      summarise(
        n_clausulas = n_distinct(clausula_slope_st),
        clausulas   = paste(unique(clausula_slope_st), collapse = ", "),
        dirs        = paste(na.omit(unique(dir_slope_st)), collapse = ", "),
        n_grupos    = n(),
        .groups     = "drop"
      ) |>
      mutate(categoria = case_when(
        n_clausulas == 1 ~ "1 clausula",
        n_clausulas == 2 ~ "2 clausulas",
        TRUE             ~ "3+ clausulas"
      ))

    cat_summary <- enun_detail |>
      count(categoria, name = "n_enun") |>
      right_join(tibble(categoria = c("1 clausula", "2 clausulas", "3+ clausulas")),
                 by = "categoria") |>
      mutate(
        n_enun    = coalesce(n_enun, 0L),
        categoria = factor(categoria, levels = c("1 clausula", "2 clausulas", "3+ clausulas")),
        pct       = n_enun / sum(n_enun) * 100,
        label     = paste0(n_enun, "\n(", round(pct, 1), "%)")
      ) |>
      arrange(categoria)

    list(summary = cat_summary, detail = enun_detail)
  })

  # NON-CONSECUTIVE chi-square: utterance × clause-direction contingency
  chi_data_noc <- reactive({
    req(resultado())
    res <- resultado()
    tbl <- res |>
      filter(!is.na(dir_slope_st)) |>
      count(enunciado, dir_slope_st) |>
      pivot_wider(names_from = dir_slope_st, values_from = n, values_fill = 0L)
    if (nrow(tbl) < 2) return(NULL)
    cont <- as.matrix(tbl[, -1, drop = FALSE])
    rownames(cont) <- tbl$enunciado
    cont
  })

  chi_test_noc <- reactive({
    m <- chi_data_noc()
    if (is.null(m) || any(dim(m) < 2)) return(NULL)
    tryCatch(chisq.test(m, simulate.p.value = TRUE, B = 5000), error = function(e) NULL)
  })

  # Reactive: ordenar clausulas numericamente para que el heatmap sea coherente
  chi_data <- reactive({
    req(resultado())
    res <- resultado()

    # Orden numerico de clausulas
    cl_all    <- unique(res$clausula_slope_st)
    cl_nums   <- suppressWarnings(as.integer(gsub("clausula", "", cl_all)))
    cl_ord    <- cl_all[order(cl_nums)]

    res <- res |>
      mutate(clausula_f = factor(clausula_slope_st, levels = cl_ord))

    cont <- table(Enunciado = res$enunciado, Clausula = res$clausula_f)
    list(res = res, cont = cont, cl_ord = cl_ord)
  })

  chi_test <- reactive({
    req(chi_data())
    ct <- chi_data()$cont
    if (any(dim(ct) < 2)) return(NULL)
    tryCatch(
      chisq.test(ct, simulate.p.value = TRUE, B = 5000),
      error = function(e) NULL
    )
  })

  # Stat cards
  mk_stat_chi <- function(val, lbl, cls = "stat-card") {
    div(class = cls,
        div(class = "stat-val", val),
        div(class = "stat-lbl", lbl))
  }

  output$chi_stat_chi2 <- renderUI({
    tst <- if (isTRUE(input$consecutivo)) chi_test() else chi_test_noc()
    req(tst)
    mk_stat_chi(round(tst$statistic, 2), "chi cuadrado")
  })
  output$chi_stat_p <- renderUI({
    tst <- if (isTRUE(input$consecutivo)) chi_test() else chi_test_noc()
    req(tst)
    p   <- tst$p.value
    lbl <- if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 4))
    cls <- if (p < 0.05) "stat-card stat-asc" else "stat-card stat-plan"
    mk_stat_chi(lbl, "p-valor (Monte Carlo)", cls)
  })
  output$chi_stat_v <- renderUI({
    tst <- if (isTRUE(input$consecutivo)) chi_test() else chi_test_noc()
    req(tst)
    if (isTRUE(input$consecutivo)) {
      ct <- chi_data()$cont; n <- sum(ct); k <- min(nrow(ct), ncol(ct)) - 1
    } else {
      ct <- chi_data_noc();  n <- sum(ct); k <- min(nrow(ct), ncol(ct)) - 1
    }
    v <- if (!is.null(ct) && k > 0) round(sqrt(tst$statistic / (n * k)), 3) else NA
    mk_stat_chi(v, "V de Cramer")
  })
  output$chi_stat_mix <- renderUI({
    if (isTRUE(input$consecutivo)) {
      req(chi_data())
      n_mix <- chi_data()$res |>
        group_by(clausula_slope_st) |>
        summarise(nd = n_distinct(enunciado), .groups = "drop") |>
        filter(nd > 1) |> nrow()
      mk_stat_chi(n_mix, "clausulas mixtas", "stat-card stat-desc")
    } else {
      req(cat_data_noc())
      det  <- cat_data_noc()$detail
      n_m  <- sum(det$n_clausulas > 1)
      mk_stat_chi(n_m, "enunciados con 2+ clausulas", "stat-card stat-desc")
    }
  })

  # Categorias de clausulas segun numero de enunciados y tipo de cruce
  cat_data <- reactive({
    req(resultado())
    res <- resultado() |> arrange(x)

    # Limites globales por enunciado
    enun_max <- res |> group_by(enunciado) |>
      summarise(max_g = max(grupo), .groups = "drop")
    enun_min <- res |> group_by(enunciado) |>
      summarise(min_g = min(grupo), .groups = "drop")

    # Clausulas donde un par consecutivo cruza exactamente la frontera:
    # grupo_i == max de su enunciado  Y  grupo_{i+1} == min del enunciado siguiente
    cls_frontera <- res |>
      mutate(
        next_enun  = lead(enunciado),
        next_grupo = lead(grupo),
        next_cl    = lead(clausula_slope_st)
      ) |>
      filter(!is.na(next_cl),
             clausula_slope_st == next_cl,
             enunciado != next_enun) |>
      left_join(enun_max, by = "enunciado") |>
      left_join(enun_min |> rename(next_enun = enunciado, min_g_next = min_g),
                by = "next_enun") |>
      filter(grupo == max_g, next_grupo == min_g_next) |>
      pull(clausula_slope_st) |>
      unique()

    res |>
      group_by(clausula_slope_st) |>
      summarise(n_enun = n_distinct(enunciado), .groups = "drop") |>
      mutate(categoria = case_when(
        n_enun == 1                                 ~ "Pura\n(1 enunciado)",
        clausula_slope_st %in% cls_frontera         ~ "Mixta\ncruce en frontera",
        TRUE                                        ~ "Mixta\nsin cruce directo"
      )) |>
      count(categoria, name = "n_clausulas") |>
      mutate(
        categoria = factor(categoria, levels = c(
          "Pura\n(1 enunciado)",
          "Mixta\ncruce en frontera",
          "Mixta\nsin cruce directo"
        )),
        pct   = n_clausulas / sum(n_clausulas) * 100,
        label = paste0(n_clausulas, "\n(", round(pct, 1), "%)")
      )
  })

  output$barplot_chi <- renderPlotly({
    if (isTRUE(input$consecutivo)) {
      # ---- consecutive: clauses by utterance-span category ----
      req(cat_data())
      df <- cat_data()
      pal <- c(
        "Pura\n(1 enunciado)"       = "#4a4a4a",
        "Mixta\ncruce en frontera"  = "#d68910",
        "Mixta\nsin cruce directo"  = "#c0392b"
      )
      n_total <- sum(df$n_clausulas)
      n_mix   <- sum(df$n_clausulas[df$categoria != "Pura\n(1 enunciado)"])
      subtt   <- paste0("Total clausulas: ", n_total,
                        "  |  Clausulas mixtas (2+ enunciados): ", n_mix,
                        " (", round(n_mix / n_total * 100, 1), "%)")
      ylab  <- "% de clausulas"
      title <- "Distribucion de clausulas segun enunciados que las componen"
      fname <- "clausulas_chi"
    } else {
      # ---- non-consecutive: utterances by clause-count category ----
      req(cat_data_noc())
      df  <- cat_data_noc()$summary
      pal <- c(
        "1 clausula"   = "#4a4a4a",
        "2 clausulas"  = "#d68910",
        "3+ clausulas" = "#c0392b"
      )
      n_total <- sum(df$n_enun)
      n_multi <- sum(df$n_enun[df$categoria != "1 clausula"])
      subtt   <- paste0("Total enunciados: ", n_total,
                        "  |  Enunciados con varias clausulas: ", n_multi,
                        " (", round(n_multi / n_total * 100, 1), "%)")
      # Rename n_enun → n_clausulas so the ggplot aes works with a single block
      df    <- df |> rename(n_clausulas = n_enun)
      ylab  <- "% de enunciados"
      title <- "Distribucion de enunciados segun clausulas en su interior"
      fname <- "enunciados_clausulas"
    }

    p <- ggplot(df, aes(x = categoria, y = pct, fill = categoria)) +
      geom_col(width = 0.55, color = "white", linewidth = 0.4) +
      geom_text(aes(y = pct / 2, label = label),
                size = 4, fontface = "bold", color = "white",
                lineheight = 0.9) +
      scale_fill_manual(values = pal, drop = FALSE) +
      scale_y_continuous(labels = function(x) paste0(x, "%"),
                         expand = expansion(mult = c(0, 0.18)),
                         limits = c(0, NA)) +
      labs(title = title, subtitle = subtt, x = NULL, y = ylab) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(color = "grey50", size = 10),
        axis.text.x   = element_text(size = 11),
        legend.position = "none",
        plot.background = element_rect(fill = "white", color = NA)
      )
    ggplotly(p, tooltip = c("x", "y")) |> pl_config(fname)
  })

  # Tabla de contingencia como DT
  output$tabla_chi <- renderDT({
    if (isTRUE(input$consecutivo)) {
      req(chi_data())
      df <- as.data.frame.matrix(chi_data()$cont)
    } else {
      m <- chi_data_noc()
      req(!is.null(m))
      df <- as.data.frame.matrix(m)
    }
    datatable(df,
      options = list(pageLength = 20, scrollX = TRUE, dom = "lrtip",
        language = list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
      class = "stripe hover")
  })

  # Clausulas mixtas (consecutive) / enunciados con varias clausulas (non-consecutive)
  output$tabla_mixtas <- renderDT({
    req(resultado())

    if (!isTRUE(input$consecutivo)) {
      req(cat_data_noc())
      det <- cat_data_noc()$detail |>
        filter(n_clausulas > 1) |>
        arrange(desc(n_clausulas), enunciado)
      if (nrow(det) == 0) {
        return(datatable(
          tibble(Nota = "Todos los enunciados tienen una sola clausula interna."),
          rownames = FALSE))
      }
      return(
        det |>
          rename("Enunciado"   = enunciado,
                 "N clausulas" = n_clausulas,
                 "Clausulas"   = clausulas,
                 "Direcciones" = dirs,
                 "N grupos"    = n_grupos) |>
          datatable(options = list(pageLength = 15, scrollX = TRUE, dom = "lrtip",
            language = list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
            rownames = FALSE, class = "stripe hover")
      )
    }

    res <- resultado() |> arrange(x)

    enun_max <- res |> group_by(enunciado) |>
      summarise(max_g = max(grupo), .groups = "drop")
    enun_min <- res |> group_by(enunciado) |>
      summarise(min_g = min(grupo), .groups = "drop")

    cls_frontera <- res |>
      mutate(next_enun  = lead(enunciado),
             next_grupo = lead(grupo),
             next_cl    = lead(clausula_slope_st)) |>
      filter(!is.na(next_cl), clausula_slope_st == next_cl, enunciado != next_enun) |>
      left_join(enun_max, by = "enunciado") |>
      left_join(enun_min |> rename(next_enun = enunciado, min_g_next = min_g),
                by = "next_enun") |>
      filter(grupo == max_g, next_grupo == min_g_next) |>
      pull(clausula_slope_st) |> unique()

    mix <- res |>
      group_by(clausula_slope_st, dir_slope_st) |>
      summarise(
        n_enunciados = n_distinct(enunciado),
        enunciados   = paste(sort(unique(enunciado)), collapse = ", "),
        n_grupos     = n(),
        .groups = "drop"
      ) |>
      filter(n_enunciados > 1) |>
      mutate(tipo = ifelse(clausula_slope_st %in% cls_frontera,
                           "cruce en frontera", "sin cruce directo")) |>
      arrange(clausula_slope_st)

    if (nrow(mix) == 0) {
      return(datatable(
        tibble(Nota = "No hay clausulas mixtas: cada clausula pertenece a un unico enunciado."),
        rownames = FALSE))
    }

    mix |>
      rename("Clausula"     = clausula_slope_st,
             "Direccion"    = dir_slope_st,
             "N enunciados" = n_enunciados,
             "Enunciados"   = enunciados,
             "N grupos"     = n_grupos,
             "Tipo"         = tipo) |>
      datatable(options = list(pageLength = 15, scrollX = TRUE, dom = "lrtip",
        language = list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames = FALSE, class = "stripe hover")
  })

}

shinyApp(ui, server)