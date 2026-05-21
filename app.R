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

# Assign clausulaN labels using vis-style logic:
# - Reajuste >= umbral  → always new clause
# - Current is plano    → join previous clause
# - Prev was plano      → check last non-plano direction
# - Both non-plano      → new clause if directions differ
asignar_clausulas_vis <- function(tend, rej_st, rej_pct, pausa_ant,
                                   umbral_pausa = 0,
                                   umbral_rej_st = 0, umbral_rej_pct = 0) {
  n <- length(tend)
  if (n == 0) return(character(0))
  out    <- character(n)
  cnt    <- 1L
  out[1] <- "clausula1"
  for (i in seq(2, n)) {
    rs <- if (is.na(rej_st[i]))   0 else abs(rej_st[i])
    rp <- if (is.na(rej_pct[i]))  0 else abs(rej_pct[i])
    pa <- if (is.na(pausa_ant[i])) 0 else pausa_ant[i]
    # Priority 1: preceding pause >= threshold
    pausa_big <- umbral_pausa > 0 & pa >= umbral_pausa
    # Priority 2: reajuste >= threshold (only if pause didn't already fire)
    rej_big <- !pausa_big & (
      (umbral_rej_st  > 0 & rs >= umbral_rej_st) |
      (umbral_rej_pct > 0 & rp >= umbral_rej_pct)
    )
    tc <- tend[i]; tp <- tend[i - 1]
    dir_nuevo <- if (is.na(tc) || tc == "plano") {
      FALSE
    } else if (is.na(tp) || tp == "plano") {
      tb <- NA_character_
      for (j in seq(i - 1, 1)) {
        if (!is.na(tend[j]) && tend[j] != "plano") { tb <- tend[j]; break }
      }
      !is.na(tb) && tb != tc
    } else {
      tp != tc
    }
    nuevo <- pausa_big || rej_big || dir_nuevo
    if (nuevo) { cnt <- cnt + 1L; out[i] <- paste0("clausula", cnt) }
    else out[i] <- out[i - 1]
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
calcular <- function(dl, umbral_st = 0.8, umbral_prog_st = 2.0,
                     consecutivo = TRUE,
                     umbral_pausa = 0.3,
                     umbral_reajuste_st = 0, umbral_reajuste_pct = 0) {

  dl <- dl |> mutate(x = row_number())

  renumber_clauses <- function(enun, local_cl) {
    keys  <- paste0(enun, "|||", local_cl)
    ukeys <- unique(keys)
    mapping <- setNames(paste0("clausula", seq_along(ukeys)), ukeys)
    mapping[keys]
  }

  grp <- if (consecutivo) dl else group_by(dl, enunciado)

  dl2 <- grp |>
    mutate(
      inflexion_pct = (f0_fin - f0_ini) / f0_ini * 100,
      inflexion_st  = hz_a_st(f0_ini, f0_fin),
      dir_grupo     = etiquetar(inflexion_st, umbral_st),
      delta_hz      = f0_ini - lag(f0_ini),
      delta_pct     = delta_hz / lag(f0_ini) * 100,
      delta_st      = hz_a_st(lag(f0_ini), f0_ini),
      dir_local     = etiquetar(delta_st, umbral_st),
      reajuste_hz   = f0_ini - lag(f0_fin),
      reajuste_pct  = (f0_ini - lag(f0_fin)) / lag(f0_fin) * 100,
      reajuste_st   = hz_a_st(lag(f0_fin), f0_ini),
      .cl = asignar_clausulas_vis(dir_grupo, reajuste_st, reajuste_pct, lag(pausa),
                                   umbral_pausa, umbral_reajuste_st, umbral_reajuste_pct)
    ) |>
    ungroup() |>
    mutate(
      clausula_slope_st = if (consecutivo) .cl else renumber_clauses(enunciado, .cl),
      clausula          = clausula_slope_st
    ) |>
    select(-.cl)

  # Per-clause interleaved regression: corrects direction label for multi-group clauses
  cl_stats <- dl2 |>
    group_by(clausula) |>
    summarise(
      n_cl        = n(),
      slope_cl_st = if (n() >= 2) slope_st_interleaved(f0_ini, f0_fin) else NA_real_,
      .groups = "drop"
    ) |>
    mutate(dir_cl_st = etiquetar(slope_cl_st, umbral_st))

  dl3 <- dl2 |>
    left_join(cl_stats, by = "clausula") |>
    mutate(
      dir_slope_st = case_when(
        n_cl >= 2 & !is.na(dir_cl_st) ~ dir_cl_st,
        TRUE                           ~ dir_grupo
      ),
      slope_st_v = coalesce(slope_cl_st, inflexion_st)
    )

  prog <- progresion_fn(dl3$reajuste_st, 0.1, umbral_prog_st)

  dl3 |> mutate(acum_st = prog$acum, progresion_st = prog$etiq)
}

# =============================================================================
# GRAFICO
# =============================================================================
hacer_grafico <- function(res, umbral_st,
                          nombre_cl   = "clausula",
                          nombre_ais  = "singular",
                          dir_labels  = NULL,
                          tend_fin    = FALSE,
                          tend_ini    = FALSE,
                          tend_med    = FALSE,
                          modo_tend   = "global",
                          mostrar_enun = TRUE) {

  if (!is.null(dir_labels))
    res <- res |> mutate(dir_slope_st = recode(dir_slope_st, !!!dir_labels))

  # Build display labels: multi-group clauses → nombre_cl+N, single → nombre_ais+N
  cl_n     <- res |> group_by(clausula) |> summarise(n_g = n(), .groups = "drop")
  cl_multi <- cl_n |> filter(n_g >= 2) |> pull(clausula)
  all_cls  <- res |> arrange(x) |> distinct(clausula) |> pull()
  m_i <- 1L; s_i <- 1L
  cl_map <- setNames(character(length(all_cls)), all_cls)
  for (cl in all_cls) {
    if (cl %in% cl_multi) { cl_map[[cl]] <- paste0(nombre_cl,  m_i); m_i <- m_i + 1L }
    else                   { cl_map[[cl]] <- paste0(nombre_ais, s_i); s_i <- s_i + 1L }
  }
  res <- res |> mutate(cl_disp = cl_map[clausula])

  cl_dirs  <- res |> filter(!is.na(dir_slope_st)) |> distinct(cl_disp) |> pull()
  cl_nodir <- setdiff(unique(res$cl_disp), cl_dirs)
  pal <- c(setNames(hue_pal()(max(length(cl_dirs), 1)), cl_dirs),
           setNames(rep("grey65", length(cl_nodir)), cl_nodir))

  ancho     <- 0.35
  y_raw_max <- max(c(res$f0_ini, res$f0_fin), na.rm = TRUE)
  y_raw_min <- min(c(res$f0_ini, res$f0_fin), na.rm = TRUE)

  conexiones <- res |>
    mutate(x_start = x, x_end = lead(x), y_start = f0_fin, y_end = lead(f0_ini)) |>
    filter(!is.na(x_end))

  flechas <- res |>
    filter(!is.na(dir_slope_st)) |>
    group_by(cl_disp, dir_slope_st) |>
    summarise(x_ini = min(x), x_fin = max(x), .groups = "drop") |>
    arrange(x_ini) |>
    mutate(fila   = row_number() %% 2,
           y_arr  = y_raw_max * 1.07 + fila * y_raw_max * 0.05,
           y_lbl  = y_arr + y_raw_max * 0.055,
           x_ini2 = ifelse(x_ini == x_fin, x_ini - 0.25, x_ini),
           x_fin2 = ifelse(x_ini == x_fin, x_fin + 0.25, x_fin))

  p <- ggplot(res)

  # Background shading for multi-group clauses
  multi_disp <- cl_map[cl_multi]
  shade_data <- res |>
    filter(cl_disp %in% multi_disp) |>
    group_by(cl_disp) |>
    summarise(x_lo = min(x) - 0.46, x_hi = max(x) + 0.46, .groups = "drop")
  if (nrow(shade_data) > 0 && length(cl_dirs) > 0) {
    shade_cols <- setNames(hue_pal()(length(cl_dirs)), cl_dirs)
    for (i in seq_len(nrow(shade_data))) {
      cd  <- shade_data$cl_disp[i]
      col <- if (cd %in% names(shade_cols)) shade_cols[[cd]] else "grey65"
      p   <- p + annotate("rect",
               xmin = shade_data$x_lo[i], xmax = shade_data$x_hi[i],
               ymin = y_raw_min * 0.97,   ymax = y_raw_max * 1.05,
               fill = col, alpha = 0.07, color = NA)
    }
  }

  # Labels for inline group annotation and x-axis
  if (mostrar_enun) {
    lbl_inline <- paste0(res$enunciado, "\ng", res$grupo)
    lbl_xaxis  <- paste0(res$enunciado, "\ng", res$grupo)
  } else {
    lbl_inline <- paste0("g", seq_len(nrow(res)))
    lbl_xaxis  <- paste0("g", seq_len(nrow(res)))
  }

  # Utterance boundary positions (gap between last group of one utterance and first of next)
  boundaries <- if (mostrar_enun) {
    res |> mutate(next_enun = lead(enunciado)) |>
      filter(!is.na(next_enun), enunciado != next_enun) |>
      mutate(xline = x + 0.5) |> pull(xline)
  } else numeric(0)

  p <- p +
    geom_segment(data = conexiones,
      aes(x = x_start + ancho, xend = x_end - ancho, y = y_start, yend = y_end),
      color = "grey72", linewidth = 0.55, linetype = "dashed") +
    geom_segment(aes(x = x - ancho, xend = x + ancho, y = f0_ini, yend = f0_fin,
                     color = cl_disp), linewidth = 2, lineend = "round") +
    geom_point(aes(x = x - ancho, y = f0_ini, fill = cl_disp),
               shape = 21, size = 2.8, color = "white", stroke = 0.7) +
    geom_point(aes(x = x + ancho, y = f0_fin, fill = cl_disp),
               shape = 24, size = 2.8, color = "white", stroke = 0.7) +
    geom_text(aes(x = x, y = pmin(f0_ini, f0_fin), label = lbl_inline),
              vjust = 2.5, size = 2.4, color = "grey40", fontface = "italic") +
    geom_segment(data = flechas,
      aes(x = x_ini2 - ancho, xend = x_fin2 + ancho, y = y_arr, yend = y_arr,
          color = cl_disp),
      linewidth = 1.05,
      arrow = arrow(length = unit(0.2, "cm"), ends = "last", type = "closed")) +
    geom_text(data = flechas,
      aes(x = (x_ini2 + x_fin2) / 2, y = y_lbl,
          label = paste0(cl_disp, "\n(", dir_slope_st, ")"), color = cl_disp),
      size = 2.4, lineheight = 0.85, show.legend = FALSE) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values  = pal) +
    scale_x_continuous(breaks = res$x, labels = lbl_xaxis,
                       expand = expansion(mult = 0.03)) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.42))) +
    labs(title = "Progresiones de F0 por clausula",
         subtitle = sprintf("o=F0 ini  ^=F0 fin  |  Umbral dir.: %.2f ST", umbral_st),
         x = NULL, y = "F0 (Hz)", color = NULL, fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey50", size = 9),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          axis.text.x = element_text(size = 8, color = "grey40"),
          legend.position = "bottom")

  # Utterance boundary lines (overlay)
  if (length(boundaries) > 0) {
    p <- p + geom_vline(xintercept = boundaries,
                        linetype = "dashed", color = "grey50", linewidth = 0.5, alpha = 0.7)
  }

  # Regression lines
  add_smooth_g <- function(plt, x_v, y_v, col, lty) {
    ok <- !is.na(y_v)
    if (sum(ok) < 2) return(plt)
    df_t <- data.frame(x = x_v[ok], y = y_v[ok])
    plt + geom_smooth(data = df_t, aes(x = x, y = y),
                      method = "lm", formula = y ~ x, color = col, fill = col,
                      linewidth = 0.9, linetype = lty, alpha = 0.08, se = TRUE,
                      inherit.aes = FALSE)
  }

  if (tend_fin || tend_ini || tend_med) {
    if (modo_tend == "global") {
      if (tend_fin) p <- add_smooth_g(p, res$x, res$f0_fin,                    "grey40",  "dashed")
      if (tend_ini) p <- add_smooth_g(p, res$x, res$f0_ini,                    "#27ae60", "dotted")
      if (tend_med) p <- add_smooth_g(p, res$x, (res$f0_ini+res$f0_fin)/2,     "#d68910", "solid")
    } else {
      for (cl in intersect(unique(res$cl_disp), multi_disp)) {
        sub <- res |> filter(cl_disp == cl)
        if (nrow(sub) < 2) next
        if (tend_fin) p <- add_smooth_g(p, sub$x, sub$f0_fin,                  "grey40",  "dashed")
        if (tend_ini) p <- add_smooth_g(p, sub$x, sub$f0_ini,                  "#27ae60", "dotted")
        if (tend_med) p <- add_smooth_g(p, sub$x, (sub$f0_ini+sub$f0_fin)/2,   "#d68910", "solid")
      }
    }
  }

  p
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body{background:#f4f3f0;font-family:'Helvetica Neue',Arial,sans-serif;color:#1a1a1a;margin:0;padding:0 24px 40px}
    .main-title{font-family:monospace;font-weight:700;font-size:22px;letter-spacing:-0.5px;padding:24px 0 2px}
    .main-subtitle{font-size:12px;color:#777;font-family:monospace}
    .btn-ayuda{background:transparent;color:#888;font-family:monospace;font-size:12px;font-weight:700;border:1.5px solid #ccc;border-radius:50%;width:22px;height:22px;padding:0;line-height:1;cursor:pointer;display:inline-flex;align-items:center;justify-content:center;vertical-align:middle}
    .btn-ayuda:hover{border-color:#555;color:#333;background:transparent}
    .params-bar{background:#1a1a1a;border-radius:10px;padding:14px 22px;margin-bottom:16px;display:flex;align-items:flex-end;gap:18px;flex-wrap:wrap}
    .params-bar label{color:#aaa!important;font-size:11px!important;font-family:monospace;text-transform:uppercase;letter-spacing:0.5px}
    .params-bar .form-control{background:#2d2d2d!important;border:1px solid #444!important;color:#f0f0f0!important;font-family:monospace;font-size:13px;border-radius:6px;width:90px!important;height:34px}
    .params-bar .form-control:focus{border-color:#8b9eff!important;box-shadow:0 0 0 2px rgba(139,158,255,0.2)!important}
    .params-bar select.form-control{width:130px!important}
    .params-bar .checkbox{margin:0!important;display:inline-block}
    .params-bar .checkbox label{color:#aaa!important;font-size:11px!important;font-family:monospace;text-transform:uppercase;letter-spacing:0.5px;padding-left:4px;white-space:nowrap}
    .params-bar input[type=checkbox]{accent-color:#8b9eff;cursor:pointer}
    .graf-sep{display:inline-block;width:1px;background:#3a3a3a;height:30px;margin:0 4px 5px;align-self:flex-end}
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

  uiOutput("main_title"),
  div(style="display:flex;align-items:center;gap:8px;margin-bottom:18px;",
    div(class = "main-subtitle", "Analisis prosodico -- slope y direccion entonativa"),
    actionButton("btn_ayuda", "?", class="btn-ayuda")
  ),

  # PARAMS BAR
  div(class = "params-bar",
    div(div(class="param-label","Umbral dir. ST"),
        numericInput("umbral_st",           NULL, value=0.8, min=0.05, step=0.05, width="90px")),
    div(div(class="param-label","Pausa (0=off)"),
        numericInput("umbral_pausa",        NULL, value=0.3, min=0,    step=0.05, width="90px")),
    div(div(class="param-label","Reajuste ST (0=off)"),
        numericInput("umbral_reajuste_st",  NULL, value=0,   min=0,    step=0.5,  width="110px")),
    div(div(class="param-label","Reajuste % (0=off)"),
        numericInput("umbral_reajuste_pct", NULL, value=0,   min=0,    step=5,    width="110px"))
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
        div(style="display:flex;gap:10px;align-items:center;margin-bottom:8px;",
          actionButton("btn_parse","Procesar datos", class="btn-parse"),
          actionButton("btn_clear","Limpiar",        class="btn-clear"),
          div(style="margin-left:6px;",
            selectInput("modo_analisis", NULL,
              choices  = c("Por enunciado"="por_enunciado",
                           "Por grupo"    ="por_grupo"),
              selected = "por_enunciado", width = "180px"))),
        tags$textarea(id="raw_text", class="data-textarea", rows="10",
          placeholder="enunciados\tF01\tF02\tpausa1\tF03\tF04\tpausa2\nenun1\t120\t135\t0.3\t128\t145\t0.5\nenun2\t138\t155\t0.4\t\t\t"),
        br(),
        div(class="sec-label","O subir archivo (.txt / .csv / .tsv)"),
        fileInput("file_upload", NULL, accept=c(".txt",".csv",".tsv"),
                  buttonLabel="Elegir archivo", placeholder="ningun archivo")
      ),

      # TAB 2: Tabla
      tabPanel("02 | Tabla",
        div(class="sec-label","Datos transformados -- un grupo tonal por fila"),
        DTOutput("tabla_larga")
      ),

      # TAB 3: Grafico
      tabPanel("03 | Grafico",
        uiOutput("tab3_label_top"),
        div(class="params-bar", style="margin-bottom:12px;",
          div(div(class="param-label","Primer grupo"),
              numericInput("graf_desde",NULL,value=1,min=1,step=1,width="70px")),
          div(div(class="param-label","Ultimo grupo"),
              numericInput("graf_hasta",NULL,value=5,min=1,step=1,width="70px")),
          div(style="padding-bottom:5px;",
              actionButton("btn_primero5","Primeros 5",class="btn-clear")),
          uiOutput("lbl_n_grupos"),
          tags$span(class="graf-sep"),
          div(div(class="param-label","Nombre clausulas"),
              tags$input(id="graf_nombre_cl", type="text", value="clausula",
                class="form-control", style="width:110px!important;")),
          div(div(class="param-label","Nombre aislados"),
              tags$input(id="graf_nombre_ais", type="text", value="singular",
                class="form-control", style="width:110px!important;")),
          tags$span(class="graf-sep"),
          div(div(class="param-label","Asc."),
              tags$input(id="dir_nombre_asc", type="text", value="ascendente",
                class="form-control", style="width:110px!important;")),
          div(div(class="param-label","Desc."),
              tags$input(id="dir_nombre_desc", type="text", value="descendente",
                class="form-control", style="width:110px!important;")),
          div(div(class="param-label","Plano"),
              tags$input(id="dir_nombre_plan", type="text", value="plano",
                class="form-control", style="width:110px!important;")),
          tags$span(class="graf-sep"),
          div(div(class="param-label","Regresion"),
              uiOutput("graf_modo_tend_ui")),
          div(div(class="param-label","Tendencias"),
              div(style="display:flex;gap:10px;padding-top:2px;",
                checkboxInput("graf_tend_fin","fin",  value=FALSE),
                checkboxInput("graf_tend_ini","ini",  value=FALSE),
                checkboxInput("graf_tend_med","media",value=FALSE))),
          uiOutput("graf_enun_ui")
        ),
        plotlyOutput("grafico", height="520px"),
        br(),
        div(class="sec-label","Valores del grafico"),
        DTOutput("tabla_graf")
      ),

      # TAB 4: Clausulas
      tabPanel("04 | Clausulas",
        uiOutput("tab4_label_dist"),
        fluidRow(
          column(3, uiOutput("stat_total")),
          column(3, uiOutput("stat_asc")),
          column(3, uiOutput("stat_desc")),
          column(3, uiOutput("stat_plan"))
        ),
        br(),
        uiOutput("tab4_label_det"),
        DTOutput("tabla_clausulas"),
        br(),
        div(class="sec-label","Comparacion slope ST vs slope %"),
        DTOutput("tabla_comp")
      ),

      # TAB 5: Descriptiva
      tabPanel("05 | Descriptiva",
        tabsetPanel(id="tabs_desc", type="pills",

          tabPanel("Resumen", br(),
            uiOutput("desc_resumen")
          ),

          tabPanel("General", br(),
            div(class="sec-label","Distribucion de variables -- todos los grupos"),
            DTOutput("desc_general")
          ),

          tabPanel("Por enunciado", br(),
            div(class="sec-label","Estadisticos por enunciado"),
            DTOutput("desc_enun")
          ),

          tabPanel("Por posicion", br(),
            div(class="sec-label","Estadisticos por posicion de grupo"),
            DTOutput("desc_pos")
          ),

          tabPanel("Por longitud", br(),
            div(class="sec-label","Enunciados segun numero de grupos"),
            DTOutput("desc_longitud")
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
      ),

      # =========================================================
      # TAB 7: Visualizador manual de progresión tonal
      # =========================================================
      tabPanel("07 | Visualizador tonal",

        div(class = "sec-label", "Configuración general"),
        div(class = "info-box", HTML(
          "Introduce manualmente los valores de F0 inicial y final de cada grupo entonativo.<br>
           El <b>umbral ST</b> de la barra superior clasifica la tendencia interna de cada GE.<br>
           Los grupos de declinación necesitan <b>mínimo 2 GE</b>; los GE aislados reciben su propia etiqueta.<br>
           En modo automático, un GE <i>plano</i> puede unirse al grupo contiguo si el reajuste no supera el umbral."
        )),

        fluidRow(
          # ── Columna de controles ──────────────────────────────
          column(4,
            div(class = "main-box", style = "padding: 18px;",

              div(class = "sec-label", "Grupos entonativos"),
              numericInput("vis_n_ge", "Número de grupos (máx. 25):",
                           value = 4, min = 1, max = 25, step = 1),

              div(style = "height:10px;"),
              div(class = "sec-label", "Etiquetas"),

              # Etiqueta para grupos de 2+ GE
              div(class = "param-label", "Nombre para grupos de declinación (≥ 2 GE):"),
              tags$input(id = "vis_nombre_grupo", type = "text",
                         class = "form-control",
                         style = "font-family:monospace;font-size:12px;
                                  background:#fafaf8;border:1.5px solid #e0dfd9;
                                  border-radius:6px;padding:7px 10px;width:100%;
                                  box-sizing:border-box;color:#1a1a1a;margin-bottom:8px;",
                         placeholder = "unidad, bloque, cláusula…",
                         value = "unidad"),

              # Etiqueta para GE aislados (sin grupo)
              div(class = "param-label", "Nombre para GE aislados (1 GE solo):"),
              tags$input(id = "vis_nombre_aislado", type = "text",
                         class = "form-control",
                         style = "font-family:monospace;font-size:12px;
                                  background:#fafaf8;border:1.5px solid #e0dfd9;
                                  border-radius:6px;padding:7px 10px;width:100%;
                                  box-sizing:border-box;color:#1a1a1a;",
                         placeholder = "singular, aislado, átono…",
                         value = "singular"),

              div(style = "height:14px;"),
              div(class = "sec-label", "Modo de agrupación"),
              radioButtons("vis_modo_grupos", label = NULL,
                           choices = list("Manual (selector por GE)" = "manual",
                                          "Automático" = "auto"),
                           selected = "manual", inline = TRUE),

              conditionalPanel("input.vis_modo_grupos == 'auto'",
                div(style = "background:#f0f4ff;border-radius:8px;padding:11px 14px;
                             margin-top:6px;border-left:3px solid #8b9eff;",
                  div(class = "sec-label", style = "margin-top:0;",
                      "Umbral de reajuste para nuevo grupo"),
                  radioButtons("vis_unidad_umbral", label = NULL,
                               choices = list("Porcentaje (%)" = "pct",
                                              "Semitonos (ST)" = "st"),
                               selected = "st", inline = TRUE),
                  numericInput("vis_umbral_valor", "Valor del umbral:",
                               value = 1.5, min = 0.1, max = 100, step = 0.1),
                  div(style = "font-size:11px;color:#888;margin-top:4px;font-family:monospace;",
                    "Reajuste F0_fin(n) → F0_ini(n+1) que, si se supera, abre nuevo grupo.
                     Los GE planos se unen al contiguo si el reajuste no supera el umbral.")
                )
              ),

              div(style = "height:14px;"),
              div(class = "sec-label", "Líneas de tendencia (regresión lineal)"),
              radioButtons("vis_modo_tend", label = NULL,
                           choices = list("Global (todos los GE)" = "global",
                                          "Por grupo de declinación" = "por_grupo"),
                           selected = "global", inline = TRUE),
              checkboxInput("vis_tend_finales",   "F0 finales   (-- gris)",   value = TRUE),
              checkboxInput("vis_tend_iniciales", "F0 iniciales (·· verde)",  value = FALSE),
              checkboxInput("vis_tend_medias",    "F0 medias    (— naranja)", value = FALSE),

              div(style = "height:10px;"),
              checkboxInput("vis_etiquetas",
                            "Mostrar etiquetas de Hz (solo si ≤ 12 GE)",
                            value = TRUE),

              hr(),
              div(class = "sec-label", "Valores de F0 por grupo entonativo"),
              uiOutput("vis_inputs_ge"),

              div(style = "height:12px;"),
              actionButton("vis_actualizar", "Generar gráfico",
                           class = "btn-parse", style = "width:100%;")
            )
          ),

          # ── Columna de gráfico y resumen ──────────────────────
          column(8,
            plotOutput("vis_grafico", height = "500px"),
            div(style = "height:10px;"),
            uiOutput("vis_resumen"),
            div(style = "height:10px;"),
            uiOutput("vis_leyenda")
          )
        )
      )

    ) # fin tabsetPanel
  )   # fin main-box
)     # fin fluidPage

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
    calcular(datos_largos(),
             umbral_st           = input$umbral_st          %||% 0.8,
             umbral_prog_st      = 2.0,
             consecutivo         = (input$modo_analisis %||% "por_enunciado") == "por_grupo",
             umbral_pausa        = input$umbral_pausa        %||% 0.3,
             umbral_reajuste_st  = input$umbral_reajuste_st  %||% 0,
             umbral_reajuste_pct = input$umbral_reajuste_pct %||% 0)
  })

  # User-chosen display names for clauses and isolated groups
  ncl  <- reactive({ v <- trimws(input$graf_nombre_cl  %||% ""); if (nchar(v) == 0) "clausula" else v })
  nais <- reactive({ v <- trimws(input$graf_nombre_ais %||% ""); if (nchar(v) == 0) "singular" else v })

  # User-chosen direction labels (internal values → display)
  dir_map <- reactive({
    asc  <- trimws(input$dir_nombre_asc  %||% ""); if (nchar(asc)  == 0) asc  <- "ascendente"
    desc <- trimws(input$dir_nombre_desc %||% ""); if (nchar(desc) == 0) desc <- "descendente"
    plan <- trimws(input$dir_nombre_plan %||% ""); if (nchar(plan) == 0) plan <- "plano"
    c(ascendente = asc, descendente = desc, plano = plan)
  })

  # Full mapping: internal clause ID → display label (same logic as hacer_grafico)
  cl_disp_map <- reactive({
    req(resultado())
    res     <- resultado()
    cl_n    <- res |> group_by(clausula) |> summarise(n_g = n(), .groups = "drop")
    cl_multi <- cl_n |> filter(n_g >= 2) |> pull(clausula)
    all_cls  <- res |> arrange(x) |> distinct(clausula) |> pull()
    m_i <- 1L; s_i <- 1L
    mp  <- setNames(character(length(all_cls)), all_cls)
    for (cl in all_cls) {
      if (cl %in% cl_multi) { mp[[cl]] <- paste0(ncl(),  m_i); m_i <- m_i + 1L }
      else                   { mp[[cl]] <- paste0(nais(), s_i); s_i <- s_i + 1L }
    }
    mp
  })

  # Dynamic labels and controls that depend on ncl()
  output$main_title <- renderUI({
    div(class = "main-title", paste0(tools::toTitleCase(ncl()), " F0"))
  })
  output$tab3_label_top <- renderUI({
    div(class = "sec-label", paste0("Progresiones de F0 por ", ncl()))
  })
  output$graf_modo_tend_ui <- renderUI({
    selectInput("graf_modo_tend", NULL,
      choices  = setNames(c("global","por_cl"), c("Global", paste0("Por ", ncl()))),
      selected = isolate(input$graf_modo_tend) %||% "global",
      width    = "130px")
  })
  output$tab4_label_dist <- renderUI({
    div(class = "sec-label", paste0("Distribucion de ", ncl(), "s entonativas"))
  })
  output$tab4_label_det <- renderUI({
    div(class = "sec-label", paste0("Detalle por ", ncl()))
  })
  # Only show "Mostrar enunciados" checkbox when mode is "por_grupo"
  output$graf_enun_ui <- renderUI({
    req(input$modo_analisis)
    if (input$modo_analisis == "por_grupo") {
      div(div(class="param-label","Enunciados"),
          div(style="padding-top:4px;",
            checkboxInput("graf_mostrar_enun", "Mostrar", value=TRUE)))
    }
  })

  # ---------------------------------------------------------------------------
  # Tab 2: Tabla
  # ---------------------------------------------------------------------------
  output$tabla_larga <- renderDT({
    req(resultado(), cl_disp_map())
    mp <- cl_disp_map(); dm <- dir_map()
    resultado() |>
      select(enunciado, grupo, f0_ini, f0_fin, pausa,
             inflexion_pct, inflexion_st,
             reajuste_hz, reajuste_pct, reajuste_st,
             dir_grupo, slope_st_v, dir_slope_st,
             clausula_slope_st) |>
      mutate(across(where(is.double), \(x) round(x, 3)),
             clausula_slope_st = mp[clausula_slope_st],
             dir_grupo    = recode(dir_grupo,    !!!dm),
             dir_slope_st = recode(dir_slope_st, !!!dm)) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,"F0 ini"=f0_ini,"F0 fin"=f0_fin,
             "Pausa"=pausa,
             "Inflexion %"=inflexion_pct,"Inflexion ST"=inflexion_st,
             "Reajuste Hz"=reajuste_hz,"Reajuste %"=reajuste_pct,"Reajuste ST"=reajuste_st,
             "Dir. propio"=dir_grupo,"Slope ST"=slope_st_v,
             "Dir. slope"=dir_slope_st) |>
      rename(!!ncl() := clausula_slope_st) |>
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
  observeEvent(input$btn_primero5, {
    n <- n_grupos_total()
    updateNumericInput(session, "graf_desde", value=1)
    updateNumericInput(session, "graf_hasta", value=min(5,n))
  })

  observeEvent(input$btn_ayuda, {
    showModal(modalDialog(
      title = "Ayuda — variables y parametros",
      size  = "l",
      easyClose = TRUE,
      footer = modalButton("Cerrar"),
      div(class="sec-label","Variables de analisis"),
      div(class="info-box", style="line-height:1.9;", HTML("
        <b>Enunciado</b>: identificador del enunciado (frase, turno u oracion). Unidad de agrupacion superior, tomada del archivo de entrada.<br>
        <b>Grupo</b>: numero de grupo tonal dentro del enunciado. Cada trio (F0 ini, F0 fin, pausa) en el archivo define un grupo.<br>
        <b>F0 ini / F0 fin</b>: frecuencia fundamental en Hz al inicio y al final del grupo. Definen el contorno interno.<br>
        <b>Pausa</b>: duracion en segundos de la pausa al final del grupo (0 = sin pausa).<br>
        <br>
        <b>Inflexion %</b>: movimiento interno del grupo [(F0_fin &minus; F0_ini) / F0_ini &times; 100]. Positivo = sube; negativo = baja.<br>
        <b>Inflexion ST</b>: mismo movimiento en semitonos [12 &times; log2(F0_fin / F0_ini)]. Escala perceptivamente uniforme. Determina <i>Dir. propio</i>.<br>
        <br>
        <b>Reajuste Hz / % / ST</b>: salto de F0 en la juntura entre grupos [F0_ini(n) &minus; F0_fin(n&minus;1)]. Mide el reset de altura al inicio de cada grupo respecto al final del anterior. El ST es perceptivamente uniforme (1 ST &asymp; umbral de percepcion tonal). Valores positivos = subida en la transicion; negativos = bajada.<br>
        <b>Acum. ST</b>: suma acumulada de Reajuste ST desde el ultimo cambio de direccion. Se reinicia cuando la direccion se invierte (reajuste supera el umbral local en sentido contrario).<br>
        <b>Progresion</b>: etiqueta derivada de Acum. ST. <i>prog_asc</i> / <i>prog_desc</i> cuando la acumulacion supera el umbral configurado; <i>neutro</i> si no se ha alcanzado aun.<br>
        <br>
        <b>Dir. propio</b>: direccion del contorno interno del grupo (ascendente / descendente / plano) segun Inflexion ST frente al Umbral dir. ST.<br>
        <b>dST</b>: diferencia en ST entre F0 ini consecutivos [12 &times; log2(F0_ini(n) / F0_ini(n&minus;1))]. Alimenta Dir. local.<br>
        <b>Dir. local</b>: direccion segun dST frente al Umbral dir. ST. Asignacion alternativa de referencia.<br>
        <br>
        <b>Slope ST</b>: para clausulas de 2+ grupos, pendiente de regresion sobre la trayectoria interleaved completa (F0 ini1, F0 fin1, F0 ini2, F0 fin2, &hellip;). Para clausulas de un solo grupo, equivale a Inflexion ST. Refleja la tendencia global de la clausula.<br>
        <b>Dir. slope</b>: direccion definitiva (ascendente / descendente / plano) segun Slope ST frente al Umbral dir. ST. Usada en el grafico, boxplots y chi cuadrado.<br>
        <b>Clausula</b>: clausula entonativa asignada al grupo. Se agrupa con el anterior si comparten direccion y el reajuste no supera el umbral configurado. Variable de segmentacion principal.<br>
        <br>
        <b>Pausa anterior (s)</b>: duracion en segundos de la pausa registrada al final del grupo anterior, antes de la transicion analizada.<br>
        <b>N</b>: numero de transiciones disponibles en la categoria (pestana 05).<br>
        <br>
        <b>Parametros de la barra superior</b><br>
        <b>Umbral dir. ST</b>: inflexion minima en ST para etiquetar un grupo como ascendente o descendente. Por debajo del umbral = plano.<br>
        <b>Reajuste ST / %</b>: si el salto de F0 entre dos grupos consecutivos supera este valor (en semitonos o porcentaje), se fuerza una nueva clausula independientemente de la direccion. Valor 0 = desactivado.
      ")),
      br(),
      div(class="sec-label","Pestana 06 — chi cuadrado"),
      div(class="info-box", style="line-height:1.9;", HTML("
        <b>Modo: Por grupo (activo)</b><br>
        Las clausulas pueden abarcar varios enunciados.<br>
        <b>Pura</b>: todos los grupos de la clausula pertenecen al mismo enunciado.<br>
        <b>Mixta — cruce en frontera</b>: la clausula contiene el grupo final de un enunciado y el inicial del siguiente.<br>
        <b>Mixta — sin cruce directo</b>: agrupa grupos de enunciados distintos sin cruzar la frontera inmediata.<br>
        El test chi cuadrado evalua la asociacion global entre enunciados y clausulas (Monte Carlo 5000 sim.).<br>
        <br>
        <b>Modo: Por enunciado (activo)</b><br>
        Las clausulas estan contenidas dentro de cada enunciado.<br>
        <b>1 clausula</b>: correspondencia plena enunciado = clausula.<br>
        <b>2 clausulas / 3+ clausulas</b>: el enunciado se segmenta internamente en varias clausulas.<br>
        El test chi cuadrado evalua la homogeneidad de la distribucion de direcciones por enunciado.
      "))
    ))
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
    modo <- input$modo_analisis %||% "por_enunciado"
    p <- hacer_grafico(resultado_graf(),
           umbral_st    = input$umbral_st %||% 0.8,
           nombre_cl    = trimws(input$graf_nombre_cl  %||% "clausula"),
           nombre_ais   = trimws(input$graf_nombre_ais %||% "singular"),
           dir_labels   = dir_map(),
           tend_fin     = isTRUE(input$graf_tend_fin),
           tend_ini     = isTRUE(input$graf_tend_ini),
           tend_med     = isTRUE(input$graf_tend_med),
           modo_tend    = input$graf_modo_tend %||% "global",
           mostrar_enun = if (modo == "por_grupo") isTRUE(input$graf_mostrar_enun) else TRUE)
    ggplotly(p, tooltip = c("x","y","colour")) |> pl_config("clausulas_f0")
  })

  output$tabla_graf <- renderDT({
    req(resultado_graf(), cl_disp_map())
    mp <- cl_disp_map(); dm <- dir_map()
    resultado_graf() |>
      select(enunciado, grupo, f0_ini, f0_fin, pausa,
             reajuste_hz, reajuste_pct, reajuste_st,
             clausula_slope_st, dir_slope_st) |>
      mutate(across(where(is.double), \(x) round(x, 3)),
             clausula_slope_st = mp[clausula_slope_st],
             dir_slope_st      = recode(dir_slope_st, !!!dm)) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,"F0 ini"=f0_ini,"F0 fin"=f0_fin,
             "Pausa"=pausa,"Reajuste Hz"=reajuste_hz,"Reajuste %"=reajuste_pct,
             "Reajuste ST"=reajuste_st,"Dir."=dir_slope_st) |>
      rename(!!ncl() := clausula_slope_st) |>
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
  output$stat_total <- renderUI({ req(resumen()); mk_stat(n_distinct(resumen()$clausula_slope_st), paste0(ncl(), " totales")) })
  output$stat_asc   <- renderUI({ req(resumen()); dm <- dir_map(); mk_stat(sum(resumen()$dir_slope_st=="ascendente",na.rm=TRUE),  dm["ascendente"],  "stat-card stat-asc") })
  output$stat_desc  <- renderUI({ req(resumen()); dm <- dir_map(); mk_stat(sum(resumen()$dir_slope_st=="descendente",na.rm=TRUE), dm["descendente"], "stat-card stat-desc") })
  output$stat_plan  <- renderUI({ req(resumen()); dm <- dir_map(); mk_stat(sum(resumen()$dir_slope_st=="plano",na.rm=TRUE),       dm["plano"],       "stat-card stat-plan") })

  output$tabla_clausulas <- renderDT({
    req(resumen(), cl_disp_map())
    mp <- cl_disp_map(); dm <- dir_map()
    resumen() |>
      mutate(clausula_slope_st = mp[clausula_slope_st],
             dir_slope_st      = recode(dir_slope_st, !!!dm)) |>
      rename("Direccion"=dir_slope_st,"N grupos"=n_grupos,
             "Enunciados"=enunciados,"Grupos"=grupos,"F0 media"=f0_media,"Slope ST"=slope_medio) |>
      rename(!!ncl() := clausula_slope_st) |>
      datatable(options=list(pageLength=15, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  output$tabla_comp <- renderDT({
    req(resultado(), cl_disp_map())
    mp <- cl_disp_map(); dm <- dir_map()
    resultado() |>
      select(enunciado, grupo, clausula_slope_st, dir_slope_st,
             inflexion_st, slope_st_v) |>
      mutate(across(where(is.double), \(x) round(x, 3)),
             clausula_slope_st = mp[clausula_slope_st],
             dir_slope_st      = recode(dir_slope_st, !!!dm)) |>
      rename("Enunciado"=enunciado,"Grupo"=grupo,
             "Dir."=dir_slope_st,
             "Inflexion ST"=inflexion_st,
             "Slope interleaved ST"=slope_st_v) |>
      rename(!!ncl() := clausula_slope_st) |>
      datatable(options=list(pageLength=15, scrollX=TRUE, dom="lrtip",
        language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames=FALSE, class="stripe hover")
  })

  # ---------------------------------------------------------------------------
  # Tab 5: Descriptiva
  # ---------------------------------------------------------------------------

  # Narrative summary helper: "media X, mediana X, min X, max X, SD X"
  fmt_nar <- function(x, dig = 3) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return("sin datos")
    paste0("media&nbsp;", round(mean(x),dig),
           ", mediana&nbsp;", round(median(x),dig),
           ", min&nbsp;", round(min(x),dig),
           ", max&nbsp;", round(max(x),dig),
           ", SD&nbsp;", round(sd(x),dig))
  }

  output$desc_resumen <- renderUI({
    req(resultado())
    res  <- resultado()
    modo <- input$modo_analisis %||% "por_enunciado"

    # --- Reajuste general (all non-NA transitions) ---
    rej_st  <- res$reajuste_st[!is.na(res$reajuste_st)]
    rej_pct <- res$reajuste_pct[!is.na(res$reajuste_pct)]
    pau_all <- res$pausa[!is.na(res$pausa)]
    pau_cero <- sum(pau_all == 0, na.rm = TRUE)

    # --- Within-utterance vs between-utterance reajuste ---
    res2 <- res |> mutate(prev_enun = lag(enunciado))
    if (modo == "por_enunciado") {
      # In por_enunciado mode all non-NA reajuste is already within-utterance
      rej_intra_st  <- rej_st
      rej_intra_pct <- rej_pct
      pau_intra     <- pau_all
      rej_inter_st  <- numeric(0)
      rej_inter_pct <- numeric(0)
      pau_inter     <- numeric(0)
    } else {
      mask_intra <- !is.na(res2$reajuste_st) & !is.na(res2$prev_enun) &
                    res2$enunciado == res2$prev_enun
      mask_inter <- !is.na(res2$reajuste_st) & !is.na(res2$prev_enun) &
                    res2$enunciado != res2$prev_enun
      rej_intra_st  <- res2$reajuste_st[mask_intra]
      rej_intra_pct <- res2$reajuste_pct[mask_intra]
      pau_intra     <- res2$pausa[mask_intra]
      rej_inter_st  <- res2$reajuste_st[mask_inter]
      rej_inter_pct <- res2$reajuste_pct[mask_inter]
      pau_inter     <- res2$pausa[mask_inter]
    }

    # --- Utterance structure ---
    enun_n  <- res |> count(enunciado, name = "ng")
    n_enun  <- nrow(enun_n)
    media_g <- round(mean(enun_n$ng), 1)
    med_g   <- median(enun_n$ng)
    rng_g   <- paste0(min(enun_n$ng), "–", max(enun_n$ng))
    modal_n <- as.integer(names(which.max(table(enun_n$ng))))
    n_modal <- sum(enun_n$ng == modal_n)
    pct_modal <- round(n_modal / n_enun * 100, 1)
    modal_enuns <- enun_n |> filter(ng == modal_n) |> pull(enunciado)
    res_modal <- res |> filter(enunciado %in% modal_enuns, !is.na(reajuste_st))

    lnk <- function(txt) paste0("<b>", txt, "</b>")

    HTML(paste0(
      # Bloque 1: reajuste general
      "<p><b>Reajuste entre grupos</b> (N&nbsp;=&nbsp;", length(rej_st), " transiciones)<br>",
      "&emsp;ST:&nbsp;&nbsp;",  fmt_nar(rej_st,  3), "<br>",
      "&emsp;%:&nbsp;&nbsp;&nbsp;", fmt_nar(rej_pct, 2), "<br>",
      "&emsp;Pausa:&nbsp;", fmt_nar(pau_all, 3), "&nbsp;s",
      "&nbsp;&mdash;&nbsp;", pau_cero, " grupos sin pausa (",
      round(pau_cero / length(pau_all) * 100, 1), "%)</p>",

      # Bloque 2: reajuste intra-enunciado
      "<p><b>Reajuste entre grupos del mismo enunciado</b> (N&nbsp;=&nbsp;", length(rej_intra_st), ")<br>",
      "&emsp;ST:&nbsp;&nbsp;",  fmt_nar(rej_intra_st,  3), "<br>",
      "&emsp;%:&nbsp;&nbsp;&nbsp;", fmt_nar(rej_intra_pct, 2), "<br>",
      "&emsp;Pausa:&nbsp;", fmt_nar(pau_intra, 3), "&nbsp;s</p>",

      # Bloque 3: reajuste inter-enunciado (solo modo por_grupo)
      if (length(rej_inter_st) > 0) paste0(
        "<p><b>Reajuste en fronteras de enunciado</b> (N&nbsp;=&nbsp;", length(rej_inter_st), ")<br>",
        "&emsp;ST:&nbsp;&nbsp;",  fmt_nar(rej_inter_st,  3), "<br>",
        "&emsp;%:&nbsp;&nbsp;&nbsp;", fmt_nar(rej_inter_pct, 2), "<br>",
        "&emsp;Pausa:&nbsp;", fmt_nar(pau_inter, 3), "&nbsp;s</p>"
      ) else "",

      # Bloque 4: estructura del corpus
      "<p><b>Estructura del corpus</b><br>",
      "&emsp;", n_enun, " enunciados &mdash; ",
      media_g, " grupos de media (mediana&nbsp;", med_g, ", rango&nbsp;", rng_g, ")<br>",
      "&emsp;El tamaño más frecuente es ", lnk(paste0(modal_n, " grupo", if (modal_n > 1) "s" else "")),
      " (N&nbsp;=&nbsp;", n_modal, " enunciados, ", pct_modal, "% del corpus)</p>",

      # Bloque 5: reajuste en enunciados del tamaño modal
      if (nrow(res_modal) > 0) paste0(
        "<p><b>Enunciados con ", modal_n, " grupo", if (modal_n > 1) "s" else "",
        "</b> &mdash; reajuste entre sus grupos:<br>",
        "&emsp;ST:&nbsp;&nbsp;",  fmt_nar(res_modal$reajuste_st,  3), "<br>",
        "&emsp;%:&nbsp;&nbsp;&nbsp;", fmt_nar(res_modal$reajuste_pct, 2), "<br>",
        "&emsp;Pausa:&nbsp;", fmt_nar(res_modal$pausa, 3), "&nbsp;s</p>"
      ) else ""
    ))
  })

  desc_vec <- function(x, dig = 3) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 1) return(tibble(N=0L, Min=NA_real_, Max=NA_real_, Media=NA_real_,
                               Mediana=NA_real_, SD=NA_real_, Asimetria=NA_real_, Curtosis=NA_real_))
    m <- mean(x); s <- sd(x)
    tibble(
      N         = as.integer(n),
      Min       = round(min(x),    dig),
      Max       = round(max(x),    dig),
      Media     = round(m,         dig),
      Mediana   = round(median(x), dig),
      SD        = round(s,         dig),
      Asimetria = round(if (!is.na(s) && s > 0 && n >= 3) mean((x-m)^3)/s^3  else NA_real_, 3),
      Curtosis  = round(if (!is.na(s) && s > 0 && n >= 4) mean((x-m)^4)/s^4 - 3 else NA_real_, 3)
    )
  }

  make_dt5 <- function(df) {
    datatable(df, options=list(pageLength=20, scrollX=TRUE, dom="lrtip",
      language=list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
      rownames=FALSE, class="stripe hover")
  }

  output$desc_general <- renderDT({
    req(resultado())
    res <- resultado()
    nm  <- c("F0 ini (Hz)", "F0 fin (Hz)", "Inflexion %", "Inflexion ST",
             "Reajuste Hz", "Reajuste %",  "Reajuste ST", "Pausa (s)")
    vl  <- list(res$f0_ini, res$f0_fin, res$inflexion_pct, res$inflexion_st,
                res$reajuste_hz, res$reajuste_pct, res$reajuste_st, res$pausa)
    map2_dfr(nm, vl, \(n, v) bind_cols(Variable = n, desc_vec(v))) |>
      make_dt5()
  })

  output$desc_enun <- renderDT({
    req(resultado())
    resultado() |>
      group_by(enunciado) |>
      summarise(
        N          = n(),
        F0ini_med  = round(mean(f0_ini,       na.rm=TRUE), 1),
        F0ini_SD   = round(sd(f0_ini,         na.rm=TRUE), 1),
        F0fin_med  = round(mean(f0_fin,       na.rm=TRUE), 1),
        F0fin_SD   = round(sd(f0_fin,         na.rm=TRUE), 1),
        InflST_med = round(mean(inflexion_st, na.rm=TRUE), 3),
        InflST_SD  = round(sd(inflexion_st,   na.rm=TRUE), 3),
        RejST_med  = round(mean(reajuste_st,  na.rm=TRUE), 3),
        RejST_SD   = round(sd(reajuste_st,    na.rm=TRUE), 3),
        Pausa_med  = round(mean(pausa,        na.rm=TRUE), 3),
        Pausa_SD   = round(sd(pausa,          na.rm=TRUE), 3),
        .groups    = "drop"
      ) |>
      rename("Enunciado"=enunciado, "N grupos"=N,
             "F0ini med"=F0ini_med, "F0ini SD"=F0ini_SD,
             "F0fin med"=F0fin_med, "F0fin SD"=F0fin_SD,
             "Infl.ST med"=InflST_med, "Infl.ST SD"=InflST_SD,
             "Rej.ST med"=RejST_med,  "Rej.ST SD"=RejST_SD,
             "Pausa med"=Pausa_med,   "Pausa SD"=Pausa_SD) |>
      make_dt5()
  })

  output$desc_pos <- renderDT({
    req(resultado())
    resultado() |>
      group_by(grupo) |>
      summarise(
        N          = n(),
        F0ini_med  = round(mean(f0_ini,       na.rm=TRUE), 1),
        F0ini_SD   = round(sd(f0_ini,         na.rm=TRUE), 1),
        F0fin_med  = round(mean(f0_fin,       na.rm=TRUE), 1),
        F0fin_SD   = round(sd(f0_fin,         na.rm=TRUE), 1),
        InflST_med = round(mean(inflexion_st, na.rm=TRUE), 3),
        InflST_SD  = round(sd(inflexion_st,   na.rm=TRUE), 3),
        RejST_med  = round(mean(reajuste_st,  na.rm=TRUE), 3),
        RejST_SD   = round(sd(reajuste_st,    na.rm=TRUE), 3),
        Pausa_med  = round(mean(pausa,        na.rm=TRUE), 3),
        Pausa_SD   = round(sd(pausa,          na.rm=TRUE), 3),
        .groups    = "drop"
      ) |>
      rename("Posicion"=grupo, "N enunciados"=N,
             "F0ini med"=F0ini_med, "F0ini SD"=F0ini_SD,
             "F0fin med"=F0fin_med, "F0fin SD"=F0fin_SD,
             "Infl.ST med"=InflST_med, "Infl.ST SD"=InflST_SD,
             "Rej.ST med"=RejST_med,  "Rej.ST SD"=RejST_SD,
             "Pausa med"=Pausa_med,   "Pausa SD"=Pausa_SD) |>
      make_dt5()
  })

  output$desc_longitud <- renderDT({
    req(resultado())
    res <- resultado()
    # Count groups per utterance, then join back to get per-group values
    enun_ng <- res |> count(enunciado, name = "n_grupos")
    res |>
      left_join(enun_ng, by = "enunciado") |>
      group_by(n_grupos) |>
      summarise(
        N_enun     = n_distinct(enunciado),
        N_grupos   = n(),
        F0ini_med  = round(mean(f0_ini,       na.rm=TRUE), 1),
        F0ini_SD   = round(sd(f0_ini,         na.rm=TRUE), 1),
        F0fin_med  = round(mean(f0_fin,       na.rm=TRUE), 1),
        F0fin_SD   = round(sd(f0_fin,         na.rm=TRUE), 1),
        InflST_med = round(mean(inflexion_st, na.rm=TRUE), 3),
        InflST_SD  = round(sd(inflexion_st,   na.rm=TRUE), 3),
        RejST_med  = round(mean(reajuste_st,  na.rm=TRUE), 3),
        RejST_SD   = round(sd(reajuste_st,    na.rm=TRUE), 3),
        Pausa_med  = round(mean(pausa,        na.rm=TRUE), 3),
        Pausa_SD   = round(sd(pausa,          na.rm=TRUE), 3),
        .groups    = "drop"
      ) |>
      arrange(n_grupos) |>
      rename("N grupos"=n_grupos, "N enunciados"=N_enun, "N obs."=N_grupos,
             "F0ini med"=F0ini_med, "F0ini SD"=F0ini_SD,
             "F0fin med"=F0fin_med, "F0fin SD"=F0fin_SD,
             "Infl.ST med"=InflST_med, "Infl.ST SD"=InflST_SD,
             "Rej.ST med"=RejST_med,  "Rej.ST SD"=RejST_SD,
             "Pausa med"=Pausa_med,   "Pausa SD"=Pausa_SD) |>
      make_dt5()
  })

  # ---------------------------------------------------------------------------
  # Tab 6: Chi cuadrado
  # ---------------------------------------------------------------------------

  # Dynamic section labels and info box
  output$tab6_infobox <- renderUI({
    cl <- ncl()
    if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") {
      div(class = "info-box", HTML(paste0(
        "Cada ", cl, " se clasifica en tres categorias:<br>
         <b>Pura</b>: todos sus grupos pertenecen al mismo enunciado.<br>
         <b>Mixta — cruce en frontera</b>: contiene el grupo final de un enunciado y el grupo inicial del siguiente.<br>
         <b>Mixta — sin cruce directo</b>: agrupa grupos de enunciados distintos sin estar en la frontera inmediata.<br>
         El test chi cuadrado evalua la asociacion global entre enunciados y ", cl, "s.")))
    } else {
      div(class = "info-box", HTML(paste0(
        "Modo no consecutivo: las ", cl, "s estan contenidas dentro de cada enunciado.<br>
         Los enunciados se clasifican segun cuantas ", cl, "s tonales integran en su interior:<br>
         <b>1 ", cl, "</b>: correspondencia plena entre enunciado y ", cl, " entonativa.<br>
         <b>2 ", cl, "s</b>: el enunciado se divide en dos ", cl, "s de distinta direccion.<br>
         <b>3+ ", cl, "s</b>: el enunciado presenta tres o mas ", cl, "s internas.<br>
         El test chi cuadrado evalua si la distribucion de direcciones por enunciado es homogenea.")))
    }
  })
  output$tab6_label_bar <- renderUI({
    cl <- ncl()
    lbl <- if ((input$modo_analisis %||% "por_enunciado") == "por_grupo")
      paste0(tools::toTitleCase(cl), "s segun numero de enunciados que las componen")
    else
      paste0("Enunciados segun numero de ", cl, "s en su interior")
    div(class = "sec-label", lbl)
  })
  output$tab6_label_table2 <- renderUI({
    cl <- ncl()
    lbl <- if ((input$modo_analisis %||% "por_enunciado") == "por_grupo")
      paste0(tools::toTitleCase(cl), "s mixtas (grupos de mas de un enunciado)")
    else
      paste0("Enunciados con varias ", cl, "s internas")
    div(class = "sec-label", lbl)
  })

  # NON-CONSECUTIVE: utterances classified by number of internal clauses
  cat_data_noc <- reactive({
    req(resultado())
    res <- resultado()

    mp <- cl_disp_map()
    enun_detail <- res |>
      group_by(enunciado) |>
      summarise(
        n_clausulas = n_distinct(clausula_slope_st),
        clausulas   = paste(mp[unique(clausula_slope_st)], collapse = ", "),
        dirs        = paste(na.omit(unique(dir_slope_st)), collapse = ", "),
        n_grupos    = n(),
        .groups     = "drop"
      ) |>
      mutate(categoria = case_when(
        n_clausulas == 1 ~ paste0("1 ", ncl()),
        n_clausulas == 2 ~ paste0("2 ", ncl(), "s"),
        TRUE             ~ paste0("3+ ", ncl(), "s")
      ))

    lvls <- c(paste0("1 ", ncl()), paste0("2 ", ncl(), "s"), paste0("3+ ", ncl(), "s"))
    cat_summary <- enun_detail |>
      count(categoria, name = "n_enun") |>
      right_join(tibble(categoria = lvls), by = "categoria") |>
      mutate(
        n_enun    = coalesce(n_enun, 0L),
        categoria = factor(categoria, levels = lvls),
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
    tst <- if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") chi_test() else chi_test_noc()
    req(tst)
    mk_stat_chi(round(tst$statistic, 2), "chi cuadrado")
  })
  output$chi_stat_p <- renderUI({
    tst <- if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") chi_test() else chi_test_noc()
    req(tst)
    p   <- tst$p.value
    lbl <- if (p < 0.001) "p < 0.001" else paste0("p = ", round(p, 4))
    cls <- if (p < 0.05) "stat-card stat-asc" else "stat-card stat-plan"
    mk_stat_chi(lbl, "p-valor (Monte Carlo)", cls)
  })
  output$chi_stat_v <- renderUI({
    tst <- if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") chi_test() else chi_test_noc()
    req(tst)
    if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") {
      ct <- chi_data()$cont; n <- sum(ct); k <- min(nrow(ct), ncol(ct)) - 1
    } else {
      ct <- chi_data_noc();  n <- sum(ct); k <- min(nrow(ct), ncol(ct)) - 1
    }
    v <- if (!is.null(ct) && k > 0) round(sqrt(tst$statistic / (n * k)), 3) else NA
    mk_stat_chi(v, "V de Cramer")
  })
  output$chi_stat_mix <- renderUI({
    if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") {
      req(chi_data())
      n_mix <- chi_data()$res |>
        group_by(clausula_slope_st) |>
        summarise(nd = n_distinct(enunciado), .groups = "drop") |>
        filter(nd > 1) |> nrow()
      mk_stat_chi(n_mix, paste0(ncl(), "s mixtas"), "stat-card stat-desc")
    } else {
      req(cat_data_noc())
      det  <- cat_data_noc()$detail
      n_m  <- sum(det$n_clausulas > 1)
      mk_stat_chi(n_m, paste0("enunciados con 2+ ", ncl(), "s"), "stat-card stat-desc")
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
    cl <- ncl()
    if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") {
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
      subtt   <- paste0("Total ", cl, "s: ", n_total,
                        "  |  ", tools::toTitleCase(cl), "s mixtas (2+ enunciados): ", n_mix,
                        " (", round(n_mix / n_total * 100, 1), "%)")
      ylab  <- paste0("% de ", cl, "s")
      title <- paste0("Distribucion de ", cl, "s segun enunciados que las componen")
      fname <- paste0(cl, "s_chi")
    } else {
      # ---- non-consecutive: utterances by clause-count category ----
      req(cat_data_noc())
      df  <- cat_data_noc()$summary
      pal <- setNames(c("#4a4a4a", "#d68910", "#c0392b"),
                      c(paste0("1 ", cl), paste0("2 ", cl, "s"), paste0("3+ ", cl, "s")))
      n_total <- sum(df$n_enun)
      n_multi <- sum(df$n_enun[df$categoria != paste0("1 ", cl)])
      subtt   <- paste0("Total enunciados: ", n_total,
                        "  |  Enunciados con varias ", cl, "s: ", n_multi,
                        " (", round(n_multi / n_total * 100, 1), "%)")
      # Rename n_enun → n_clausulas so the ggplot aes works with a single block
      df    <- df |> rename(n_clausulas = n_enun)
      ylab  <- "% de enunciados"
      title <- paste0("Distribucion de enunciados segun ", cl, "s en su interior")
      fname <- paste0("enunciados_", cl, "s")
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
    if ((input$modo_analisis %||% "por_enunciado") == "por_grupo") {
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

    if (!(input$modo_analisis %||% "por_enunciado") == "por_grupo") {
      req(cat_data_noc())
      det <- cat_data_noc()$detail |>
        filter(n_clausulas > 1) |>
        arrange(desc(n_clausulas), enunciado)
      if (nrow(det) == 0) {
        return(datatable(
          tibble(Nota = paste0("Todos los enunciados tienen una sola ", ncl(), " interna.")),
          rownames = FALSE))
      }
      return(
        det |>
          mutate(dirs = sapply(strsplit(dirs, ", "), function(v) {
            dm <- dir_map(); paste(ifelse(v %in% names(dm), dm[v], v), collapse = ", ")
          })) |>
          rename("Enunciado"   = enunciado,
                 "Direcciones" = dirs,
                 "N grupos"    = n_grupos) |>
          rename(!!paste0("N ", ncl(), "s") := n_clausulas,
                 !!paste0(tools::toTitleCase(ncl()), "s") := clausulas) |>
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
        tibble(Nota = paste0("No hay ", ncl(), "s mixtas: cada ", ncl(), " pertenece a un unico enunciado.")),
        rownames = FALSE))
    }

    mix |>
      mutate(clausula_slope_st = cl_disp_map()[clausula_slope_st],
             dir_slope_st      = recode(dir_slope_st, !!!dir_map())) |>
      rename("Direccion"    = dir_slope_st,
             "N enunciados" = n_enunciados,
             "Enunciados"   = enunciados,
             "N grupos"     = n_grupos,
             "Tipo"         = tipo) |>
      rename(!!tools::toTitleCase(ncl()) := clausula_slope_st) |>
      datatable(options = list(pageLength = 15, scrollX = TRUE, dom = "lrtip",
        language = list(url="//cdn.datatables.net/plug-ins/1.10.11/i18n/Spanish.json")),
        rownames = FALSE, class = "stripe hover")
  })

  # ===========================================================================
  # TAB 7: Visualizador manual de progresión tonal
  # ===========================================================================

  vis_colores_ge <- c(
    "#2ecc71","#e67e22","#3498db","#e74c3c","#9b59b6",
    "#f1c40f","#1abc9c","#e91e63","#8bc34a","#d35400",
    "#00bcd4","#ff5722","#4caf50","#ffc107","#7e57c2",
    "#26c6da","#ef5350","#66bb6a","#ffa726","#ab47bc",
    "#26a69a","#ffca28","#42a5f5","#ec407a","#8d6e63"
  )
  vis_colores_grupo <- c(
    "#c0392b","#d68910","#27ae60","#2980b9","#8e44ad",
    "#16a085","#d35400","#7f8c8d","#1a5276","#145a32"
  )

  # ── Detección automática de grupos ──────────────────────────────────────────
  # Lógica:
  #   1. Se abre nuevo grupo si cambia la tendencia NO-plana (asc↔desc).
  #   2. Un GE "plano" puede unirse al grupo contiguo si el reajuste
  #      F0_fin(prev) → F0_ini(plano) no supera el umbral.
  #   3. Si el reajuste SÍ supera el umbral, el GE plano inicia un grupo nuevo.
  #   4. Grupos de 1 solo GE serán etiquetados como "aislados".
  vis_detectar_grupos <- function(df, unidad, umbral) {
    n <- nrow(df)
    if (n == 0) return(df)
    grupo    <- integer(n)
    grupo[1] <- 1L
    g        <- 1L

    for (i in seq(2, n)) {
      reajuste <- if (unidad == "pct") {
        abs(df$f0_ini[i] - df$f0_fin[i-1]) / df$f0_fin[i-1] * 100
      } else {
        abs(hz_a_st(df$f0_fin[i-1], df$f0_ini[i]))
      }

      tend_prev <- df$tendencia[i-1]
      tend_curr <- df$tendencia[i]

      nuevo_grupo <- if (reajuste >= umbral) {
        # Reajuste brusco → siempre nuevo grupo
        TRUE
      } else if (tend_curr == "plano") {
        # GE plano con reajuste suave → se une al grupo anterior
        FALSE
      } else if (tend_prev == "plano") {
        # GE no-plano tras plano → mismo grupo si dirección compatible
        # buscar la última tendencia no-plana anterior
        tend_base <- NA_character_
        for (j in seq(i-1, 1)) {
          if (df$tendencia[j] != "plano") { tend_base <- df$tendencia[j]; break }
        }
        !is.na(tend_base) && tend_base != tend_curr
      } else {
        # Ambos no-planos: nuevo grupo si cambia la dirección
        tend_prev != tend_curr
      }

      if (nuevo_grupo) g <- g + 1L
      grupo[i] <- g
    }

    df$grupo_num <- grupo
    df
  }

  # ── Inputs dinámicos por GE ─────────────────────────────────────────────────
  output$vis_inputs_ge <- renderUI({
    n       <- min(input$vis_n_ge %||% 4, 25L)
    opciones <- c("Sin grupo", paste("Grupo", seq_len(n)))
    lapply(seq_len(n), function(i) {
      col_i <- vis_colores_ge[((i-1) %% length(vis_colores_ge)) + 1]
      div(style = paste0("background:#f8f7f4;border-left:3px solid ", col_i,
                         ";border-radius:5px;padding:10px 12px;margin-bottom:8px;"),
        div(style = paste0("font-family:monospace;font-size:11px;font-weight:700;
                            color:", col_i, ";margin-bottom:6px;"),
            paste0("GE", i)),
        fluidRow(
          column(6,
            div(class = "param-label", "F0 inicial (Hz)"),
            numericInput(paste0("vis_ini_", i), NULL,
                         value = 130 - (i-1)*3, min = 50, max = 500, step = 1,
                         width = "100%")
          ),
          column(6,
            div(class = "param-label", "F0 final (Hz)"),
            numericInput(paste0("vis_fin_", i), NULL,
                         value = 116 - (i-1)*3, min = 50, max = 500, step = 1,
                         width = "100%")
          )
        ),
        if (isolate(input$vis_modo_grupos) == "manual")
          div(
            div(class = "param-label", style = "margin-top:6px;", "Grupo declinación"),
            selectInput(paste0("vis_grupo_", i), NULL,
                        choices = opciones, selected = "Sin grupo", width = "100%")
          )
      )
    })
  })

  # ── Datos reactivos ──────────────────────────────────────────────────────────
  vis_datos <- eventReactive(input$vis_actualizar, {
    n              <- min(input$vis_n_ge %||% 4, 25L)
    nombre         <- trimws(input$vis_nombre_grupo   %||% "unidad")
    nombre_aislado <- trimws(input$vis_nombre_aislado %||% "singular")
    if (nchar(nombre)         == 0) nombre         <- "unidad"
    if (nchar(nombre_aislado) == 0) nombre_aislado <- "singular"
    umbral_st_app <- input$umbral_st %||% 0.8

    registros <- lapply(seq_len(n), function(i) {
      f0_ini    <- input[[paste0("vis_ini_", i)]]   %||% 120
      f0_fin    <- input[[paste0("vis_fin_", i)]]   %||% 110
      grupo_raw <- input[[paste0("vis_grupo_", i)]] %||% "Sin grupo"
      if (is.na(f0_ini)) f0_ini <- 120
      if (is.na(f0_fin)) f0_fin <- 110

      inf_st  <- hz_a_st(f0_ini, f0_fin)
      inf_pct <- (f0_fin - f0_ini) / f0_ini * 100
      tend <- if (inf_st >  umbral_st_app) "ascendente"
              else if (inf_st < -umbral_st_app) "descendente"
              else "plano"

      grupo_num <- if (input$vis_modo_grupos == "auto" || grupo_raw == "Sin grupo")
                     NA_integer_
                   else as.integer(sub("Grupo ", "", grupo_raw))

      data.frame(
        ge           = i,
        ge_label     = paste0("GE", i),
        x_ini        = i - 0.38,
        x_fin        = i + 0.38,
        f0_ini       = f0_ini,
        f0_fin       = f0_fin,
        f0_media     = (f0_ini + f0_fin) / 2,
        delta_hz     = f0_fin - f0_ini,
        inflexion_st  = round(inf_st,  3),
        inflexion_pct = round(inf_pct, 2),
        grupo_num    = grupo_num,
        tendencia    = tend,
        stringsAsFactors = FALSE
      )
    })

    df <- bind_rows(registros)

    # Añadir reajuste respecto al GE anterior
    df$reajuste_hz  <- c(NA, df$f0_ini[-1] - df$f0_fin[-nrow(df)])
    df$reajuste_pct <- c(NA, (df$f0_ini[-1] - df$f0_fin[-nrow(df)]) /
                               df$f0_fin[-nrow(df)] * 100)
    df$reajuste_st  <- c(NA, sapply(seq(2, nrow(df)), function(i)
                               hz_a_st(df$f0_fin[i-1], df$f0_ini[i])))
    df$reajuste_hz  <- round(df$reajuste_hz,  1)
    df$reajuste_pct <- round(df$reajuste_pct, 2)
    df$reajuste_st  <- round(df$reajuste_st,  3)

    # Modo automático
    if (input$vis_modo_grupos == "auto") {
      unidad <- input$vis_unidad_umbral %||% "st"
      umbral <- input$vis_umbral_valor  %||% 1.5
      if (is.na(umbral) || umbral <= 0) umbral <- 1.5
      df <- vis_detectar_grupos(df, unidad, umbral)
    }

    # ── Asignar etiquetas ────────────────────────────────────────────────────
    # Contar cuántos GE hay en cada grupo_num
    if (any(!is.na(df$grupo_num))) {
      conteo_grupos <- table(df$grupo_num[!is.na(df$grupo_num)])

      # Grupos con ≥ 2 GE → reciben nombre + número correlativo
      grupos_multi  <- as.integer(names(conteo_grupos[conteo_grupos >= 2]))
      grupos_solo   <- as.integer(names(conteo_grupos[conteo_grupos == 1]))

      idx_multi <- setNames(seq_along(grupos_multi), as.character(grupos_multi))
      cnt_solo  <- 0L

      df$grupo_label <- NA_character_
      for (i in seq_len(nrow(df))) {
        gn <- df$grupo_num[i]
        if (is.na(gn)) next
        if (gn %in% grupos_multi) {
          df$grupo_label[i] <- paste0(nombre, idx_multi[as.character(gn)])
        } else {
          cnt_solo <- cnt_solo + 1L
          df$grupo_label[i] <- paste0(nombre_aislado, cnt_solo)
        }
      }
    } else {
      df$grupo_label <- NA_character_
    }

    df
  }, ignoreNULL = FALSE)

  # ── Gráfico ──────────────────────────────────────────────────────────────────
  output$vis_grafico <- renderPlot({
    df <- vis_datos()
    if (is.null(df)) return(NULL)

    n       <- nrow(df)
    col_vec <- vis_colores_ge[((seq_len(n)-1) %% length(vis_colores_ge)) + 1]
    names(col_vec) <- df$ge_label

    y_vals <- c(df$f0_ini, df$f0_fin)
    y_min  <- min(y_vals) * 0.89
    y_max  <- max(y_vals) * 1.10

    p <- ggplot(df)

    # Sombreado: solo para grupos con ≥ 2 GE (grupos_multi)
    grupos_v      <- df[!is.na(df$grupo_label), ]
    nombre_ais    <- trimws(input$vis_nombre_aislado %||% "singular")
    if (nchar(nombre_ais) == 0) nombre_ais <- "singular"
    grupos_multi_v <- grupos_v[!grepl(paste0("^", nombre_ais), grupos_v$grupo_label), ]
    etiq_multi     <- unique(grupos_multi_v$grupo_label)

    if (length(etiq_multi) > 0) {
      for (gi in seq_along(etiq_multi)) {
        gdf    <- grupos_multi_v[grupos_multi_v$grupo_label == etiq_multi[gi], ]
        cfondo <- vis_colores_grupo[((gi-1) %% length(vis_colores_grupo)) + 1]
        p <- p + annotate("rect",
          xmin = min(gdf$ge) - 0.46, xmax = max(gdf$ge) + 0.46,
          ymin = y_min, ymax = y_max, fill = cfondo, alpha = 0.08)
      }
    }

    p <- p +
      geom_segment(aes(x = x_ini, xend = x_fin,
                       y = f0_ini, yend = f0_fin, color = ge_label),
                   linewidth = 2.2, lineend = "round") +
      geom_point(aes(x = x_ini, y = f0_ini, fill = ge_label),
                 shape = 21, size = 3.2, color = "white", stroke = 0.8) +
      geom_point(aes(x = x_fin, y = f0_fin, fill = ge_label),
                 shape = 24, size = 3.2, color = "white", stroke = 0.8) +
      scale_color_manual(values = col_vec, name = "GE") +
      scale_fill_manual(values  = col_vec, name = "GE") +
      scale_x_continuous(breaks = seq_len(n), labels = df$ge_label,
                         limits = c(0.3, n + 0.7),
                         expand = expansion(mult = 0.03)) +
      scale_y_continuous(limits = c(y_min, y_max),
                         expand = expansion(mult = c(0.04, 0.20)),
                         labels = function(x) paste0(x, " Hz")) +
      labs(
        title    = "Progresión tonal manual — grupos entonativos",
        subtitle = sprintf(
          "○ = F0 ini  △ = F0 fin  |  Umbral ST: %.2f  |  Modo: %s",
          input$umbral_st %||% 0.8,
          if (input$vis_modo_grupos == "auto")
            paste0("auto (", input$vis_unidad_umbral %||% "st", " ",
                   input$vis_umbral_valor %||% 1.5, ")")
          else "manual"
        ),
        x = NULL, y = "F0 (Hz)", color = NULL, fill = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title         = element_text(face = "bold", size = 13),
        plot.subtitle      = element_text(color = "grey50", size = 9),
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x        = element_text(size = 9, color = "grey40"),
        legend.position    = if (n <= 10) "right" else "bottom",
        plot.background    = element_rect(fill = "white", color = NA)
      )

    # Etiquetas Hz
    if (isTRUE(input$vis_etiquetas) && n <= 12) {
      p <- p +
        geom_text(aes(x = x_ini - 0.05, y = f0_ini,
                      label = paste0(round(f0_ini), " Hz"), color = ge_label),
                  hjust = 1, size = 2.8, show.legend = FALSE) +
        geom_text(aes(x = x_fin + 0.05, y = f0_fin,
                      label = paste0(round(f0_fin), " Hz"), color = ge_label),
                  hjust = 0, size = 2.8, show.legend = FALSE)
    }

    # Etiquetas de grupo (solo grupos ≥ 2 GE)
    if (nrow(grupos_multi_v) > 0) {
      eg <- grupos_multi_v |>
        group_by(grupo_label) |>
        summarise(xc = mean(ge), .groups = "drop")
      p <- p + geom_text(data = eg,
                         aes(x = xc, y = y_max * 0.984, label = grupo_label),
                         color = "#d68910", size = 3.2,
                         fontface = "bold", vjust = 1, inherit.aes = FALSE)
    }

    # Etiquetas de GE aislados (nombre_ais + número)
    aislados_v <- grupos_v[grepl(paste0("^", nombre_ais), grupos_v$grupo_label), ]
    if (nrow(aislados_v) > 0) {
      p <- p + geom_text(data = aislados_v,
                         aes(x = ge, y = y_max * 0.975, label = grupo_label),
                         color = "grey55", size = 2.8,
                         fontface = "italic", vjust = 1, inherit.aes = FALSE)
    }

    # Líneas de tendencia
    add_smooth_vis <- function(plt, x_v, y_v, col, lty) {
      if (length(x_v) < 2) return(plt)
      df_t <- data.frame(x = x_v, y = y_v)
      plt + geom_smooth(data = df_t, aes(x = x, y = y),
                        method = "lm", formula = y ~ x,
                        color = col, fill = col,
                        linewidth = 0.9, linetype = lty,
                        alpha = 0.08, se = TRUE, inherit.aes = FALSE)
    }

    t_fin  <- isTRUE(input$vis_tend_finales)
    t_ini  <- isTRUE(input$vis_tend_iniciales)
    t_med  <- isTRUE(input$vis_tend_medias)
    modo_t <- input$vis_modo_tend %||% "global"

    if (modo_t == "global") {
      xp <- seq_len(n)
      if (t_fin && n >= 2) p <- add_smooth_vis(p, xp, df$f0_fin,   "grey40",  "dashed")
      if (t_ini && n >= 2) p <- add_smooth_vis(p, xp, df$f0_ini,   "#27ae60", "dotted")
      if (t_med && n >= 2) p <- add_smooth_vis(p, xp, df$f0_media, "#d68910", "solid")
    } else {
      if (nrow(grupos_multi_v) > 0) {
        for (etiq in etiq_multi) {
          sub_df <- grupos_multi_v[grupos_multi_v$grupo_label == etiq, ]
          if (nrow(sub_df) < 2) next
          xs <- sub_df$ge
          if (t_fin) p <- add_smooth_vis(p, xs, sub_df$f0_fin,   "grey40",  "dashed")
          if (t_ini) p <- add_smooth_vis(p, xs, sub_df$f0_ini,   "#27ae60", "dotted")
          if (t_med) p <- add_smooth_vis(p, xs, sub_df$f0_media, "#d68910", "solid")
        }
      }
    }

    print(p)
  }, bg = "white")

  # ── Resumen numérico ─────────────────────────────────────────────────────────
  output$vis_resumen <- renderUI({
    df <- vis_datos()
    if (is.null(df)) return(NULL)
    n_asc  <- sum(df$tendencia == "ascendente")
    n_desc <- sum(df$tendencia == "descendente")
    n_plan <- sum(df$tendencia == "plano")
    nombre_ais <- trimws(input$vis_nombre_aislado %||% "singular")
    grupos_v   <- df[!is.na(df$grupo_label), ]
    n_grup_multi <- length(unique(
      grupos_v$grupo_label[!grepl(paste0("^", nombre_ais), grupos_v$grupo_label)]))
    n_ais <- length(unique(
      grupos_v$grupo_label[grepl(paste0("^", nombre_ais), grupos_v$grupo_label)]))
    delta  <- round(df$f0_fin[nrow(df)] - df$f0_ini[1], 1)
    tend_g <- if (delta > 0) "↑ asc" else if (delta < 0) "↓ desc" else "→ plana"
    modo_txt <- if (input$vis_modo_grupos == "auto")
      paste0("auto · ", input$vis_unidad_umbral %||% "st",
             " ", input$vis_umbral_valor %||% 1.5)
    else "manual"
    div(class = "stat-card",
      style = "display:flex;gap:22px;flex-wrap:wrap;align-items:center;",
      div(div(class="stat-val",style="font-size:20px;",nrow(df)),
          div(class="stat-lbl","GE totales")),
      div(div(class="stat-val stat-asc",style="font-size:20px;",n_asc),
          div(class="stat-lbl","ascendentes")),
      div(div(class="stat-val stat-desc",style="font-size:20px;",n_desc),
          div(class="stat-lbl","descendentes")),
      div(div(class="stat-val",style="font-size:20px;",n_plan),
          div(class="stat-lbl","planos")),
      div(div(class="stat-val",style="font-size:20px;color:#d68910;",n_grup_multi),
          div(class="stat-lbl",paste0("grupos (", modo_txt, ")"))),
      div(div(class="stat-val",style="font-size:20px;color:#888;",n_ais),
          div(class="stat-lbl","aislados")),
      div(div(class="stat-val",style="font-size:17px;color:#8b9eff;",
              paste0(tend_g, " (", ifelse(delta >= 0, "+", ""), delta, " Hz)")),
          div(class="stat-lbl","tendencia conjunta"))
    )
  })

  # ── Leyenda detallada por GE ─────────────────────────────────────────────────
  output$vis_leyenda <- renderUI({
    df <- vis_datos()
    if (is.null(df)) return(NULL)

    # Cabecera
    cabecera <- div(
      style = "display:flex;gap:10px;padding:4px 0 8px;border-bottom:2px solid #e8e7e4;
               font-family:monospace;font-size:10px;color:#999;text-transform:uppercase;
               letter-spacing:0.5px;font-weight:700;",
      div(style="width:12px;flex-shrink:0;"),
      div(style="min-width:30px;","GE"),
      div(style="min-width:90px;","Dirección"),
      div(style="min-width:65px;","Δ Hz"),
      div(style="min-width:65px;","Δ %"),
      div(style="min-width:65px;","ST"),
      div(style="min-width:90px;","Reajuste Hz"),
      div(style="min-width:70px;","Reajuste %"),
      div(style="min-width:65px;","Reaj. ST"),
      div(style="margin-left:auto;","Grupo")
    )

    filas <- lapply(seq_len(nrow(df)), function(i) {
      row   <- df[i, ]
      color <- vis_colores_ge[((i-1) %% length(vis_colores_ge)) + 1]

      col_dir <- switch(row$tendencia,
        ascendente  = "#27ae60",
        descendente = "#e74c3c",
        "#f39c12")
      lbl_dir <- switch(row$tendencia,
        ascendente  = "↑ ascendente",
        descendente = "↓ descendente",
        "→ plano")

      delta_hz  <- row$f0_fin - row$f0_ini
      delta_txt <- paste0(ifelse(delta_hz >= 0, "+", ""), round(delta_hz, 1), " Hz")
      pct_txt   <- paste0(ifelse(row$inflexion_pct >= 0, "+", ""),
                          round(row$inflexion_pct, 1), "%")
      st_txt    <- paste0(ifelse(row$inflexion_st >= 0, "+", ""),
                          row$inflexion_st, " ST")

      rea_hz  <- if (is.na(row$reajuste_hz))  "—"
                 else paste0(ifelse(row$reajuste_hz  >= 0, "+", ""),
                              row$reajuste_hz,  " Hz")
      rea_pct <- if (is.na(row$reajuste_pct)) "—"
                 else paste0(ifelse(row$reajuste_pct >= 0, "+", ""),
                              row$reajuste_pct, "%")
      rea_st  <- if (is.na(row$reajuste_st))  "—"
                 else paste0(ifelse(row$reajuste_st  >= 0, "+", ""),
                              row$reajuste_st,  " ST")

      grupo_txt <- if (!is.na(row$grupo_label)) row$grupo_label else "—"
      col_grupo <- if (!is.na(row$grupo_label)) "#d68910" else "#bbb"

      div(style = "display:flex;align-items:center;gap:10px;
                   padding:5px 0;border-bottom:1px solid #f0efe9;
                   font-family:monospace;font-size:11px;",
        div(style = paste0("width:12px;height:12px;border-radius:50%;
                            background:", color, ";flex-shrink:0;")),
        div(style = "min-width:30px;color:#1a1a1a;font-weight:700;", row$ge_label),
        div(style = paste0("color:", col_dir, ";min-width:90px;"),   lbl_dir),
        div(style = "color:#555;min-width:65px;",  delta_txt),
        div(style = "color:#555;min-width:65px;",  pct_txt),
        div(style = "color:#555;min-width:65px;",  st_txt),
        div(style = "color:#888;min-width:90px;",  rea_hz),
        div(style = "color:#888;min-width:70px;",  rea_pct),
        div(style = "color:#888;min-width:65px;",  rea_st),
        div(style = paste0("color:", col_grupo, ";margin-left:auto;font-weight:600;"),
            grupo_txt)
      )
    })

    div(class = "main-box", style = "padding:14px 18px;margin-top:4px;",
      div(class = "sec-label", style = "margin-top:0;",
          "Detalle por grupo entonativo"),
      cabecera,
      filas
    )
  })


}

shinyApp(ui, server)