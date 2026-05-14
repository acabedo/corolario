# Corolario

**Análisis prosódico de progresiones de F0 por cláusula entonativa**

Aplicación Shiny para el estudio de la entonación a partir de datos de frecuencia fundamental (F0). Permite detectar cláusulas entonativas, analizar transiciones entre grupos tonales, calcular correlaciones entre reajuste y pausa, e identificar patrones de coincidencia entre enunciados y cláusulas.

**Demo:** [acabedo.github.io/corolario](https://acabedo.github.io/corolario)

---

## Descripción

Corolario toma como entrada una tabla de grupos tonales (F0 inicial, F0 final, pausa) y calcula automáticamente:

- La **dirección entonativa** de cada grupo mediante regresión de slope sobre una ventana deslizante
- La **cláusula entonativa** a la que pertenece cada grupo (agrupación de runs consecutivos con la misma dirección), con corrección por trayectoria interleaved y detección de patrones valle
- El **reajuste de F0** en cada juntura entre grupos: diferencia entre F0 inicial del grupo actual y F0 final del grupo anterior
- La **inflexión tonal** interna de cada grupo
- La **progresión acumulada**: suma de reajustes en la dirección actual hasta superar un umbral configurable
- Estadísticos de transición, correlaciones reajuste–pausa y análisis de coincidencia cláusula–enunciado

---

## Instalación

### Requisitos

- R ≥ 4.2
- Paquetes:

```r
install.packages(c("shiny", "tidyverse", "zoo", "ggplot2", "scales", "DT", "readr"))
```

### Ejecutar localmente

```r
# Clonar el repositorio
# git clone https://github.com/acabedo/corolario.git

# Desde R:
shiny::runApp("ruta/a/corolario")
```

---

## Formato de datos

La entrada es una tabla separada por tabulaciones (`.txt`, `.csv` o `.tsv`) con la siguiente estructura:

```
enunciados   F01   F02   pausa1   F03   F04   pausa2   ...
enun1        120   135   0.3      128   145   0.5
enun2        138   155   0.4
```

Cada **trío de columnas** `(F0_ini, F0_fin, pausa)` corresponde a un grupo tonal. Las celdas vacías indican que el grupo no existe para ese enunciado.

- **F0_ini / F0_fin**: frecuencia fundamental en Hz al inicio y al final del grupo
- **pausa**: duración en segundos de la pausa al final del grupo (`0` = sin pausa)

Los datos pueden pegarse directamente en el campo de texto de la pestaña `01 | Datos` o cargarse mediante el botón de subida de archivo.

---

## Uso

### Barra de parámetros (global)

| Parámetro | Descripción | Valor por defecto |
|---|---|---|
| Ventana slope | Número de grupos en la ventana deslizante de regresión | 3 |
| Umbral % / enun. | Umbral mínimo de slope (%) para etiquetar dirección | 5 |
| Umbral ST / enun. | Umbral mínimo de slope (ST) para etiquetar dirección | 0.8 |
| Umbral local ST | Umbral para dirección entre grupos consecutivos y para detección de valles | 0.5 |
| Nivel slope | Unidad de análisis del slope: grupo o enunciado | Grupo |
| Prog. umbral ST | Acumulación mínima en ST para declarar progresión | 2.0 |
| Prog. umbral % | Acumulación mínima en % para declarar progresión | 15 |

---

### Pestañas

#### 01 \| Datos
Entrada de datos por pegado de texto o subida de archivo. Incluye un ejemplo de formato en el campo de texto.

#### 02 \| Tabla
Tabla completa con un grupo tonal por fila. Incluye todas las variables calculadas y un glosario descriptivo en la parte inferior.

Variables principales:

| Variable | Descripción |
|---|---|
| F0 ini / F0 fin | Frecuencia fundamental en Hz al inicio y fin del grupo |
| Inflexion % / ST | Movimiento interno del grupo: `(F0_fin − F0_ini) / F0_ini` |
| Reajuste % / ST | Salto en la juntura entre grupos: `F0_ini(n) − F0_fin(n−1)` |
| Acum. ST | Suma acumulada de reajuste ST desde el último cambio de dirección |
| Progresion | `prog_asc` / `prog_desc` / `neutro` según supere el umbral configurado |
| dST | Delta entre F0_ini consecutivos; solo para Dir. local |
| Slope ST | Pendiente de regresión sobre la trayectoria interleaved de la cláusula |
| Dir. slope | Dirección definitiva del grupo (ascendente / descendente / plano) |
| Clausula | Cláusula entonativa asignada |

#### 03 \| Gráfico
Visualización de las progresiones de F0 por cláusula. Muestra F0 ini (círculo) y F0 fin (triángulo) de cada grupo, conectados por líneas de transición. Las flechas y etiquetas indican la dirección de cada cláusula.

- Control de rango de grupos visualizados
- Exportación a PNG
- Tabla de valores correspondientes a los grupos del gráfico

#### 04 \| Cláusulas
Distribución de cláusulas entonativas: totales, ascendentes, descendentes y planas. Detalle por cláusula y comparación entre slope ST y slope %.

#### 05 \| Transiciones
Análisis de reajustes y pausas entre grupos consecutivos.

- **Por cláusula**: boxplots de reajuste (Hz, %, ST) y pausa agrupados por dirección entonativa
- **Global**: medias generales, conteo de pausas = 0 y correlación entre reajuste y pausa (solo sobre transiciones con pausa > 0), con gráfico de dispersión y coeficientes de Pearson y Spearman

#### 06 \| Chi cuadrado
Análisis de coincidencia entre enunciados y cláusulas.

- **Gráfico de barras**: porcentaje de cláusulas según el número de enunciados que las componen (pura / mixta con cruce en frontera / mixta sin cruce directo)
- **Test chi cuadrado** (Monte Carlo, 5000 simulaciones) con V de Cramer
- **Tabla de contingencia** enunciado × cláusula
- **Tabla de cláusulas mixtas** con clasificación del tipo de cruce

---

## Algoritmo de detección de cláusulas

La asignación de cláusulas sigue un proceso de cuatro pasos:

1. **Slope deslizante** sobre F0_ini con ventana configurable (regresión lineal normalizada en ST)
2. **Asignación provisional** de cláusulas por runs consecutivos de igual dirección, con detección de **patrones valle**: si un grupo sube internamente y el anterior bajaba (ambos superando el umbral local), se fuerza una ruptura de cláusula
3. **Corrección por cláusula**: para cláusulas de 2+ grupos, el slope se recalcula sobre la trayectoria interleaved completa `[F0_ini₁, F0_fin₁, F0_ini₂, F0_fin₂, …]`, capturando tanto el contorno interno de cada grupo como la tendencia global
4. **Fallback**: las cláusulas de un único grupo usan la dirección del contorno interno (inflexión) en lugar del slope

La ruptura por patrón valle es **asimétrica**: los arcos (sube→baja dentro de la misma cláusula) no se parten, ya que representan unidades fraseales completas.

---

## Tecnología

- [R](https://www.r-project.org/) + [Shiny](https://shiny.posit.co/)
- [tidyverse](https://www.tidyverse.org/) · [ggplot2](https://ggplot2.tidyverse.org/) · [DT](https://rstudio.github.io/DT/) · [zoo](https://cran.r-project.org/package=zoo)
- Demo desplegada con [Shinylive](https://posit-dev.github.io/r-shinylive/) vía GitHub Pages

---

## Licencia

MIT © [acabedo](https://github.com/acabedo)
