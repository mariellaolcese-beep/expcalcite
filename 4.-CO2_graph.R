library(tidyverse)

# -----------------------------------------------------------------------------------------------------
#  0. AJUSTES
# -----------------------------------------------------------------------------------------------------
archivo      <- "C:/Experimentcalcite/data/planilla_total.csv"  # planilla total (Partes A, B y C)
carpeta_out  <- "C:/Experimentcalcite/resultados"
codificacion <- "UTF-8"

tiempos_mostrar <- NULL          # NULL = todos los tiempos; o por ejemplo c(0, 1) para los dos primeros
etiq_tiempo <- c(`0` = "T0", `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4")

col_D <- "#222222"   # Dry (negro)
col_W <- "#E06D9B"   # Wet (rosa medio)
col_trt <- c(Dry = col_D, Wet = col_W)

etiq_site  <- c(A = "Low OM", B = "High OM")
etiq_stage <- c(initial = "Initial (I)", final = "Final (F)")

tema_poster <- theme_minimal(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 15),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", size = 14),
        legend.title = element_text(face = "bold", size = 14),
        legend.text = element_text(size = 13),
        legend.position = "top",
        panel.spacing = unit(1, "lines"),
        panel.grid.minor = element_blank())

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------------------------------
#  1. LEER Y LIMPIAR LA PLANILLA TOTAL 
# -----------------------------------------------------------------------------------------------------
lineas <- readLines(archivo, n = 20, warn = FALSE, encoding = codificacion)
fila_header <- grep("sample_id", lineas)[1]
if (is.na(fila_header)) stop("No se encontró 'sample_id' en las primeras 20 líneas del archivo.")

linea_hdr <- lineas[fila_header]
sep <- if (grepl(";", linea_hdr)) ";" else if (grepl("\t", linea_hdr)) "\t" else ","

pl <- read_delim(archivo, delim = sep, skip = fila_header - 1,
                 col_types = cols(.default = col_character()),
                 locale = locale(encoding = codificacion), na = c("", "NA"), show_col_types = FALSE)
pl <- pl[, !grepl("^\\.\\.\\.[0-9]+$", names(pl))]

# Detectar nombres posibles de las columnas de CO2 e isótopos
col_ppm   <- intersect(c("CO2_corrected_ppm", "CO2_ppm"), names(pl))[1]
col_d_co2 <- intersect(c("13_CO2_δ", "13_CO2_?", "13_CO2"), names(pl))[1]

if (is.na(col_ppm))   stop("No encuentro la columna de CO2 en ppm (CO2_corrected_ppm o CO2_ppm).")
if (is.na(col_d_co2)) stop("No encuentro la columna del δ13C del CO2 (13_CO2_δ o 13_CO2).")

faltan <- setdiff(c("sample_id", "extraction_time", "stage", "treatment", "site", "mmol de CO2"), names(pl))
if (length(faltan) > 0) stop("Faltan columnas en la planilla: ", paste(faltan, collapse = ", "))

datos_co2 <- pl %>%
  transmute(
    sample_id,
    tiempo   = as.integer(gsub("\\D", "", extraction_time)),
    stage    = factor(tolower(stage), levels = c("initial", "final")),
    trt      = factor(treatment, levels = c("Dry", "Wet")),
    site     = factor(site, levels = c("A", "B")),
    CO2_ppm  = num(.data[[col_ppm]]),
    mmol_CO2 = num(`mmol de CO2`),
    d_CO2    = num(.data[[col_d_co2]])
  ) %>%
  filter(!is.na(site), !is.na(trt), !is.na(tiempo), !is.na(stage))

if (!is.null(tiempos_mostrar)) {
  datos_co2 <- datos_co2 %>% filter(tiempo %in% tiempos_mostrar)
}

# -----------------------------------------------------------------------------------------------------
#  2. RESUMEN Y GRÁFICOS DE BARRAS (A: CO2 ppm , B: mmol CO2 y δ13C-CO2)
# -----------------------------------------------------------------------------------------------------

###DENTRO DE LA FUNCION BARRA CAMBIO LOS TITULOS, NOMBRES Y ESAS COSAS 
resumen_co2 <- function(var) {
  datos_co2 %>%
    filter(!is.na(.data[[var]])) %>%
    group_by(tiempo, site, trt, stage) %>%
    summarise(m = mean(.data[[var]]), sd = sd(.data[[var]]), n = n(), .groups = "drop") %>%
    mutate(
      site   = factor(site, levels = c("A", "B"), labels = unname(etiq_site)),
      stage  = factor(stage, levels = c("initial", "final"), labels = unname(etiq_stage)),
      tiempo = factor(tiempo, levels = as.integer(names(etiq_tiempo)), labels = unname(etiq_tiempo)) %>% droplevels()
    )
}

barras <- function(df, ylab, titulo, linea_cero = TRUE) {
  g <- ggplot(df, aes(x = site, y = m, fill = trt))
  if (linea_cero) g <- g + geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3)
  g +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd),
                  position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6) +
    facet_grid(stage ~ tiempo) +
    scale_fill_manual(values = col_trt, name = "Treatment") +
    labs(title = titulo, x = "Site", y = ylab) +
    tema_poster
}

# A. CO2 en ppm 
p_co2  <- barras(resumen_co2("CO2_ppm"),  expression(CO[2]*" (ppm)"),
                 "CO2_HS_concentration_temporallity", linea_cero = FALSE)

# B. mmol de CO2 y δ13C-CO2 (‰)
p_mmol <- barras(resumen_co2("mmol_CO2"), expression(CO[2]*" (mmol)"),
                 "CO2 in headspace (mmol)")
p_d13  <- barras(resumen_co2("d_CO2"),    expression(delta^{13}*"C-CO"[2]*" (\u2030)"),
                 "δ13C of the CO2 in each Site, Treatment and Time")

cat("\nFrascos con CO2 (ppm/mmol) o δ13C-CO2 vacío (revisar):\n")
datos_co2 %>% filter(is.na(CO2_ppm) | is.na(mmol_CO2) | is.na(d_CO2)) %>% print(n = 30)

# -----------------------------------------------------------------------------------------------------
#  3. DISPERSIÓN: CO2 (mmol) vs δ13C-CO2 (‰), frasco a frasco
# -----------------------------------------------------------------------------------------------------
p_disp <- datos_co2 %>%
  filter(!is.na(mmol_CO2), !is.na(d_CO2)) %>%
  ggplot(aes(mmol_CO2, d_CO2, color = trt, shape = stage)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ site, labeller = labeller(site = etiq_site)) +
  scale_color_manual(values = col_trt, name = "Treatment") +
  scale_shape_manual(values = c(initial = 1, final = 16), labels = etiq_stage, name = "Stage") +
  labs(x = expression(CO[2]*" (mmol)"), y = expression(delta^{13}*"C-CO"[2]*" (\u2030)"),
       title = "CO2 vs δ13C-CO2 (one point per bottle)") +
  tema_poster

# -----------------------------------------------------------------------------------------------------
#  4. MOSTRAR EN PANTALLA Y GUARDAR
# -----------------------------------------------------------------------------------------------------
graphics.off()  # limpia cualquier ventana gráfica oculta previa
for (p in list(p_co2, p_mmol, p_d13, p_disp)) {
  print(p)
  Sys.sleep(0.2)
}

ggsave(file.path(carpeta_out, "CO2_fig1_ppm.png"),        p_co2,  width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig2_mmol.png"),       p_mmol, width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig3_d13C.png"),       p_d13,  width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig4_dispersion.png"), p_disp, width = 11, height = 6, dpi = 300)

write_csv(datos_co2, file.path(carpeta_out, "CO2_isotopos_por_frasco.csv"))
