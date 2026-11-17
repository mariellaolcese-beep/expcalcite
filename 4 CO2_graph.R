library(ggplot2)
library(dplyr)

# Paleta consistente
col_D <- "#B08968"  # Dry
col_W <- "#4A7C82"  # Wet

# 1. Cargar datos
co2_data <- read.csv("co2_data.csv", fileEncoding = "UTF-8-BOM")

# 2. Asignar nombres limpios a las columnas
names(co2_data) <- c("SAMPLE", "TIME", "SITE", "STAGE", "TREATMENT", "CO2_ppm")

# Identificar los dos primeros tiempos únicos que existen en tus datos
tiempos_unicos <- sort(unique(co2_data$TIME))[1:2]

# 3. Resumen, filtrado automático de los 2 primeros tiempos y Mutate
co2_summary <- co2_data %>%
  filter(TIME %in% tiempos_unicos) %>%  # Toma automáticamente los dos primeros tiempos
  group_by(TIME, SITE, TREATMENT, STAGE) %>%
  summarise(
    mean_CO2 = mean(CO2_ppm, na.rm = TRUE),
    sd_CO2 = sd(CO2_ppm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SITE = factor(SITE, levels = c("A", "B"), labels = c("site low OM", "site high OM")),
    TIME = factor(TIME, levels = tiempos_unicos, labels = c("Day 0 (T0)", "Day 14 (T1)")),
    STAGE = factor(STAGE, levels = c("F", "I"), labels = c("Final (F)", "Initial (I)"))
  )

# 4. Gráfico ajustado con ambos tiempos y texto legible
p_co2 <- ggplot(co2_summary, aes(x = SITE, y = mean_CO2, fill = TREATMENT)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_CO2 - sd_CO2, ymax = mean_CO2 + sd_CO2),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6
  ) +
  facet_grid(STAGE ~ TIME) +
  scale_fill_manual(values = c(D = col_D, W = col_W), name = "Treatment") +
  labs(
    title = "CO2 concentration in each Site, Treatment and Time",
    x = "Site",
    y = expression(CO[2]*" (ppm)")
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # Título centrado
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13, color = "black"),
    strip.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 14),
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# Guardar la imagen
ggsave("poster_panels_co2.png", p_co2, width = 11, height = 8, dpi = 300)

library(tidyverse)

# -----------------------------------------------------------------------------------------------------
#  0. AJUSTES
# -----------------------------------------------------------------------------------------------------
archivo_co2  <- "co2_data.csv"                                  # tu archivo de ppm (Parte A)
archivo      <- "C:/Experimentcalcite/data/planilla_total.csv"  # planilla total (Partes B y C)
carpeta_out  <- "C:/Experimentcalcite/resultados"
codificacion <- "UTF-8"

tiempos_mostrar <- NULL          # NULL = todos los tiempos; o por ejemplo c(0, 1) para los dos primeros
etiq_tiempo <- c(`0` = "T0", `1` = "T1", `2` = "T2", `3` = "T3", `4` = "T4")   # cámbialas por "Day 0 (T0)", etc., si quieres

col_D <- "#B08968"   # Dry
col_W <- "#4A7C82"   # Wet
col_trt <- c(Dry = col_D, Wet = col_W)
etiq_site  <- c(A = "site low OM", B = "site high OM")
etiq_stage <- c(initial = "Initial (I)", final = "Final (F)")

tema_poster <- theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        axis.title = element_text(face = "bold", size = 16),
        axis.text = element_text(size = 13, color = "black"),
        strip.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.text = element_text(size = 14),
        legend.position = "top", panel.grid.minor = element_blank())

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------------------------------
#  A. CO2 (ppm): tu figura
# -----------------------------------------------------------------------------------------------------
planilla_tota <- read.csv(archivo_co2, fileEncoding = "UTF-8-BOM")
names(co2_data) <- c("SAMPLE", "TIME", "SITE", "STAGE", "TREATMENT", "CO2_ppm")

tiempos_A <- if (is.null(tiempos_mostrar)) sort(unique(co2_data$TIME)) else tiempos_mostrar

co2_summary <- co2_data %>%
  filter(TIME %in% tiempos_A) %>%
  group_by(TIME, SITE, TREATMENT, STAGE) %>%
  summarise(mean_CO2 = mean(CO2_ppm, na.rm = TRUE), sd_CO2 = sd(CO2_ppm, na.rm = TRUE), .groups = "drop") %>%
  mutate(SITE  = factor(SITE, levels = c("A", "B"), labels = unname(etiq_site)),
         TIME  = factor(TIME, levels = tiempos_A, labels = etiq_tiempo[as.character(tiempos_A)]),
         STAGE = factor(STAGE, levels = c("I", "F"), labels = c("Initial (I)", "Final (F)")),
         TREATMENT = factor(TREATMENT, levels = c("D", "W")))

p_co2 <- ggplot(co2_summary, aes(x = SITE, y = mean_CO2, fill = TREATMENT)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_CO2 - sd_CO2, ymax = mean_CO2 + sd_CO2),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6) +
  facet_grid(STAGE ~ TIME) +
  scale_fill_manual(values = c(D = col_D, W = col_W), name = "Treatment") +
  labs(title = "CO2 concentration in each Site, Treatment and Time", x = "Site", y = expression(CO[2]*" (ppm)")) +
  tema_poster

# -----------------------------------------------------------------------------------------------------
#  B. ISÓTOPOS: δ13C-CO2 (‰) y mmol de CO2, desde la planilla total
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

col_d_co2 <- intersect(c("13_CO2_δ", "13_CO2"), names(pl))[1]       # el δ del CO2 puede tener dos nombres
if (is.na(col_d_co2)) stop("No encuentro la columna del δ13C del CO2 (13_CO2_δ o 13_CO2).")
faltan <- setdiff(c("sample_id", "extraction_time", "stage", "treatment", "site", "mmol de CO2"), names(pl))
if (length(faltan) > 0) stop("Faltan columnas en la planilla: ", paste(faltan, collapse = ", "))

iso <- pl %>%
  transmute(sample_id,
            tiempo = as.integer(gsub("\\D", "", extraction_time)),
            stage  = factor(tolower(stage), levels = c("initial", "final")),
            trt    = factor(treatment, levels = c("Dry", "Wet")),
            site   = factor(site, levels = c("A", "B")),
            mmol_CO2 = num(`mmol de CO2`),
            d_CO2    = num(.data[[col_d_co2]])) %>%
  filter(!is.na(site), !is.na(trt), !is.na(tiempo), !is.na(stage))

if (!is.null(tiempos_mostrar)) iso <- iso %>% filter(tiempo %in% tiempos_mostrar)

resumen_iso <- function(var) {
  iso %>% filter(!is.na(.data[[var]])) %>%
    group_by(tiempo, site, trt, stage) %>%
    summarise(m = mean(.data[[var]]), sd = sd(.data[[var]]), n = n(), .groups = "drop") %>%
    mutate(site  = factor(site, levels = c("A", "B"), labels = unname(etiq_site)),
           stage = factor(stage, levels = c("initial", "final"), labels = unname(etiq_stage)),
           tiempo = factor(tiempo, levels = as.integer(names(etiq_tiempo)), labels = unname(etiq_tiempo)) %>% droplevels())
}

barras <- function(df, ylab, titulo) {
  ggplot(df, aes(x = site, y = m, fill = trt)) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6) +
    facet_grid(stage ~ tiempo) +
    scale_fill_manual(values = col_trt, name = "Treatment") +
    labs(title = titulo, x = "Site", y = ylab) + tema_poster
}

p_mmol <- barras(resumen_iso("mmol_CO2"), expression(CO[2]*" (mmol)"), "CO2 in the headspace (mmol)")
p_d13  <- barras(resumen_iso("d_CO2"), expression(delta^{13}*"C-CO"[2]*" (\u2030)"),
                 "δ13C of the CO2 in each Site, Treatment and Time")

cat("\nFrascos con CO2 o δ13C-CO2 vacío (revisar):\n")
iso %>% filter(is.na(mmol_CO2) | is.na(d_CO2)) %>% print(n = 30)

# -----------------------------------------------------------------------------------------------------
#  C. DISPERSIÓN: CO2 (mmol) vs δ13C-CO2 (‰), frasco a frasco
#     Un δ más negativo = más carbono de respiración (≈ -25 ‰); más positivo = carbono de calcita disuelta
#     o aire de la atmósfera (≈ -8 ‰).
# -----------------------------------------------------------------------------------------------------
p_disp <- iso %>% filter(!is.na(mmol_CO2), !is.na(d_CO2)) %>%
  ggplot(aes(mmol_CO2, d_CO2, color = trt, shape = stage)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ site, labeller = labeller(site = etiq_site)) +
  scale_color_manual(values = col_trt, name = "Treatment") +
  scale_shape_manual(values = c(initial = 1, final = 16), labels = etiq_stage, name = "Stage") +
  labs(x = expression(CO[2]*" (mmol)"), y = expression(delta^{13}*"C-CO"[2]*" (\u2030)"),
       title = "CO2 vs δ13C-CO2 (one point per bottle)") + tema_poster

# -----------------------------------------------------------------------------------------------------
#  D. MOSTRAR Y GUARDAR
# -----------------------------------------------------------------------------------------------------
print(p_co2); print(p_mmol); print(p_d13); print(p_disp)

ggsave(file.path(carpeta_out, "CO2_fig1_ppm.png"),         p_co2,  width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig2_mmol.png"),        p_mmol, width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig3_d13C.png"),        p_d13,  width = 11, height = 8, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig4_dispersion.png"),  p_disp, width = 11, height = 6, dpi = 300)
write_csv(iso, file.path(carpeta_out, "CO2_isotopos_por_frasco.csv"))


