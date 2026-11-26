# ======================================================================================================================
#  CO2 DEL HEADSPACE Y d13C-CO2: LARGO PLAZO (initial), CAMBIO EN 24 h (final - initial), TABLAS, ESTADISTICA Y KEELING (SOLO PARA EXPLORAR, NO ESTA VERIFICADO POR MI)
#  (Dry y Wet; sitio A = MO baja, sitio B = MO alta)
# ======================================================================================================================
#
#  DE QUE SE TRATA
#  Mira el gas (headspace) de las botellas: cuanto CO2 hay y de que "tipo" de carbono es (d13C-CO2), y como cambia
#  (a) a lo largo del experimento y (b) dentro de las 24 h cerradas. Compara Dry con Wet y A con B.
#
#  QUE NECESITA
#  planilla_total.csv con: sample_id, extraction_time, stage, treatment, site, "mmol de CO2", una columna de CO2 en ppm
#  (CO2_corrected_ppm o CO2_ppm) y una del d13C del CO2 (13_CO2_δ, 13_CO2_? o 13_CO2). Paquetes: tidyverse, car, emmeans.
#  Origen de los datos (confirma en tu planilla): CO2_corrected_ppm = CO2 del headspace ya corregido; "mmol de CO2" = cantidad
#  de CO2 del headspace; 13_CO2 = d13C medido por el laboratorio de isotopos.
#  No usa constantes fisicas. Ajustes: tiempos_mostrar (NULL = todos) y x_keeling ("CO2_ppm" o "mmol_CO2").
#
#  PASOS
#    0. Ajustes       rutas, tiempos, colores
#    1. Leer          lee la planilla y detecta las columnas; formato largo (una fila por botella y variable)
#    2. Delta         DELTA = cada botella final - promedio de las initial de su tiempo, sitio y tratamiento
#    3. Figuras       Fig. 1-2 CO2 ppm (largo plazo initial / delta 24 h); Fig. 3-4 CO2 mmol; Fig. 5-6 d13C-CO2;
#                    
#    4. Tablas        resumen por grupo (initial, final, delta, SD combinada, "distinguible"), Dry vs Wet en el delta 
#    5. Estadistica   Shapiro, Levene, ANOVA sitio x tiempo y Tukey entre tiempos, por variable, tratamiento y etapa
#                     (initial, final y delta)
#   #    7. Guardar       figuras y CSV en carpeta_out
#
#  COMO LEER EL KEELING (ESTO ME FALTA CHEQUEAR SI ES REAL O NO AQUI LA IA SE FLIPO)
#  Si el CO2 de una botella es una mezcla de un CO2 de fondo (aire, constante) y un CO2 añadido (respiracion y/o calcita), los
#  puntos (1/CO2, d13C) caen en una recta. Su INTERCEPTO (donde 1/CO2 = 0, es decir CO2 infinito) es el d13C del CO2 añadido:
#    - cerca de -25 permil: la fuente es sobre todo respiracion
#    - cerca de -5 permil o mas alto: la fuente es sobre todo calcita (nativa y/o marcada)
#  Se calcula con botellas initial + final juntas (initial = poco CO2, final = mas CO2). Se da por grupo (tratamiento x sitio) y
#  por tiempo. Mira el R2 y el intervalo de confianza: si el rango de CO2 es pequeño, el intercepto sale muy incierto.
#
#  LIMITES A TENER EN CUENTA
#    - El Keeling supone UNA fuente añadida y un fondo constante. Aqui la fuente mezcla respiracion + calcita nativa + calcita
#      marcada, asi que el intercepto es el d13C PROMEDIO del CO2 añadido, no el de la respiracion sola. Para el fondo "limpio"
#      del CO2 habria que hacerlo con sedimento SIN marcar.
#    - El delta del d13C (final - initial) mezcla cantidad y composicion; leerlo junto al CO2.
#    - ppm y mmol son equivalentes solo si el headspace tiene el mismo volumen y temperatura en todas las botellas.
#    - Con n = 4 por grupo, un solo dato raro mueve el promedio.
#    - Aqui no hay calcita marcada disuelta: para eso usa CO2_total_y_calcita.R.
# ======================================================================================================================

library(tidyverse)
library(car)        # leveneTest
library(emmeans)    # Tukey

# -----------------------------------------------------------------------------------------------------
#  0. AJUSTES
# -----------------------------------------------------------------------------------------------------
archivo      <- "C:/Experimentcalcite/data/planilla_total.csv"
carpeta_out  <- "C:/Experimentcalcite/resultados"
codificacion <- "UTF-8"

tiempos_mostrar <- NULL          # NULL = todos los tiempos; o por ejemplo c(0, 1)
x_keeling       <- "CO2_ppm"     # variable para el Keeling: "CO2_ppm" o "mmol_CO2"

col_site   <- c(A = "#1A1A1A", B = "#D81B60")       # Negro = Sitio A, Rosa = Sitio B
etiq_site  <- c(A = "A (Low OM)", B = "B (High OM)")
col_trt    <- c(Dry = "#222222", Wet = "#E06D9B")
etiq_stage <- c(initial = "Initial (I)", final = "Final (F)")

tema <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 12),
    axis.title       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 10),
    plot.caption     = element_text(size = 9),
    panel.spacing    = unit(1.2, "lines")
  )

num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
dir.create(carpeta_out, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------------------------------
#  1. LEER Y LIMPIAR LA PLANILLA
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
  filter(!is.na(site), !is.na(trt), !is.na(tiempo), !is.na(stage)) %>%
  mutate(TIME = factor(tiempo, levels = 0:4))

if (!is.null(tiempos_mostrar)) datos_co2 <- datos_co2 %>% filter(tiempo %in% tiempos_mostrar)

cat("\nFrascos con CO2 (ppm/mmol) o δ13C-CO2 vacío (revisar):\n")
datos_co2 %>% filter(is.na(CO2_ppm) | is.na(mmol_CO2) | is.na(d_CO2)) %>% print(n = 30)

# Formato largo: una fila por frasco y variable
nombres_var <- c(CO2_ppm = "CO2 (ppm)", mmol_CO2 = "CO2 (mmol)", d_CO2 = "δ13C-CO2 (‰)")
datos_l <- datos_co2 %>%
  pivot_longer(c(CO2_ppm, mmol_CO2, d_CO2), names_to = "var", values_to = "v") %>%
  mutate(var = factor(nombres_var[var], levels = unname(nombres_var)))

cat("\nFrascos con dato por variable, tratamiento y etapa:\n")
datos_l %>% filter(!is.na(v)) %>% count(var, trt, stage) %>% print(n = Inf)

# -----------------------------------------------------------------------------------------------------
#  2. DELTA 24 h: cada frasco final - promedio de los initial de su tiempo, sitio, tratamiento y variable
# -----------------------------------------------------------------------------------------------------
base_ini <- datos_l %>%
  filter(stage == "initial") %>%
  group_by(var, tiempo, site, trt) %>%
  summarise(m_ini = mean(v, na.rm = TRUE), .groups = "drop")

datos_delta <- datos_l %>%
  filter(stage == "final") %>%
  left_join(base_ini, by = c("var", "tiempo", "site", "trt")) %>%
  mutate(delta = v - m_ini)

# -----------------------------------------------------------------------------------------------------
#  3. FIGURAS (ROMBO = PROMEDIO + PUNTOS INDIVIDUALES; columnas = Dry y Wet)
# -----------------------------------------------------------------------------------------------------
eje_t <- scale_x_continuous(breaks = 0:4, labels = paste0("t", 0:4))

fig_largo <- function(var_sel, ylab, titulo) {
  pts  <- datos_l %>% filter(stage == "initial", var == var_sel, !is.na(v))
  prom <- pts %>% group_by(tiempo, site, trt) %>% summarise(m = mean(v), .groups = "drop")
  ggplot() +
    geom_line(data = prom, aes(tiempo, m, color = site, group = site), linewidth = 0.65) +
    geom_point(data = pts, aes(tiempo, v, color = site), size = 1.9, alpha = 0.55) +
    geom_point(data = prom, aes(tiempo, m, color = site), shape = 18, size = 4.2) +
    facet_wrap(~ trt, nrow = 1, scales = "free_y") +
    eje_t +
    scale_color_manual(values = col_site, labels = etiq_site, name = "Sitio") +
    labs(x = NULL, y = ylab, title = titulo,
         subtitle = "Frascos initial (headspace al empezar la ventana de 24 h). Puntos pequeños = cada frasco; rombo = promedio",
         caption = "Si los puntos de un tiempo se agrupan alrededor del rombo, el cambio es consistente entre réplicas") +
    tema
}

fig_delta <- function(var_sel, ylab, titulo) {
  pts  <- datos_delta %>% filter(var == var_sel, !is.na(delta))
  prom <- pts %>% group_by(tiempo, site, trt) %>% summarise(m = mean(delta), .groups = "drop")
  ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.45) +
    geom_line(data = prom, aes(tiempo, m, color = site, group = site), linewidth = 0.65) +
    geom_point(data = pts, aes(tiempo, delta, color = site), size = 1.9, alpha = 0.55) +
    geom_point(data = prom, aes(tiempo, m, color = site), shape = 18, size = 4.2) +
    facet_wrap(~ trt, nrow = 1, scales = "free_y") +
    eje_t +
    scale_color_manual(values = col_site, labels = etiq_site, name = "Sitio") +
    labs(x = NULL, y = ylab, title = titulo,
         subtitle = "Cada punto pequeño = un frasco final menos el promedio de los initial",
         caption = "Puntos a ambos lados de la línea de 0 = el signo no es consistente entre réplicas") +
    tema
}

p_ppm_largo  <- fig_largo(nombres_var["CO2_ppm"],  "CO2 (ppm)",  "CO2 del headspace a largo plazo (ppm)")
p_ppm_delta  <- fig_delta(nombres_var["CO2_ppm"],  "Δ CO2 final - initial (ppm en 24 h)",  "CO2 acumulado en 24 h (ppm)")
p_mmol_largo <- fig_largo(nombres_var["mmol_CO2"], "CO2 (mmol)", "CO2 del headspace a largo plazo (mmol)")
p_mmol_delta <- fig_delta(nombres_var["mmol_CO2"], "Δ CO2 final - initial (mmol en 24 h)", "CO2 acumulado en 24 h (mmol)")
p_d13_largo  <- fig_largo(nombres_var["d_CO2"],    "δ13C-CO2 (‰)", "δ13C del CO2 a largo plazo")
p_d13_delta  <- fig_delta(nombres_var["d_CO2"],    "Δ δ13C-CO2 final - initial (‰ en 24 h)", "Cambio del δ13C del CO2 en 24 h")

# Dispersión: CO2 (mmol) vs δ13C-CO2, un punto por frasco
p_disp <- datos_co2 %>%
  filter(!is.na(mmol_CO2), !is.na(d_CO2)) %>%
  ggplot(aes(mmol_CO2, d_CO2, color = trt, shape = stage)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ site, labeller = labeller(site = etiq_site)) +
  scale_color_manual(values = col_trt, name = "Treatment") +
  scale_shape_manual(values = c(initial = 1, final = 16), labels = etiq_stage, name = "Stage") +
  labs(x = expression(CO[2]*" (mmol)"), y = expression(delta^{13}*"C-CO"[2]*" (\u2030)"),
       title = "CO2 vs δ13C-CO2 (un punto por frasco)") +
  tema

# -----------------------------------------------------------------------------------------------------
#  4. TABLAS: resumen por grupo y Dry vs Wet en el delta
# -----------------------------------------------------------------------------------------------------
# 4a. Por variable, tratamiento, sitio y tiempo: media ± SD de initial y final, delta, SD combinada y "distinguible"
est <- datos_l %>%
  filter(!is.na(v)) %>%
  group_by(var, trt, site, tiempo, stage) %>%
  summarise(m = mean(v), sd = sd(v), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = stage, values_from = c(m, sd, n)) %>%
  mutate(dif = m_final - m_initial,
         dif_sd = sqrt(sd_final^2 + sd_initial^2),
         distinguible = abs(dif) > 2 * dif_sd) %>%
  filter(!is.na(dif)) %>%
  arrange(var, trt, site, tiempo)

cat("\n=== Resumen por grupo: initial, final y delta de 24 h (distinguible = |delta| > 2 x SD combinada) ===\n")
print(est %>% mutate(across(c(m_initial, m_final, sd_initial, sd_final, dif, dif_sd), ~ signif(.x, 3))), n = Inf, width = Inf)

# 4b. Dry vs Wet en el delta (t de Welch por variable, sitio y tiempo; p ajustada con Holm dentro de cada variable)
dw <- datos_delta %>%
  filter(!is.na(delta)) %>%
  group_by(var, site, tiempo) %>%
  group_modify(~ {
    a <- .x$delta[.x$trt == "Dry"]; b <- .x$delta[.x$trt == "Wet"]
    if (length(a) < 2 || length(b) < 2) return(tibble(media_Dry = NA_real_, media_Wet = NA_real_, p = NA_real_))
    tibble(media_Dry = mean(a), media_Wet = mean(b),
           p = tryCatch(t.test(a, b)$p.value, error = function(e) NA_real_))
  }) %>%
  ungroup() %>%
  group_by(var) %>% mutate(p_Holm = p.adjust(p, method = "holm")) %>% ungroup()

cat("\n=== Dry vs Wet en el delta de 24 h (t de Welch; p_Holm = p ajustada por comparaciones multiples) ===\n")
print(dw %>% mutate(across(c(media_Dry, media_Wet), ~ signif(.x, 3)), p = round(p, 3), p_Holm = round(p_Holm, 3)), n = Inf)

# -----------------------------------------------------------------------------------------------------
#  5. ESTADISTICA (por variable, tratamiento y etapa)
# -----------------------------------------------------------------------------------------------------
cat("\n================ Shapiro-Wilk por grupo (n >= 3) ================\n")
datos_l %>% filter(!is.na(v)) %>%
  group_by(var, trt, site, stage, tiempo) %>%
  summarise(n = n(), p_shapiro = if (n() >= 3 && length(unique(v)) > 1) shapiro.test(v)$p.value else NA_real_,
            .groups = "drop") %>%
  mutate(normal = if_else(p_shapiro > 0.05, "SÍ", "NO")) %>% print(n = Inf)

stats_var <- function(d_all, etiqueta) {
  for (var_sel in levels(d_all$var)) for (tr in c("Dry", "Wet")) {
    d <- d_all %>% filter(var == var_sel, trt == tr, !is.na(v)) %>% droplevels()
    cat("\n>>>", var_sel, "|", etiqueta, "|", tr, "<<<\n")
    if (nrow(d) == 0 || nlevels(d$TIME) < 2) { cat("  (sin datos o menos de 2 tiempos)\n"); next }
    lev <- tryCatch(leveneTest(v ~ interaction(site, TIME), data = d)$`Pr(>F)`[1], error = function(e) NA)
    cat("  Levene p =", signif(lev, 3), "\n")
    mod <- tryCatch(aov(v ~ site * TIME, data = d), error = function(e) NULL)
    if (is.null(mod)) { cat("  (no se pudo ajustar el ANOVA)\n"); next }
    print(summary(mod))
    tk <- tryCatch(emmeans(mod, pairwise ~ TIME, adjust = "tukey")$contrasts, error = function(e) NULL)
    if (!is.null(tk)) { cat("  Tukey entre tiempos:\n"); print(tk) }
  }
}
stats_var(datos_l %>% filter(stage == "initial"), "initial (largo plazo)")
stats_var(datos_l %>% filter(stage == "final"),   "final")
stats_var(datos_delta %>% mutate(v = delta),      "delta 24 h (final - initial)")

# -----------------------------------------------------------------------------------------------------
#  6. KEELING: δ13C-CO2 frente a 1/CO2. Intercepto = δ13C del CO2 que se añadió (fuente)
# -----------------------------------------------------------------------------------------------------
keel <- datos_co2 %>%
  filter(!is.na(.data[[x_keeling]]), .data[[x_keeling]] > 0, !is.na(d_CO2)) %>%
  mutate(inv = 1 / .data[[x_keeling]])

keeling_fit <- function(d) {
  if (nrow(d) < 4 || length(unique(d$inv)) < 3)
    return(tibble(n = nrow(d), intercepto = NA_real_, ee = NA_real_, ic95_min = NA_real_, ic95_max = NA_real_,
                  R2 = NA_real_, pendiente = NA_real_))
  m  <- lm(d_CO2 ~ inv, data = d)
  s  <- summary(m)$coefficients
  ci <- confint(m)["(Intercept)", ]
  tibble(n = nrow(d), intercepto = s[1, 1], ee = s[1, 2], ic95_min = ci[1], ic95_max = ci[2],
         R2 = summary(m)$r.squared, pendiente = s[2, 1])
}

k_global <- keel %>% group_by(trt, site) %>% group_modify(~ keeling_fit(.x)) %>% ungroup()
k_tiempo <- keel %>% group_by(trt, site, tiempo) %>% group_modify(~ keeling_fit(.x)) %>% ungroup()

cat("\n=== Keeling por tratamiento y sitio (todos los tiempos; initial + final): intercepto = δ13C del CO2 añadido (‰) ===\n")
print(k_global %>% mutate(across(c(intercepto, ee, ic95_min, ic95_max, R2, pendiente), ~ signif(.x, 3))), n = Inf, width = Inf)
cat("\n=== Keeling por tratamiento, sitio y tiempo ===\n")
print(k_tiempo %>% mutate(across(c(intercepto, ee, ic95_min, ic95_max, R2, pendiente), ~ signif(.x, 3))), n = Inf, width = Inf)

p_keel <- ggplot(keel, aes(inv, d_CO2, color = site)) +
  geom_point(aes(shape = stage), size = 2.2, alpha = 0.8) +
  geom_smooth(aes(group = site), method = "lm", se = FALSE, fullrange = TRUE, linewidth = 0.6) +
  facet_wrap(~ trt, nrow = 1) +
  expand_limits(x = 0) +
  scale_color_manual(values = col_site, labels = etiq_site, name = "Sitio") +
  scale_shape_manual(values = c(initial = 1, final = 16), labels = etiq_stage, name = "Stage") +
  labs(x = paste0("1 / ", x_keeling), y = "δ13C-CO2 (‰)",
       title = "Gráfico de Keeling: δ13C del CO2 frente a 1/CO2",
       subtitle = "El punto donde la recta corta el eje y (1/CO2 = 0) es el δ13C del CO2 que se añadió en las 24 h",
       caption = "Initial (círculo abierto) y final (círculo lleno) juntos; una recta por sitio") +
  tema

p_keel_t <- ggplot(keel, aes(inv, d_CO2, color = site)) +
  geom_point(aes(shape = stage), size = 1.8, alpha = 0.8) +
  geom_smooth(aes(group = site), method = "lm", se = FALSE, fullrange = TRUE, linewidth = 0.5) +
  facet_grid(trt ~ tiempo, scales = "free_x", labeller = labeller(tiempo = function(x) paste0("t", x))) +
  expand_limits(x = 0) +
  scale_color_manual(values = col_site, labels = etiq_site, name = "Sitio") +
  scale_shape_manual(values = c(initial = 1, final = 16), labels = etiq_stage, name = "Stage") +
  labs(x = paste0("1 / ", x_keeling), y = "δ13C-CO2 (‰)",
       title = "Keeling por tiempo",
       caption = "Con pocos puntos por panel la recta es muy incierta: mira la tabla (R2 e intervalo)") +
  tema

# -----------------------------------------------------------------------------------------------------
#  7. MOSTRAR EN PANTALLA Y GUARDAR
# -----------------------------------------------------------------------------------------------------
graphics.off()
for (p in list(p_ppm_largo, p_ppm_delta, p_mmol_largo, p_mmol_delta, p_d13_largo, p_d13_delta, p_disp, p_keel, p_keel_t)) {
  print(p)
  Sys.sleep(0.2)
}

ggsave(file.path(carpeta_out, "CO2_fig1_ppm_largo.png"),   p_ppm_largo,  width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig2_ppm_delta.png"),   p_ppm_delta,  width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig3_mmol_largo.png"),  p_mmol_largo, width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig4_mmol_delta.png"),  p_mmol_delta, width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig5_d13C_largo.png"),  p_d13_largo,  width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig6_d13C_delta.png"),  p_d13_delta,  width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig7_dispersion.png"),  p_disp,       width = 11,  height = 6,   dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig8_Keeling.png"),     p_keel,       width = 9,   height = 5,   dpi = 300)
ggsave(file.path(carpeta_out, "CO2_fig9_Keeling_tiempo.png"), p_keel_t,  width = 13,  height = 6,   dpi = 300)

write_csv(datos_co2,   file.path(carpeta_out, "CO2_isotopos_por_frasco.csv"))
write_csv(datos_delta, file.path(carpeta_out, "CO2_delta_por_frasco.csv"))
write_csv(est,         file.path(carpeta_out, "CO2_tabla_resumen_delta.csv"))
write_csv(dw,          file.path(carpeta_out, "CO2_DryWet_delta_Welch.csv"))
write_csv(k_global,    file.path(carpeta_out, "CO2_Keeling_global.csv"))
write_csv(k_tiempo,    file.path(carpeta_out, "CO2_Keeling_por_tiempo.csv"))
