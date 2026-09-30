# =========================================================
# 04_graficos.R
# =========================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)

col_D <- "#B08968"
col_W <- "#4A7C82"

library(ggplot2)
library(dplyr)
library(tidyr)

# Paleta consistente
col_D <- "#B08968"  # Dry
col_W <- "#4A7C82"  # Wet

# 1. Cargar y procesar datos
ph_raw <- read.csv("ph_data.csv", sep=",", dec=".", fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)

names(ph_raw)[names(ph_raw) == "site..A.B."] <- "site"
names(ph_raw)[names(ph_raw) == "treatment..D.W."] <- "treatment"
names(ph_raw)[names(ph_raw) == "tiempo..t0.etc."] <- "time"
names(ph_raw)[names(ph_raw) == "sacrificio.I.F"] <- "sacrificio"

ph_long <- ph_raw %>%
  select(site, treatment, time, sacrificio, Replicate, pH_OL, pH_PW) %>%
  pivot_longer(cols = c(pH_OL, pH_PW), names_to = "layer", values_to = "pH") %>%
  mutate(layer = ifelse(layer == "pH_OL", "OL", "PW"))

####  OPCIONAL SI QUIERO SACAR EL DELTA FINAL VS INICIAL
ph_delta <- ph_long %>%
  group_by(site, treatment, time, layer, Replicate) %>%
  summarise(
    pH_I = pH[sacrificio == "I"][1],
    pH_F = pH[sacrificio == "F"][1],
    .groups = "drop"
  ) %>%
  mutate(delta_pH = pH_F - pH_I) %>%
  group_by(site, treatment, time, layer) %>%
  summarise(
    mean_delta = mean(delta_pH, na.rm=TRUE),
    sd_delta = sd(delta_pH, na.rm=TRUE), 
    .groups = "drop"
  ) %>%
  filter(time == "t0") %>%
  mutate(
    site = factor(site, levels = c("A", "B"), labels = c("site low OM", "site high OM"))
  )

# 2. Gráfico con formato homogéneo
p_ph <- ggplot(ph_delta, aes(x = site, y = mean_delta, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_delta - sd_delta, ymax = mean_delta + sd_delta),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6
  ) +
  facet_wrap(~layer) +
  scale_fill_manual(values = c(D = col_D, W = col_W), name = "Treatment") +
  geom_hline(yintercept = 0, linewidth = 0.5, color = "black") +
  labs(
    x = "Site", 
    y = expression(Delta*"pH (Final - Initial)"), 
    title = "pH change during 24h incubation (t0)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # Título centrado idéntico
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13, color = "black"),
    strip.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 14),
    legend.position = "top",                                          # Leyenda arriba
    panel.grid.minor = element_blank()
  )

# 3. Guardar la imagen
ggsave("poster_panels_ph.png", p_ph, width = 11, height = 8, dpi = 300)

cat("\nListo. Imagen guardada como poster_panels_ph.png\n")

# =========================================================
# Verificación del Supuesto de Normalidad para el pH
# =========================================================

library(dplyr)
library(tidyr)

# 1. Cargar y preparar datos de pH en formato largo
ph_raw <- read.csv("ph_data.csv", sep=",", dec=".", fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)

names(ph_raw)[names(ph_raw) == "site..A.B."] <- "site"
names(ph_raw)[names(ph_raw) == "treatment..D.W."] <- "treatment"
names(ph_raw)[names(ph_raw) == "tiempo..t0.etc."] <- "time"
names(ph_raw)[names(ph_raw) == "sacrificio.I.F"] <- "sacrificio"

ph_long_stats <- ph_raw %>%
  select(site, treatment, time, sacrificio, Replicate, pH_OL, pH_PW) %>%
  pivot_longer(cols = c(pH_OL, pH_PW), names_to = "layer", values_to = "pH") %>%
  mutate(layer = ifelse(layer == "pH_OL", "OL", "PW")) %>%
  filter(!is.na(pH))

# 2. Ajustar el modelo ANOVA sobre la variable pH
modelo_ph_abs <- aov(pH ~ site * treatment * layer * time * sacrificio, data = ph_long_stats)

# 3. Extraer residuos y prueba de Shapiro-Wilk
residuos_ph_abs <- residuals(modelo_ph_abs)
shapiro_ph_abs <- shapiro.test(residuos_ph_abs)

cat("\n--- RESULTADO DEL TEST DE SHAPIRO-WILK (pH) ---\n")
print(shapiro_ph_abs)

if (shapiro_ph_abs$p.value > 0.05) {
  cat("\nConclusión: Se CUMPLE el supuesto de normalidad (p > 0.05).\n")
} else {
  cat("\nConclusión: NO se cumple el supuesto de normalidad (p < 0.05). Conviene usar pruebas no paramétricas.\n")
}

# 4. Diagnóstico visual
par(mfrow = c(1, 2))

qqnorm(residuos_ph_abs, main = "Q-Q Plot de Residuos (pH)", pch = 19, col = "darkblue")
qqline(residuos_ph_abs, col = "red", lwd = 2)

hist(residuos_ph_abs, 
     main = "Histograma de Residuos", 
     xlab = "Residuos de pH", 
     col = "steelblue", 
     border = "white",
     probability = TRUE)
lines(density(residuos_ph_abs), col = "red", lwd = 2)

par(mfrow = c(1, 1))

# ==============================================================================
# SCRIPT DE ANÁLISIS NO PARAMÉTRICO PARA pH
# Pruebas: Wilcoxon (Mann-Whitney) y Scheirer-Ray-Hare
# ==============================================================================

library(dplyr)
if(!require(rcompanion)) install.packages("rcompanion")
library(rcompanion)

# ------------------------------------------------------------------------------
# 1. COMPARACIÓN ENTRE CAPAS (OL vs PW)
# ------------------------------------------------------------------------------

cat("=== 1A. Diferencia Global entre OL y PW ===\n")
test_capas_global_ph <- wilcox.test(pH ~ layer, data = ph_long_stats)
print(test_capas_global_ph)

cat("\n=== 1B. Diferencia entre OL y PW dentro de cada Sitio ===\n")
capas_por_sitio_ph <- ph_long_stats %>%
  group_by(site) %>%
  summarise(
    W_stat  = wilcox.test(pH ~ layer)$statistic,
    p_value = wilcox.test(pH ~ layer)$p.value,
    .groups = "drop"
  )
print(capas_por_sitio_ph)


# ------------------------------------------------------------------------------
# 2. COMPARACIÓN ENTRE SITIOS (Site A vs Site B)
# ------------------------------------------------------------------------------

cat("\n=== 2A. Diferencia Global entre Sitios ===\n")
test_sitios_global_ph <- wilcox.test(pH ~ site, data = ph_long_stats)
print(test_sitios_global_ph)

cat("\n=== 2B. Diferencia entre Sitios desglosada por Capa (OL y PW) ===\n")
sitios_por_capa_ph <- ph_long_stats %>%
  group_by(layer) %>%
  summarise(
    W_stat  = wilcox.test(pH ~ site)$statistic,
    p_value = wilcox.test(pH ~ site)$p.value,
    .groups = "drop"
  )
print(sitios_por_capa_ph)


# ------------------------------------------------------------------------------
# 3. COMPARACIÓN ENTRE TRATAMIENTOS (Dry vs Wet)
# ------------------------------------------------------------------------------

cat("\n=== 3A. Diferencia Global entre Tratamientos (D vs W) ===\n")
test_tratamiento_global_ph <- wilcox.test(pH ~ treatment, data = ph_long_stats)
print(test_tratamiento_global_ph)

cat("\n=== 3B. Diferencia entre Tratamientos por Capa ===\n")
tratamiento_por_capa_ph <- ph_long_stats %>%
  group_by(layer) %>%
  summarise(
    W_stat  = wilcox.test(pH ~ treatment)$statistic,
    p_value = wilcox.test(pH ~ treatment)$p.value,
    .groups = "drop"
  )
print(tratamiento_por_capa_ph)


# ------------------------------------------------------------------------------
# 4. PRUEBA GLOBAL MULTIFACTORIAL (Scheirer-Ray-Hare)
# ------------------------------------------------------------------------------

cat("\n=== 4. ANOVA No Paramétrico (Scheirer-Ray-Hare) ===\n")
scheirer_ph <- scheirerRayHare(pH ~ site * treatment * layer * time, data = ph_long_stats)
print(scheirer_ph)

# Crear tabla de resumen de resultados para pH
tabla_resumen_ph <- data.frame(
  Comparacion = c("Capa (Global: OL vs PW)", "Capa en Sitio A", "Capa en Sitio B",
                  "Sitio (Global: A vs B)", "Sitio A vs B (en OL)", "Sitio A vs B (en PW)",
                  "Tratamiento (Global: D vs W)", "Tratamiento D vs W (en PW)", "Tratamiento D vs W (en OL)"),
  Prueba = c("Wilcoxon", "Wilcoxon", "Wilcoxon", 
             "Wilcoxon / Scheirer-Ray-Hare", "Wilcoxon", "Wilcoxon", 
             "Wilcoxon / Scheirer-Ray-Hare", "Wilcoxon", "Wilcoxon"),
  Estadistico = c("W = 5183", "W = 1215", "W = 1359", 
                  "W = 3690 / H = 2.12", "W = 880.5", "W = 916.5", 
                  "W = 2397 / H = 1.02", "W = 826", "W = 390"),
  p_value = c("< 0.0001", "< 0.0001", "< 0.0001", "0.166", "0.447", "0.457", "0.365", "0.005", "0.146"),
  Significativo = c("SÍ (p < 0.001)", "SÍ (p < 0.001)", "SÍ (p < 0.001)", 
                   "No", "No", "No", 
                   "No", "SÍ (p < 0.01)", "No")
)

# Guardar en archivo CSV
write.csv(tabla_resumen_ph, "Resumen_Resultados_pH.csv", row.names = FALSE)



# =========================================================
# 2. DIC
# =========================================================
library(ggplot2)
library(dplyr)

# Paleta consistente
col_D <- "#B08968"  # Dry
col_W <- "#4A7C82"  # Wet

# 1. Cargar datos
dic <- read.csv("DIC.csv", sep = ",", dec = ".", fileEncoding = "UTF-8-BOM")
names(dic) <- c("SAMPLE", "TIME", "LAYER", "SITE", "TREATMENT", "DIC")

# 2. Resumen y ajuste de factores/etiquetas
tiempos_unicos_dic <- sort(unique(dic$TIME))[1:2]

dic_summary <- dic %>%
  filter(TIME %in% tiempos_unicos_dic) %>%
  group_by(TIME, SITE, TREATMENT, LAYER) %>%
  summarise(
    mean_DIC = mean(DIC, na.rm = TRUE),
    sd_DIC = sd(DIC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    SITE = factor(SITE, levels = c("A", "B"), labels = c("site low OM", "site high OM")),
    TIME = factor(TIME, levels = tiempos_unicos_dic, labels = c("Day 0 (T0)", "Day 14 (T1)")),
    LAYER = factor(LAYER)
  )

# 3. Gráfico 
p_dic <- ggplot(dic_summary, aes(x = SITE, y = mean_DIC, fill = TREATMENT)) +
  geom_col(position = position_dodge2(preserve = "single", width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_DIC - sd_DIC, ymax = mean_DIC + sd_DIC),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.6
  ) +
  facet_grid(LAYER ~ TIME) +
  scale_fill_manual(values = c(D = col_D, W = col_W), name = "Treatment") +
  labs(
    title = "DIC concentration in each Site, Treatment and Time",
    x = "Site",
    y = "DIC (mM)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5), # Título centrado igual al anterior
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 13, color = "black"),
    strip.text = element_text(face = "bold", size = 15),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 14),
    legend.position = "top",                                        # Leyenda en la parte superior
    panel.grid.minor = element_blank()
  )

# 4. Guardar la imagen
ggsave("poster_panels_dic.png", p_dic, width = 11, height = 8, dpi = 300)

cat("\nListo. Imagen guardada como poster_panels_dic.png\n")

# ==============================================================================
# SCRIPT DE ANÁLISIS NO PARAMÉTRICO PARA DIC (Tratamiento WET)
# Pruebas: Wilcoxon (Mann-Whitney) y Scheirer-Ray-Hare
# ==============================================================================

#En caso de que los datos no tengan normalidad.
# Cargar librerías necesarias
library(dplyr)

if(!require(rcompanion)) install.packages("rcompanion")
library(rcompanion)

# ------------------------------------------------------------------------------
# 1. COMPARACIÓN ENTRE CAPAS (OL vs PW)
# ------------------------------------------------------------------------------

cat("=== 1A. Diferencia Global entre OL y PW ===\n")
test_capas_global <- wilcox.test(DIC ~ LAYER, data = dic_wet_stats)
print(test_capas_global)

cat("\n=== 1B. Diferencia entre OL y PW dentro de cada Sitio ===\n")
capas_por_sitio <- dic_wet_stats %>%
  group_by(SITE) %>%
  summarise(
    W_stat  = wilcox.test(DIC ~ LAYER)$statistic,
    p_value = wilcox.test(DIC ~ LAYER)$p.value,
    .groups = "drop"
  )
print(capas_por_sitio)

cat("\n=== 1C. Diferencia entre OL y PW por Sitio y Tiempo ===\n")
capas_sitio_tiempo <- dic_wet_stats %>%
  group_by(SITE, TIME) %>%
  summarise(
    W_stat  = wilcox.test(DIC ~ LAYER)$statistic,
    p_value = wilcox.test(DIC ~ LAYER)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    significativo = ifelse(p_adj < 0.05, "SÍ", "No")
  )
print(capas_sitio_tiempo)


# ------------------------------------------------------------------------------
# 2. COMPARACIÓN ENTRE SITIOS (Site A vs Site B), Simil a ANOVA.
# ------------------------------------------------------------------------------

cat("\n=== 2A. Diferencia Global entre Sitios ===\n")
test_sitios_global <- wilcox.test(DIC ~ SITE, data = dic_wet_stats)
print(test_sitios_global)

cat("\n=== 2B. Diferencia entre Sitios desglosada por Capa (OL y PW) ===\n")
sitios_por_capa <- dic_wet_stats %>%
  group_by(LAYER) %>%
  summarise(
    W_stat  = wilcox.test(DIC ~ SITE)$statistic,
    p_value = wilcox.test(DIC ~ SITE)$p.value,
    .groups = "drop"
  )
print(sitios_por_capa)


# ------------------------------------------------------------------------------
# 3. PRUEBA GLOBAL MULTIFACTORIAL (Scheirer-Ray-Hare) SIMIL A UN ANOVAMULTIFAC
# ------------------------------------------------------------------------------

cat("\n=== 3. ANOVA No Paramétrico (Scheirer-Ray-Hare) ===\n")
scheirer_global <- scheirerRayHare(DIC ~ LAYER * SITE * TIME, data = dic_wet_stats)
print(scheirer_global)



# ===# =========================================================
# 3. CaCO3 y OM (Con filtro temporal para el dato raro de OM) ESTE FILTRO DEBO SACARLO LUEGO!!!!
# =========================================================
sed <- read.csv("sedimento_clean.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(pct_OM) & pct_OM > 0) %>% 
  filter(pct_OM < 15) # <-- FILTRO TEMPORAL: Ajusta este número según tu dato raro

p_caco3 <- ggplot(sed, aes(x=site, y=pct_CaCO3_corrected, fill=site)) +
  geom_boxplot(width=0.5, alpha=0.7) +
  geom_jitter(width=0.08, alpha=0.5) +
  stat_compare_means(method="wilcox.test", label="p.format") +
  labs(x="Site", y="%CaCO3", title="CaCO3 content") +
  theme_minimal(base_size=13) + theme(legend.position="none")

p_om <- ggplot(sed, aes(x=site, y=pct_OM, fill=site)) +
  geom_boxplot(width=0.5, alpha=0.7) +
  geom_jitter(width=0.08, alpha=0.5) +
  stat_compare_means(method="wilcox.test", label="p.format") +
  labs(x="Site", y="%OM", title="Organic matter content") +
  theme_minimal(base_size=13) + theme(legend.position="none")


# Combinar y exportar
# =========================================================
final_plot <- (p_ph) 
ggsave("poster_panels_ph.png", final_plot, width=12, height=9, dpi=300)

cat("\nListo. Imagen guardada como poster_panels.png\n")

final_plot <- (p_dic) 
ggsave("poster_panels_dic.png", final_plot, width=12, height=9, dpi=300)

cat("\nListo. Imagen guardada como poster_panels.png\n")

final_plot <- (p_caco3) 
ggsave("poster_panels_caco3.png", final_plot, width=12, height=9, dpi=300)

final_plot <- (p_om) 
ggsave("poster_panels_om.png", final_plot, width=12, height=9, dpi=300)

