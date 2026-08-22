# =========================================================
# 01_instalar_paquetes.R
# Ejecutar UNA SOLA VEZ si no los tienes instalados.
# =========================================================
# install.packages(c("ggplot2", "dplyr", "tidyr", "patchwork", "ggpubr", "car", "stringr", "broom"))

# =========================================================
# 02_procesar_sedimento.R
# =========================================================
library(dplyr)
library(stringr)

# --- 1. Leer los datos crudos ---
sed <- read.csv("SED.csv", stringsAsFactors = FALSE)

# --- 2. Extraer metadatos desde lab_code ---
sed <- sed %>%
  mutate(
    time       = str_extract(lab_code, "^t[0-9]+"),
    sacrificio = str_extract(lab_code, "(?<=t[0-9])[IF]"),
    site       = str_extract(lab_code, "(?<=-)[AB]"),
    treatment  = str_extract(lab_code, "(?<=[AB])[DW]"),
    replicate  = as.numeric(str_extract(lab_code, "[0-9]+$"))
  )

# --- 3. Calcular %OM y %CaCO3 ---
FACTOR_CACO3 <- 1 / 0.4397   # = 2.274

sed <- sed %>%
  mutate(
    pct_OM              = (pre450_g - post450_g) / pre450_g * 100,
    post900_neto        = post900_g - crisol_g,
    weight_loss_900      = post450_g - post900_neto,
    pct_CaCO3_raw       = weight_loss_900 / post450_g * 100,
    pct_CaCO3_corrected = pct_CaCO3_raw * FACTOR_CACO3
  )

# --- 4. Revision de valores atipicos / imposibles ---
sospechosos <- sed %>% filter(pct_CaCO3_raw < 0 | pct_CaCO3_raw > 100)
if (nrow(sospechosos) > 0) {
  cat("ATENCION -- revisa estas filas, tienen valores imposibles:\n")
  print(sospechosos$lab_code)
}

# --- 5. Guardar el archivo limpio ---
write.csv(sed, "sedimento_clean.csv", row.names = FALSE)
cat("\nListo. Archivo guardado en sedimento_clean.csv\n")
cat("Filas procesadas:", nrow(sed), "\n")


# =========================================================
# 03_analisis_y_graficos.R
# =========================================================
library(car)
library(ggplot2)
library(broom)

# --- 1. Filtrar ceros y NAs para el analisis estadistico ---
sed_clean <- read.csv("sedimento_clean.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(pct_OM) & pct_OM > 0 & pct_OM < 100) %>%
  mutate(
    site = as.factor(site),
    treatment = as.factor(treatment)
  )

sed_caco3 <- sed_clean %>% filter(!is.na(pct_CaCO3_corrected) & pct_CaCO3_corrected > 0)

# --- 2. Diagnostico de Normalidad y Varianzas ---
modelo_OM <- aov(pct_OM ~ site * treatment, data = sed_clean)
modelo_CaCO3 <- aov(pct_CaCO3_corrected ~ site * treatment, data = sed_caco3)

cat("\n===== PRUEBAS DE NORMALIDAD (Shapiro-Wilk) =====\n")
print(shapiro.test(residuals(modelo_OM)))
print(shapiro.test(residuals(modelo_CaCO3)))

# --- 3. Pruebas No Parametricas (Wilcoxon / Kruskal-Wallis) ---
cat("\n===== PRUEBAS DE WILCOXON POR SITIO =====\n")
print(wilcox.test(pct_OM ~ site, data = sed_clean))
print(wilcox.test(pct_CaCO3_corrected ~ site, data = sed_caco3))

cat("\n===== PRUEBAS DE WILCOXON POR TRATAMIENTO =====\n")
print(wilcox.test(pct_OM ~ treatment, data = sed_clean))
print(wilcox.test(pct_CaCO3_corrected ~ treatment, data = sed_caco3))

# --- 4. Exportar Tablas ANOVA a CSV ---
write.csv(tidy(modelo_OM), "resultados_ANOVA_OM.csv", row.names = FALSE)
write.csv(tidy(modelo_CaCO3), "resultados_ANOVA_CaCO3.csv", row.names = FALSE)

# --- 5. Graficos Boxplot ---
g1 <- ggplot(sed_clean, aes(x = site, y = pct_OM, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Materia Organica (%OM) por Sitio y Tratamiento",
    x = "Sitio",
    y = "% Materia Organica",
    fill = "Tratamiento"
  )

g2 <- ggplot(sed_caco3, aes(x = site, y = pct_CaCO3_corrected, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Carbonato de Calcio (%CaCO3) por Sitio y Tratamiento",
    x = "Sitio",
    y = "% CaCO3 (Corregido)",
    fill = "Tratamiento"
  )

ggsave("grafico_OM.png", plot = g1, width = 6, height = 4, dpi = 300)
ggsave("grafico_CaCO3.png", plot = g2, width = 6, height = 4, dpi = 300)

cat("\n Script completado con exito. Graficos y tablas guardadas.\n")



# =========================================================
#  ESTA PARTE ES PARA CUANDO QUIERA VER LA ESTRUCTURA PARA ANÁLISIS MULTIFACTORIAL, OSEA CUANDO QUIERA VER LOS DISTINTOS TIEMPOS (Sitio x Tratamiento x Tiempo)
# =========================================================
library(dplyr)
library(car)
library(ggplot2)
library(broom)

# --- 1. Cargar y asegurar que 'time' sea factor ---
sed_multi <- read.csv("sedimento_clean.csv", stringsAsFactors = FALSE) %>%
  filter(!is.na(pct_OM) & pct_OM > 0 & pct_OM < 100) %>%
  mutate(
    site      = as.factor(site),
    treatment = as.factor(treatment),
    time      = as.factor(time) # Se define tiempo como factor categórico
  )

sed_multi_caco3 <- sed_multi %>% 
  filter(!is.na(pct_CaCO3_corrected) & pct_CaCO3_corrected > 0)

# --- 2. ANOVA de 3 Vías (Sitio x Tratamiento x Tiempo) ---
mod_3way_OM    <- aov(pct_OM ~ site * treatment * time, data = sed_multi)
mod_3way_CaCO3 <- aov(pct_CaCO3_corrected ~ site * treatment * time, data = sed_multi_caco3)

cat("\n===== ANOVA 3 VÍAS: %OM =====\n")
summary(mod_3way_OM)

cat("\n===== ANOVA 3 VÍAS: %CaCO3 =====\n")
summary(mod_3way_CaCO3)

# --- 3. Verificación de Normalidad de Residuos ---
cat("\n===== PRUEBAS DE NORMALIDAD (MODELO 3 VÍAS) =====\n")
print(shapiro.test(residuals(mod_3way_OM)))
print(shapiro.test(residuals(mod_3way_CaCO3)))

# --- 4. Alternativa No Paramétrica si la Normalidad Falla ---
# Evaluar el efecto del tiempo de forma independiente por variable
cat("\n===== KRUSKAL-WALLIS PARA FACTOR TIEMPO =====\n")
print(kruskal.test(pct_OM ~ time, data = sed_multi))
print(kruskal.test(pct_CaCO3_corrected ~ time, data = sed_multi_caco3))

# --- 5. Gráfico de Tendencia Temporal por Sitio y Tratamiento ---
# Muestra la evolución a lo largo del tiempo separada por sitio
g_tiempo_OM <- ggplot(sed_multi, aes(x = time, y = pct_OM, fill = treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(0.2), alpha = 0.4) +
  facet_wrap(~ site, labeller = label_both) + # Separa el gráfico en paneles (Sitio A y Sitio B)
  theme_classic() +
  labs(
    title = "Evolución Temporal de Materia Orgánica (%OM)",
    x = "Tiempo",
    y = "% OM",
    fill = "Tratamiento"
  )

ggsave("grafico_tiempo_OM.png", plot = g_tiempo_OM, width = 8, height = 4, dpi = 300)