rm(list = ls())

library(ggplot2)
library(dplyr)
library(lubridate)
library(tidyr)

setwd("C:/Experimentcalcite/data")

# 1. Leer archivo
od_raw <- read.csv("SENSORS_OD.csv", stringsAsFactors = FALSE)
names(od_raw) <- c("sensor", "date", "treatment", "site", "replicate", "OD")

# 2. Conversión limpia de tipos
od <- od_raw %>%
  mutate(
    OD = as.numeric(gsub(",", ".", as.character(OD))),
    date = parse_date_time(date, orders = c("dmy", "ymd", "mdy", "dmy HMS", "ymd HMS")),
    date = as.Date(date)
  )

# 3. Filtrar Wet
od_wet <- od %>% 
  filter(treatment == "W") %>% 
  drop_na(OD, date, site)

col_A <- "#C97B4A"; col_B <- "#3E7C8C"

# 4. Boxplot con stat_summary CORREGIDO (Retorna data.frame)
g_od_wet <- ggplot(od_wet, aes(x = factor(site), y = OD, fill = factor(site))) +
  geom_boxplot(width = 0.5, alpha = 0.75, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = c(A = col_A, B = col_B)) +
  stat_summary(
    fun.data = function(x) {
      data.frame(
        y = max(x, na.rm = TRUE) + (max(x, na.rm = TRUE) * 0.03), # Coordenada Y numérica pura
        label = paste0("n = ", length(x))                          # Etiqueta en texto
      )
    }, 
    geom = "text", size = 4, fontface = "bold"
  ) +
  labs(title = "OD — Wet only, A vs B", x = "Site", y = "OD (umol/L)") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", size = 12))

# 5. Evolución temporal
resumen_wet <- od_wet %>%
  group_by(date, site) %>%
  summarise(mean_OD = mean(OD, na.rm = TRUE), 
            sd_OD = sd(OD, na.rm = TRUE), 
            .groups = "drop")

g_od_time <- ggplot(resumen_wet, aes(x = date, y = mean_OD, color = factor(site))) +
  geom_line(linewidth = 0.8) + 
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_OD - sd_OD, ymax = mean_OD + sd_OD, fill = factor(site)), 
              alpha = 0.15, color = NA) +
  scale_color_manual(values = c(A = col_A, B = col_B)) +
  scale_fill_manual(values = c(A = col_A, B = col_B), guide = "none") +
  labs(title = "OD over time — Wet only", x = "Date", y = "OD (umol/L)") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 12))

# 6. Desplegar y Guardar
print(g_od_wet)
print(g_od_time)

ggsave("grafico_OD_wet_boxplot.png", g_od_wet, width = 5, height = 5, dpi = 300)
ggsave("grafico_OD_wet_tiempo.png", g_od_time, width = 8, height = 5, dpi = 300)



library(dplyr)

# 1. Test de Normalidad por Sitio (A y B por separado)
norm_A <- shapiro.test(od_wet$OD[od_wet$site == "A"])
norm_B <- shapiro.test(od_wet$OD[od_wet$site == "B"])

cat("=========================================\n")
cat("1. PRUEBAS DE NORMALIDAD (Shapiro-Wilk)\n")
cat("=========================================\n")
cat("Sitio A: p-value =", norm_A$p.value, "\n")
cat("Sitio B: p-value =", norm_B$p.value, "\n")

# 2. Decisión automática e informe de resultados
cat("\n=========================================\n")
cat("2. COMPARACIÓN ESTADÍSTICA (A vs B)\n")
cat("=========================================\n")

if (norm_A$p.value > 0.05 & norm_B$p.value > 0.05) {
  cat("-> Ambos sitios son normales. Usando Student's t-test:\n\n")
  test_res <- t.test(OD ~ site, data = od_wet)
  print(test_res)
} else {
  cat("-> Al menos un sitio NO es normal. Usando Wilcoxon rank sum test:\n\n")
  test_res <- wilcox.test(OD ~ site, data = od_wet)
  print(test_res)
}

