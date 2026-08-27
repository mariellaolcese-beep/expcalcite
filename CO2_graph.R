library(ggplot2)
library(dplyr)

# Paleta consistente
col_D <- "#B08968"  # Dry
col_W <- "#4A7C82"  # Wet

# 1. Cargar datos (si tu CSV usa coma como separador y punto como decimal)
co2_data <- read.csv("co2_data.csv", fileEncoding = "UTF-8-BOM")

# NOTA: Si lo anterior sigue dando 1 sola columna, prueba con sep="," o sep="\t":
# co2_data <- read.csv("co2_data.csv", sep = ",", dec = ".", fileEncoding = "UTF-8-BOM")

# Verificar cuántas columnas se leyeron realmente antes de renombrar:
ncol(co2_data)  # Debería dar 6
head(co2_data)  # Te muestra las primeras filas

# 2. Asignar nombres limpios a las 6 columnas
names(co2_data) <- c("SAMPLE", "TIME", "SITE", "STAGE", "TREATMENT", "CO2_ppm")

# 3. Resumen y Mutate
co2_summary <- co2_data %>%
  group_by(TIME, SITE, TREATMENT, STAGE) %>%
  summarise(
    mean_CO2 = mean(CO2_ppm, na.rm = TRUE),
    sd_CO2 = sd(CO2_ppm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    TIME = factor(TIME, labels = c("Day 0 (T0)", "Day 14 (T1)", "Day 28 (T2)")),
    STAGE = factor(STAGE, levels = c("F", "I"), labels = c("Final (F)", "Initial (I)"))
  )

# 4. Gráfico
p_co2 <- ggplot(co2_summary, aes(x = SITE, y = mean_CO2, fill = TREATMENT)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean_CO2 - sd_CO2, ymax = mean_CO2 + sd_CO2),
    position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5
  ) +
  facet_grid(STAGE ~ TIME) +
  scale_fill_manual(values = c(D = col_D, W = col_W), name = "Treatment") +
  labs(
    title = "CO2 concentration in each Site, Treatment and Time",
    x = "Site",
    y = expression(CO[2]*" (ppm)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave("poster_panels_co2.png", p_co2, width = 12, height = 9, dpi = 300)

