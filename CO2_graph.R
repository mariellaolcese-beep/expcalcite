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
