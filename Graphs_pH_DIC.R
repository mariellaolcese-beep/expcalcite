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

# =========================================================
# 1. pH -- Delta I->F, para OL Y PW (por eso antes faltaba PW)
# =========================================================
ph_raw <- read.csv("ph_data.csv", sep=",", dec=".", fileEncoding="UTF-8-BOM", stringsAsFactors = FALSE)

# Ajusta estos nombres si en tu csv se llaman distinto:
names(ph_raw)[names(ph_raw) == "site..A.B."] <- "site"
names(ph_raw)[names(ph_raw) == "treatment..D.W."] <- "treatment"
names(ph_raw)[names(ph_raw) == "tiempo..t0.etc."] <- "time"
names(ph_raw)[names(ph_raw) == "sacrificio.I.F"] <- "sacrificio"

ph_long <- ph_raw %>%
  select(site, treatment, time, sacrificio, Replicate, pH_OL, pH_PW) %>%
  pivot_longer(cols = c(pH_OL, pH_PW), names_to = "layer", values_to = "pH") %>%
  mutate(layer = ifelse(layer == "pH_OL", "OL", "PW"))

ph_delta <- ph_long %>%
  group_by(site, treatment, time, layer, Replicate) %>%
  summarise(pH_I = pH[sacrificio == "I"][1],
            pH_F = pH[sacrificio == "F"][1],
            .groups = "drop") %>%
  mutate(delta_pH = pH_F - pH_I) %>%
  group_by(site, treatment, time, layer) %>%
  summarise(mean_delta = mean(delta_pH, na.rm=TRUE),
            sd_delta = sd(delta_pH, na.rm=TRUE), .groups="drop") %>%
  filter(time == "t0")   # cambia o quita este filtro si quieres t1 tambien

p_ph <- ggplot(ph_delta, aes(x=site, y=mean_delta, fill=treatment)) +
  geom_col(position=position_dodge(width=0.7), width=0.6) +
  geom_errorbar(aes(ymin=mean_delta-sd_delta, ymax=mean_delta+sd_delta),
                position=position_dodge(width=0.7), width=0.15) +
  facet_wrap(~layer) +
  scale_fill_manual(values=c(D=col_D, W=col_W), name="Treatment") +
  geom_hline(yintercept=0, linewidth=0.4) +
  labs(x="Site", y=expression(Delta*"pH (Final - Initial)"), title="pH change during 24h incubation (t0)") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", size=12))

# =========================================================
# 2. DIC
# =========================================================
dic <- read.csv("DIC.csv", sep=",", dec=".", fileEncoding="UTF-8-BOM")
names(dic) <- c("SAMPLE","TIME","LAYER","SITE","TREATMENT","DIC")

dic_summary <- dic %>%
  group_by(TIME, SITE, TREATMENT, LAYER) %>%
  summarise(mean_DIC = mean(DIC, na.rm=TRUE),
            sd_DIC = sd(DIC, na.rm=TRUE), .groups="drop") %>%
  mutate(TIME = factor(TIME, labels=c("Day 0","Day 14")))

p_dic <- ggplot(dic_summary, aes(x=SITE, y=mean_DIC, fill=TREATMENT)) +
  geom_col(position=position_dodge(width=0.7), width=0.6) +
  geom_errorbar(aes(ymin=mean_DIC-sd_DIC, ymax=mean_DIC+sd_DIC),
                position=position_dodge(width=0.7), width=0.2) +
  facet_grid(LAYER ~ TIME) +
  scale_fill_manual(values=c(D=col_D, W=col_W), name="Treatment") +
  labs(x="Site", y="DIC (mM)", title="DIC by site, treatment and time") +
  theme_minimal(base_size=13) +
  theme(plot.title=element_text(face="bold", size=12))

# =========================================================
# 3. CaCO3 y OM (desde el archivo ya procesado) ESTA PARTE ESTA HORRENDA 
# =========================================================
SED <- read.csv("sedimento_clean.csv", stringsAsFactors = FALSE)

p_caco3 <- ggplot(sed, aes(x=site, y=pct_CaCO3_corrected, fill=site)) +
  geom_boxplot(width=0.5, alpha=0.7) +
  geom_jitter(width=0.08, alpha=0.5) +
  stat_compare_means(method="t.test", label="p.format") +
  labs(x="Site", y="%CaCO3", title="CaCO3 content") +
  theme_minimal(base_size=13) + theme(legend.position="none")

p_om <- ggplot(sed, aes(x=site, y=pct_OM, fill=site)) +
  geom_boxplot(width=0.5, alpha=0.7) +
  geom_jitter(width=0.08, alpha=0.5) +
  stat_compare_means(method="t.test", label="p.format") +
  labs(x="Site", y="%OM", title="Organic matter content") +
  theme_minimal(base_size=13) + theme(legend.position="none")

# =========================================================
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
