##

library(data.table)
theme_set(ggpubr::theme_classic2())
ribflv_data <- fread("../../ElentaFBA/ribflv_alone_biomass_vary.txt")
ribflv_plot <- ggplot(ribflv_data[seq(1, nrow(ribflv_data), by = 10)], aes(x = -1*Var1, y = growth_rate_results)) + geom_line() + geom_point(size = 0.5) + xlim(0, 0.002) + 
  xlab("Biomass coefficient for riboflavin") + ylab("FBA max growth rate") + 
  geom_vline(xintercept = 0.00005, linetype = 2, color = "blue")+ theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10))

fad_data <- fread("../../ElentaFBA/fad_biomass_vary.txt")
fad_plot <- ggplot(fad_data[seq(1, nrow(fad_data), by = 20)], aes(x = -1*Var1, y = growth_rate_results)) + geom_line() + geom_point(size = 0.5) + xlim(0, 0.005) + 
  xlab("Biomass coefficient for FAD") + ylab("FBA max growth rate") + 
  geom_vline(xintercept = 0.001, linetype = 2, color = "blue") + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10))
library(patchwork)
ribflv_plot + fad_plot

atp_data <- fread("../../ElentaFBA/atp_biomass_vary.txt")
atp_plot <- ggplot(atp_data, aes(x = -1*Var1, y = growth_rate_results)) + geom_line() + geom_point(size = 0.5) + xlim(35, 90) + 
  xlab("Biomass coefficient for ATP") + ylab("FBA max growth rate") + 
  geom_vline(xintercept = 40.1102, linetype = 2, color = "blue") + theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10))

combined_plot <- ribflv_plot + fad_plot + atp_plot
ggsave(combined_plot, file = "figures/biomass_coefficient_sensitivity.pdf", width = 8, height = 2.8)
ggsave(combined_plot, file = "figures/biomass_coefficient_sensitivity.png", width = 8, height = 2.8)

combined_dat <- rbind(data.table(ribflv_data[seq(1, nrow(ribflv_data), by = 10)], Compound = "Riboflavin", 
                                 ModelValue = 0.00005), 
                      data.table(fad_data[seq(1, nrow(fad_data), by = 20)], Compound = "FAD",
                                 ModelValue = 0.001),
                      data.table(atp_data, Compound = "ATP", ModelValue = 40.1102))
combined_plot2 <- ggplot(combined_dat, aes(x = -1*Var1, y = growth_rate_results)) + geom_line() + geom_point(size = 0.5) + 
  xlab("Biomass coefficient") + ylab("FBA max growth rate") + facet_wrap(~Compound, scales = "free") + 
  geom_vline(aes(xintercept = ModelValue, color = "Value in iEL2243_2"), linetype = 2) +
  scale_color_manual(values = "blue", name = "") +
  theme(axis.text = element_text(size = 9), axis.title = element_text(size = 10), 
        strip.background = element_blank(), legend.position = "bottom")
ggsave(combined_plot2, file = "figures/biomass_coefficient_sensitivity_v2.pdf", width = 7.5, height = 3.3)
ggsave(combined_plot2, file = "figures/biomass_coefficient_sensitivity_v2.png", width = 7.5, height = 3.3)

all_comps <- fread("../../ElentaFBA/all_biomass_components_vary.txt")
all_comps[growth_rate_results_2x < 0.958]
