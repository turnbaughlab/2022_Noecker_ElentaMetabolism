library(data.table)
library(readxl)
#library(ggplot2)
#library(cowplot)
#library(ggrepel)
#library(sva)
#library(UpSetR)
#library(ggpubr)
library(tidyverse)
library(lmerTest)
library(broom.mixed)
#theme_set(theme_pubr())

load("mouse_data_models.rda")
feature_list <- met_data_sub[Duplicate == 0,sort(unique(IDUniq))]
metadata2 <- data.table(read_xlsx("origData/ForCN_elen_mono.xlsx", skip = 2))
metadata2 <- metadata2[!is.na(Metabolomics)]
#metadata2[,table(Cage, Group)]
met_data_sub[,Animal2:=as.numeric(gsub(".*_", "", Animal))]
met_data_sub <- merge(met_data_sub, metadata2[,list(`Mouse ID`, Cage)], by.x = "Animal2", by.y = "Mouse ID", all.x = T)
mods_all = vector("list", length(feature_list))
effects_all_all =vector("list", length(feature_list))
contrasts_all <- vector("list", length(feature_list))
constrasts_all_strains <- vector("list", length(feature_list))
for(i in 1:length(feature_list)){
  x = feature_list[i]
  mod_data = met_data_sub[IDUniq==x]
  mod_results = lmer(log10value~SampleType+Group+SampleType:Group + (1|Cage/Animal), data = mod_data)
  mod_coefs = data.table(tidy(mod_results), Model = "Full")
  site_contrasts_full <- difflsmeans(mod_results) %>% as.data.frame() %>% rownames_to_column() %>% data.table()
  site_contrasts_full <- site_contrasts_full[rowname %in% c("SampleTypeSI:GroupGF - SampleTypeSI:Group2243",
                                                            "SampleTypeSI:GroupGF - SampleTypeSI:Group15644",
                                                            "SampleTypeSI:GroupGF - SampleTypeSI:GroupAB12n2",
                                                            "SampleTypeSI:Group2243 - SampleTypeSI:Group15644",
                                                            "SampleTypeSI:Group2243 - SampleTypeSI:GroupAB12n2",
                                                            "SampleTypeSI:Group15644 - SampleTypeSI:GroupAB12n2",
                                                            "SampleTypeCecal:GroupGF - SampleTypeCecal:Group2243",
                                                            "SampleTypeCecal:GroupGF - SampleTypeCecal:Group15644",
                                                            "SampleTypeCecal:GroupGF - SampleTypeCecal:GroupAB12n2",
                                                            "SampleTypeCecal:Group2243 - SampleTypeCecal:Group15644",
                                                            "SampleTypeCecal:Group2243 - SampleTypeCecal:GroupAB12n2",
                                                            "SampleTypeCecal:Group15644 - SampleTypeCecal:GroupAB12n2",
                                                            "SampleTypeLI:GroupGF - SampleTypeLI:Group2243",
                                                            "SampleTypeLI:GroupGF - SampleTypeLI:Group15644",
                                                            "SampleTypeLI:GroupGF - SampleTypeLI:GroupAB12n2",
                                                            "SampleTypeLI:Group2243 - SampleTypeLI:Group15644",
                                                            "SampleTypeLI:Group2243 - SampleTypeLI:GroupAB12n2",
                                                            "SampleTypeLI:Group15644 - SampleTypeLI:GroupAB12n2")]
  
  #type_mod = lmer(value~SampleType + (1|Animal), data = mod_data)
  site_contrasts_full[,IDUniq:=x]
  nostrain_mod = lmer(log10value~SampleType+Colonized+SampleType:Colonized + (1|Animal), data = mod_data)
  ns_coefs = data.table(tidy(nostrain_mod), Model = "NoStrains")
  # e_full = effects::allEffects(mod_results)
  # e_tab = as.data.table(e_full[[1]])
  # e_tab[,Model:="Full"]
  # e_nostrain = effects::allEffects(nostrain_mod)
  # e_ns_tab = as.data.table(e_nostrain[[1]])
  # e_ns_tab[,Model:="NoStrains"]
  site_contrasts <- difflsmeans(nostrain_mod) %>% as.data.frame() %>% rownames_to_column() %>% data.table()
  site_contrasts <- site_contrasts[rowname %in% c("SampleTypeSI:ColonizedGF - SampleTypeSI:ColonizedEl", "SampleTypeLI:ColonizedGF - SampleTypeLI:ColonizedEl",
                                                  "SampleTypeCecal:ColonizedGF - SampleTypeCecal:ColonizedEl")]
  site_contrasts[,IDUniq:=x]
  chisq_pval = anova(mod_results, nostrain_mod)["Pr(>Chisq)"][2,1]
  mod_all = rbind(mod_coefs, ns_coefs, fill = T)
#  effects_all = rbind(e_tab, e_ns_tab, fill = T)
  if(chisq_pval < 0.05){
    mod_all[,Best:=ifelse(Model == "Full", 1, 0)]
#    effects_all[,Best:=ifelse(Model == "Full", 1, 0)]
  } else {
    mod_all[,Best:=ifelse(Model == "NoStrains", 1, 0)]
#    effects_all[,Best:=ifelse(Model == "NoStrains", 1, 0)]
  }
  mod_all[,IDUniq:=x]
#  effects_all[,ID2:=x]
  mods_all[[i]] = mod_all
#  effects_all_all[[i]] = effects_all
  contrasts_all[[i]] <- site_contrasts
  constrasts_all_strains[[i]] <- site_contrasts_full
}
mods_all = rbindlist(mods_all)
#effects_all = rbindlist(effects_all_all)
#effects_all[,Group2:=ifelse(Model == "Full", as.character(Group), as.character(Colonized))]
contrasts_all <- rbindlist(contrasts_all)
constrasts_all_strains <- rbindlist(constrasts_all_strains)
#save(mods_all, effects_all, feature_list, contrasts_all, contrasts_all_strains, file = "mouse_lmer_results_log10.rda")
save(mods_all, feature_list, contrasts_all, constrasts_all_strains, file = "mouse_lmer_results_log10_wCage.rda")
