### Make standardized processed versions of all datasets
library(data.table)
library(tidyverse)
library(readxl)
library(webchem)
library(omu)
library(broom)
library(lmerTest)
library(fuzzyjoin)
library(broom.mixed)

################################# Pilot
pos_dataset <- data.table(read_xlsx("MetabolomicsElenta5MediasDeFelice2019-11/All Aligned_updatedextra2020-5-18.xlsx", sheet = 1, skip = 1))
pos_dataset <- pos_dataset[,1:81,with=F]
pos_dataset[,IonMode:="Pos"]
samp_table <- data.table(Sample = names(pos_dataset)[38:67])
samp_table[,Condition:=case_when(
  grepl("[A-E]1", Sample) ~ "BHI+Arg",
  grepl("[A-E]2", Sample) ~ "EDM2",
  grepl("[A-E]3", Sample) ~ "EDM1",
  grepl("[A-E]4", Sample) ~ "EDM3",
  grepl("[A-E]5", Sample) ~ "ISP-2+Arg"
)]
samp_table[,Strain:=ifelse(grepl("_[1-3]", Sample), "c", "2243")]
pos_dataset[,IDUniq:=paste0(identifier, "_Positive")]
ids_keep <- c("IDUniq", "Alignment ID", "Average Mz", "Average Rt(min)", "Adduct type", "MSI", "Metabolite name", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode")
samps_all <- c(samp_table[,Sample], names(pos_dataset)[grepl("Pool", names(pos_dataset))], names(pos_dataset)[grepl("Blank", names(pos_dataset))])
pos_dataset <- melt(pos_dataset, id.vars = ids_keep, measure.vars = samps_all, variable.name = "Sample")

neg_dataset <- data.table(read_xlsx("MetabolomicsElenta5MediasDeFelice2019-11/Neg_submit_extraProcessing2020-5-18.xlsx", sheet = 1))
data_with_ms2_neg <- data.table(read_xlsx("DataWithMS2/Pilot spectra.xlsx", sheet = 1))
neg_dataset[,IonMode:="Neg"]
neg_dataset[,IDUniq:=paste0(identifier, "_Negative")]
setnames(neg_dataset, gsub("HILIC_NEG_", "", names(neg_dataset)))
samps_neg <- c(samp_table[,Sample], names(neg_dataset)[grepl("Pool", names(neg_dataset))], names(neg_dataset)[grepl("Blank", names(neg_dataset))])
neg_dataset <- melt(neg_dataset, id.vars = ids_keep[ids_keep %in% names(neg_dataset)], measure.vars = samps_neg, variable.name = "Sample")

add_cols <- ids_keep[ids_keep %in% names(data_with_ms2_neg) & !ids_keep %in% names(neg_dataset)]
neg_dataset <- merge(neg_dataset,
                     data_with_ms2_neg[,names(data_with_ms2_neg) %in% c(add_cols, "Alignment ID"), with = F], by = "Alignment ID", all.x = T)
data_all <- rbind(pos_dataset, neg_dataset, fill = T)
data_all[,MinFeatValue:=min(value[value != 0]), by=IDUniq]
data_all[,log10value:=log10(value + 0.25*MinFeatValue)]
data_all <- merge(data_all, samp_table, by = "Sample", all.x = T)

### Check combining features
setnames(data_all, c("Alignment ID", "Average Mz", "Average Rt(min)", "Adduct type", "Metabolite name"),c("AlignmentID", "AvgMZ", "AvgRt", "Adduct", "MetName"))
ids_keep <- c("IDUniq", "AlignmentID", "AvgMZ", "AvgRt", "Adduct", "MetName", "MSI", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode")
pos_mat <- data_all[IonMode=="Pos",list(IDUniq, AvgMZ, AvgRt, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value")
neg_mat <- pos_mat %>% filter(grepl("Neg", IDUniq))
pos_mat <- pos_mat %>% filter(grepl("Pos", IDUniq))
pos_adducts <- data_all[grepl("Pos", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
neg_adducts <- data_all[grepl("Neg", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
adduct_table <- data.table(expand_grid(pos_adducts, neg_adducts))
adduct_table[pos_adducts == "[M+H]+" & neg_adducts == "[M-H]-", MassDiff:=2*1.007276]
adduct_table[pos_adducts == "[M+H-H2O]+" & neg_adducts == "[M-H]-", MassDiff:=-2*1.007276+18.010565]

pos_mat <- data_all[,list(IDUniq, AvgMZ, AvgRt, Adduct, MetName, MSI, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("Pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M-H]-" ~ as.numeric(AvgMZ)+ 1.007276
)]
pos_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
)]
test_join <- difference_inner_join(neg_mat_sub, pos_mat_sub, by = "ExpectedM", max_dist = 0.002) %>% data.table()
test_join[,RTDiff:=as.numeric(AvgRt.x)-as.numeric(AvgRt.y)]
test_join <- test_join[abs(RTDiff) < 0.1]
### Get correlations too but this seems to have worked well
pos_samps <- data_all[IDUniq %in% test_join$IDUniq.x]
neg_samps <- data_all[IDUniq %in% test_join$IDUniq.y]
pos_samps <- merge(pos_samps, test_join[,list(IDUniq.x, IDUniq.y)], by.x = "IDUniq", by.y = "IDUniq.x", allow.cartesian = T)
all_samps <- merge(pos_samps, neg_samps, by.x = c("IDUniq.y", "Sample"), by.y = c("IDUniq", "Sample"))
corre_results <- all_samps[,cor(value.x, value.y), by = list(IDUniq, IDUniq.y)]
test_join <- merge(test_join, corre_results, by.x = c("IDUniq.x", "IDUniq.y"), by.y = c("IDUniq", "IDUniq.y"), all = T)
test_join[,hist(V1, breaks = 50)]
test_join <- test_join[V1 > 0.7]
test_join[grepl("nknown", MetName.x), MSI.x:=7]
test_join[grepl("nknown", MetName.y), MSI.y:=7]
test_join[MSI.y=="1,4", MSI.y:=1]
exclude_feats <- sapply(1:nrow(test_join), function(x){
  if(is.na(test_join[x, as.numeric(MSI.x)]) & is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y])
  } else if(test_join[x, as.numeric(MSI.x)] < test_join[x, as.numeric(MSI.y)] | is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y]) } else {
      return(test_join[x, IDUniq.x])
    }})

data_all[,Duplicate:=ifelse(IDUniq %in% exclude_feats, 1, 0)]

mean_data <- data_all[,list(mean(value), sd(value), min(value), max(value), median(value),
                            mean(log10value), sd(log10value), min(log10value),
                            max(log10value), median(log10value)), by=c(ids_keep, "Condition", "Strain", "MinFeatValue", "Duplicate")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                      "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

ctrl_means <- mean_data[Strain == "c"]
el_means <- mean_data[Strain == "2243"]
el_means <- merge(el_means, ctrl_means, by = c(ids_keep, "Condition", "MinFeatValue"))
setnames(el_means, gsub(".x", "_Elenta", names(el_means), fixed = T))
setnames(el_means, gsub(".y", "_Ctrl", names(el_means), fixed = T))
el_means[,log2FC:=log2((meanValue_Elenta+0.25*MinFeatValue)/(meanValue_Ctrl+0.25*MinFeatValue))]

ttest_results <- tibble(data_all) %>% filter(!is.na(Condition) & Duplicate == 0) %>%
  group_by(IDUniq, Condition) %>%
  group_modify(
    ~ {
      tryCatch(t.test(log10value~Strain, data = .x) %>% tidy() %>% select(estimate, estimate1, estimate2, statistic, p.value), error = function(e){return(data.frame(NA))})#, otherwise = data.frame(NA)
    }
    #tryCatch(t.test(log10value[Strain == "2243"], log10value[Strain == "c"]) %>% tidy() %>%
    #            select(estimate, estimate1, estimate2, statistic, p.value), error = function(e){return(NA)})

  ) %>% data.table()
ttest_results[,`NA.`:=NULL]
ttest_results[,FDRCorrect:=p.adjust(p.value, method = "BH"), by = Condition]

el_means <- merge(el_means, ttest_results, by = c("IDUniq", "Condition"), all.x = T)
el_means[,Duplicate:=Duplicate_Elenta]
el_means[,Duplicate_Elenta:=NULL]
el_means[,Duplicate_Ctrl:=NULL]

##### Add GNPS annot
molnet_data_pos <- fread("DataWithMS2/GNPS_results/Pilot_Pos/DB_result/5585604187e843e993cef991950926da.tsv")
molnet_data_pos[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Pilot_Pos/output_network/ClassyFireResults_Network.txt")
molnet_data_pos <- merge(molnet_data_pos, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_neg <- fread("DataWithMS2/GNPS_results/Pilot_Neg/DB_result/8167765c239b4f47a0e5fe72b4cfd47f.tsv")
molnet_data_neg[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Pilot_Neg/output_network/ClassyFireResults_Network.txt")
molnet_data_neg <- merge(molnet_data_neg, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_pos[,IonMode:="Pos"]
molnet_data_neg[,IonMode:="Neg"]
molnet_data <- rbind(molnet_data_pos, molnet_data_neg, fill = T)
molnet_data[,AlignmentID:=as.character(AlignmentID)]
el_means <- merge(el_means, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
data_all <- merge(data_all, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
## Leaving out all the carefully merged IDs here

fwrite(el_means, file = "processedDatasets/finalAllDatasets/pilot_mean_summaries.csv")
fwrite(data_all, file = "processedDatasets/finalAllDatasets/pilot_all_data.csv")

inchikeys_all <- el_means[,sort(unique(INCHIKEY))]
fwrite(data.table(INCHIKEY = inchikeys_all), file = "processedDatasets/pilot_inchikeys.txt")
rm(list = ls())

########################## Time course
ids_keep <- c("IDUniq", "AlignmentID", "AvgMZ", "AvgRt", "Adduct", "MSI", "MetName", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode", "BK")
growth_file = "../TimeCourse_take2_2020-8-17/ODdata.xlsx"
growth_data = read_xlsx(growth_file, sheet = 1)
time_pts = read_xlsx(growth_file, sheet = 2)
time_pts = data.table(time_pts %>% mutate(TimePoint = as.numeric(gsub("TP", "", TimePoint))))
growth_data = data.table(growth_data) %>% melt(id.var = c("Strain", "AcConc", "Replicate"))
growth_data[,Time:=as.numeric(gsub(".*_", "", variable))]
growth2 = growth_data
setnames(growth2, "value", "OD600")
growth2[Time==41.5, Time:=42]
growth2[Time==26.5, Time:=26]
growth2[Time==30.5, Time:=30]
mean_growth2 = growth2[, .(meanOD =mean(OD600), sdOD = sd(OD600)), by=list(Strain, AcConc, Time)]
time_course2 = fread("processedDatasets/timeCourse2_allData.txt")

ms2_dat <- data.table(read_xlsx("../DataWithMS2/Neg Combined for Cecilia.xlsx", skip = 4, sheet = "Timecourse"))
ms2_dat <- ms2_dat[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`, BK)]
ms2_dat2 <- data.table(read_xlsx("../DataWithMS2/Pos Combined for Cecilia.xlsx", skip = 4, sheet = "TimeCourse"))
ms2_dat2 <- ms2_dat2[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`, BK)]
ms2_dat[,IonMode:="Neg"]
ms2_dat2[,IonMode:="Pos"]
ms2_dat <- rbind(ms2_dat, ms2_dat2, fill = T)
setnames(ms2_dat, "Alignment ID", "AlignmentID")
time_course2 <- merge(time_course2, ms2_dat, by = c("IonMode", "AlignmentID"))
bad_samps <- time_course2[!grepl("iSTD", MetName),sum(value), by=list(Sample, IonMode)] %>% dcast(Sample~IonMode, value.var= "V1") %>% filter(Neg < 1.7e9 & Pos < 7e9) %>% pull(Sample)
#bad_samps = time_course2[!grepl("iSTD", MetName),sum(value), by=Sample][V1 < 1.5e10, Sample]
time_course2 = time_course2[!Sample %in% bad_samps]
time_course2[,MinFeatValue:=min(value[value != 0 & !is.na(value)]), by = IDUniq]
time_course2[,log10value:=log10(value + MinFeatValue*0.25)]

## Get rid of pools and NAT samples here
time_course2 <- time_course2[Strain %in% c("2243", "AB8", "Val", "Ctrl") & AcConc %in% c("0", "1", "10", "10Add")]

setnames(time_course2, "Adduct type", "Adduct")
pos_mat <- time_course2[IonMode=="Pos",list(IDUniq, AvgMZ, AvgRt, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value")
neg_mat <- pos_mat %>% filter(grepl("Neg", IDUniq))
pos_mat <- pos_mat %>% filter(grepl("Pos", IDUniq))
pos_adducts <- time_course2[grepl("Pos", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
neg_adducts <- time_course2[grepl("Neg", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
adduct_table <- data.table(expand_grid(pos_adducts, neg_adducts))
adduct_table[pos_adducts == "[M+H]+" & neg_adducts == "[M-H]-", MassDiff:=2*1.007276]
adduct_table[pos_adducts == "[M+H-H2O]+" & neg_adducts == "[M-H]-", MassDiff:=-2*1.007276+18.010565]

pos_mat <- time_course2[,list(IDUniq, AvgMZ, AvgRt, Adduct, MetName, MSI, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("Pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M-H]-" ~ as.numeric(AvgMZ)+ 1.007276
)]
pos_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
)]
test_join <- difference_inner_join(neg_mat_sub, pos_mat_sub, by = "ExpectedM", max_dist = 0.002) %>% data.table()
test_join[,RTDiff:=as.numeric(AvgRt.x)-as.numeric(AvgRt.y)]
test_join <- test_join[abs(RTDiff) < 0.1]
### Get correlations too but this seems to have worked well
pos_samps <- time_course2[IDUniq %in% test_join$IDUniq.x]
neg_samps <- time_course2[IDUniq %in% test_join$IDUniq.y]
pos_samps <- merge(pos_samps, test_join[,list(IDUniq.x, IDUniq.y)], by.x = "IDUniq", by.y = "IDUniq.x", allow.cartesian = T)
all_samps <- merge(pos_samps, neg_samps, by.x = c("IDUniq.y", "Sample"), by.y = c("IDUniq", "Sample"))
corre_results <- all_samps[,cor(value.x, value.y), by = list(IDUniq, IDUniq.y)]
test_join <- merge(test_join, corre_results, by.x = c("IDUniq.x", "IDUniq.y"), by.y = c("IDUniq", "IDUniq.y"), all = T)
test_join[,hist(V1, breaks = 50)]
test_join <- test_join[V1 > 0.7]
test_join[grepl("nknown", MetName.x), MSI.x:=7]
test_join[grepl("nknown", MetName.y), MSI.y:=7]
test_join[MSI.y=="1,4", MSI.y:=1]
exclude_feats <- sapply(1:nrow(test_join), function(x){
  if(is.na(test_join[x, as.numeric(MSI.x)]) & is.na(test_join[x, as.numeric(MSI.y)])){ ## Return the one we want to get rid of
    return(test_join[x, IDUniq.y])
  } else if(test_join[x, as.numeric(MSI.x)] < test_join[x, as.numeric(MSI.y)] | is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y]) } else {
      return(test_join[x, IDUniq.x])
    }})

time_course2[,Duplicate:=ifelse(IDUniq %in% exclude_feats, 1, 0)]

## Alternate version with filtered features
#time_course2 <- fread("processedDatasets/time_course_all_data_clustered_features_test.csv")

# Sample != "c_10_TP1_1_R1"
mean_data <- time_course2[,list(mean(value), sd(value), min(value), max(value), median(value),
                               mean(log10value), sd(log10value), min(log10value),
                               max(log10value), median(log10value)), by=c(ids_keep, "Strain", "AcConc", "Time", "TimePoint", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

mean_data[,Strain:=factor(Strain, levels = c("2243", "AB8", "Val", "Ctrl"))]

ctrl_data2 = mean_data[Strain == "Ctrl", list(TimePoint, AcConc, meanValue, sdValue, IDUniq)]
setnames(ctrl_data2, c("meanValue", "sdValue"), c("ctrl_meanValue", "ctrl_sdValue"))
mean_data <- merge(mean_data, ctrl_data2, by = c("IDUniq", "TimePoint", "AcConc"), all.x = T)
mean_data = mean_data[!grepl("iSTD", MetName)]
mean_data[,log2FC:=log2((meanValue+0.25*MinFeatValue)/(ctrl_meanValue+0.25*MinFeatValue))] # Min nonzero value is 133

## Due to late-time point contamination of ctrls, calculating log2fc based on ctrl means from first 4 time points
ctrl_data3 <- time_course2[Strain == "Ctrl" & TimePoint < 5,list(meanValue = mean(value), sdValue = sd(value)), by=list(AcConc, Strain, IonMode, IDUniq, AlignmentID, MetName, MSI)]
setnames(ctrl_data3, c("meanValue", "sdValue"), c("ctrlMeanTP1To4", "ctrlSDTP1To4"))
mean_data = merge(mean_data, ctrl_data3[,list(IDUniq, AcConc, ctrlMeanTP1To4, ctrlSDTP1To4)], by = c("IDUniq", "AcConc"), all.x = T)
mean_data[,log2FC_adj:=log2((meanValue+0.25*MinFeatValue)/(ctrlMeanTP1To4+0.25*MinFeatValue))] # Min nonzero value is 133

time_course2 <- time_course2[!grepl("iSTD", MetName)]
  #data_sub <- time_course2[Strain %in% c("2243", "Ctrl") & AcConc==1 & !grepl("iSTD", MetName)]

## We need to compare with ctrl data from initial time points bc these are clearly contaminated
ttest_results4 <- rbindlist(lapply(c("2243", "AB8", "Val"), function(x){
  data_all <- time_course2[Strain %in% c(x, "Ctrl") & TimePoint <= 4 & Duplicate==0]
  foo <- data_all %>%
  group_by(IDUniq, AcConc, TimePoint) %>%
  group_modify(
    ~ tryCatch(t.test(log10value~Strain, data=.x) %>% tidy() %>% select(estimate, estimate1, estimate2, statistic, p.value)
             , error = function(e){return(data.frame(NA))})
  ) %>% data.table()
  foo[,Strain:=x]
}))
ttest_results4[,NA.:=NULL]
#diff_test_results_4 <- time_course2[TimePoint <= 4,t.test(log10value[Strain != "Ctrl"], log10value[Strain == "Ctrl"], use = "complete.obs")[c("statistic", "p.value")], by = list(IDUniq, TimePoint, Strain)]

ctrl_data4 <- time_course2[TimePoint==4 & Strain == "Ctrl"]

data_sub56 <- merge(time_course2[TimePoint > 4 & Strain %in% c("2243", "Val", "AB8")], ctrl_data4[,list(IDUniq, AcConc, Replicate, value, log10value)], by = c("IDUniq", "AcConc", "Replicate"), all = T)
ttest_results5_6 <- rbindlist(lapply(c("2243", "AB8", "Val"), function(x){
  data_all <- data_sub56[Strain == x & Duplicate == 0] #time_course2[Strain %in% c(x, "Ctrl") & TimePoint <= 4]
  foo <- data_all %>%
    group_by(IDUniq, AcConc, TimePoint) %>%
    group_modify(
      ~ tryCatch(t.test(.x$log10value.x, .x$log10value.y) %>% tidy() %>% select(estimate, estimate1, estimate2, statistic, p.value)
               , error = function(e){return(data.frame(NA))})
    ) %>% data.table()
  foo[,Strain:=x]
}))
ttest_results5_6[,NA.:=NULL]
#diff_test_results_5_6 <- data_sub56[,t.test(log10value.x, log10value.y, use = "complete.obs")[c("statistic", "p.value")], by = list(IDUniq, TimePoint, Strain)]
diff_test_results <- rbind(ttest_results4, ttest_results5_6, fill = T)
diff_test_results[,FDRCorrect:=p.adjust(p.value, method = "BH"), by = list(Strain, AcConc, TimePoint)]

## Going to leave out santaR analysis here for now
mean_data <- merge(mean_data, diff_test_results, by = c("IDUniq", "Strain", "AcConc", "TimePoint"), all.x = T)
#mean_data <- mean_data[!is.na(Strain)]

##### Add GNPS annot
molnet_data_pos <- fread("../DataWithMS2/GNPS_results/TimeCourse_Pos/ProteoSAFe-MOLNETENHANCER-56259b0f-download_data/DB_result/50d72a7f29974149976a7f6a976b66f4.tsv")
molnet_data_pos[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("../DataWithMS2/GNPS_results/TimeCourse_Pos/ProteoSAFe-MOLNETENHANCER-56259b0f-download_data/output_network/ClassyFireResults_Network.txt")
molnet_data_pos <- merge(molnet_data_pos, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_neg <- fread("../DataWithMS2/GNPS_results/TimeCourse_Neg/ProteoSAFe-MOLNETENHANCER-c21d49a4-download_data/DB_result/ec6e29b73d464a8b884d0b972207aa4b.tsv")
molnet_data_neg[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("../DataWithMS2/GNPS_results/TimeCourse_Neg/ProteoSAFe-MOLNETENHANCER-c21d49a4-download_data/output_network/ClassyFireResults_Network.txt")
neg_molnet3 <- fread("../DataWithMS2/GNPS_results/TimeCourse_Neg/ProteoSAFe-MOLNETENHANCER-c21d49a4-download_data/clusterinfo_summary/ee907712595a43b0933ead4b685e833e.tsv")
molnet_data_neg <- merge(molnet_data_neg, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_pos[,IonMode:="Pos"]
molnet_data_neg[,IonMode:="Neg"]
molnet_data <- rbind(molnet_data_pos, molnet_data_neg, fill = T)
mean_data <- merge(mean_data, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
time_course2 <- merge(time_course2, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)

#Alternative version
fwrite(mean_data, file = "processedDatasets/timecourse_mean_summaries_clustered_features_test.csv")

fwrite(mean_data, file = "processedDatasets/finalAllDatasets/timecourse_mean_summaries.csv")

fwrite(time_course2, file = "processedDatasets/finalAllDatasets/timecourse_all_data.csv")

inchikeys_all <- mean_data[,sort(unique(INCHIKEY))]
fwrite(data.table(INCHIKEY = inchikeys_all), file = "processedDatasets/timecourse_inchikeys.txt")

rm(list = ls())
########################## Strains GL
met_data = readRDS(file = "MetabolomicsStrains2020-1-25/strainDataProcessedGrowthPhyloExtraAnnot.rds")
presence_cutoff = 100000

ids_keep <- c("IDUniq", "AlignmentID", "AverageMZ", "AverageRt", "MetName", "MSI", "INCHIKEY", "IonMode", "Formula",
              "Adduct", "MS/MS spectrum")

ms2_dat_neg <- data.table(read_xlsx("DataWithMS2/Neg Combined for Cecilia.xlsx", sheet = "GL", skip = 4))
ms2_dat_neg <- ms2_dat_neg[,list(`Alignment ID`, `Adduct type`, Formula, `MS/MS spectrum`)]
ms2_dat_neg[,IonMode:="Negative"]
ms2_dat_pos <- data.table(read_xlsx("DataWithMS2/Pos Combined for Cecilia.xlsx", sheet = "GL", skip = 4))
ms2_dat_pos <- ms2_dat_pos[,list(`Alignment ID`, `Adduct type`, Formula, `MS/MS spectrum`)]
ms2_dat_pos[,IonMode:="Positive"]
ms2_dat_neg <- rbind(ms2_dat_neg, ms2_dat_pos, fill = T)
met_data <- merge(met_data, ms2_dat_neg, by.x = c("IonMode", "AlignmentID"), by.y = c("IonMode", "Alignment ID"), all.x = T)
setnames(met_data, "Adduct type", "Adduct")

good_samps = met_data[BadSampOD=="None" & !is.na(BadSampOD), sort(unique(SampleID))]
good_features = met_data[value > presence_cutoff & SampleID %in% good_samps, unique(IDUniq)] #OK almost all of them are above in something that passes (and no identified ones are not)
met_data = met_data[IDUniq %in% good_features]

met_data[grepl("Pool", SampleID), Strain:="Pool"]
met_data[grepl("Blank", SampleID), Strain:="Blank"]
met_data[grepl("Ctrl", SampleID), Strain:="Ctrl"]
met_data[,SampleType:=ifelse(Strain=="Pool", "Pool", "Culture")]
met_data[,SampleType:=ifelse(Strain=="Blank", "Blank", SampleType)]
met_data[,SampleType:=ifelse(Strain=="Ctrl", "Ctrl", SampleType)]
met_data[,Present:=ifelse(value > presence_cutoff, 1, 0)]

#Fix met names
met_data[,MetName:=gsub("_\\[.*$", "", MetName)]
met_data[,MetName:=gsub(";.*", "", MetName)]
met_data[,MetName:=gsub("Unknown_", "", MetName)]
met_data[,MetName:=gsub("\\.", "", MetName)]
met_data[,MetName:=gsub("-L-", "-", MetName)]
met_data[,MetName:=gsub(" [A-B]$", "", MetName)]
met_data[,MetName:=gsub("_[A-Za-z]+$", "", MetName)]
met_data[,MetName:=gsub("__.*$", "", MetName)]

strain_names = unique(met_data[!Strain %in% c("Blank", "Ctrl", "Pool"),list(Strain, Name2)])
#Remove features that have inconsistent non-zero values in pooled samples and >50% of non-zero strain samples
#Remove features with cv > 50% in controls and >50% of cultures
met_data_filt = met_data[!BadSampOD %in% c("Ctrl","Sample")]
met_data_means = met_data_filt[,list(mean(value), sd(value), min(value), max(value), 
                                     sum(Present, na.rm = T)/length(Present), sum(value==0, na.rm = T)/length(value)), by = list(IDUniq, MSI, MetName, Strain, SampleType)]
setnames(met_data_means, c("V1", "V2", "V3", "V4", "V5", "V6"), c("meanValue", "sdValue", "minValue", "maxValue", "PresentFraction", "ZeroFraction"))

#Do this including blanks and pools
bad_features1 = met_data_means[(SampleType=="Pool" & ZeroFraction > 0.2), unique(IDUniq)]
bad_features1a = met_data_means[SampleType %in% c("Culture", "Ctrl") & maxValue > presence_cutoff, sum(ZeroFraction > 0)/length(ZeroFraction), by=IDUniq][V1 > 0.5, IDUniq]
bad_features_zeros = intersect(bad_features1, bad_features1a)

#Also those with median cv across nonzero culture samples > 100% or cv in controls > 100%
bad_features_cv = dcast(met_data_means[meanValue > presence_cutoff, median(sdValue/meanValue), by=list(IDUniq, SampleType)], IDUniq~SampleType, value.var = "V1")[Culture > 1 | Ctrl > 1, IDUniq]
#ggplot(bad_features_cv, aes(x=Ctrl, y = Culture, color = ifelse(IDUniq %in% met_data_means[!is.na(MSI), IDUniq], 1, 0))) + geom_point()
met_data_means[!is.na(MSI) & IDUniq %in% bad_features_cv]

met_data_good = met_data_filt[!IDUniq %in% c(bad_features_zeros, bad_features_cv) & !SampleType %in% c("Pool", "Blank")]

all_features = met_data_good[,unique(IDUniq)]

## Less clear that transformation is needed here
ggplot(met_data_good[IDUniq %in% sample(all_features, 30)], aes(x=value)) + geom_histogram() + facet_wrap(~IDUniq, scales = "free")
ggplot(met_data_good[IDUniq %in% sample(all_features, 30)], aes(x=log10(value+100))) + geom_histogram() + facet_wrap(~IDUniq, scales = "free")

met_data_good[,MinFeatValue:=min(value[value != 0 & !is.na(value)]), by = IDUniq]
met_data_good[,log10value:=log10(value + MinFeatValue*0.25)]
met_data_good = met_data_good[SampleID != "GL_PTJ0088_1"]

e_data_good <- met_data_good[grepl("Eggerthella", Name2)|grepl("Ctrl", Name2)]
setnames(e_data_good, c("AverageRt", "AverageMZ"), c("AvgRt", "AvgMZ"))
ids_keep[3:4] <- c("AvgRt", "AvgMZ")

pos_mat <- e_data_good[,list(IDUniq, AvgMZ, AvgRt, SampleID, value)] %>% pivot_wider(names_from = "SampleID", values_from = "value")
neg_mat <- pos_mat %>% filter(grepl("Neg", IDUniq))
pos_mat <- pos_mat %>% filter(grepl("Pos", IDUniq))
pos_adducts <- e_data_good[grepl("Pos", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
neg_adducts <- e_data_good[grepl("Neg", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
adduct_table <- data.table(expand_grid(pos_adducts, neg_adducts))
adduct_table[pos_adducts == "[M+H]+" & neg_adducts == "[M-H]-", MassDiff:=2*1.007276]
adduct_table[pos_adducts == "[M+H-H2O]+" & neg_adducts == "[M-H]-", MassDiff:=-2*1.007276+18.010565]

pos_mat <- e_data_good[,list(IDUniq, AvgMZ, AvgRt, Adduct, MetName, MSI, SampleID, value)] %>% pivot_wider(names_from = "SampleID", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("Pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M-H]-" ~ as.numeric(AvgMZ)+ 1.007276
)]
pos_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
)]
test_join <- difference_inner_join(neg_mat_sub, pos_mat_sub, by = "ExpectedM", max_dist = 0.002) %>% data.table()
test_join[,RTDiff:=as.numeric(AvgRt.x)-as.numeric(AvgRt.y)]
test_join <- test_join[abs(RTDiff) < 0.1]
### Get correlations too but this seems to have worked well
pos_samps <- e_data_good[IDUniq %in% test_join$IDUniq.x]
neg_samps <- e_data_good[IDUniq %in% test_join$IDUniq.y]
pos_samps <- merge(pos_samps, test_join[,list(IDUniq.x, IDUniq.y)], by.x = "IDUniq", by.y = "IDUniq.x", allow.cartesian = T)
all_samps <- merge(pos_samps, neg_samps, by.x = c("IDUniq.y", "SampleID"), by.y = c("IDUniq", "SampleID"))
corre_results <- all_samps[,cor(value.x, value.y), by = list(IDUniq, IDUniq.y)]
test_join <- merge(test_join, corre_results, by.x = c("IDUniq.x", "IDUniq.y"), by.y = c("IDUniq", "IDUniq.y"), all = T)
test_join[,hist(V1, breaks = 50)]
test_join <- test_join[V1 > 0.7]
test_join[grepl("nknown", MetName.x), MSI.x:=7]
test_join[grepl("nknown", MetName.y), MSI.y:=7]
test_join[MSI.y=="1,4", MSI.y:=1]
test_join[MSI.x=="1,4", MSI.x:=1]
exclude_feats <- sapply(1:nrow(test_join), function(x){
  if(is.na(test_join[x, as.numeric(MSI.x)]) & is.na(test_join[x, as.numeric(MSI.y)])){ 
    return(test_join[x, IDUniq.y]) 
  } else if(test_join[x, as.numeric(MSI.x)] < test_join[x, as.numeric(MSI.y)] | is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y]) } else {
      return(test_join[x, IDUniq.x])
    }})

e_data_good[,Duplicate:=ifelse(IDUniq %in% exclude_feats, 1, 0)]

print(names(e_data_good))
mean_data <- e_data_good[,list(mean(value), sd(value), min(value), max(value), median(value), 
                            mean(log10value), sd(log10value), min(log10value), 
                            max(log10value), median(log10value)), by=c(ids_keep, "Name2", "Strain", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue", 
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

all_features = e_data_good[Duplicate == 0,unique(IDUniq)]

met_diff_abund <- rbindlist(lapply(all_features, function(x){
  if(e_data_good[IDUniq==x, length(value[value != 0]) > 31]){
    foo = tidy(e_data_good[IDUniq==x, lm(value~Strain)]) %>% mutate(IDUniq = x) %>% as.data.table #%>% dcast(IDUniq~term, value.var = c("estimate", "p.value"))
    return(foo)
  }
}))
met_diff_abund <- met_diff_abund[term != "(Intercept)"]
met_diff_abund[,Strain:=gsub("Strain", "", term)]

ctrl_means = e_data_good[grepl("Ctrl", Condition),list(mean(value), sd(value), min(value), max(value)), by=list(IDUniq, MSI, MetName)]
setnames(ctrl_means, c("V1", "V2", "V3", "V4"), c("ctrl_meanValue", "ctrl_sdValue", "ctrl_minValue", "ctrl_maxValue"))

mean_data = merge(mean_data, met_diff_abund[,list(Strain, IDUniq, estimate, p.value)], by = c("IDUniq", "Strain"), 
                  all.x = T)
mean_data[,FDRCorrect:=p.adjust(p.value, method = "BH"), by = IDUniq]
mean_data <- merge(mean_data, ctrl_means, by = c("IDUniq", "MSI", "MetName"), all.x = T)
mean_data[,sum(FDRCorrect < 0.1, na.rm = T), by=IDUniq][,table(V1)]
mean_data[,log2FC:=log2((meanValue+0.25*MinFeatValue)/(ctrl_meanValue+0.25*MinFeatValue))]


##### Add GNPS annot
molnet_data_pos <- fread("DataWithMS2/GNPS_results/Strains_Pos/DB_result/2a749710faa4449ebd07416e74e62128.tsv")
molnet_data_pos[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Strains_Pos/output_network/ClassyFireResults_Network.txt")
molnet_data_pos <- merge(molnet_data_pos, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_neg <- fread("DataWithMS2/GNPS_results/Strains_Neg/DB_result/cd095eece21046458f2e8ed83df3b40c.tsv")
molnet_data_neg[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Strains_Neg/output_network/ClassyFireResults_Network.txt")
molnet_data_neg <- merge(molnet_data_neg, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_pos[,IonMode:="Positive"]
molnet_data_neg[,IonMode:="Negative"]
molnet_data <- rbind(molnet_data_pos, molnet_data_neg, fill = T)
molnet_data[,AlignmentID:=as.character(AlignmentID)]
mean_data <- merge(mean_data, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
e_data_good <- merge(e_data_good, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)

inchikeys_all <- mean_data[,sort(unique(INCHIKEY))]
fwrite(data.table(INCHIKEY = inchikeys_all), file = "processedDatasets/strains_gl_inchikeys.txt")

####### Since this is for the paper, only including the Eggerthella strains here
fwrite(mean_data, file = "processedDatasets/finalAllDatasets/strains_edm_mean_summaries.csv")
e_data_good[,KEGG:=NULL]
fwrite(e_data_good, file = "processedDatasets/finalAllDatasets/strain_edm_all_data.csv")

rm(list = ls())

########################## Strains BHI
ids_keep <- c("IDUniq", "Alignment ID", "Average Mz", "Average Rt(min)", "Metabolite name", "MSI", "INCHIKEY", "IonMode", "Formula",
              "Adduct type", "MS/MS spectrum")
ms2_dat_neg <- data.table(read_xlsx("DataWithMS2/Neg Combined for Cecilia.xlsx", sheet = "BHI", skip = 4))
ms2_dat_neg <- ms2_dat_neg[,1:185,with=F]
ms2_dat_neg[,IonMode:="Negative"]
setnames(ms2_dat_neg, gsub("_NEG", "", names(ms2_dat_neg)))
ms2_dat_pos <- data.table(read_xlsx("DataWithMS2/Pos Combined for Cecilia.xlsx", sheet = "BHI", skip = 4))
ms2_dat_pos <- ms2_dat_pos[,1:187,with=F]
ms2_dat_pos[,IonMode:="Positive"]
setnames(ms2_dat_pos, gsub("_POS", "", names(ms2_dat_pos)))
ms2_dat_neg <- rbind(ms2_dat_neg, ms2_dat_pos, fill = T)
met_data <- ms2_dat_neg[,1:186]
met_data[,IDUniq:=paste0(`Average Mz`, "_", `Average Rt(min)`, "_", `Alignment ID`, "_", IonMode)]
met_data <- melt(met_data, id.vars = ids_keep, measure.vars = names(met_data)[grepl("BHI_", names(met_data))], variable.name = "Sample")

samp_table <- data.table(Sample = met_data[,sort(unique(Sample))])
samp_table[,Strain:=gsub("BHI_", "", gsub("_[1-9]$", "", Sample))]
samp_table[,Replicate:=gsub(".*_", "", Sample)]
samp_table[grepl("Pool",Strain), Strain:="Pool"]
samp_table[grepl("Blank",Strain), Strain:="Blank"]
samp_table[,SampleType:=ifelse(Strain=="Pool", "Pool", "Culture")]
samp_table[,SampleType:=ifelse(Strain=="Blank", "Blank", SampleType)]
samp_table[,SampleType:=ifelse(grepl("Ctrl", Strain), "Ctrl", SampleType)]


met_data[grepl("Ctrl", Strain), Strain:="Ctrl"]
met_data <- merge(met_data, samp_table, by = "Sample", all.x = T)
met_data_good <- met_data[!SampleType %in% c("Blank", "Pool")]

## Less clear that transformation is needed here
ggplot(met_data_good[IDUniq %in% sample(all_features, 30)], aes(x=value)) + geom_histogram() + facet_wrap(~IDUniq, scales = "free")
ggplot(met_data_good[IDUniq %in% sample(all_features, 30)], aes(x=log10(value+100))) + geom_histogram() + facet_wrap(~IDUniq, scales = "free")

met_data_good[,MinFeatValue:=min(value[value != 0 & !is.na(value)]), by = IDUniq]
met_data_good[,log10value:=log10(value + MinFeatValue*0.25)]

#e_data_good <- met_data_good[grepl("Eggerthella", Name2)|grepl("Ctrl", Name2)]

mean_data <- met_data_good[,list(mean(value), sd(value), min(value), max(value), median(value), 
                               mean(log10value), sd(log10value), min(log10value), 
                               max(log10value), median(log10value)), by=c(ids_keep, "Strain", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue", 
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

all_features = met_data_good[,unique(IDUniq)]

met_diff_abund <- rbindlist(lapply(all_features, function(x){
  if(met_data_good[IDUniq==x, length(value[value != 0]) > 31]){
    foo = tidy(met_data_good[IDUniq==x, lm(value~Strain)]) %>% mutate(IDUniq = x) %>% as.data.table #%>% dcast(IDUniq~term, value.var = c("estimate", "p.value"))
    return(foo)
  }
}))
met_diff_abund <- met_diff_abund[term != "(Intercept)"]
met_diff_abund[,Strain:=gsub("Strain", "", term)]

ctrl_means = met_data_good[grepl("Ctrl", Strain),list(mean(value), sd(value), min(value), max(value)), by=IDUniq]
setnames(ctrl_means, c("V1", "V2", "V3", "V4"), c("ctrl_meanValue", "ctrl_sdValue", "ctrl_minValue", "ctrl_maxValue"))

mean_data = merge(mean_data, met_diff_abund[,list(Strain, IDUniq, estimate, p.value)], by = c("IDUniq", "Strain"), all.x = T)
mean_data[,FDRCorrect:=p.adjust(p.value, method = "BH"), by = IDUniq]
mean_data <- merge(mean_data, ctrl_means, by = "IDUniq", all.x = T)
mean_data[,sum(FDRCorrect < 0.1, na.rm = T), by=IDUniq][,table(V1)]
mean_data[,log2FC:=log2((meanValue+0.25*MinFeatValue)/(ctrl_meanValue+0.25*MinFeatValue))]

####### Since this is for the paper, only including the Eggerthella strains here
fwrite(mean_data, file = "processedDatasets/finalAllDatasets/strains_bhi_mean_summaries.csv")
fwrite(met_data_good, file = "processedDatasets/finalAllDatasets/strain_bhi_all_data.csv")

rm(list = ls())

########################## Mice Intestinal
ids_keep <- c("IDUniq", "AlignmentID", "AvgMZ", "AvgRt", "Adduct type", "MSI", "MetName", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode")

data_file = "GnotoDataMaggie/origData/Intestinal Contents Submit Dec2020.xlsx"

all_met_data = data.table(read_xlsx(data_file, sheet = 2))
all_met_data <- all_met_data[,1:103]
all_met_data[,IonMode:="Negative"]
pos_data = data.table(read_xlsx(data_file, sheet = 3))
pos_data <- pos_data[,1:104]
pos_data[,IonMode:="Positive"]
pos_data[,Polarity:=NULL]
all_met_data = rbind(all_met_data, pos_data, fill = T)
setnames(all_met_data, names(all_met_data)[1:6], c("MSI", "AlignmentID", "AvgRt", "AvgMZ", "MetName", "INCHIKEY"))

all_met_data[,IDUniq:=paste0(AvgMZ, "_", AvgRt, "_", AlignmentID, "_", IonMode)]

ms2_dat <- data.table(read_xlsx("DataWithMS2/Neg Combined for Cecilia.xlsx", skip = 4, sheet = "Intestinal contents"))
ms2_dat <- ms2_dat[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`)]
ms2_dat2 <- data.table(read_xlsx("DataWithMS2/Pos Combined for Cecilia.xlsx", skip = 4, sheet = "Intestinal Contents"))
ms2_dat2 <- ms2_dat2[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`)]
ms2_dat[,IonMode:="Negative"]
ms2_dat2[,IonMode:="Positive"]
ms2_dat <- rbind(ms2_dat, ms2_dat2, fill = T)
setnames(ms2_dat, "Alignment ID", "AlignmentID")
all_met_data <- merge(all_met_data, ms2_dat, by = c("IonMode", "AlignmentID"))
all_met_data = melt(all_met_data, id.var = ids_keep, variable.name = "Sample", variable.factor = F, 
                    measure.vars = names(all_met_data)[grepl("_El", names(all_met_data))|grepl("_GF", names(all_met_data))|grepl("pool", names(all_met_data))|grepl("Bk", names(all_met_data))])
all_met_data[,MinFeatValue:=min(value[value != 0]), by=IDUniq]
all_met_data[,log10value:=log10(value+0.25*MinFeatValue)]

sample_tab = data.table(Sample =all_met_data[,unique(Sample)])
sample_tab[,c("SampleType", "Group", "Replicate"):=tstrsplit(Sample, split = "_")]
sample_tab[SampleType=="Bk", Replicate:=Group]
sample_tab[SampleType=="Bk", Group:="Bk"]
all_met_data = merge(all_met_data, sample_tab, by = "Sample", all.x = T)
# Initial data checks
all_met_data = all_met_data[!is.na(Group)]
good_samps = all_met_data[!Group %in% c("Bk", "pool"), sort(unique(Sample))]
presence_cutoff = 5*10^4
good_features = all_met_data[value > presence_cutoff & Sample %in% good_samps, unique(IDUniq)] #All
all_met_data = all_met_data[IDUniq %in% good_features]

#Fix met names
all_met_data[,MetName:=gsub("_\\[.*$", "", MetName)]
all_met_data[,MetName:=gsub(";.*", "", MetName)]
all_met_data[,MetName:=gsub("Unknown_", "", MetName)]
all_met_data[,MetName:=gsub("\\.", "", MetName)]
all_met_data[,MetName:=gsub("-L-", "-", MetName)]
all_met_data[,MetName:=gsub(" [A-B]$", "", MetName)]
all_met_data[,MetName:=gsub("_[A-Za-z]+$", "", MetName)]
all_met_data[,MetName:=gsub("__.*$", "", MetName)]
all_met_data[,Present:=ifelse(value > presence_cutoff, 1, 0)]

met_data_means = all_met_data[,list(mean(value), sd(value), sd(value)/sqrt(length(value)), min(value), max(value), sum(Present, na.rm = T)/length(Present), sum(value==0, na.rm = T)/length(value)), by = list(IDUniq, MSI, MetName,  AlignmentID, Group, SampleType)]
setnames(met_data_means, c("V1", "V2", "V3", "V4", "V5", "V6", "V7"), c("meanValue", "sdValue", "seValue", "minValue", "maxValue", "PresentFraction", "ZeroFraction"))
# #Feature filtering
# #Do this including blanks and pools
bad_features1 = met_data_means[(Group=="pool" & ZeroFraction > 0.5), unique(IDUniq)]
bad_features1a = met_data_means[!Group %in% c("pool", "bk") & maxValue > presence_cutoff, sum(ZeroFraction > 0)/length(ZeroFraction), by=IDUniq][V1 > 0.5, IDUniq]
bad_features_zeros = intersect(bad_features1, bad_features1a)

tot_counts =all_met_data[!grepl("iSTD", MetName),sum(value), by = list(Sample, SampleType, IonMode)][order(V1)]
tot_counts = dcast(tot_counts, Sample+SampleType~IonMode, value.var = "V1")
#ggplot(tot_counts, aes(x=Positive, y = Negative, label = Sample)) + geom_point() + geom_text_repel()
#ggplot(tot_counts, aes(x=pos, y = Neg, label = Sample)) + geom_point() + geom_text_repel() + xlim(2.5e10, 4.2e10) + ylim(7.5e9, 1.3e10)
#Clearly we should analyze each sample type separately?
bad_samps = tot_counts[SampleType %in% c("SI", "Cecal") & Positive < 3.5e10, Sample]

istd_data = all_met_data[grepl("iSTD", MetName)]
istd_wide = dcast(istd_data, Sample+IonMode~MetName, value.var = "value")
setnames(istd_wide, gsub(" ", "_", gsub("-", "", names(istd_wide))))
good_samps <- good_samps[!good_samps %in% bad_samps]

met_data_sub = all_met_data[Sample %in% good_samps & !IDUniq %in% bad_features_zeros]
met_data_sub = met_data_sub[!grepl("iSTD", MetName)]
met_data_sub[,Group:=gsub("El", "", Group)]

met_data_sub[,SampleType:=factor(SampleType, levels = c("SI", "Cecal", "LI"))]

met_data_sub[,Present:=ifelse(value > presence_cutoff, 1, 0)]

setnames(met_data_sub, "Adduct type", "Adduct")
ids_keep[ids_keep == "Adduct type"] <- "Adduct"
pos_mat <- met_data_sub[,list(IDUniq, AvgMZ, AvgRt, Adduct, MetName, MSI, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("Pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M-H]-" ~ as.numeric(AvgMZ)+ 1.007276
)]
pos_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
)]
test_join <- difference_inner_join(neg_mat_sub, pos_mat_sub, by = "ExpectedM", max_dist = 0.002) %>% data.table()
test_join[,RTDiff:=as.numeric(AvgRt.x)-as.numeric(AvgRt.y)]
test_join <- test_join[abs(RTDiff) < 0.1]
### Get correlations too but this seems to have worked well
pos_samps <- met_data_sub[IDUniq %in% test_join$IDUniq.x]
neg_samps <- met_data_sub[IDUniq %in% test_join$IDUniq.y]
pos_samps <- merge(pos_samps, test_join[,list(IDUniq.x, IDUniq.y)], by.x = "IDUniq", by.y = "IDUniq.x", allow.cartesian = T)
all_samps <- merge(pos_samps, neg_samps, by.x = c("IDUniq.y", "Sample"), by.y = c("IDUniq", "Sample"))
corre_results <- all_samps[,cor(value.x, value.y), by = list(IDUniq, IDUniq.y)]
test_join <- merge(test_join, corre_results, by.x = c("IDUniq.x", "IDUniq.y"), by.y = c("IDUniq", "IDUniq.y"), all = T)
test_join[,hist(V1, breaks = 50)]
test_join <- test_join[V1 > 0.7]
test_join[grepl("nknown", MetName.x), MSI.x:=7]
test_join[grepl("nknown", MetName.y), MSI.y:=7]
test_join[MSI.y=="1,4", MSI.y:=1]
test_join[MSI.x=="1,4", MSI.x:=1]
exclude_feats <- sapply(1:nrow(test_join), function(x){
  if(is.na(test_join[x, as.numeric(MSI.x)]) & is.na(test_join[x, as.numeric(MSI.y)])){ 
    return(test_join[x, IDUniq.y]) 
  } else if(test_join[x, as.numeric(MSI.x)] < test_join[x, as.numeric(MSI.y)] | is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y]) } else {
      return(test_join[x, IDUniq.x])
    }})

met_data_sub[,Duplicate:=ifelse(IDUniq %in% exclude_feats, 1, 0)]


mean_data <- met_data_sub[,list(mean(value), sd(value), min(value), max(value), median(value), 
                                mean(log10value), sd(log10value), min(log10value), 
                                max(log10value), median(log10value)), by=c(ids_keep, "SampleType", "Group", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue", 
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

gf_data_means = mean_data[Group=="GF", list(IDUniq, meanValue, sdValue, meanLog10Value, sdLog10Value, SampleType)]
setnames(gf_data_means, c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value"), paste0("GF_", c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value")))

mean_data = merge(mean_data, gf_data_means, by = c("IDUniq", "SampleType"), all.x = T)
mean_data[,log2FC_GF:=log2((meanValue + 0.25*MinFeatValue)/(GF_meanValue + 0.25*MinFeatValue))]

met_data_sub[,Colonized:=factor(ifelse(Group=="GF", "GF", "El"), levels = c("GF", "El"))]
met_data_sub[,Animal:=gsub("^[A-Za-z]+_", "", Sample)]
mean_data[,Colonized:=factor(ifelse(Group=="GF", "GF", "El"), levels = c("GF", "El"))]

feature_list <- met_data_sub[Duplicate == 0,sort(unique(IDUniq))]
met_data_sub[,SampleType:=factor(SampleType, levels = c("SI", "Cecal", "LI"))]
met_data_sub[,Group:=factor(Group, levels = c("GF", "2243", "15644", "AB12n2"))]

#### Add cage metadata
cage_metadata <- data.table(read_xlsx("GnotoDataMaggie/origData/ForCN_elen_mono.xlsx", skip = 2))
cage_metadata <- cage_metadata[!is.na(Metabolomics)]
met_data_sub[,Animal2:=as.numeric(gsub(".*_", "", Animal))]
met_data_sub <- merge(met_data_sub, cage_metadata[,list(`Mouse ID`, Cage)], by.x = "Animal2", by.y = "Mouse ID", all.x = T)

#save(met_data_sub, feature_list, file = "mouse_data_models.rda")

# See run_mouse_lmer.R
# mods_all = vector("list", length(feature_list))
# effects_all_all =vector("list", length(feature_list))
# contrasts_all <- vector("list", length(feature_list))
# constrasts_all_strains <- vector("list", length(feature_list))
# for(i in 1:length(feature_list)){
#   x = feature_list[i]
#   mod_data = met_data_sub[IDUniq==x]
#   mod_results = lmer(log10value~SampleType+Group+SampleType:Group + (1|Animal), data = mod_data)
#   mod_coefs = data.table(tidy(mod_results), Model = "Full")
#   site_contrasts_full <- difflsmeans(mod_results) %>% as.data.frame() %>% rownames_to_column() %>% data.table()
#   site_contrasts_full <- site_contrasts_full[rowname %in% c("SampleTypeSI:GroupGF - SampleTypeSI:Group2243",
#                                                             "SampleTypeSI:GroupGF - SampleTypeSI:Group15644",
#                                                             "SampleTypeSI:GroupGF - SampleTypeSI:GroupAB12n2",
#                                                             "SampleTypeSI:Group2243 - SampleTypeSI:Group15644",
#                                                             "SampleTypeSI:Group2243 - SampleTypeSI:GroupAB12n2",
#                                                             "SampleTypeSI:Group15644 - SampleTypeSI:GroupAB12n2",
#                                                             "SampleTypeCecal:GroupGF - SampleTypeCecal:Group2243",
#                                                             "SampleTypeCecal:GroupGF - SampleTypeCecal:Group15644",
#                                                             "SampleTypeCecal:GroupGF - SampleTypeCecal:GroupAB12n2",
#                                                             "SampleTypeCecal:Group2243 - SampleTypeCecal:Group15644",
#                                                             "SampleTypeCecal:Group2243 - SampleTypeCecal:GroupAB12n2",
#                                                             "SampleTypeCecal:Group15644 - SampleTypeCecal:GroupAB12n2",
#                                                             "SampleTypeLI:GroupGF - SampleTypeLI:Group2243",
#                                                             "SampleTypeLI:GroupGF - SampleTypeLI:Group15644",
#                                                             "SampleTypeLI:GroupGF - SampleTypeLI:GroupAB12n2",
#                                                             "SampleTypeLI:Group2243 - SampleTypeLI:Group15644",
#                                                             "SampleTypeLI:Group2243 - SampleTypeLI:GroupAB12n2",
#                                                             "SampleTypeLI:Group15644 - SampleTypeLI:GroupAB12n2")]
#   
#   #type_mod = lmer(value~SampleType + (1|Animal), data = mod_data)
#   site_contrasts_full[,IDUniq:=x]
#   nostrain_mod = lmer(log10value~SampleType+Colonized+SampleType:Colonized + (1|Animal), data = mod_data)
#   ns_coefs = data.table(tidy(nostrain_mod), Model = "NoStrains")
#   # e_full = effects::allEffects(mod_results)
#   # e_tab = as.data.table(e_full[[1]])
#   # e_tab[,Model:="Full"]
#   # e_nostrain = effects::allEffects(nostrain_mod)
#   # e_ns_tab = as.data.table(e_nostrain[[1]])
#   # e_ns_tab[,Model:="NoStrains"]
#   site_contrasts <- difflsmeans(nostrain_mod) %>% as.data.frame() %>% rownames_to_column() %>% data.table()
#   site_contrasts <- site_contrasts[rowname %in% c("SampleTypeSI:ColonizedGF - SampleTypeSI:ColonizedEl", "SampleTypeLI:ColonizedGF - SampleTypeLI:ColonizedEl",
#                                                   "SampleTypeCecal:ColonizedGF - SampleTypeCecal:ColonizedEl")]
#   site_contrasts[,IDUniq:=x]
#   chisq_pval = anova(mod_results, nostrain_mod)["Pr(>Chisq)"][2,1]
#   mod_all = rbind(mod_coefs, ns_coefs, fill = T)
#   #  effects_all = rbind(e_tab, e_ns_tab, fill = T)
#   if(chisq_pval < 0.05){
#     mod_all[,Best:=ifelse(Model == "Full", 1, 0)]
#     #    effects_all[,Best:=ifelse(Model == "Full", 1, 0)]
#   } else {
#     mod_all[,Best:=ifelse(Model == "NoStrains", 1, 0)]
#     #    effects_all[,Best:=ifelse(Model == "NoStrains", 1, 0)]
#   }
#   mod_all[,IDUniq:=x]
#   #  effects_all[,ID2:=x]
#   mods_all[[i]] = mod_all
#   #  effects_all_all[[i]] = effects_all
#   contrasts_all[[i]] <- site_contrasts
#   constrasts_all_strains[[i]] <- site_contrasts_full
# }
# mods_all = rbindlist(mods_all)
# #effects_all = rbindlist(effects_all_all)
# #effects_all[,Group2:=ifelse(Model == "Full", as.character(Group), as.character(Colonized))]
# contrasts_all <- rbindlist(contrasts_all)
# constrasts_all_strains <- rbindlist(constrasts_all_strains)
# #save(mods_all, effects_all, feature_list, contrasts_all, contrasts_all_strains, file = "mouse_lmer_results_log10.rda")
# save(mods_all, feature_list, contrasts_all, constrasts_all_strains, file = "mouse_lmer_results_log10.rda")

#load("GnotoDataMaggie/mouse_lmer_results_log10.rda")
load("GnotoDataMaggie/mouse_lmer_results_log10_wCage.rda")

constrasts_all_strains <- constrasts_all_strains[!IDUniq %in% exclude_feats]
constrasts_all_strains[,FDRCorrect:=p.adjust(`Pr(>|t|)`, method = "BH"), by = rowname]
constrasts_all_strains[,contrast:=gsub(" - [A-Za-z]+:", "-",gsub("Group", "", gsub("SampleType", "", rowname)))]

constrasts_keep <- constrasts_all_strains[grepl("GF", contrast)]
constrasts_keep[,Site:=gsub(":.*", "", contrast)]
constrasts_keep[,Site:=factor(Site, levels = c("SI", "Cecal", "LI"))]
levels(constrasts_keep$Site) <- c("Ileum", "Cecum", "Colon")
constrasts_keep[,Comparison:=gsub(".*:", "", contrast)]
constrasts_keep[,c("Group1", "Group2"):=tstrsplit(Comparison, split = "-")]
constrasts_keep[,Comparison:=paste0(Group2, "-", Group1)]
constrasts_keep[,Estimate:=-1*Estimate]
levels(constrasts_keep$Site) <- c("Ileum", "Cecum", "Colon")

constrasts_all_strains[,Site:=gsub(":.*", "", contrast)]
constrasts_all_strains[,Site:=factor(Site, levels = c("SI", "Cecal", "LI"))]
levels(constrasts_all_strains$Site) <- c("Ileum", "Cecum", "Colon")
constrasts_all_strains[,Comparison:=gsub(".*:", "", contrast)]

constrasts_keep[,SampleType:=case_when(
  Site == "Ileum" ~ "SI",
  Site == "Cecum" ~ "Cecal",
  Site == "Colon" ~ "LI"
)]
## Combine with mean values
mean_data <- merge(mean_data, constrasts_keep[,list(IDUniq, contrast, Estimate, `Std. Error`, df, `t value`, lower, upper, `Pr(>|t|)`, FDRCorrect, Site, Comparison, Group2, SampleType)], 
                   by.x = c("IDUniq", "Group", "SampleType"), by.y = c("IDUniq", "Group2", "SampleType"), all.x = T)


##### Add GNPS annot
molnet_data_pos <- fread("DataWithMS2/GNPS_results/Intestinal_Pos/DB_result/0ff00d1334f84d2d92861ba47cc9c29b.tsv")
molnet_data_pos[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Intestinal_Pos/output_network/ClassyFireResults_Network.txt")
molnet_data_pos <- merge(molnet_data_pos, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_neg <- fread("DataWithMS2/GNPS_results/Intestinal_Neg/DB_result/d910f577fb6b41daa92d7ca849c98af7.tsv")
molnet_data_neg[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Intestinal_Neg/output_network/ClassyFireResults_Network.txt")
molnet_data_neg <- merge(molnet_data_neg, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_pos[,IonMode:="Positive"]
molnet_data_neg[,IonMode:="Negative"]
molnet_data <- rbind(molnet_data_pos, molnet_data_neg, fill = T)
#molnet_data[,AlignmentID:=as.character(AlignmentID)]
mean_data <- merge(mean_data, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
met_data_sub <- merge(met_data_sub, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)

inchikeys_all <- mean_data[,sort(unique(INCHIKEY))]
fwrite(data.table(INCHIKEY = inchikeys_all), file = "processedDatasets/mouse_contents_inchikeys.txt")


fwrite(mean_data, file = "processedDatasets/finalAllDatasets/mouse_intestinal_mean_summaries.csv")
fwrite(met_data_sub, file = "processedDatasets/finalAllDatasets/mouse_intestinal_all_data.csv")

rm(list = ls())

########################## Mice Serum
ids_keep <- c("IDUniq", "Alignment ID", "Average Mz", "Average Rt(min)", "Adduct type", "MSI", "Metabolite name", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode")

data_file = "GnotoDataMaggie/origData/EL Serum Combined Submit Nov 2020.xlsx"

all_met_data = data.table(read_xlsx(data_file, sheet = 2))#,
neg_data = data.table(read_xlsx(data_file, sheet = 3))
all_met_data = rbind(all_met_data, neg_data, fill = T)

setnames(all_met_data, names(all_met_data)[1:7], c("MSI", "AlignmentID", "AvgRt", "AvgMZ", "MetName", "INCHIKEY", "IonMode"))
all_met_data[,IDUniq:=paste0(AvgMZ, "_", AvgRt, "_", AlignmentID, "_", IonMode)]
ids_keep <- c("IDUniq", "AlignmentID", "AvgMZ", "AvgRt", "Adduct type", "MSI", "MetName", "Formula", "INCHIKEY", "MS/MS spectrum", "IonMode")


ms2_dat <- data.table(read_xlsx("DataWithMS2/Neg Combined for Cecilia.xlsx", skip = 4, sheet = "Serum"))
ms2_dat <- ms2_dat[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`)]
ms2_dat2 <- data.table(read_xlsx("DataWithMS2/Pos Combined for Cecilia.xlsx", skip = 4, sheet = "Serum"))
ms2_dat2 <- ms2_dat2[,list(`Alignment ID`, `Adduct type`, `Formula`, `MS/MS spectrum`)]
ms2_dat[,IonMode:="Neg"]
ms2_dat2[,IonMode:="pos"]
ms2_dat <- rbind(ms2_dat, ms2_dat2, fill = T)
setnames(ms2_dat, "Alignment ID", "AlignmentID")
all_met_data <- merge(all_met_data, ms2_dat, by = c("IonMode", "AlignmentID"))
all_met_data = melt(all_met_data, id.var = ids_keep, variable.name = "Sample", variable.factor = F, 
                    measure.vars = names(all_met_data)[grepl("Serum", names(all_met_data))|grepl("pool", names(all_met_data))|grepl("bk", names(all_met_data))])
all_met_data[,MinFeatValue:=min(value[value != 0]), by=IDUniq]
all_met_data[,log10value:=log10(value+0.25*MinFeatValue)]

cage_metadata <- data.table(read_xlsx("GnotoDataMaggie/origData/ForCN_elen_mono.xlsx", skip = 2))
cage_metadata <- cage_metadata[!is.na(Metabolomics)]
#metadata2[,table(Cage, Group)]
all_met_data[,Animal:=as.numeric(gsub("^Serum_.*_", "", Sample))]
all_met_data <- merge(all_met_data, cage_metadata[,list(`Mouse ID`, Cage)], by.x = "Animal", by.y = "Mouse ID", all.x = T)

sample_tab = data.table(Sample =all_met_data[,unique(Sample)])
sample_tab[,c("SampleType", "Group", "Replicate"):=tstrsplit(Sample, split = "_")]
sample_tab[SampleType %in% c("bk", "pool"), Replicate:=Group]
sample_tab[SampleType %in% c("bk", "pool"), Group:=SampleType]

all_met_data = merge(all_met_data, sample_tab, by = "Sample", all.x = T)
# Initial data checks
good_samps = all_met_data[!Group %in% c("bk", "pool"), sort(unique(Sample))]
presence_cutoff = 5*10^4
good_features = all_met_data[value > presence_cutoff & Sample %in% good_samps, unique(IDUniq)] #OK almost all of them are above in something that passes (and no identified ones are not)
all_met_data = all_met_data[IDUniq %in% good_features]

#Fix met names
all_met_data[,MetName:=gsub("_\\[.*$", "", MetName)]
all_met_data[,MetName:=gsub(";.*", "", MetName)]
all_met_data[,MetName:=gsub("Unknown_", "", MetName)]
all_met_data[,MetName:=gsub("\\.", "", MetName)]
all_met_data[,MetName:=gsub("-L-", "-", MetName)]
all_met_data[,MetName:=gsub(" [A-B]$", "", MetName)]
all_met_data[,MetName:=gsub("_[A-Za-z]+$", "", MetName)]
all_met_data[,MetName:=gsub("__.*$", "", MetName)]

met_data_sub = all_met_data[!Group %in% c("bk", "pool")]
gf_data = all_met_data[Group == "GF"]


mean_data <- met_data_sub[,list(mean(value), sd(value), min(value), max(value), median(value), 
                            mean(log10value), sd(log10value), min(log10value), 
                            max(log10value), median(log10value)), by=c(ids_keep, "SampleType", "Group", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue", 
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

gf_data_means = mean_data[Group=="GF", list(IDUniq, meanValue, sdValue, meanLog10Value, sdLog10Value)]
setnames(gf_data_means, c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value"), paste0("GF_", c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value")))

mean_data = merge(mean_data, gf_data_means, by = "IDUniq", all.x = T)
mean_data[,log2FC_GF:=log2((meanValue + 0.25*MinFeatValue)/(GF_meanValue + 0.25*MinFeatValue))]

#Remove features that have inconsistent non-zero values in pooled samples and >50% of non-zero strain sample groups
#Remove features with cv > 50% in controls and >50% of cultures
all_met_data[,Present:=ifelse(value > presence_cutoff, 1, 0)]
met_data_means = all_met_data[,list(mean(value), sd(value), min(value), max(value), sum(Present, na.rm = T)/length(Present), sum(value==0, na.rm = T)/length(value)), by = list(IDUniq, MSI, MetName,  AlignmentID,  Group)]
setnames(met_data_means, c("V1", "V2", "V3", "V4", "V5", "V6"), c("meanValue", "sdValue", "minValue", "maxValue", "PresentFraction", "ZeroFraction"))

bad_features1 = met_data_means[(Group=="pool" & ZeroFraction > 0.2), unique(IDUniq)]
bad_features1a = met_data_means[!Group %in% c("pool", "bk") & maxValue > presence_cutoff, sum(ZeroFraction > 0)/length(ZeroFraction), by=IDUniq][V1 > 0.5, IDUniq]
bad_features_zeros = intersect(bad_features1, bad_features1a)

bad_features2 = met_data_sub[!Group %in% c("pool", "Bk", "bk"), sum(value < presence_cutoff), by = IDUniq][V1 > 18, IDUniq]

met_data_sub = met_data_sub[!IDUniq %in% bad_features_zeros & !IDUniq %in% bad_features2]
met_data_sub = met_data_sub[!grepl("iSTD", MetName)]
met_data_sub = met_data_sub[!Sample %in% c("Serum_El15644_6")] #We'll just exclude this outlier for now

## Adduct checking
setnames(met_data_sub, "Adduct type", "Adduct")
ids_keep[ids_keep == "Adduct type"] <- "Adduct"
pos_mat <- met_data_sub[,list(IDUniq, AvgMZ, AvgRt, Adduct, MetName, MSI, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M-H]-" ~ as.numeric(AvgMZ)+ 1.007276
)]
pos_mat_sub[,ExpectedM:=case_when(
  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
)]
test_join <- difference_inner_join(neg_mat_sub, pos_mat_sub, by = "ExpectedM", max_dist = 0.002) %>% data.table()
test_join[,RTDiff:=as.numeric(AvgRt.x)-as.numeric(AvgRt.y)]
test_join <- test_join[abs(RTDiff) < 0.1]
### Get correlations too but this seems to have worked well
pos_samps <- met_data_sub[IDUniq %in% test_join$IDUniq.x]
neg_samps <- met_data_sub[IDUniq %in% test_join$IDUniq.y]
pos_samps <- merge(pos_samps, test_join[,list(IDUniq.x, IDUniq.y)], by.x = "IDUniq", by.y = "IDUniq.x", allow.cartesian = T)
all_samps <- merge(pos_samps, neg_samps, by.x = c("IDUniq.y", "Sample"), by.y = c("IDUniq", "Sample"))
corre_results <- all_samps[,cor(value.x, value.y), by = list(IDUniq, IDUniq.y)]
test_join <- merge(test_join, corre_results, by.x = c("IDUniq.x", "IDUniq.y"), by.y = c("IDUniq", "IDUniq.y"), all = T)
test_join[,hist(V1, breaks = 50)]
test_join <- test_join[V1 > 0.7]
test_join[grepl("nknown", MetName.x), MSI.x:=7]
test_join[grepl("nknown", MetName.y), MSI.y:=7]
test_join[MSI.y=="1,4", MSI.y:=1]
test_join[MSI.x=="1,4", MSI.x:=1]
test_join[MSI.y=="2,4", MSI.y:=2]
test_join[MSI.x=="2,4", MSI.x:=2]

exclude_feats <- sapply(1:nrow(test_join), function(x){
  if(is.na(test_join[x, as.numeric(MSI.x)]) & is.na(test_join[x, as.numeric(MSI.y)])){ 
    return(test_join[x, IDUniq.y]) 
  } else if(test_join[x, as.numeric(MSI.x)] < test_join[x, as.numeric(MSI.y)] | is.na(test_join[x, as.numeric(MSI.y)])){
    return(test_join[x, IDUniq.y]) } else {
      return(test_join[x, IDUniq.x])
    }})

met_data_sub[,Duplicate:=ifelse(IDUniq %in% exclude_feats, 1, 0)]


mean_data <- met_data_sub[,list(mean(value), sd(value), min(value), max(value), median(value), 
                                mean(log10value), sd(log10value), min(log10value), 
                                max(log10value), median(log10value)), by=c(ids_keep, "SampleType", "Group", "MinFeatValue")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue", 
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

gf_data_means = mean_data[Group=="GF", list(IDUniq, meanValue, sdValue, meanLog10Value, sdLog10Value)]
setnames(gf_data_means, c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value"), paste0("GF_", c("meanValue", "sdValue", "meanLog10Value", "sdLog10Value")))

mean_data = merge(mean_data, gf_data_means, by = "IDUniq", all.x = T)
mean_data[,log2FC_GF:=log2((meanValue + 0.25*MinFeatValue)/(GF_meanValue + 0.25*MinFeatValue))]

########## Define increased/decresaed, Combine all the datasets

good_features = met_data_sub[Duplicate == 0,sort(unique(IDUniq))]
## This one seems to be better with log10values
met_data_sub[,Group:=factor(Group, levels = c("GF", "El2243", "El15644", "ElAB12n2"))]
# met_diff_abund2 = rbindlist(lapply(good_features, function(x){
#   ctrl_vals = met_data_sub[IDUniq==x & Group=="GF", value]
#   foo = data.table(broom::tidy(met_data_sub[IDUniq==x,lm(log10value~Group)]))
#   foo[,IDUniq:=x]
#   foo[,CtrlMean:=mean(ctrl_vals, na.rm = T)]
#   return(foo[term != "(Intercept)"])
# }))

met_diff_abund2 = rbindlist(lapply(good_features, function(x){
  ctrl_vals = met_data_sub[IDUniq==x & Group=="GF", value]
  foo = data.table(broom.mixed::tidy(met_data_sub[IDUniq==x,lmer(log10value~Group + (1|Cage))]))
  foo[,IDUniq:=x]
  foo[,CtrlMean:=mean(ctrl_vals, na.rm = T)]
  return(foo[effect == "fixed" & term != "(Intercept)"])
}))
met_diff_abund2[,group:=NULL]
met_diff_abund2[,effect:=NULL]
met_diff_abund2[,Group:=gsub("Group", "", term)]
met_diff_abund2[,FDRCorrect:=p.adjust(p.value, method = "BH"), by=Group]
mean_data = merge(mean_data, met_diff_abund2[,list(Group, IDUniq, estimate, statistic, p.value,FDRCorrect)], by = c("IDUniq", "Group"), all.x = T)

molnet_data_pos <- fread("DataWithMS2/GNPS_results/Serum_Pos/DB_result/9ac4218781ce49aba78d34457e910cb9.tsv")
molnet_data_pos[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Serum_Pos/output_network/ClassyFireResults_Network.txt")
molnet_data_pos <- merge(molnet_data_pos, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_neg <- fread("DataWithMS2/GNPS_results/Serum_Neg/DB_result/d8ce7865c8ba49ccb8b5d66985147de6.tsv")
molnet_data_neg[,AlignmentID:=as.numeric(gsub("spectra_filtered/specs_ms.mgf","", FileScanUniqueID ))]
molnet_class_data <- fread("DataWithMS2/GNPS_results/Serum_Neg/output_network/ClassyFireResults_Network.txt")
molnet_data_neg <- merge(molnet_data_neg, molnet_class_data, by.x = "AlignmentID", by.y = "cluster index", all = T)

molnet_data_pos[,IonMode:="pos"]
molnet_data_neg[,IonMode:="Neg"]
molnet_data <- rbind(molnet_data_pos, molnet_data_neg, fill = T)
#molnet_data[,AlignmentID:=as.character(AlignmentID)]
mean_data <- merge(mean_data, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)
met_data_sub <- merge(met_data_sub, molnet_data[,list(AlignmentID, IonMode, Compound_Name, CF_kingdom, CF_superclass, CF_class, CF_subclass, CF_Dparent)], by = c("IonMode", "AlignmentID"), all.x = T)

inchikeys_all <- mean_data[,sort(unique(INCHIKEY))]
fwrite(data.table(INCHIKEY = inchikeys_all), file = "processedDatasets/mouse_serum_inchikeys.txt")


fwrite(mean_data, file = "../processedDatasets/finalAllDatasets/mouse_serum_mean_summaries.csv")
fwrite(met_data_sub, file = "../processedDatasets/finalAllDatasets/mouse_serum_all_data.csv")


################ Read in all files, get all unique INCHIKEYs, map to KEGG, 
inchikeys_pilot <- fread("processedDatasets/pilot_inchikeys.txt")
inchikeys_time <- fread("processedDatasets/timecourse_inchikeys.txt")
inchikeys_strains <- fread("processedDatasets/strains_gl_inchikeys.txt")
inchikeys_mouse <- fread("processedDatasets/mouse_contents_inchikeys.txt")
inchikeys_serum <- fread("processedDatasets/mouse_serum_inchikeys.txt")
inchikeys_all <- data.table(INCHIKEY = sort(unique(c(inchikeys_pilot$INCHIKEY, inchikeys_time$INCHIKEY, inchikeys_strains$INCHIKEY, 
                               inchikeys_mouse$INCHIKEY, inchikeys_serum$INCHIKEY))))
inchikeys_all[,KEGG:=""]
inchikeys_all[,HMDB:=""]
for(j in 1:nrow(inchikeys_all)){
  inchikeys_all[j,KEGG:=cts_convert(INCHIKEY, from = "InChIKey", to = "KEGG", match = "first")]
  inchikeys_all[j, HMDB:=cts_convert(INCHIKEY, from = "InChiKey", to = "human metabolome database", match = "first")]
}
kegg_class <- assign_hierarchy(inchikeys_all, keep_unknowns = T, identifier = "KEGG")
kegg_class <- merge(kegg_class, inchikeys_all[,list(INCHIKEY, HMDB)], by = "INCHIKEY", all.x = T)
fwrite(kegg_class, file = "processedDatasets/finalAllDatasets/inchikeys_kegg_all.csv")

# hmdb_categories <- fread("metaboanalystR_master_compound_db.tsv")
# hmdb_categories[,HMDB2:=gsub("HMDB", "HMDB00", HMDB2)]
# hmdb_categories <- hmdb_categories[HMDB2 %in% inchikeys_all[,HMDB]]
# hmdb_categories <- hmdb_categories[-hmdb_categories[, which(duplicated(HMDB2))],]

#hmdb_categories2 <- nest(hmdb_categories[,list(kegg_id, exactmass, formula, super_class, main_class, sub_class, synonym, HMDB2)],data = c(kegg_id, exactmass, formula, super_class, main_class, sub_class, synonym))
#inchikeys_all <- merge(inchikeys_all, hmdb_categories[,list(kegg_id, exactmass, formula, super_class, main_class, sub_class, synonym, HMDB2)], 
#                       by.x = "HMDB", by.y = "HMDB2", all.x = T)
# library(omu)
# class_assign <- assign_hierarchy(data.table(KEGG = comp_set[!is.na(kegg_id) & kegg_id != "",kegg_id], IDUniq = comp_set[!is.na(kegg_id) & kegg_id != "",IDUniq]), keep_unknowns = T, identifier = "KEGG")

################################### Targeted acetate dataset
rm(list = ls())
datadir <- "../SILAcetateTimeCourse/"
  
layout_file <- paste0(datadir, "Dodd_AcetatePlateSubset2021-05-18.xlsx")
targeted_file <- paste0(datadir, "Dodd_20210616 UCSF_Cecilia_invitro_sample_NPH.xlsx")

layout_tab <- data.table(read_xlsx(layout_file))
layout_tab <- layout_tab[1:8, 1:13, with=F]
setnames(layout_tab, c("Row", 1:12))
layout_tab <- melt(layout_tab, id.var = "Row", variable.name = "Column", value.name = "Sample")
layout_tab[,well:=paste0(Row, Column)]
layout_tab <- layout_tab[!is.na(Sample)]
layout_tab[,c("Strain", "AcetateConc", "TimePoint", "Replicate"):=tstrsplit(Sample, split = "_")]
time_pt_tab <- data.table(TimePoint=paste0("TP", 0:8), Time = c(0, 15, 18.5, 23,28.5, 39, 43, 47.5, 64))

met_data <- read_xlsx(targeted_file)
loqs <- data.table(met_data[3:4, seq(from = 3, to = ncol(met_data), by = 4)])
setnames(loqs, unlist(met_data[5, seq(from = 3, to = ncol(met_data), by = 4)]))
loqs[,Limit:=c("ULOQ", "LLOQ")]
loqs <- dcast(melt(loqs, id.var = "Limit", variable.name = "Compound"), Compound~Limit)

#met_data <- met_data[5:nrow(met_data),]
met_list <- loqs[,unique(Compound)]
setnames(met_data, c(c("Sample", "Location"), sapply(met_list, function(x){ return(c(as.character(x), paste0(x, "_Accuracy"), paste0(x, "_S_N"), paste0(x, "_iSTD_Area")))})))
met_data <- met_data[7:nrow(met_data),]
met_data <- data.table(met_data)
met_data <- melt(met_data, id.var = c("Sample", "Location"))
met_data[,Compound:=gsub("_.*", "", variable)]
met_data[,Type:=ifelse(variable == Compound, "value", gsub("^[A-Za-z0-9\\ \\-]+_", "", variable))]
met_data <- dcast(met_data, Sample+Location+Compound~Type, value.var = "value")
setnames(met_data, "value", "mM")
met_data[,well:=gsub("S.", "", Sample)]
setnames(met_data, "Sample", "MS_Sample")
met_data <- merge(met_data, loqs, by = "Compound", all.x = T)
met_data <- merge(met_data, layout_tab, by = "well", all.x = T)
met_data[,mM:=as.numeric(mM)]
met_data[,HypAcetateConc:=ifelse((TimePoint == "TP4" | (TimePoint=="TP7" & Replicate < 3))& AcetateConc == 1, "10", AcetateConc)]
met_data[,HypAcetateConc:=ifelse((TimePoint == "TP4" | TimePoint=="TP7" & Replicate < 3) & AcetateConc == 10, "1", HypAcetateConc)]
met_data[,HypAcetateConc2:=case_when(
  TimePoint %in% c("TP4", "TP7") & AcetateConc == 1 ~ "10",
  TimePoint %in% c("TP4", "TP7") & AcetateConc == 10 ~ "1",
  TRUE ~ AcetateConc
)]

met_data[,LLOQ:=as.numeric(LLOQ)]
met_data[,ULOQ:=as.numeric(ULOQ)]


growth_data <- data.table(read_xlsx(paste0(datadir, "ODdata.xlsx"), col_types = c("text", "text", "text", 
                                                                                  rep("numeric", 9))))
growth_data <- melt(growth_data, id.var = c("Strain", "AcConc", "Replicate"), variable.name = "TimePoint", value.name = "OD600")
growth_data[,Time:=gsub(".*_", "", TimePoint)]
growth_data[,TimePoint:=gsub("_.*", "", TimePoint)]
growth_data[,TimePoint:=paste0("TP", as.numeric(gsub("TP", "", TimePoint))-1)]
setnames(growth_data, "AcConc", "AcetateConc")
growth_data[Strain == "Valencia", Strain:="Val"]
growth_data[Strain == "AB8n2", Strain:="AB8"]
met_data <- merge(met_data, growth_data, by = c("Strain", "AcetateConc", "TimePoint", "Replicate"), all.x = T)

met_data[,AcetateConc:=HypAcetateConc]
met_data[,HypAcetateConc:=NULL]
met_data[,HypAcetateConc2:=NULL]

fwrite(met_data, file = "processedDatasets/targeted_met_data.csv", sep = ",")
