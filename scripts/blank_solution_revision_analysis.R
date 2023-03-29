### Process and analyze trace solutions
library(data.table)
library(tidyverse)
library(readxl)
library(webchem)
library(omu)
library(broom)
library(lmerTest)
library(fuzzyjoin)
library(broom.mixed)
theme_set(ggpubr::theme_classic2())

pos_dataset <- data.table(read_xlsx("../Solutions_TraceCheck_Revisions/PETU006_CombineAlign_Submit.xlsx", sheet = "Pos All Features", skip = 3))
pos_dataset <- pos_dataset[,1:28,with=F]
pos_dataset[,IonMode:="Pos"]
samp_table <- data.table(Sample = names(pos_dataset)[8:28])
samp_table[,Solution:=case_when(
  grepl("AA_", Sample) ~ "Amino acids",
  grepl("Arg_", Sample) ~ "Arginine",
  grepl("Ura_", Sample) ~ "Uracil"
)]
samp_table[,Conc:=gsub("_rep.*", "", gsub("PETU006_Pos_[A-Za-z]+_", "", Sample))]
samp_table[grepl("BK_", Sample), Conc:=NA]
samp_table[,Rep:=ifelse(!grepl("BK", Sample), gsub(".*x_rep", "", gsub("_QE1.*", "", Sample)), NA)]
samp_table[,Strain:=ifelse(grepl("_[1-3]", Sample), "c", "2243")]
pos_dataset[,IDUniq:=paste0(`Average Mz`, "_", `Average Rt(min)`, "_", `Alignment ID`, "_Positive")]

ids_keep <- c("IDUniq", "Alignment ID", "Average Mz", "Average Rt(min)",  "MSI", "Metabolite name",  "INCHIKEY", "IonMode")
pos_dataset <- melt(pos_dataset, id.vars = ids_keep, measure.vars = samp_table[,Sample], variable.name = "Sample")

neg_dataset <- data.table(read_xlsx("../Solutions_TraceCheck_Revisions/PETU006_CombineAlign_Submit.xlsx", sheet = "Neg All Features", skip = 3))
neg_dataset[,IonMode:="Neg"]
neg_dataset[,IDUniq:=paste0(`Average Mz`, "_", `Average Rt(min)`, "_", `Alignment ID`, "_Negative")]

samp_table[,NegSample:=names(neg_dataset)[8:28]] #Same alphabetical order
# Fix samp table
samp_table_neg <- samp_table[,2:ncol(samp_table)]
setnames(samp_table_neg, "NegSample", "Sample")
samp_table_pos <- samp_table[,1:(ncol(samp_table)-1)]
samp_table <- rbind(samp_table_pos, samp_table_neg, fill = T)
neg_dataset <- melt(neg_dataset, id.vars = ids_keep[ids_keep %in% names(neg_dataset)], measure.vars = samp_table[grepl("Neg", Sample),Sample], variable.name = "Sample")
data_all <- rbind(pos_dataset, neg_dataset, fill = T)
data_all[,MinFeatValue:=min(value[value != 0]), by=IDUniq]
data_all[,log10value:=log10(value + 0.25*MinFeatValue)]
data_all <- merge(data_all, samp_table, by = "Sample", all.x = T)

setnames(data_all, c("Alignment ID", "Average Mz", "Average Rt(min)", "Metabolite name"),c("AlignmentID", "AvgMZ", "AvgRt", "MetName"))

# pos_mat <- data_all[IonMode=="Pos",list(IDUniq, AvgMZ, AvgRt, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value")
# neg_mat <- pos_mat %>% filter(grepl("Neg", IDUniq))
# pos_mat <- pos_mat %>% filter(grepl("Pos", IDUniq))
# pos_adducts <- data_all[grepl("Pos", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
# neg_adducts <- data_all[grepl("Neg", IDUniq), length(unique(IDUniq)), by = Adduct][V1 > 2, Adduct]
# adduct_table <- data.table(expand_grid(pos_adducts, neg_adducts))
# adduct_table[pos_adducts == "[M+H]+" & neg_adducts == "[M-H]-", MassDiff:=2*1.007276]
# adduct_table[pos_adducts == "[M+H-H2O]+" & neg_adducts == "[M-H]-", MassDiff:=-2*1.007276+18.010565]
data_all[,Sample:=gsub("Neg_", "", gsub("Pos_", "", Sample))]
data_all[,Sample:=gsub("_QE1.*", "", Sample)]
samp_table[,SampleAll:=gsub("Neg_", "", gsub("Pos_", "", Sample))]
samp_table[,SampleAll:=gsub("_QE1.*", "", SampleAll)]

pos_mat <- data_all[,list(IDUniq, AvgMZ, AvgRt, MetName, MSI, Sample, value)] %>% pivot_wider(names_from = "Sample", values_from = "value") %>% data.table()
neg_mat_sub <- data.table(pos_mat[grepl("Neg", IDUniq),1:6, with = F])
pos_mat_sub <- data.table(pos_mat[grepl("Pos", IDUniq),1:6,with=F])
neg_mat_sub[,ExpectedM:=as.numeric(AvgMZ)+ 1.007276]
pos_mat_sub[,ExpectedM:=as.numeric(AvgMZ)- 1.007276]
#  Adduct=="[M+H]+" ~ as.numeric(AvgMZ)- 1.007276,
#  Adduct=="[M+H-H2O]+" ~ as.numeric(AvgMZ) -1.007276+18.010565,
#)]
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
                            max(log10value), median(log10value)), by=c("IDUniq", "AlignmentID", "AvgMZ", "AvgRt", "MSI", "MetName", "INCHIKEY", "IonMode", "Solution", "Conc", "MinFeatValue", "Duplicate")]
setnames(mean_data, paste0("V", 1:10), c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                                         "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"))

blank_means <- mean_data[is.na(Solution)]
setnames(blank_means, c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                        "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value"), paste0("Blank_", c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                                                                                                                                    "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value")))
mean_data <- merge(mean_data, blank_means[,c("IDUniq", paste0("Blank_", c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                                                                          "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value")))], by = "IDUniq", all.x = T)

mean_data[,Present1:=ifelse(meanValue > 3*Blank_meanValue, 1, 0)]

data_all <- merge(data_all, blank_means[,c("IDUniq", paste0("Blank_", c("meanValue", "sdValue", "minValue", "maxValue", "medianValue",
                                                                        "meanLog10Value", "sdLog10Value", "minLog10Value", "maxLog10Value", "medianLog10Value")))], by = "IDUniq", all.x = T)
data_all[,Present1:=ifelse(value > 3*Blank_meanValue, 1, 0)]


comp_order <- mean_data[Solution == "Amino acids" & Conc == "2x"][order(meanLog10Value), paste0(MetName, " (", MSI, ")_", IonMode, AlignmentID)]
ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][3589]", MetName) & Duplicate == 0], aes(x = Sample, y = factor(paste0(MetName, "_", IonMode, AlignmentID), levels = comp_order), fill = log10value)) + 
  geom_tile() + scale_fill_distiller(palette = "Blues", direction = 1) + ylab("") #+ facet_wrap(~IonMode, scales = "free_y")

ggplot(mean_data[!is.na(MSI) & !grepl("[_\\-][Dd][3589]", MetName) & Duplicate == 0], aes(x = paste0(Solution, Conc), y = factor(paste0(MetName, "_", IonMode, AlignmentID), levels = comp_order), fill = meanLog10Value)) + 
  geom_tile() + scale_fill_distiller(palette = "Blues", direction = 1) + ylab("") + facet_wrap(~Present1, scales = "free") #+ facet_wrap(~IonMode, scales = "free_y")

ggplot(mean_data[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & Present1 == 1], aes(x = paste0(Solution, Conc), y = factor(paste0(MetName, "_", IonMode, AlignmentID), levels = comp_order), fill = meanLog10Value)) + 
  geom_tile() + scale_fill_distiller(palette = "Blues", direction = 1) + ylab("") #+ facet_wrap(~Present1, scales = "free") #+ facet_wrap(~IonMode, scales = "free_y")

ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & Present1 == 1], aes(x = gsub("PETU006_", "", Sample), y = factor(paste0(MetName, "_", IonMode, AlignmentID), levels = comp_order), fill = log10value)) + 
  geom_tile(color = "black") + scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey90") + ylab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_line(color = "black")) + 
  xlab("Sample")#+ facet_wrap(~IonMode, scales = "free_y") 

ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & Present1 == 1 & grepl("1x", Sample)], aes(x = gsub("PETU006_", "", Sample), y = factor(paste0(MetName, "_", IonMode, AlignmentID), levels = comp_order), fill = log10value)) + 
  geom_tile(color = "black") + scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey90") + ylab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_line(color = "black")) + 
  xlab("Sample")#+ facet_wrap(~IonMode, scales = "free_y") 

met_names_all <- data.table(MetName = data_all[!grepl("[_\\-][Dd][35789]", MetName),sort(unique(MetName))])
met_names_all[,Intended:=c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 
                           0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1,1, 0, 1,1)]

data_all <- merge(data_all, met_names_all, by = "MetName", all.x = T)
mean_data <- merge(mean_data, met_names_all, by = "MetName", all.x = T)
data_all[Solution == "Arginine" & !grepl("Arginine", MetName), Intended:=0]
data_all[Solution == "Uracil" & !grepl("Uracil", MetName), Intended:=0]
simple_heatmap_plot <- ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & 
                                         Present1 == 1 & grepl("1x", Sample)], 
                              aes(x = gsub("PETU006_", "", Sample), 
                                  y = factor(paste0(MetName," (", MSI, ")_", IonMode, AlignmentID), levels = comp_order), fill = log10value)) + 
  geom_tile(color = "black") + scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey90") + ylab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_line(color = "black")) + 
  xlab("Sample") + geom_point(data = data_all[Intended==1 & Duplicate == 0 & grepl("1x", Sample)], size = 0.5)#+ facet_wrap(~IonMode, scales = "free_y") 
ggsave(simple_heatmap_plot, file = "solutions_simpleHeatmap.pdf", width = 8, height = 9)

time_series_dat <- fread("processedDatasets/timecourse_mean_summaries.csv")

#Check for dipeptide-like features -nope
dipeptides <- time_series_dat[grepl("Arg", MetName) & grepl("\\-", MetName)]

dipeptide_info <- unique(dipeptides[,list(IonMode, AlignmentID, IDUniq, AvgMZ, AvgRt, Adduct, MSI, MetName, Formula, INCHIKEY)])
new_dat <- unique(mean_data[abs(AvgMZ-300) < 35 & is.na(MSI) & AvgRt > 8, list(IDUniq, AlignmentID, AvgMZ, AvgRt, MSI, MetName, INCHIKEY, IonMode)])
compare_test <- difference_inner_join(new_dat, dipeptide_info, by = "AvgMZ", max_dist = 0.02)

dipeptides[TimePoint==1, mean(ctrlMeanTP1To4), by=list(MetName, IDUniq, AcConc)]
fwrite(dipeptide_info[order(IonMode, AvgMZ)], file = "time_course_dipeptide_features.txt", quote = F, sep = "\t")

## Check how many detected features are also in the time course dataset
met_names_all <- merge(met_names_all, unique(data_all[,list(MetName, INCHIKEY)]), by = "MetName", all.x = T) #[,INCHIKEY:=sapply(MetName, function(x){ return(data_all[MetName==x, unique(INCHIKEY)])})]
met_names_all[,PresentInTimeCourseCtrl:=ifelse(INCHIKEY %in% time_series_dat[Strain == "Ctrl" & TimePoint == 1, INCHIKEY], 1, 0)]
met_names_all[,PresentInTimeCourseCtrl2:=ifelse(tolower(MetName) %in% time_series_dat[Strain == "Ctrl" & TimePoint == 1, tolower(MetName)], 1, 0)]
library(webchem)
cts_convert(query = "Valine", from = "Chemical Name", to= "InChiKey")
met_names_all[,OstensibleName:=cts_convert(query = INCHIKEY, from = "InChiKey", to = "Chemical Name")]
## Brian's Inchikeys are different in this dataset!!

new_dat_feature_tab <- unique(data_all[,list(IDUniq, AvgMZ, AvgRt, IonMode, MSI, MetName, INCHIKEY)])
time_series_features <- unique(time_series_dat[,list(IDUniq, AvgMZ, AvgRt, IonMode, MSI, MetName, INCHIKEY)])
library(fuzzyjoin)
compare_dat <- data.table(difference_left_join(new_dat_feature_tab, time_series_features, by = "AvgMZ", max_dist = 0.001))
compare_dat[,MzDiff:=abs(AvgMZ.x-AvgMZ.y)]
compare_dat[abs(AvgRt.x-AvgRt.y) > 0.5 | IonMode.x != IonMode.y, MzDiff:=NA]
#compare_dat <- compare_dat[abs(AvgRt.x-AvgRt.y) < 0.5 | is.na(IDUniq.y)]
compare_dat[,PairRank:=rank(MzDiff, ties = "first"), by = IDUniq.x]

compare_dat <- compare_dat[PairRank==1|is.na(IDUniq.y)]
#compare_dat[,GoodPair:=ifelse((abs(AvgRt.x-AvgRt.y) < 0.2 & IonMode.x == IonMode.y), 1, 0)]

View(compare_dat[,list(IDUniq.x, IDUniq.y, MetName.x, MetName.y, MSI.x, MSI.y)])
compare_dat[,length(unique(IDUniq.y[!is.na(IDUniq.y)])), by=IDUniq.x][,table(V1)]

data_all[,HasComparisonInTimeSeries:=ifelse(IDUniq %in% compare_dat[!is.na(IDUniq.y) & !is.na(MzDiff), IDUniq.x], 1, 0)]
fwrite(data_all[,list(Sample, Solution, Conc, Rep, Strain, IDUniq, MetName, AlignmentID, AvgMZ, AvgRt, MSI, INCHIKEY, IonMode, Duplicate, Present1, Intended, HasComparisonInTimeSeries, value, log10value)], file = "figure_source_data/figureS1e_blankSolutions.tsv", quote=F, row.names = F, sep = "\t")
simple_heatmap_plot <- ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & Present1 == 1 & !grepl("BK", Sample)], 
                              aes(x = gsub("PETU006_", "", Sample), y = factor(paste0(MetName, " (", MSI, ")_", IonMode, AlignmentID), levels = comp_order), 
                                  fill = log10value)) + 
  geom_tile(color = "black") + scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey90") + ylab("") + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.background = element_blank(), panel.grid = element_line(color = "black"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title.y = element_blank()) + facet_grid(.~paste0(Solution, " ", Conc), scales = "free_x") +
  xlab("Sample") + scale_x_discrete(labels = c("", "Amino acids 1x", "", "", "Amino acids 2x", "", "", "Arginine 1x", "","", "Arginine working stock (6x)", "", "", "Uracil 1x", "", "", "Uracil working stock (9x)", "")) #+ geom_point(data = data_all[Intended==1 & Duplicate == 0& !grepl("BK", Sample)], size = 1.5) + 
  #geom_point(data = data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0 & HasComparisonInTimeSeries==1& Present1 == 1& !grepl("BK", Sample)], color = "white", size = 0.7)#+ facet_wrap(~IonMode, scales = "free_y") 

comp_order2 <- comp_order[comp_order %in% data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0, paste0(MetName, " (", MSI, ")_", IonMode, AlignmentID)]]
prev_info_panel <- ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0  & !grepl("BK", Sample) & HasComparisonInTimeSeries==1], 
                          aes(x = 1, y = factor(paste0(MetName, " (", MSI, ")_", IonMode, AlignmentID), levels = comp_order))) + geom_point(color = "darkred") + 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank()) + scale_y_discrete(limits = comp_order2)
intended_panel <- ggplot(data_all[!is.na(MSI) & !grepl("[_\\-][Dd][35789]", MetName) & Duplicate == 0  & !grepl("BK", Sample) & Intended==1], 
                         aes(x = 1, y = factor(paste0(MetName, " (", MSI, ")_", IonMode, AlignmentID), levels = comp_order))) + geom_point( color = "darkorange") +
  scale_y_discrete(limits = comp_order2, labels = gsub("_.*$", "", comp_order2)) + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title = element_blank())
library(patchwork)
combined_plot <- intended_panel+prev_info_panel+simple_heatmap_plot + plot_layout(nrow = 1, widths = c(1, 1, 32))
ggsave(combined_plot, file = "solutions_simpleHeatmap.pdf", width = 11, height = 7)

ggplot(data_all[Solution == "Arginine" & Present1 == 1 & !is.na(MSI) & !grepl("iSTD", MetName) & !grepl("[_\\-][Dd][35789]", MetName)], aes(x = Conc, y = log10value, color = MetName)) + geom_point() + facet_wrap(~paste0(MetName, AlignmentID, IonMode)) + guides(color ="none")
ggplot(data_all[Solution == "Arginine" & Present1 == 1 & !is.na(MSI) & !grepl("iSTD", MetName) & !grepl("[_\\-][Dd][35789]", MetName)], aes(x = Conc, y = value, color = MetName)) + geom_point() + facet_wrap(~paste0(MetName, AlignmentID, IonMode), scales = "free_y") + guides(color ="none")

data_all[,PresentInArginine:=ifelse(IDUniq %in% data_all[Solution == "Arginine" & Present1 == 1 & !grepl("iSTD", MetName) & !grepl("[_\\-][Dd][35789]", MetName), IDUniq], 1, 0)]
data_all[,PresentInAAStock:=ifelse(IDUniq %in% data_all[Solution == "Amino acids" & Present1 == 1 & !grepl("iSTD", MetName) & !grepl("[_\\-][Dd][35789]", MetName), IDUniq], 1, 0)]
comp_data <- unique(data_all[!grepl("iSTD", MetName) & !grepl("[_\\-][Dd][35789]", MetName),
                             list(IDUniq, MetName, MSI, HasComparisonInTimeSeries, Intended, PresentInArginine, PresentInAAStock)])
comp_data[,fisher.test(table(PresentInArginine, HasComparisonInTimeSeries))]
comp_data[,fisher.test(table(PresentInAAStock & PresentInArginine==0, HasComparisonInTimeSeries))]
comp_data[,fisher.test(table(PresentInAAStock & PresentInArginine==1, HasComparisonInTimeSeries))]


### Comparison of retention times
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
time_course2 = fread("../processedDatasets/timeCourse2_allData.txt")

ms2_dat <- data.table(read_xlsx("../DataWithMS2/Neg Combined for Cecilia.xlsx", skip = 4, sheet = "Timecourse"))
ms2_dat <- ms2_dat[,list(`Alignment ID`, `Reference RT`, `Reference m/z`, `Adduct type`, `Formula`, `MS/MS spectrum`, BK)]
ms2_dat2 <- data.table(read_xlsx("../DataWithMS2/Pos Combined for Cecilia.xlsx", skip = 4, sheet = "TimeCourse"))
ms2_dat2 <- ms2_dat2[,list(`Alignment ID`, `Reference RT`,`Reference m/z`,  `Adduct type`, `Formula`, `MS/MS spectrum`, BK)]
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

annot_features <- unique(time_course2[,list(IonMode, AlignmentID, IDUniq, MSI, AvgRt, AvgMZ, MetName, INCHIKEY, `Reference RT`, `Reference m/z`, Adduct, Duplicate)])
annot_features[Duplicate == 0][`Reference RT` != "null",hist(abs(AvgRt-as.numeric(`Reference RT`)), breaks = 50)]
compare_rts <- annot_features[ MSI != ""]
compare_rts[,RTDiff:=abs(AvgRt-as.numeric(`Reference RT`))]
compare_rts[,MZDiff:=abs(AvgMZ - as.numeric(`Reference m/z`))]
compare_rts <- compare_rts[MSI %in% c(1,2,3,4)]
fwrite(compare_rts, file = "figure_source_data/figureS12_timeSeries_compareAnnots.tsv", quote=F, row.names = F, sep = "\t")

ggplot(compare_rts, aes(x = MZDiff)) + facet_wrap(~MSI) + geom_histogram(bins = 26) #+ xlim(0, 0.25)
rt_compare <- ggplot(compare_rts[MSI %in% c(1,2)], aes(x = RTDiff)) + facet_wrap(~paste0("MSI ", MSI)) + geom_histogram(binwidth = 0.02) +
  xlab("Difference between observed and\nlibrary standard retention time (minutes)") + ylab("Number of features") +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 11))
mz_compare <- ggplot(compare_rts[MSI %in% c(1,2, 3, 4)], aes(x = MZDiff)) + facet_wrap(~paste0("MSI ", MSI)) + geom_histogram(binwidth = 0.0002) +
  xlab("Difference between observed and\nlibrary standard MS1 m/z") + ylab("Number of features") +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 11))
library(patchwork)
final_plot_comp <- mz_compare+rt_compare + plot_annotation(tag_levels = "A")
ggsave(final_plot_comp, file = "fig_sX_libraryRefDiffs.pdf", width = 7.5, height = 3.8)
ggsave(final_plot_comp, file = "fig_sX_libraryRefDiffs.png", width = 7.5, height = 3.8)
compare_rts[MSI %in% c(1,2),table(MZDiff < 0.002 & RTDiff < 0.25)]
compare_rts[MSI %in% c(1,2) & (MZDiff > 0.002 | RTDiff > 0.25)]
compare_rts[(MZDiff > 0.003 | RTDiff > 0.25)]
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

time_series_orig <- read_xl()
