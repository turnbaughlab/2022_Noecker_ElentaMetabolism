library(data.table)
library(readxl)
library(ggplot2)
library(cowplot)
library(igraph)
#library(MSnbase)
library(MsCoreUtils)
library(Spectra)
### Processing for negative and positive data
# 
datadir <- "../../DataWithMS2/"
data_all_neg <- "Neg Combined for Cecilia.xlsx"
data_all_pos <- "Pos Combined for Cecilia.xlsx"
additional_file <- "Pilot spectra.xlsx"
theme_set(theme_cowplot())

#### Read in all datasets
datasets_pos <- c("BHI", "GL", "TimeCourse", "Intestinal Contents", "Serum", "HuMix", "SIL PreTest")
all_pos_data = rbindlist(lapply(datasets_pos, function(x){
  foo <- data.table(read_xlsx(paste0(datadir, data_all_pos), sheet = x, skip = 4))
  foo2 <- data.table(read_xlsx(paste0(datadir, data_all_pos), sheet = x, skip = 0, n_max = 4, trim_ws = F, col_names = F))
  num_id_vars <- which(names(foo) == "MS/MS spectrum")
  idvars <- names(foo)[1:num_id_vars]
  final_column <- which(foo2[4,] == "Average")[1] + num_id_vars - 2
  foo <- foo[,1:final_column, with = F]
  foo <- melt(foo, id.var = idvars, variable.name = "Sample", variable.factor = F)
  foo[,Dataset:=x]
  #foo = foo[!grepl("...", Sample, fixed = T)]
  foo[,IDUniq:=paste0(Dataset, "_", `Alignment ID`, "_Pos")]
}), fill = T)

pilot_pos <- data.table(read_xlsx(paste0(datadir, additional_file), sheet = 2, skip = 4))
foo2 <- data.table(read_xlsx(paste0(datadir, additional_file), sheet = 2, skip = 0, n_max = 4, trim_ws = F, col_names = F))
num_id_vars <- which(names(pilot_pos) == "MS/MS spectrum")
idvars <- names(pilot_pos)[1:num_id_vars]
final_column <- which(foo2[4,] == "Average")[1] + num_id_vars - 2
pilot_pos <- pilot_pos[,1:final_column, with = F]
pilot_pos <- melt(pilot_pos, id.var = idvars, variable.name = "Sample", variable.factor = F)
pilot_pos[,Dataset:="Pilot"]
#foo = foo[!grepl("...", Sample, fixed = T)]
pilot_pos[,IDUniq:=paste0(Dataset, "_", `Alignment ID`, "_Pos")]
rm(foo2)
all_pos_data <- rbind(all_pos_data, pilot_pos, fill = T)

all_pos_data[,AvgMZ:=lapply(`Average Mz`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
all_pos_data[,AvgRt:=lapply(`Average Rt(min)`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]

### Get feature data
feature_data <- unique(all_pos_data[,list(Dataset, IDUniq, MSI, `Metabolite name`, INCHIKEY, `Average Mz`, `Average Rt(min)`, `MS/MS spectrum`, `Adduct type`)])
feature_data[,AvgMZ:=lapply(`Average Mz`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
feature_data[,AvgRt:=lapply(`Average Rt(min)`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
features_all <- feature_data[,sort(unique(IDUniq))]

spectra_data <- feature_data[!is.na(`MS/MS spectrum`)]
feature_data[,`MS/MS spectrum`:=NULL]

### Make tables of MS2 spectra
bad_spectra <-  all_pos_data[grepl("bk", `Spectrum reference file name`, ignore.case = T) & `RT matched`==T, unique(IDUniq)]
spectra_data <- spectra_data[!IDUniq %in% bad_spectra]
spectra_all <- list()
spec_mats <- list()
#spectra_all <- fread(text = spectra_data[,`MS/MS spectrum`], sep = " ", sep2 = ":")
for(j in 1:nrow(spectra_data)){
  spec <- spectra_data[j]
  peaks <- data.table(tstrsplit(spec[,`MS/MS spectrum`], split = " |_"))
  peaks[,c("MZ","Height"):=tstrsplit(V1, split = ":")]
  peaks[,V1:=NULL]
  peaks[,MZ:=as.numeric(MZ)]
  peaks[,Height:=as.numeric(Height)]
  #peaks <- peaks[Height > 2000]
  for(k in 1:(ncol(spec)-1)){
    peaks[,names(spec)[k]:=spec[,get(names(spec)[k])]]
  }
  spectra_all[[j]] <- peaks[order(MZ)]
  spec_mats[[j]] <- as.matrix(peaks[order(MZ), list(MZ, Height)])
}
names(spec_mats) = spectra_data[,IDUniq]
spectra_all <- rbindlist(spectra_all)

## Functions for efficient comparison
compare_vals <- function(val1, val2){
  foo <- expand.grid(val1, val2)
  foo$diff <- abs(foo$Var1 - foo$Var2)
  return(min(foo$diff, na.rm = T))
}

compare_vals_character <- function(val1, val2){
  return(length(intersect(as.character(val1), as.character(val2))) > 0)
}

compare_ms2 <- function(id1, id2){
  mat1 <- spec_mats[[id1]]
  mat2 <- spec_mats[[id2]]
  foo <- joinPeaks(mat1, mat2, tolerance = 0.01)
  return(ndotproduct(foo$x, foo$y))
}

### Get all close data pairs, calculate MS2 similairty between them
## Do this in chunks for memory
mz_diff_cutoff <- 0.01
rt_diff_cutoff <- 1
close_data_all <- data.table()
for(j in 1:floor(length(features_all)/1000)){
  mz_rt_compare1 <- CJ(features_all, features_all[((j-1)*1000+1):(j*1000)], unique = T)[features_all < V2]
  setnames(mz_rt_compare1, c("Var1", "Var2"))
  mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var1", by.y = "IDUniq")
  mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var2", by.y = "IDUniq")
  mz_rt_compare1[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
  close_pairs <- mz_rt_compare1[MZDiff < mz_diff_cutoff]
  rm(mz_rt_compare1)
  close_pairs[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
  close_pairs <- close_pairs[RtDiff < rt_diff_cutoff]
  close_data_all <- rbind(close_data_all, close_pairs, fill = T)
  saveRDS(close_data_all, file = paste0("posDataClose", j, ".rds"))
}
mz_rt_compare1 <- CJ(features_all, features_all[(floor(length(features_all)/1000)*1000+1):length(features_all)], unique = T)[features_all < V2]
setnames(mz_rt_compare1, c("Var1", "Var2"))
mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var1", by.y = "IDUniq")
mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var2", by.y = "IDUniq")
mz_rt_compare1[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
close_pairs <- mz_rt_compare1[MZDiff < mz_diff_cutoff]
rm(mz_rt_compare1)
close_pairs[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
close_pairs <- close_pairs[RtDiff < rt_diff_cutoff]
close_data_all <- rbind(close_data_all, close_pairs, fill = T)
close_data_all[,HasMS2.x:=ifelse(Var1 %in% spectra_all[,IDUniq], 1, 0)]
close_data_all[,HasMS2.y:=ifelse(Var2 %in% spectra_all[,IDUniq], 1, 0)]

saveRDS(close_data_all, file = "posDataAllClose.rds")

# Random sample of not close pairs - FDR
not_close <- CJ(sample(features_all, 200), sample(features_all, 200))[V1 < V2]
setnames(not_close, c("Var1", "Var2"))
not_close <- merge(not_close, feature_data, by.x = "Var1", by.y = "IDUniq")
not_close <- merge(not_close, feature_data, by.x = "Var2", by.y = "IDUniq")
not_close[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
not_close[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
not_close[,HasMS2.x:=ifelse(Var1 %in% spectra_all[,IDUniq], 1, 0)]
not_close[,HasMS2.y:=ifelse(Var2 %in% spectra_all[,IDUniq], 1, 0)]
not_close <- not_close[MZDiff > mz_diff_cutoff | RtDiff > rt_diff_cutoff]
saveRDS(not_close, file = "posDataAllSampleNotClose.rds")

mz_diff_cutoff <- 0.007
rt_diff_cutoff <- 0.5

closer_mz_cutoff <- 0.001
closer_rt_cutoff <- 0.2

close_data_all <- close_data_all[MZDiff < mz_diff_cutoff & RtDiff < rt_diff_cutoff]
close_data_all[HasMS2.x & HasMS2.y, MZ_sim:=mapply(compare_ms2, Var1, Var2)]
not_close[HasMS2.x & HasMS2.y, MZ_sim:=mapply(compare_ms2, Var1, Var2)]

sim_cutoff <- not_close[,quantile(MZ_sim, 0.995, na.rm = T)] #0.205

similar_features <- close_data_all[MZ_sim > sim_cutoff | ((HasMS2.x ==0| HasMS2.y==0) & MZDiff < closer_mz_cutoff & RtDiff < closer_rt_cutoff)]

close_data_all[,CloseMZRt:=T]
not_close[,CloseMZRt:=F]
data_all <- rbind(close_data_all[,list(Var1, Var2, MZ_sim, CloseMZRt)], not_close[,list(Var1, Var2, MZ_sim, CloseMZRt)])
similarity_dist_plot <- ggplot(data_all[!is.na(MZ_sim)], aes(x = MZ_sim, color = CloseMZRt, group = CloseMZRt)) + geom_density(na.rm = T) + ylim(0, 12) + xlab("MS2 cosine similarity")
#similarity_dist_plot
ggsave(similarity_dist_plot, file = "ms2_similarity_distribution_pos.pdf", width = 5, height = 4.3)

## Add name info to similar features table, make sure they have the same adduct
met_names_pos <- unique(all_pos_data[,list(Dataset, IDUniq, MSI, `Metabolite name`, INCHIKEY, `Adduct type`)])
met_names_pos[,`Metabolite name`:=gsub("_\\[M\\+H.*$", "", `Metabolite name`)]

similar_features <- merge(similar_features, met_names_pos[,list(IDUniq, INCHIKEY, `Adduct type`)], by.x = "Var1", by.y = "IDUniq", all.x = T)
similar_features <- merge(similar_features, met_names_pos[,list(IDUniq, INCHIKEY, `Adduct type`)], by.x = "Var2", by.y = "IDUniq", all.x = T)

similar_features[,`Adduct type.x`:=sapply(`Adduct type.x`, function(x) { as.character(strsplit(x, split = "_")[[1]])})]
similar_features[,`Adduct type.y`:=sapply(`Adduct type.y`, function(x) { as.character(strsplit(x, split = "_")[[1]])})]
similar_features[,SharedAdduct:=mapply(compare_vals_character, `Adduct type.x`, `Adduct type.y`)]
similar_features <- similar_features[SharedAdduct == T]

overall_pair_counts <- ggplot(similar_features, aes(x = Dataset.y, y = Dataset.x)) + geom_bin2d() + scale_fill_distiller(palette = "Blues", direction = 1)
ggsave(overall_pair_counts, file = "numSimilarFeaturesPos.pdf", width = 5, height = 4)

### Merge all linked components into mega-features
feature_merge_graph <- graph_from_edgelist(as.matrix(similar_features[,list(Var1, Var2)]), directed = F)
new_features <- components(feature_merge_graph)
feature_assignment <- data.table(V(feature_merge_graph)$name, new_features$membership)
feature_assignment[,NewID:=paste0("Merged_", V2)]
similar_features <- similar_features[order(Var1)]
feature_assignment[,AvgMZ:=sapply(V1, function(x){
  unique(unlist(similar_features[Var1 == x, AvgMZ.x]))
})]
feature_assignment[sapply(AvgMZ, length)==0, AvgMZ:=sapply(V1, function(x){
  unique(unlist(similar_features[Var2 == x, AvgMZ.y]))
})]
feature_assignment[,AvgRt:=sapply(V1, function(x){
  unique(unlist(similar_features[Var1 == x, AvgRt.x]))
})]
feature_assignment[sapply(AvgRt, length)==0, AvgRt:=sapply(V1, function(x){
  unique(unlist(similar_features[Var2 == x, AvgRt.y]))
})]


feature_assignment[,Dataset:=gsub("_.*", "", V1)]
feature_assignment <- merge(feature_assignment, met_names_pos, by.x = c("Dataset", "V1"), by.y = c("Dataset", "IDUniq"), all.x = T)
feature_assignment[,`Metabolite name`:=gsub("_\\[M\\+H.*$", "", `Metabolite name`)]
feature_assignment[!is.na(MSI),length(unique(INCHIKEY)), by=NewID][,table(V1)]
feature_assignment[MSI==1,length(unique(INCHIKEY)), by=NewID][,table(V1)]

similar_features <- merge(similar_features, feature_assignment[,list(NewID, V1)], by.x = "Var1", by.y = "V1", all.x = T)
all_similar_pairs_plot <- ggplot(similar_features, aes(x = MZDiff, y = RtDiff, color = MZ_sim)) + geom_point() + scale_color_distiller(palette = "YlGnBu")
ggsave(all_similar_pairs_plot, file = "allSimilarPairwiseFeatures.pdf")

## Get merged names
similar_features[INCHIKEY.x != INCHIKEY.y & MSI.x==1 & MSI.y==1 & INCHIKEY.x != "null" & INCHIKEY.y != "null"] #Only 61
feature_assignment[,`Metabolite name`:=gsub(" \\[.*", "", `Metabolite name`)]
feature_assignment[,`Metabolite name`:=gsub("__.*", "", `Metabolite name`)]
feature_assignment[!is.na(MSI) & MSI != "NA",NameMSI:=paste0(`Metabolite name`, MSI)]
feature_assignment[,MergedName:=paste0(unique(NameMSI[!`Metabolite name` %in% c("Unknown", "unknown") & !is.na(`Metabolite name`) & `Metabolite name` != "NA"]), collapse = "_"), by = NewID]
feature_assignment[,MergedName:=gsub("_NA", "", MergedName)]
feature_assignment[,MergedName:=gsub("NA_", "", MergedName)]
feature_assignment[MergedName==""|MergedName=="NA", MergedName:="Unknown"]

saveRDS(feature_assignment, file = "similar_features_merge_pos_final.rds")
fwrite(feature_assignment, file = "similar_features_merge_pos.tsv", quote=F, row.names = F, sep = "\t")
fwrite(similar_features, file = "similar_features_pairwise_similarity_pos.tsv", quote=F, row.names = F, sep = "\t")

# Make version of full dataset with merged IDs
all_pos_data <- merge(all_pos_data, feature_assignment[,list(V1, NewID)], by.x = "IDUniq", by.y = "V1", all.x = T)
saveRDS(all_pos_data, file = "allPosDataWithMergedIDs.rds")

##### Repeat everything with negative data
rm(list = ls())
datadir <- "../DataWithMS2/"
data_all_neg <- "Neg Combined for Cecilia.xlsx"
data_all_pos <- "Pos Combined for Cecilia.xlsx"
additional_file <- "Pilot spectra.xlsx"

#### Read in negative data
datasets_neg <- c("BHI", "GL", "Timecourse", "Intestinal contents", "Serum", "HuMix", "SIL pretest")
all_neg_data = rbindlist(lapply(datasets_neg, function(x){
  foo <- data.table(read_xlsx(paste0(datadir, data_all_neg), sheet = x, skip = 4))
  foo2 <- data.table(read_xlsx(paste0(datadir, data_all_neg), sheet = x, skip = 0, n_max = 4, trim_ws = F, col_names = F))
  foo2_firstcol <- which(foo2[1,] == "Class")
  foo2 <- foo2[,foo2_firstcol:ncol(foo2), with = F]
  num_id_vars <- which(names(foo) == "MS/MS spectrum")
  idvars <- names(foo)[1:num_id_vars]
  final_column <- which(foo2[4,] == "Average")[1] + num_id_vars - 2
  foo <- foo[,1:final_column, with = F]
  foo <- melt(foo, id.var = idvars, variable.name = "Sample", variable.factor = F)
  foo[,Dataset:=x]
  #foo = foo[!grepl("...", Sample, fixed = T)]
  foo[,IDUniq:=paste0(Dataset, "_", `Alignment ID`, "_Neg")]
}), fill = T)

pilot_neg <- data.table(read_xlsx(paste0(datadir, additional_file), sheet = 1))
num_id_vars <- which(names(pilot_neg) == "MS/MS spectrum")
idvars <- names(pilot_neg)[1:num_id_vars]
pilot_neg <- melt(pilot_neg, id.var = idvars, variable.name = "Sample", variable.factor = F)
pilot_neg[,Dataset:="Pilot"]
pilot_neg[,IDUniq:=paste0(Dataset, "_", `Alignment ID`, "_Neg")]
all_neg_data <- rbind(all_neg_data, pilot_neg, fill = T)

all_neg_data[,AvgMZ:=lapply(`Average Mz`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
all_neg_data[,AvgRt:=lapply(`Average Rt(min)`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]

feature_data <- unique(all_neg_data[,list(Dataset, IDUniq, MSI, `Metabolite name`, INCHIKEY, `Average Mz`, `Average Rt(min)`, `MS/MS spectrum`, `Adduct type`)])
feature_data[,AvgMZ:=lapply(`Average Mz`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
feature_data[,AvgRt:=lapply(`Average Rt(min)`, function(x) { return(as.numeric(strsplit(x, split ="_")[[1]]))} )]
features_all <- feature_data[,sort(unique(IDUniq))]

##### Make table of MS2 spectra
spectra_data <- feature_data[!is.na(`MS/MS spectrum`)]
feature_data[,`MS/MS spectrum`:=NULL]

bad_spectra <-  all_neg_data[grepl("bk", `Spectrum reference file name`, ignore.case = T) & `RT matched`==T, unique(IDUniq)] 
spectra_data <- spectra_data[!IDUniq %in% bad_spectra]
spectra_all <- list()
spec_mats <- list()
#spectra_all <- fread(text = spectra_data[,`MS/MS spectrum`], sep = " ", sep2 = ":")
for(j in 1:nrow(spectra_data)){
  spec <- spectra_data[j]
  peaks <- data.table(tstrsplit(spec[,`MS/MS spectrum`], split = " |_"))
  peaks[,c("MZ","Height"):=tstrsplit(V1, split = ":")]
  peaks[,V1:=NULL]
  peaks[,MZ:=as.numeric(MZ)]
  peaks[,Height:=as.numeric(Height)]
  #peaks <- peaks[Height > 2000]
  for(k in 1:(ncol(spec)-1)){
    peaks[,names(spec)[k]:=spec[,get(names(spec)[k])]]
  }
  spectra_all[[j]] <- peaks[order(MZ)]
  spec_mats[[j]] <- as.matrix(peaks[order(MZ), list(MZ, Height)])
}
names(spec_mats) = spectra_data[,IDUniq]
spectra_all <- rbindlist(spectra_all)

compare_vals <- function(val1, val2){
  foo <- expand.grid(val1, val2)
  foo$diff <- abs(foo$Var1 - foo$Var2)
  return(min(foo$diff, na.rm = T))
}

compare_vals_character <- function(val1, val2){
  return(length(intersect(as.character(val1), as.character(val2))) > 0)
}

compare_ms2 <- function(id1, id2){
  mat1 <- spec_mats[[id1]]
  mat2 <- spec_mats[[id2]]
  foo <- joinPeaks(mat1, mat2, tolerance = 0.01)
  return(ndotproduct(foo$x, foo$y))
}

##################### Run through pairwise comparison for negative features
mz_diff_cutoff <- 0.01
rt_diff_cutoff <- 1
close_data_all <- data.table()
for(j in 1:floor(length(features_all)/1000)){
  mz_rt_compare1 <- CJ(features_all, features_all[((j-1)*1000+1):(j*1000)], unique = T)[features_all < V2]
  setnames(mz_rt_compare1, c("Var1", "Var2"))
  mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var1", by.y = "IDUniq")
  mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var2", by.y = "IDUniq")
  mz_rt_compare1[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
  close_pairs <- mz_rt_compare1[MZDiff < mz_diff_cutoff]
  rm(mz_rt_compare1)
  close_pairs[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
  close_pairs <- close_pairs[RtDiff < rt_diff_cutoff]
  close_data_all <- rbind(close_data_all, close_pairs, fill = T)
  saveRDS(close_data_all, file = paste0("negDataClose", j, ".rds"))
}
mz_rt_compare1 <- CJ(features_all, features_all[(floor(length(features_all)/1000)*1000+1):length(features_all)], unique = T)[features_all < V2]
setnames(mz_rt_compare1, c("Var1", "Var2"))
mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var1", by.y = "IDUniq")
mz_rt_compare1 <- merge(mz_rt_compare1, feature_data, by.x = "Var2", by.y = "IDUniq")
mz_rt_compare1[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
close_pairs <- mz_rt_compare1[MZDiff < mz_diff_cutoff]
rm(mz_rt_compare1)
close_pairs[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
close_pairs <- close_pairs[RtDiff < rt_diff_cutoff]
close_data_all <- rbind(close_data_all, close_pairs, fill = T)
close_data_all[,HasMS2.x:=ifelse(Var1 %in% spectra_all[,IDUniq], 1, 0)]
close_data_all[,HasMS2.y:=ifelse(Var2 %in% spectra_all[,IDUniq], 1, 0)]

saveRDS(close_data_all, file = "negDataAllClose.rds")

# Random sample of not close pairs
not_close <- CJ(sample(features_all, 200), sample(features_all, 200))[V1 < V2]
setnames(not_close, c("Var1", "Var2"))
not_close <- merge(not_close, feature_data, by.x = "Var1", by.y = "IDUniq")
not_close <- merge(not_close, feature_data, by.x = "Var2", by.y = "IDUniq")
not_close[,MZDiff:=mapply(compare_vals, AvgMZ.x, AvgMZ.y)]
not_close[,RtDiff:=mapply(compare_vals, AvgRt.x, AvgRt.y)]
not_close[,HasMS2.x:=ifelse(Var1 %in% spectra_all[,IDUniq], 1, 0)]
not_close[,HasMS2.y:=ifelse(Var2 %in% spectra_all[,IDUniq], 1, 0)]
not_close <- not_close[MZDiff > mz_diff_cutoff | RtDiff > rt_diff_cutoff]
saveRDS(not_close, file = "negDataAllSampleNotClose.rds")
############################################

close_data_all[,Var1:=gsub("Pos", "Neg" ,Var1)]
close_data_all[,Var2:=gsub("Pos", "Neg" ,Var2)]
not_close[,Var1:=gsub("Pos", "Neg" ,Var1)]
not_close[,Var2:=gsub("Pos", "Neg" ,Var2)]

mz_diff_cutoff <- 0.007
rt_diff_cutoff <- 0.5

closer_mz_cutoff <- 0.001
closer_rt_cutoff <- 0.5
closer_rt_cutoff <- 0.2

close_data_all <- close_data_all[MZDiff < mz_diff_cutoff & RtDiff < rt_diff_cutoff]
close_data_all[HasMS2.x & HasMS2.y, MZ_sim:=mapply(compare_ms2, Var1, Var2)]
not_close[HasMS2.x & HasMS2.y, MZ_sim:=mapply(compare_ms2, Var1, Var2)]

sim_cutoff <- not_close[,quantile(MZ_sim, 0.995, na.rm = T)] #0.205

similar_features <- close_data_all[MZ_sim > sim_cutoff | ((HasMS2.x ==0| HasMS2.y==0) & MZDiff < closer_mz_cutoff & RtDiff < closer_rt_cutoff)]

close_data_all[,CloseMZRt:=T]
not_close[,CloseMZRt:=F]
data_all <- rbind(close_data_all[,list(Var1, Var2, MZ_sim, CloseMZRt)], not_close[,list(Var1, Var2, MZ_sim, CloseMZRt)])
similarity_dist_plot <- ggplot(data_all[!is.na(MZ_sim)], aes(x = MZ_sim, color = CloseMZRt, group = CloseMZRt)) + geom_density(na.rm = T) + ylim(0, 12) + xlab("MS2 cosine similarity")
ggsave(similarity_dist_plot, file = "ms2_similarity_distribution_neg.pdf", width = 5, height = 4.3)

similar_features[`Metabolite name.x` != `Metabolite name.y` & !grepl("known", `Metabolite name.x`) & !grepl("known", `Metabolite name.y`)][,list(`Metabolite name.x`, `Metabolite name.y`, MSI.x, MSI.y, HasMS2.x, HasMS2.y)]

met_names_neg <- unique(all_neg_data[,list(Dataset, IDUniq, MSI, `Metabolite name`, INCHIKEY, `Adduct type`)])
met_names_neg[,`Metabolite name`:=gsub("_\\[M\\+H.*$", "", `Metabolite name`)]

similar_features <- merge(similar_features, met_names_neg[,list(IDUniq, INCHIKEY, `Adduct type`)], by.x = "Var1", by.y = "IDUniq", all.x = T)
similar_features <- merge(similar_features, met_names_neg[,list(IDUniq, INCHIKEY, `Adduct type`)], by.x = "Var2", by.y = "IDUniq", all.x = T)

similar_features[,`Adduct type.x`:=sapply(`Adduct type.x`, function(x) { as.character(strsplit(x, split = "_")[[1]])})]
similar_features[,`Adduct type.y`:=sapply(`Adduct type.y`, function(x) { as.character(strsplit(x, split = "_")[[1]])})]
similar_features[,SharedAdduct:=mapply(compare_vals_character, `Adduct type.x`, `Adduct type.y`)]
similar_features <- similar_features[SharedAdduct == T]
overall_pair_counts <- ggplot(similar_features, aes(x = Dataset.y, y = Dataset.x)) + geom_bin2d() + scale_fill_distiller(palette = "Blues", direction = 1)
ggsave(overall_pair_counts, file = "numSimilarFeaturesNeg.pdf", width = 5, height = 4)

#### Merge connected features into new mega features
feature_merge_graph <- graph_from_edgelist(as.matrix(similar_features[,list(Var1, Var2)]), directed = F)
new_features <- components(feature_merge_graph)
feature_assignment <- data.table(V(feature_merge_graph)$name, new_features$membership)
feature_assignment[,NewID:=paste0("Merged_", V2)]
similar_features <- similar_features[order(Var1)]

feature_assignment[,Dataset:=gsub("_.*", "", V1)]
feature_assignment <- merge(feature_assignment, met_names_neg, by.x = c("Dataset", "V1"), by.y = c("Dataset", "IDUniq"), all.x = T)
feature_assignment[,`Metabolite name`:=gsub("_\\[M\\+H.*$", "", `Metabolite name`)]

similar_features <- merge(similar_features, feature_assignment[,list(NewID, V1)], by.x = "Var1", by.y = "V1", all.x = T)
all_similar_pairs_plot <- ggplot(similar_features, aes(x = MZDiff, y = RtDiff, color = MZ_sim)) + geom_point() + scale_color_distiller(palette = "YlGnBu")
ggsave(all_similar_pairs_plot, file = "allSimilarPairwiseFeaturesNeg.pdf")

similar_features[INCHIKEY.x != INCHIKEY.y & MSI.x==1 & MSI.y==1 & INCHIKEY.x != "null" & INCHIKEY.y != "null"] #Only 61
feature_assignment[,`Metabolite name`:=gsub(" \\[.*", "", `Metabolite name`)]
feature_assignment[,`Metabolite name`:=gsub("__.*", "", `Metabolite name`)]
feature_assignment[!is.na(MSI) & MSI != "NA",NameMSI:=paste0(`Metabolite name`, MSI)]
feature_assignment[,MergedName:=paste0(unique(NameMSI[!`Metabolite name` %in% c("Unknown", "unknown") & !is.na(`Metabolite name`) & `Metabolite name` != "NA"]), collapse = "_"), by = NewID]
feature_assignment[,MergedName:=gsub("_NA", "", MergedName)]
feature_assignment[,MergedName:=gsub("NA_", "", MergedName)]
feature_assignment[MergedName==""|MergedName=="NA", MergedName:="Unknown"]

saveRDS(feature_assignment, file = "similar_features_merge_neg_final.rds")
fwrite(feature_assignment, file = "similar_features_merge_neg.tsv", quote=F, row.names = F, sep = "\t")
fwrite(similar_features, file = "similar_features_pairwise_similarity_neg.tsv", quote=F, row.names = F, sep = "\t")

all_neg_data <- merge(all_neg_data, feature_assignment[,list(V1, NewID)], by.x = "IDUniq", by.y = "V1", all.x = T)
saveRDS(all_neg_data, file = "allNegDataWithMergedIDs.rds")
