##### Make Metabolomics Workbench tables
library(data.table)

## Time course
time_course_processed <- fread("processedDatasets/timecourse_all_data.csv")
metadata <- unique(time_course_processed[,list(IonMode, Sample, Strain, AcConc, Replicate, Time, SampleType)])
raw_files <- list.files("MetabolomicsWorkbench/TimeCourse/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/TimeCourse/raw/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^[A-Z][a-z]+_", "", Sample)]
files_tab[!SampleOrig %in% metadata[,Sample]]
files_tab[,IonMode:=gsub("_.*$", "", Sample)]
setnames(files_tab, "SampleOrig", "Sample name")
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("Sample", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Pos", "Neg")]

study_final <- dcast(files_tab, `Sample name`+Strain+AcConc+Replicate+Time+SampleType~IonMode, value.var = c("RawFile", "MzMLFile"))
study_final[,SubjectID:=gsub("_TP[0-9]", "", `Sample name`)]
setnames(study_final, "AcConc", "AcetateConcentration_mM")
fwrite(study_final, file = "MetabolomicsWorkbench/TimeCourse/metadata_formatted.tsv", sep = "\t")

named_metabolite_table_pos <- dcast(time_course_processed[IonMode=="Pos" & MetName != "unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
named_metabolite_table_neg <- dcast(time_course_processed[IonMode=="Neg" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/TimeCourse/named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/TimeCourse/named_data_pos.tsv", sep = "\t")

pos_names <- fread("MetabolomicsWorkbench/TimeCourse/refmet_results_pos.txt")
neg_names <- fread("MetabolomicsWorkbench/TimeCourse/refmet_results_neg.txt")
named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-", MetName, `Standardized name`)]
named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-", MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 2:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 2:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 7:(ncol(named_metabolite_table_neg)-2)), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 7:(ncol(named_metabolite_table_pos)-2)), with=F]


fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/TimeCourse/named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/TimeCourse/named_data_pos.tsv", sep = "\t")
fwrite(pos_metadata, file = "MetabolomicsWorkbench/TimeCourse/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(neg_metadata, file = "MetabolomicsWorkbench/TimeCourse/metabolite_metadata_neg.tsv", sep = "\t")


## Strains
strains_processed <- fread("processedDatasets/strain_edm_all_data.csv")
metadata <- unique(strains_processed[,list(SampleID, IonMode, Condition, Strain, Name2, PlateNum, Row, Column, Medium, ODNormValue2, IsCtrl, BadSampOD)])
raw_files <- list.files("MetabolomicsWorkbench/Strains/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/Strains/raw/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^GL_[A-Z]+_", "", Sample)]
metadata[,SampleID:=gsub("GL_", "", SampleID)]
files_tab[!SampleOrig %in% metadata[,SampleID]]
files_tab[,IonMode:=gsub("_.*$", "", gsub("GL_", "", Sample))]
setnames(files_tab, "SampleOrig", "Sample name")
files_tab[,IonMode:=ifelse(IonMode == "POS", "Positive", "Negative")]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("SampleID", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]

study_final <- dcast(files_tab, `Sample name`+Strain+Name2+Condition+PlateNum+Row+Column+ODNormValue2+IsCtrl+BadSampOD~IonMode, value.var = c("RawFile", "MzMLFile"))
study_final[,SubjectID:=gsub("_[1-3]$", "", `Sample name`)]
setnames(study_final, c("Name2", "ODNormValue2"), c("Strain name", "Normalized final OD600"))
fwrite(study_final, file = "MetabolomicsWorkbench/Strains/metadata_formatted.tsv", sep = "\t")

strains_processed[,SampleID:=gsub("GL_", "", SampleID)]
named_metabolite_table_pos <- dcast(strains_processed[IonMode=="Positive" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~SampleID, value.var = "value")
named_metabolite_table_neg <- dcast(strains_processed[IonMode=="Negative" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~SampleID, value.var = "value")
named_metabolite_table_pos[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
named_metabolite_table_neg[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/Strains/named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/Strains/named_data_pos.tsv", sep = "\t")

pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]

named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 2:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 2:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 7:(ncol(named_metabolite_table_neg)-2)), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 7:(ncol(named_metabolite_table_pos)-2)), with=F]


fwrite(unique(named_metabolite_table_neg), file = "MetabolomicsWorkbench/Strains/named_data_neg.tsv", sep = "\t")
fwrite(unique(named_metabolite_table_pos), file = "MetabolomicsWorkbench/Strains/named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/Strains/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/Strains/metabolite_metadata_neg.tsv", sep = "\t")



## Mice
mouse_processed <- fread("processedDatasets/mouse_intestinal_all_data.csv")
metadata <- unique(mouse_processed[,list(Sample, IonMode, Animal2, SampleType, Group, Replicate, Colonized, Cage)])
raw_files <- list.files("MetabolomicsWorkbench/GnotoData/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/GnotoData/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^[A-Z][a-z]+_", "", Sample)]
files_tab[!SampleOrig %in% metadata[,Sample]]
files_tab[,IonMode:=gsub("_.*$", "", Sample)]
setnames(files_tab, "SampleOrig", "Sample name")
mouse_processed[,IonMode:=ifelse(grepl("eg", IonMode), "Negative", "Positive")]
files_tab[,IonMode:=ifelse(IonMode == "Pos", "Positive", "Negative")]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("Sample", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]
setnames(files_tab, "Animal2", "AnimalID")

study_final <- dcast(files_tab, `Sample name`+AnimalID+SampleType+Group+Replicate+Colonized+Cage~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/GnotoData//metadata_formatted.tsv", sep = "\t")

mouse_processed[,table(Sample %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(mouse_processed[IonMode=="Positive" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
named_metabolite_table_neg <- dcast(mouse_processed[IonMode=="Negative" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
named_metabolite_table_pos[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
named_metabolite_table_neg[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/GnotoData/named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/GnotoData//named_data_pos.tsv", sep = "\t")

pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]

named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 2:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 2:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 7:(ncol(named_metabolite_table_neg)-2)), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 7:(ncol(named_metabolite_table_pos)-2)), with=F]


fwrite(unique(named_metabolite_table_neg), file = "MetabolomicsWorkbench/GnotoData/named_data_neg.tsv", sep = "\t")
fwrite(unique(named_metabolite_table_pos), file = "MetabolomicsWorkbench/GnotoData//named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/GnotoData/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/GnotoData/metabolite_metadata_neg.tsv", sep = "\t")



# Serum
mouse_serum <- fread("processedDatasets/mouse_serum_all_data.csv")
metadata <- unique(mouse_serum[,list(Sample, IonMode, Animal, SampleType, Group, Replicate,Cage)])
raw_files <- list.files("MetabolomicsWorkbench/GnotoData/raw/Serum", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/GnotoData/mzML/Serum", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^[A-Z][a-z]+_", "", Sample)]
files_tab[!SampleOrig %in% metadata[,Sample]]
files_tab[,IonMode:=gsub("_.*$", "", Sample)]
setnames(files_tab, "SampleOrig", "Sample name")
mouse_serum[,IonMode:=ifelse(grepl("eg", IonMode), "Negative", "Positive")]
files_tab[,IonMode:=ifelse(IonMode == "Pos", "Positive", "Negative")]
metadata[,IonMode:=ifelse(IonMode == "pos", "Positive", "Negative")]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("Sample", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]
setnames(files_tab, "Animal", "AnimalID")
files_tab <- unique(files_tab)

study_final <- dcast(files_tab, `Sample name`+AnimalID+SampleType+Group+Replicate+Cage~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/GnotoData//serum_metadata_formatted.tsv", sep = "\t")

mouse_serum[,table(Sample %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(mouse_serum[IonMode=="Positive" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
named_metabolite_table_neg <- dcast(mouse_serum[IonMode=="Negative" & MetName != "Unknown"], MetName+AlignmentID+MSI+AvgMZ+AvgRt+INCHIKEY~Sample, value.var = "value")
named_metabolite_table_pos[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
named_metabolite_table_neg[duplicated(MetName), MetName:=paste0(MetName, rank(AlignmentID)), by = MetName]
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/GnotoData/serum_named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/GnotoData//serum_named_data_pos.tsv", sep = "\t")

pos_names <- fread("MetabolomicsWorkbench/GnotoData/serum_refmet_results_pos.txt")
neg_names <- fread("MetabolomicsWorkbench/GnotoData/serum_refmet_results_neg.txt")
named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]

named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

named_metabolite_table_pos[duplicated(`metabolite name`), `metabolite name`:=paste0(`metabolite name`, rank(AlignmentID)), by = `metabolite name`]
named_metabolite_table_neg[duplicated(`metabolite name`), `metabolite name`:=paste0(`metabolite name`, rank(AlignmentID)), by = `metabolite name`]

pos_metadata <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 2:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 2:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(ncol(named_metabolite_table_neg), 7:(ncol(named_metabolite_table_neg)-2)), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 7:(ncol(named_metabolite_table_pos)-2)), with=F]


fwrite(unique(named_metabolite_table_neg), file = "MetabolomicsWorkbench/GnotoData/serum_named_data_neg.tsv", sep = "\t")
fwrite(unique(named_metabolite_table_pos), file = "MetabolomicsWorkbench/GnotoData//serum_named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/GnotoData/serum_metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/GnotoData/serum_metabolite_metadata_neg.tsv", sep = "\t")


## Targeted
targeted <- fread("processedDatasets/targeted_met_data.csv")
metadata <- unique(targeted[,list(MS_Sample, Location, Strain, AcetateConc, TimePoint, Replicate, well, OD600, Time, Row, Column)])
#metadata <- metadata[MS_Sample != ""]
raw_files <- list.files("MetabolomicsWorkbench/TargetedAcetate/raw/", pattern = ".d$")
mzml_files <- list.files("MetabolomicsWorkbench/TargetedAcetate/raw/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub("\\.d$", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
#files_tab[,Sample:=gsub("P1-", "S.", Sample)]
#files_tab[grepl("S\\.[A-H]0[1-9]", Sample),Sample:=gsub("0", "", Sample)]
setnames(files_tab, "Sample", "Sample name")
files_tab <- merge(files_tab, metadata, by.x = c("Sample name"), by.y = c("Location"), all.x = T)
files_tab <- unique(files_tab)

fwrite(files_tab, file = "MetabolomicsWorkbench/TargetedAcetate/metadata_formatted.tsv", sep = "\t")

targeted[,table(Location %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(targeted, Compound+LLOQ+ULOQ~Location, value.var = "mM")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/TargetedAcetate/named_data.tsv", sep = "\t")

pos_names <- fread("MetabolomicsWorkbench/TargetedAcetate/refmet_results.txt")
named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "Compound", by.y = "Input name", all.x = T)
named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), Compound, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 2:3),with=F]

named_metabolite_table_pos <- named_metabolite_table_pos[,c(ncol(named_metabolite_table_pos), 4:(ncol(named_metabolite_table_pos)-2)), with=F]


fwrite(unique(named_metabolite_table_pos), file = "MetabolomicsWorkbench/TargetedAcetate/named_data.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/TargetedAcetate//metabolite_metadata.tsv", sep = "\t")

## SIL Acetate extracellular
########## Extracellular
datadir <- "../SILAcetateTimeCourse/"
datadir2 <- "../"

low_level_data <- readRDS(paste0(datadir2, "LowLevelDataAllProcessedFixedAcetate.rds"))
top_data <- readRDS(paste0(datadir2, "SILprocessedTopLevelDataReplicatesFixedAcetate.rds"))
top_data_mean <- readRDS(paste0(datadir2, "SILprocessedTopLevelDataMeansFixedAcetate.rds"))

low_level_data[,Strain2:=factor(Strain, levels = c("c", "2243", "AB8", "Val"))]
levels(low_level_data$Strain2) <- c("Ctrl", "2243", "AB8n2", "Valencia")

bad_feat2 <- low_level_data[!grepl("L", AcetateGroup) & AreaFrac > 1e6 & NumLabeledC > 0, list(CompID, ExLabel, NumLabeledC, FeatID)][,length(NumLabeledC), by=list(CompID, ExLabel, FeatID)][V1 > 1]
bad_feats <- low_level_data[BadFeat==1,list(CompID, ExLabel, NumLabeledC)][,length(NumLabeledC), by=list(CompID, ExLabel)][V1 > 1]
bad_feats[,BadFeatAll:=1]

metadata <- unique(low_level_data[,list(Sample, Strain2, AcetateGroup, TimePoint, SampleFull, Sample2, F_ID, Replicate, Time)])
raw_files <- list.files("MetabolomicsWorkbench/SILAcetate_extracellular/Dec2020_SIL_Extracellular/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/SILAcetate_extracellular/Dec2020_SIL_Extracellular/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^[A-Z][a-z]+_", "", Sample)]
files_tab[!SampleOrig %in% metadata[,Sample2]]
files_tab[,IonMode:=gsub("_.*$", "", Sample)]
setnames(files_tab, "SampleOrig", "Sample name")
low_level_data[,IonMode:=ifelse(grepl("eg", Mode), "Negative", "Positive")]
metadata[,IonMode:=ifelse(grepl("Neg", Sample), "Negative", "Positive")]
files_tab[,IonMode:=ifelse(IonMode == "Pos", "Positive", "Negative")]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode", "Sample"), by.y = c("Sample2", "IonMode", "Sample"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]

study_final <- dcast(files_tab, `Sample name`+Strain2+AcetateGroup+TimePoint+Time+Replicate+F_ID~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/SILAcetate_extracellular///metadata_formatted.tsv", sep = "\t")

low_level_data[,table(Sample2 %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(low_level_data[IonMode=="Positive" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~Sample2, value.var = "AreaFrac")
named_metabolite_table_neg <- dcast(low_level_data[IonMode=="Negative" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~Sample2, value.var = "AreaFrac")
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILAcetate_extracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILAcetate_extracellular///named_data_pos.tsv", sep = "\t")

# pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
# neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
# named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]
# 
# named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(1:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(1:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(1, 7:(ncol(named_metabolite_table_neg))), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(1, 7:(ncol(named_metabolite_table_pos))), with=F]

fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILAcetate_extracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILAcetate_extracellular///named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/SILAcetate_extracellular/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/SILAcetate_extracellular/metabolite_metadata_neg.tsv", sep = "\t")


########## Intracellular
datadir <- "../SILAcetateTimeCourse/"
lower_level_intra_data <- readRDS(paste0(datadir, "IntracellularExtracts/LowLevelDataAllProcessed.rds"))
bad_feat3 <- lower_level_intra_data[AreaFrac > 1000 & NumLabeledC > 0, sum(Status == "Contaminating mass"), by = FeatID] ## 
lower_level_intra_data[,BadFeat3:=ifelse(FeatID %in% bad_feat3[V1 > 11, FeatID], 1, 0)] ## contaminating in > 5% of 
lower_level_intra_data <- lower_level_intra_data[BadFeat3==0]
lower_level_intra_data <- lower_level_intra_data[is.na(AlwaysZero)]
low_level_data <- lower_level_intra_data

metadata <- unique(low_level_data[,list(Sample, Strain, AcetateGroup, TimePoint, SampleFull, Sample2, F_ID, Replicate, Time)])
raw_files <- list.files("MetabolomicsWorkbench/SILAcetate_intracellular/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/SILAcetate_intracellular/raw/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub("^[A-Z][a-z]+_", "", Sample)]
files_tab[!SampleOrig %in% metadata[,Sample2]]
files_tab[,IonMode:=gsub("_.*$", "", Sample)]
setnames(files_tab, "SampleOrig", "Sample name")
low_level_data[,IonMode:=ifelse(grepl("eg", Mode), "Negative", "Positive")]
metadata[,IonMode:=ifelse(grepl("Neg", Sample), "Negative", "Positive")]
files_tab[,IonMode:=ifelse(IonMode == "Pos", "Positive", "Negative")]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode", "Sample"), by.y = c("Sample2", "IonMode", "Sample"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]

files_tab[grepl("BK", `Sample name`), Strain:=NA]
files_tab[grepl("BK", `Sample name`), AcetateGroup:=NA]
files_tab[grepl("BK", `Sample name`), Replicate:=NA]
files_tab[grepl("BK", `Sample name`), F_ID:=NA]
study_final <- dcast(files_tab, `Sample name`+Strain+AcetateGroup+TimePoint+Time+Replicate+F_ID~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/SILAcetate_intracellular///metadata_formatted.tsv", sep = "\t")

low_level_data[,table(Sample2 %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(low_level_data[IonMode=="Positive" & !grepl("Unknown", CompID)], FeatID+CompID+Tags.y+Formula+NumLabeledC~Sample2, value.var = "AreaFrac", fill = 0)
named_metabolite_table_neg <- dcast(low_level_data[IonMode=="Negative" & !grepl("Unknown", CompID)], FeatID+CompID+Tags.y+Formula+NumLabeledC~Sample2, value.var = "AreaFrac", fill = 0)
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILAcetate_intracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILAcetate_intracellular///named_data_pos.tsv", sep = "\t")

# pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
# neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
# named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]
# 
# named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(1:5),with=F]
neg_metadata <- named_metabolite_table_neg[,c(1:5),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(1, 6:(ncol(named_metabolite_table_neg))), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(1, 6:(ncol(named_metabolite_table_pos))), with=F]

fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILAcetate_intracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILAcetate_intracellular///named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/SILAcetate_intracellular/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/SILAcetate_intracellular/metabolite_metadata_neg.tsv", sep = "\t")

## Arginine SIL extracellular
datadir <- "../SILArginine2022-03/"
low_level_data <- readRDS(paste0(datadir, "SIL_Extra_HILIC_LowLevelDataAllProcessed.rds"))
metadata <- unique(low_level_data[,list(Sample, Strain, ArgGroup, TimePoint, Time, MS_ID, F_ID, Mode)])
raw_files <- list.files("MetabolomicsWorkbench/SILArginine_extracellular/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/SILArginine_extracellular/raw//mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub(".*_", "", Sample)]
metadata[,SampleOrig:=gsub(".*_", "", MS_ID)]
files_tab[!SampleOrig %in% metadata[,SampleOrig]]
files_tab[,IonMode:=ifelse(grepl("Neg", Sample), "Negative", "Positive")]
setnames(files_tab, "SampleOrig", "Sample name")
low_level_data[,IonMode:=ifelse(grepl("eg", Mode), "Negative", "Positive")]
metadata[,IonMode:=ifelse(grepl("Neg", MS_ID), "Negative", "Positive")]
metadata[,table(IonMode)]
metadata[,Mode:=NULL]
files_tab[,Sample:=NULL]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("SampleOrig", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]
#setnames(files_tab, c("Sample name", "Sample"), c("Sample name", "MSName"))
study_final <- dcast(files_tab, `Sample name`+Strain+ArgGroup+TimePoint+Time+Sample+F_ID~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/SILArginine_extracellular///metadata_formatted.tsv", sep = "\t")

low_level_data <- merge(low_level_data, files_tab[,list(IonMode, Sample, `Sample name`)], by = c("Sample", "IonMode"), all.x = T)
low_level_data[,table(`Sample name` %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(low_level_data[IonMode=="Positive" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~`Sample name`, value.var = "AreaFrac")
named_metabolite_table_neg <- dcast(low_level_data[IonMode=="Negative" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~`Sample name`, value.var = "AreaFrac")
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILArginine_extracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILArginine_extracellular///named_data_pos.tsv", sep = "\t")

# pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
# neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
# named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]
# 
# named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(1:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(1:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(1, 7:(ncol(named_metabolite_table_neg))), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(1, 7:(ncol(named_metabolite_table_pos))), with=F]

fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILArginine_extracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILArginine_extracellular///named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/SILArginine_extracellular/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/SILArginine_extracellular/metabolite_metadata_neg.tsv", sep = "\t")

## Arginine intracellular
datadir <- "../SILArginine2022-03/"
low_level_data <- readRDS(paste0(datadir, "SIL_Intra_HILIC_LowLevelDataAllProcessed.rds"))
metadata <- unique(low_level_data[,list(Sample, Strain, ArgGroup, TimePoint, Time, Replicate, MS_ID, F_ID)])
raw_files <- list.files("MetabolomicsWorkbench/SILArginine_intracellular/raw/", pattern = ".raw$")
mzml_files <- list.files("MetabolomicsWorkbench/SILArginine_intracellular/raw/mzML/", pattern = ".mzML")
files_tab <- merge(data.table(RawFile = raw_files, Sample = gsub(".raw", "", raw_files)), data.table(MzMLFile = mzml_files, Sample = gsub(".mzML", "", mzml_files)), 
                   by = "Sample", all= T)
files_tab[,SampleOrig:=gsub(".*_", "", Sample)]
metadata[,SampleOrig:=gsub(".*_", "", MS_ID)]
files_tab[!SampleOrig %in% metadata[,SampleOrig]]
files_tab[,IonMode:=ifelse(grepl("Neg", Sample), "Negative", "Positive")]
setnames(files_tab, "SampleOrig", "Sample name")
files_tab[,Sample:=NULL]
low_level_data[,IonMode:=ifelse(grepl("eg", Mode), "Negative", "Positive")]
metadata[,IonMode:=ifelse(grepl("Neg", MS_ID), "Negative", "Positive")]
metadata[,table(IonMode)]
files_tab <- merge(files_tab, metadata, by.x = c("Sample name", "IonMode"), by.y = c("SampleOrig", "IonMode"), all.x = T)
files_tab <- files_tab[IonMode %in% c("Positive", "Negative")]
files_tab <- files_tab[!grepl("shutdown", Sample)]

study_final <- dcast(files_tab, `Sample name`+Strain+ArgGroup+TimePoint+Time+Replicate+Sample~IonMode, value.var = c("RawFile", "MzMLFile"))
fwrite(study_final, file = "MetabolomicsWorkbench/SILArginine_intracellular///metadata_formatted.tsv", sep = "\t")

low_level_data <- merge(low_level_data, files_tab[,list(IonMode, Sample, `Sample name`)], by = c("Sample", "IonMode"), all.x = T)
low_level_data[,table(`Sample name` %in% files_tab[,`Sample name`])]
named_metabolite_table_pos <- dcast(low_level_data[IonMode=="Positive" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~`Sample name`, value.var = "AreaFrac", fill = 0)
named_metabolite_table_neg <- dcast(low_level_data[IonMode=="Negative" & !grepl("Unknown", CompID)], FeatID+CompID+MolWt+MSI+Formula+NumLabeledC~`Sample name`, value.var = "AreaFrac", fill = 0)
fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILArginine_intracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILArginine_intracellular///named_data_pos.tsv", sep = "\t")

# pos_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_pos.txt")
# neg_names <- fread("MetabolomicsWorkbench/Strains/refmet_results_neg.txt")
# named_metabolite_table_pos <- merge(named_metabolite_table_pos, pos_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_pos[,`metabolite name`:=ifelse(`Standardized name` == "-"| is.na(`Standardized name`), MetName, `Standardized name`)]
# 
# named_metabolite_table_neg <- merge(named_metabolite_table_neg, neg_names[,c(1,2)], by.x = "MetName", by.y = "Input name", all.x = T)
# named_metabolite_table_neg[,`metabolite name`:=ifelse(`Standardized name` == "-" | is.na(`Standardized name`), MetName, `Standardized name`)]

pos_metadata <- named_metabolite_table_pos[,c(1:6),with=F]
neg_metadata <- named_metabolite_table_neg[,c(1:6),with=F]

named_metabolite_table_neg <- named_metabolite_table_neg[,c(1, 7:(ncol(named_metabolite_table_neg))), with=F]
named_metabolite_table_pos <- named_metabolite_table_pos[,c(1, 7:(ncol(named_metabolite_table_pos))), with=F]

fwrite(named_metabolite_table_neg, file = "MetabolomicsWorkbench/SILArginine_intracellular//named_data_neg.tsv", sep = "\t")
fwrite(named_metabolite_table_pos, file = "MetabolomicsWorkbench/SILArginine_intracellular///named_data_pos.tsv", sep = "\t")
fwrite(unique(pos_metadata), file = "MetabolomicsWorkbench/SILArginine_intracellular/metabolite_metadata_pos.tsv", sep = "\t")
fwrite(unique(neg_metadata), file = "MetabolomicsWorkbench/SILArginine_intracellular/metabolite_metadata_neg.tsv", sep = "\t")

