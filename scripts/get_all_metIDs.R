library(data.table)
library(tidyverse)
library(readxl)
library(webchem)
library(omu)
library(broom)
library(lmerTest)
library(fuzzyjoin)
library(broom.mixed)

inchikeys_pilot <- fread("../processedDatasets/pilot_inchikeys.txt")
inchikeys_time <- fread("../processedDatasets/timecourse_inchikeys.txt")
inchikeys_strains <- fread("../processedDatasets/strains_gl_inchikeys.txt")
inchikeys_mouse <- fread("../processedDatasets/mouse_contents_inchikeys.txt")
inchikeys_serum <- fread("../processedDatasets/mouse_serum_inchikeys.txt")
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
fwrite(kegg_class, file = "processedDatasets/inchikeys_kegg_all.csv")