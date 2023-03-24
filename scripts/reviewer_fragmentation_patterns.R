### feature autocorrelation

time_course_dataset <- fread("processedDatasets/timecourse_all_data.csv")

time_course_mat <- dcast(time_course_dataset[Duplicate == 0], IDUniq~Sample, value.var = "value") %>% column_to_rownames("IDUniq") %>% as.matrix()
cor_mat <- cor(t(time_course_mat)) %>% as.data.frame() %>% rownames_to_column("IDUniq") %>% as.data.table() %>% 
  melt(id.var = "IDUniq", variable.name = "IDUniq2", variable.factor= F) %>% filter(IDUniq < IDUniq2)

cor_mat_spearman <- cor(t(time_course_mat), method = "spearman") %>% as.data.frame() %>% rownames_to_column("IDUniq") %>% as.data.table() %>% 
  melt(id.var = "IDUniq", variable.name = "IDUniq2", variable.factor= F) %>% filter(IDUniq < IDUniq2)

cor_mat <- merge(cor_mat, cor_mat_spearman, by = c("IDUniq", "IDUniq2"))
ggplot(cor_mat, aes(x = value.x)) + geom_histogram(binwidth = 0.02) + xlab("Pearson correlation between features") + ylab("Number of feature pairs")
cor_mat[value.x > 0.995 & value.y > 0.99]

adduct_tab <- unique(time_course_dataset[,list(IDUniq, Adduct, MetName, MSI, CF_kingdom, CF_superclass, 
                                               CF_class, CF_subclass, CF_Dparent)])
cor_mat <- merge(cor_mat, adduct_tab, by = "IDUniq", all.x = T)
cor_mat <- merge(cor_mat, adduct_tab, by.x = "IDUniq2", by.y = "IDUniq", all.x = T)

## What are the mass differences between correlated features? Pretend everything is M+H or M-H
cor_mat[,c("MZ", "RT", "AlignmentID", "IonMode"):=tstrsplit(IDUniq, split = "_")]
cor_mat[,c("MZ2", "RT2", "AlignmentID2", "IonMode2"):=tstrsplit(IDUniq2, split = "_")]
cor_mat[,RTDiff:=abs(as.numeric(RT2)-as.numeric(RT))]
# cor_mat[,MassDiff:=case_when(
#   IonMode == IonMode2 ~ abs(as.numeric(MZ2)-as.numeric(MZ))
#   IonMode == "Pos" & IonMode2 == "Neg" ~ abs(as.numeric(MZ2)-as.numeric(MZ)) 
# )]
cor_mat[,Comparison:=case_when(
  IonMode == IonMode2 & IonMode=="Pos" ~ "BothPos",
  IonMode == IonMode2 & IonMode=="Neg" ~ "BothNeg",
  IonMode != IonMode2 ~ "PosNeg"
)]
cor_mat[,MZDiff:=abs(as.numeric(MZ2)-as.numeric(MZ))]

ggplot(cor_mat[value.x > 0.99 & value.y > 0.99], aes(x = RTDiff, y = MZDiff)) + geom_bin2d(color = "black", bins= 50) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(~Comparison)
ggplot(cor_mat[value.x > 0.95 & value.y > 0.95], aes(x = RTDiff, y = MZDiff)) + geom_bin2d(color = "black", bins= 50) + scale_fill_distiller(palette = "Blues", direction = 1) + facet_wrap(~Comparison)

ggplot(cor_mat[value.x > 0.99 & value.y > 0.99 & RTDiff < 0.2], aes(x = MZDiff)) + geom_histogram(binwidth = 2) + 
  geom_vline(xintercept = 18, color = "red") + facet_wrap(~Comparison)

ggplot(cor_mat[MZDiff < 0.002 & RTDiff < 0.2], aes(x = value.x, y = value.y)) + geom_bin2d(binwidth = 0.01, color = "black") + scale_fill_distiller(palette = "Blues", direction = 1) + 
  xlim(0.5, 1) + ylim(0.5, 1)

ggplot(cor_mat[RTDiff < 0.2 & MZDiff < 0.002], aes(x = value.x)) + geom_histogram(binwidth = 0.02) + xlab("Pearson correlation between features\nwith similar m/z and RT") + ylab("Number of feature pairs")
ggplot(cor_mat[RTDiff < 0.2 & MZDiff < 0.002], aes(x = value.y)) + geom_histogram(binwidth = 0.02) + xlab("Spearman correlation between features\nwith similar m/z and RT") + ylab("Number of feature pairs")
ggplot(cor_mat[RTDiff < 0.2], aes(x = value.x)) + geom_histogram(binwidth = 0.02) + xlab("Pearson correlation between features\nwith similar RT") + ylab("Number of feature pairs")

cor_mat_long <- melt(cor_mat, id.var = c("IDUniq", "IDUniq2", "RTDiff", "MZDiff", "Comparison"), variable.name = 'CorrType', measure.vars = c("value.x", "value.y"))
cor_mat_long[,CorrType:=ifelse(CorrType=="value.x", "Pearson", "Spearman")]
hist1 <- ggplot(cor_mat_long, aes(x = value)) + geom_histogram(binwidth = 0.02) + facet_wrap(~CorrType) + xlab("Correlation between features") + ylab("Number of feature pairs")
hist2 <- ggplot(cor_mat_long[RTDiff < 0.2], aes(x = value)) + geom_histogram(binwidth = 0.02) + facet_wrap(~CorrType) + xlab("Correlation between features\nwith similar RT") + ylab("Number of feature pairs")
hist3 <- ggplot(cor_mat_long[MZDiff < 0.002 & RTDiff < 0.2], aes(x = value)) + geom_histogram(binwidth = 0.02) + facet_wrap(~CorrType)+ xlab("Correlation between features\nwith similar RT and m/z") + ylab("Number of feature pairs")


cor_mat[value.x > 0.97 & RTDiff < 0.2 & MSI.x != "" & MSI.y != ""]
cor_mat[value.x > 0.95 & value.y > 0.95 & RTDiff < 0.2 & MSI.x != "" & MSI.y != ""]
cor_mat[value.x > 0.95 & value.y > 0.95 & RTDiff < 0.2 & (MSI.x != "" | MSI.y != "")]

cor_mat[,BadPair:=ifelse(CF_class.x != CF_class.y & CF_class.x != "no matches" & CF_class.y != "no matches", 1, 0)]
cor_mat_sub <- cor_mat[value.x > 0.95 & value.y > 0.95 & RTDiff < 0.2 & BadPair == 0] 
adduct_tab[,HasClosePair:=ifelse(IDUniq %in% cor_mat_sub[,IDUniq]|IDUniq %in% cor_mat_sub[,IDUniq2], 1, 0)]
adduct_tab[,ClosePairs:=sapply(IDUniq, function(x){
  return(nrow(cor_mat_sub[IDUniq == x|IDUniq2==x]))
})]
adduct_tab[,table(HasClosePair)]
adduct_tab[order(ClosePairs, decreasing = T)]

adduct_tab[,hist(ClosePairs, breaks = 35)]
hist4 <- ggplot(adduct_tab, aes(x = ClosePairs)) + geom_histogram(binwidth = 1) + 
  xlab("Number of linked/correlated features") + ylab("Number of features") #+ facet_wrap(~ifelse(MSI=="", "Unknown", "ID"), scales = "free_y")

plot_all <- (hist1 + hist2 + hist3) / hist4

ggsave(plot_all, file = "figures/histograms_featureCorr.pdf", width = 8, height = 5)

#### Define connected components
library(igraph)
library(ggraph)
foo1 <- graph_from_data_frame(cor_mat_sub[order(IDUniq),list(IDUniq, IDUniq2, value.x, value.y, RTDiff)], directed = F, 
                              vertices = adduct_tab[IDUniq %in% cor_mat_sub[,IDUniq]|IDUniq %in% cor_mat_sub[,IDUniq2],list(IDUniq, MSI, MetName)])

clust_def <- components(foo1)
clust_def_tab <- data.table(IDUniq = names(clust_def$membership), ClustID = clust_def$membership)
adduct_tab <- merge(adduct_tab, clust_def_tab, by = "IDUniq", all.x = T)
adduct_tab[,table(ClustID)]
V(foo1)$MetName <- ifelse(grepl("known", V(foo1)$MetName), "", V(foo1)$MetName)
correlation_network <- ggraph(foo1) + geom_edge_link() + geom_node_point(color = "cornflowerblue") + geom_node_text(aes(label = MetName), size = 3, color = "darkred") +
  theme(panel.background = element_blank())
ggsave(correlation_network, file = "figures/time_course_feature_RT_corr_network.pdf", width = 9, height = 7)
## How to pick a representative feature? Just pick the largest m/z?
adduct_tab[,c("MZ", "RT", "AlignmentID", "IonMode"):=tstrsplit(IDUniq, split = "_")]

adduct_tab[!is.na(ClustID),ClustMaxMZ:=max(as.numeric(MZ)), by = ClustID]
adduct_tab[!is.na(ClustID) & MSI != "",RepID:=1]
no_rep_yet <- adduct_tab[!is.na(ClustID), length(RepID[!is.na(RepID)]), by=ClustID][V1 ==0, ClustID]
adduct_tab[ClustID %in% no_rep_yet & MZ==ClustMaxMZ, RepID:=1]
adduct_tab[ClustID %in% adduct_tab[,sum(RepID==1 & MSI == "") > 1, by = ClustID][V1==T, ClustID] & RepID==1, 
           RepID:=rep(c(1,0), 2)] ## for m/z ties, just pick the first one randomly
rep_features <- adduct_tab[RepID==1 | is.na(ClustID)]

data_sub <- time_course_dataset[IDUniq %in% rep_features[,IDUniq]]
fwrite(data_sub, file = "processedDatasets/time_course_all_data_clustered_features_test.csv")
