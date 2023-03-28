## Replot of 2019 curves
library(data.table)
datadir <- "../GrowthExperiments/2019_10_21_GMMd_LAB_Elen2243/"
datafile <- paste0(datadir, "Experiment_elen2243_gmm.txt")
metadata_file <- paste0(datadir, "layout.xlsx")
source("scripts/growth_curve_functions.R")
# Set up metadata
layout_info = data.table(read_xlsx(metadata_file, sheet = 1))
layout_info[,`6`:=NULL]
setnames(layout_info, c("Row", 1:5))
layout_info[,Type:=ifelse(is.na(Row), "Strain", "Media")]
layout_info <- layout_info[1:8]
layout_info[,Row:=sapply(1:nrow(layout_info), function(x){
  if(!is.na(layout_info[x, Row])){
    return(layout_info[x,Row])
  } else {
    return(layout_info[x-1, Row])
  }
})]
layout_info = melt(layout_info, id.var = c("Row", "Type"), variable.name = "Column")
#layout_info[value=="" & Type=="Strain", value:="Control"]
layout_info = dcast(layout_info, Row+Column~Type, value.var = "value")
layout_info[,well:=paste0(Row, Column)]
layout_info[Strain == "c", Strain:="Control"]
#For normalization, the metadata needs to have a column named "Condition" to specify which blank wells/controls 
# should be used for normalization with different experimental groups
setnames(layout_info, "Media", "Condition") #In this case we have sterile blanks for each media group
layout_info[,Condition:=case_when(
  Condition == "GMM+Arg" ~ "EDM2", 
  Condition == "GMM+LABsub+Arg" ~ "EDM1",
  Condition == "GMM+LABall+Arg" ~ "EDMall", 
  Condition == "ISP-2" ~ "ISP-2+Arg",
  Condition == "BHI+Arg" ~ "BHI+Arg"
)]
# Get set of blank wells for normalization from metadata
blank_wells1 = layout_info[Strain=="Control", well]
names(blank_wells1) = layout_info[Strain=="Control", Condition]

# Read in file
growth_data1 = read_growthcurve_file(datafile)

# Formatting and initial plot
growth_data_long = melt(growth_data1, id.var = c("time", "temp"), variable.name = "well")
growth_data_long = merge(growth_data_long, layout_info, by = "well", all.x = T)
growth_data_long = growth_data_long[!is.na(Condition)] #wells without metadata were empty
ggplot(growth_data_long, aes(x=time, y= value)) + geom_line(aes(group = well, color = well %in% blank_wells1))
ggplot(growth_data_long[!well %in% blank_wells1], aes(x=time, y= value)) + geom_line() + facet_wrap(~well) + 
  geom_vline(xintercept = 44, linetype = 2)
growth_data_long <- growth_data_long[!is.na(value)]

# Normalize based on blank wells
corrected_data = correct_growth_data_custom(growth_data_long, blank_wells = blank_wells1, 
                                            method = "blankTimeMatch", wide_form = F)
# Change to look for a Treatment or Condition variable
growth_data_norm = dcast(corrected_data, time~well, value.var = "NewValue", fill = 0)

# Fit curves with growthcurver
model_fits = fit_growthcurver_models(growth_data_norm, indiv_fits = F, trim_after = 44)

# Fit curves with different trimming parameters for different wells
all_samps = corrected_data[,unique(well)]
trim_after_wells = rep(44, times = length(all_samps))
#trim_after_wells[grepl("10", all_samps)] =

trim_before_wells = rep(0, length(all_samps))
#trim_before_wells[all_samps=="C10"] = 5
#trim_before_wells[all_samps=="D8"] = 4

model_fits2 = fit_growthcurver_models(growth_data_norm, indiv_fits = T, 
                                      all_sample_mins = trim_before_wells, 
                                      all_sample_maxes = trim_after_wells)
model_fits2 = merge(model_fits2, layout_info, by.x = "Sample", by.y = "well", all.x = T)

# Check bad model fits
bad_fits = model_fits2[(k < 0.01 & r > 0.5) | k > 1 | t_mid < 0, Sample]
ggplot(corrected_data[well %in% bad_fits], aes(x=time, y = NewValue, group = well, color = Strain)) + geom_line() + facet_wrap(~well)

layout_info[,Condition:=factor(Condition, levels = c("BHI+Arg", "ISP-2+Arg", "EDM1", "EDM2", "EDMall"))]
# Make nice plot
plot1 = plot_growth_curves(growth_data_norm = growth_data_norm, 
                           group_metadata = layout_info, 
                           color_var = "Condition", facet_var = "Strain")
plot2 = plot_growth_curves(growth_data_norm = growth_data_norm, 
                           group_metadata = layout_info, color_var = "Condition", 
                           facet_var = "Condition", mean_se = T) + ylab("E. lenta DSM 2243 abundance (Normalized OD600)") + theme_classic2()

growth_data_norm1 = melt(growth_data_norm, id.var = "time", variable.name = "well", variable.factor = F)
growth_data_norm1 = merge(growth_data_norm1, layout_info, by = "well", all.x = T)
plot_data = growth_data_norm1[,list(mean(value), sd(value)/sqrt(length(value))), by=c("time", "Condition")]
setnames(plot_data, c("V1", "V2"), c("MeanValue", "SEValue"))
plot1 = ggplot(plot_data, aes(x = time, y = MeanValue)) + geom_ribbon(aes(ymin = MeanValue - SEValue, ymax = MeanValue + SEValue, fill = Condition)) + 
  scale_fill_brewer(palette = "Set1", guide = "none") + facet_grid(.~Condition) + theme_classic2()+ ylab("E. lenta DSM 2243 abundance (Normalized OD600)") 
plot1 = ggplot(growth_data_norm1, aes(x = time, y = value)) + geom_line(aes(group = well, color = Condition)) + 
  scale_color_brewer(palette = "Set1", guide = "none") + facet_grid(.~Condition) + theme_classic2()+ ylab("E. lenta DSM 2243 abundance (Normalized OD600)") 

ggsave(plot1, file = "orig_edm_experiment_curves.pdf", width = 6, height = 1.7)

fwrite(growth_data_norm1, file = "figure_source_data/figureS1B_growthCurves.tsv", quote=F, row.names = F, sep = "\t")

