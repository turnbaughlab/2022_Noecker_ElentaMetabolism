## Growth curve processing and visualization functions
# 9/21/2020

library(growthcurver)
library(data.table)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggrepel)
theme_set(theme_cowplot())

# Read in file in the format exported by Gen5 for repeated measurements
# Uses fread(cmd) to pull the key part of the file
read_growthcurve_file = function(filename){
  #Get everything between the Time and temp column names and the results column name
  growth_data1 = fread(cmd = paste0("awk '/^Time\tT/,/^Results/' '", path.expand(filename), "'"), blank.lines.skip = T, header = T)  
  growth_data1[,Time:=as.numeric(lubridate::seconds(lubridate::hms(Time)))/3600]
  setnames(growth_data1, "Time", "time")
  setnames(growth_data1, names(growth_data1)[2], "temp")
  growth_data1 = growth_data1[!is.na(temp)] # Remove Missing/stopped early data
  return(growth_data1)
}

# Different options for normalizing OD readings
# gd = growth data table
# blank_wells specified wells to normalize by
# skip_pts wells to skip
# method argument can be "medianBlank" or "earlyData" or "blankTimeMatch"
# correct_val provide a single background value to use for normalization
# return_long whether normalized table should be returned in long or wide form
# wide_form Can provide table in long or wide form (expecting wide by default)
# ctrlMaxTime Don't normalize after this time (default 1000 hours)
# ctrlMinTime Don't normalize before this time (default 0)
correct_growth_data_custom = function(gd, blank_wells = NULL, skip_pts = c(), method = "medianBlank", 
                                      correct_val = NULL, return_long = T, wide_form = T, 
                                      ctrlMaxTime = 1000, ctrlMinTime = 0){
  gd = as.data.table(gd) #Convert to data.table
  if(wide_form) gd = melt(gd, id.var = c("time", "temp"), variable.name = "well")
  wells_to_correct = gd[!well %in% blank_wells, unique(well)]
  
  if(!is.null(correct_val)){
    cat("Correcting based on provided value\n")
    gd[,NewValue:=value-correct_val]
    gd_corrected = gd
  } else if(method == "medianBlank"){ #Use median reading of all blank wells as the correction value
    if(!is.null(names(blank_wells))){ # If a vector of blank wells is provided
      #A different set of blanks for each media group
      cat("Correcting based on provided blank wells for each condition\n")
      media_tab = data.table(Condition = names(blank_wells), well = blank_wells)
      correct_vars = gd[well %in% blank_wells & time < ctrlMaxTime & time >= ctrlMinTime, median(value), by=Condition]
      gd_corrected = merge(gd[!well %in% blank_wells], correct_vars[,list(Condition, V1)], by = "Condition", all.x = T)
      
    } else {
      # Or else there are no groups
      cat("Correcting based on all blank wells\n")
      correct_vars = gd[well %in% blank_wells, median(value), by = well]
      if(by_column){
        correct_vars[,Column:=gsub("[A-Z]", "", well)]
        gd[,Column:=gsub("[A-Z]", "", well)]
        gd_corrected = merge(gd[!well %in% blank_wells], correct_vars[,list(Column, V1)], by = "Column", all.x = T)
      } 
    }
    gd_corrected[,NewValue:=ifelse(time %in% skip_pts, value, value-V1)] #Skip requested points
  } else if(method == "earlyData"){ 
    # Correct each well based on its early time points
    cat("Correcting based on early time points\n")
    correct_vars = gd[time < ctrlMaxTime & time >= ctrlMinTime, median(value), by=well]
    gd_corrected = merge(gd, correct_vars, by = "well", all.x = T)
    gd_corrected[,NewValue:=ifelse(time %in% skip_pts, value, value - V1)]
  } else if(method == "blankTimeMatch"){ 
    # Correct based on blank readings at the same time point
    if(!is.null(names(blank_wells))){
      media_tab = data.table(Condition = names(blank_wells), well = blank_wells)
    }
    cat("Correcting based on provided blank wells for each condition, time-matched\n")
    correct_vars = gd[well %in% blank_wells & time < ctrlMaxTime, median(as.numeric(value), na.rm = T), by=list(Condition, time)]
    gd_corrected = merge(gd[!well %in% blank_wells], correct_vars[,list(Condition, time, V1)], by = c("Condition", "time"), all.x = T)
    gd_corrected[,NewValue:=ifelse(time %in% skip_pts, as.numeric(value), as.numeric(value)-V1)]  
  }
  if(return_long) return(gd_corrected) else return(dcast(gd_corrected, time~well, value.var = "NewValue", fill = 0))
}

# Fit growthcurver models, assumes data is already normalized
# growth_data_norm Normalized data table
# indiv_fits Whether to use growthcurver function for all wells simultaneously (SummarizeGrowthByPlate) or fit curves for each well separately (SummarizeGrowth)
# trim_after Trim after this time point for all wells
# all_sample_mins Trim before these time points (one for each well)
# all_sample_maxes Trim after these time points (one for each well)
fit_growthcurver_models = function(growth_data_norm, indiv_fits = F, trim_after = 42, 
                                   all_sample_mins = NULL, all_sample_maxes = NULL){
  ## Get rid of negative values
  growth_data_norm = as.data.table(growth_data_norm)
  all_samps = names(growth_data_norm)[names(growth_data_norm) != "time"]
  for(j in 1:length(all_samps)){
    if(any(growth_data_norm[!is.na(get(all_samps[j])),get(all_samps[j])] < 0)){
      growth_data_norm[,all_samps[j] := get(all_samps[j]) - min(get(all_samps[j]), na.rm = T)]
    }
  }
  if(!indiv_fits){
    results_table = data.table(SummarizeGrowthByPlate(growth_data_norm, t_trim = trim_after, bg_correct = "none", plot_fit = F))
    setnames(results_table, c("sample", "n0"), c("Sample", "N0")) #For consistency
    return(results_table)
  } else {
    if(is.null(all_sample_mins)){
      all_sample_mins = rep(0, length(all_samps))
    }
    if(is.null(all_sample_maxes)){
      all_sample_maxes = rep(1000, length(all_samps))
    }
    if(length(all_sample_mins) != length(all_samps)|length(all_sample_maxes) != length(all_samps)){
      stop("Time cutoffs not specified for all wells")
    }
    all_samps = names(growth_data_norm)[names(growth_data_norm) != "time"]
    all_samp_fits = lapply(1:length(all_samps), function(x){
      values_fit <- growth_data_norm[growth_data_norm$time > all_sample_mins[x],get(all_samps[x])]
      good_vals <- which(!is.na(values_fit))
      time_fit <- growth_data_norm$time[growth_data_norm$time > all_sample_mins[x]][good_vals]
      values_fit <- values_fit[good_vals]
      if(length(values_fit) > 0) {
        return(SummarizeGrowth(time_fit, values_fit, 
                               t_trim = all_sample_maxes[x], bg_correct = "none"))
        
      } else { return(NULL)}
    })
    results_table = rbindlist(lapply(1:length(all_samp_fits), function(x){
      return(data.table(Sample = all_samps[x], r = all_samp_fits[[x]]$vals$r, k = all_samp_fits[[x]]$vals$k, 
                        N0 = all_samp_fits[[x]]$vals$n0,r_se = all_samp_fits[[x]]$vals$r_se, r_p = all_samp_fits[[x]]$vals$r_p,
                        k_se = all_samp_fits[[x]]$vals$k_se, k_p = all_samp_fits[[x]]$vals$k_p,
                        N0_se = all_samp_fits[[x]]$vals$N0_se, t_mid = all_samp_fits[[x]]$vals$t_mid, 
                        t_gen = all_samp_fits[[x]]$vals$t_gen, sigma = all_samp_fits[[x]]$vals$sigma, 
                        auc_l = all_samp_fits[[x]]$vals$auc_l, auc_e = all_samp_fits[[x]]$vals$auc_e))
    }), fill = T)
    return(results_table)
  }
}

#calculate the harmonic mean (estimated growth rate across all wells in a group)
get_group_growth_rates = function(results_table, grouping_vars){
  harmonic_mean_growth_rate = results_table[,length(r)/mean(sapply(r, function(x){ 1/x})), by=grouping_vars]
  setnames(harmonic_mean_growth_rate, "V1", "HarmonicMean")
  return(harmonic_mean_growth_rate)
}


# Plot growth curves
# growth_data_norm : Normalized growth data
# group_metadata : Metadata table
# color_var and facet_var : one or more column names in group_metadata, to use for color and faceting
# mean_se : Whether to show one line for each well or a ribbon plot of summarized growth across replicates, default F
# add_fitted_curves : Whether to add logistic model fits (must provide table of model fit results), deafult F
# model_fits : Table of model fit results (if add_fitted_curves = T)
plot_growth_curves = function(growth_data_norm, group_metadata, color_var = NULL, facet_var = NULL, mean_se = F, add_fitted_curves = F, model_fits = NULL){
  growth_data_norm = as.data.table(growth_data_norm)
  group_metadata = as.data.table(group_metadata)
  growth_data_norm1 = melt(growth_data_norm, id.var = "time", variable.name = "well", variable.factor = F)
  growth_data_norm1 = merge(growth_data_norm1, group_metadata, by = "well", all.x = T)
  if(mean_se){
    plot_data = growth_data_norm1[,list(mean(value), sd(value)/sqrt(length(value))), by=c("time", color_var, facet_var)]
    setnames(plot_data, c("V1", "V2"), c("MeanValue", "SEValue"))
    plot1 = ggplot(plot_data, aes(x = time, y = MeanValue)) + geom_ribbon(aes(ymin = MeanValue - SEValue, ymax = MeanValue + SEValue, fill = get(color_var)), alpha = 0.7) + scale_fill_brewer(palette = "Set1", name = color_var)
    if(!is.null(facet_var)) plot1 = plot1 + facet_wrap(formula(paste0("~", paste(facet_var, collapse = "+"))))
  } else {
    plot_data = growth_data_norm1
    plot1 = ggplot(plot_data, aes(x = time, y = value)) + geom_line(aes(color = get(color_var), group = well)) + scale_color_brewer(palette = "Set1", name = color_var)
    if(!is.null(facet_var)) plot1 = plot1 + facet_wrap(formula(paste0("~", paste(facet_var, collapse = "+"))))
  }
  if(add_fitted_curves & !is.null(model_fits)){
    fitted_data = get_fitted_data(model_fits, time_max = growth_data_norm[,max(time)])
    fitted_data = merge(fitted_data, group_metadata, by.x = "Sample", by.y = "well", all.x = T)
    plot1 = ggplot(fitted_data[Sample %in% growth_data_norm1[,unique(well)]], aes(color = get(color_var))) + geom_line(aes(x = time, y = PredictedValue, group = Sample, color = get(color_var))) + 
      geom_point(data = growth_data_norm1, aes(x = time, y = value), size = 0.5, alpha = 0.5) 
    if(!is.null(facet_var)) plot1 = plot1 + facet_wrap(formula(paste0("~", paste(facet_var, collapse = "+"))))+ scale_color_brewer(palette = "Set1", name = color_var)
  }
  return(plot1 + ylab("Normalized OD600"))
}

#logistic function
logistic_curve = function(t, r, k, n0){
  return(n0 * exp(r*t)/(1 + n0 * (exp(r*t) - 1)/k))
}

#get fitted values for a set of curves
get_fitted_data = function(results_table, time_min = 0, time_max = 48, time_step = 0.25){
  pred_data = rbindlist(lapply(1:nrow(results_table), function(x){
    n0 = results_table[x, N0]
    r = results_table[x, r]
    k = results_table[x, k]
    pred_table1 = data.table(Sample = results_table[x, Sample], time = seq(time_min, time_max, by = time_step))
    pred_table1[,PredictedValue:=logistic_curve(time, r = r, k = k, n0 = n0)]
    return(pred_table1)
  }))
  return(pred_data)
}

