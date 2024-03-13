library(pwr)
library(lme4)
library(tidyr)
library(dplyr)

#################################### AIM 1 #####################################

# Define basic inputs for experimental simulations 
sig_value = 0.05 # Alpha 
power_value = 0.8 # 1 - beta 
synthesized_data_k = 500 # Number of simulated experiments per N value (higher means more reliable power estimation)
synthesized_data_n_array = 15 # Maximum group N value assessed (should not be too high anyway)
margin = 0.8 # Conservative safeguard for final N
bio_var = 0.15 # Expected between-subject variation 

# Initialize array to hold power calculations in 
power_values = data.frame(matrix(nrow = synthesized_data_n_array, ncol = 3)) 

# Calculate statistical power for models fit to experiments simulated at every N value
for (kk in 3:synthesized_data_n_array){
print(kk);
synthesized_data_n = kk 

# Initialize array to hold p values for each fit in this round of simulated experiments  
p_values = data.frame(matrix(nrow = synthesized_data_k, ncol = 3)) 

# It's 1H-MRS o'clock
for (mm in 1:synthesized_data_k){
  
### WILDTYPE 
# Generate the core data vector based on preliminary empirical findings 
preliminary_experiment = as.data.frame(matrix(nrow=7,ncol=7))
colnames(preliminary_experiment) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
preliminary_experiment$id = c(1, 1, 1, 1, 1, 1, 1)
preliminary_experiment$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
preliminary_experiment$metab = c(6.036, 5.814, 5.676, 5.158, 4.843, 4.727, 4.833)
preliminary_experiment$crlb = c(5, 4, 5, 5, 5, 6, 3)
preliminary_experiment$metab_kx = c(4.843, 4.843, 5.158, 4.833, 4.833, 4.833, 4.833)
preliminary_experiment$crlb_kx = c(5, 5, 5, 3, 3, 3, 3)
preliminary_experiment$strain = c('wt', 'wt', 'wt', 'wt', 'wt', 'wt', 'wt')

nrow_prelimsynth = length(preliminary_experiment$time)*synthesized_data_n
prelimsynth = as.data.frame(matrix(nrow=0,ncol=7))
colnames(prelimsynth) = colnames(preliminary_experiment)

# Generate further synthesized data vectors according to defined biological variability and measured CRLB
for (ii in 1:synthesized_data_n) {
  
  if (ii==1){
    prelimsynthnew = preliminary_experiment
  }
  else{
    prelimsynthnew = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(prelimsynthnew) = colnames(prelimsynth)
    prelimsynthnew$id = c(ii, ii, ii, ii, ii, ii, ii) 
    prelimsynthnew$time = preliminary_experiment$time
    prelimsynthnew$crlb = preliminary_experiment$crlb
    prelimsynthnew$crlb_kx = preliminary_experiment$crlb_kx
    prelimsynthnew$strain = preliminary_experiment$strain
  
    metabsynth = preliminary_experiment$metab
    metab_kx_synth = preliminary_experiment$metab_kx
    
    for (jj in 1:length(preliminary_experiment$metab)){
      metab_orig = metabsynth[jj]; 
      metab_kx_orig = metab_kx_synth[jj]; 
      
      #Add biological variability 
      sd_bio_variability = bio_var*metab_orig 
      metab_orig_mean = rnorm(1, mean = metab_orig, sd = sd_bio_variability)
      sd_bio_kx_variability = bio_var*metab_kx_orig 
      metab_kx_orig_mean = rnorm(1, mean = metab_kx_orig, sd = sd_bio_kx_variability)
      
      # Add measurement variability according to CRLB
      metab_orig_crlb = preliminary_experiment$crlb[jj]
      metab_orig_sd = metab_orig_mean*metab_orig_crlb/100;
      metabsynth[jj] = rnorm(1, mean = metab_orig_mean, sd = metab_orig_sd)
      
      metab_kx_orig_crlb = preliminary_experiment$crlb_kx[jj]
      metab_kx_orig_sd = metab_kx_orig_mean*metab_kx_orig_crlb/100;
      metab_kx_synth[jj] = rnorm(1, mean = metab_kx_orig_mean, sd = metab_kx_orig_sd)
      
      # Control for NA metabolite values
      if (is.na(metabsynth[jj])) {
        metabsynth[jj]=0
      }
      
      # Control for negative metabolite values, which are not output by LCModel
      if (metabsynth[jj]<0) {
        metabsynth[jj]=0
      }
      
      # Control for NA metabolite values
      if (is.na(metab_kx_synth[jj])) {
        metab_kx_synth[jj]=0
      }
      
      # Control for negative metabolite values, which are not output by LCModel
      if (metab_kx_synth[jj]<0) {
        metab_kx_synth[jj]=0
      }
    }
    
    prelimsynthnew$metab = metabsynth
    prelimsynthnew$metab_kx = metab_kx_synth
  }
  
  prelimsynth = rbind(prelimsynth, prelimsynthnew)
}

## AGED
# Generate the core data vector based on preliminary empirical findings and hypothesis regarding age effects
preliminary_experiment_path = as.data.frame(matrix(nrow=7,ncol=7))
colnames(preliminary_experiment_path) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
idnum_path = synthesized_data_n + 1
preliminary_experiment_path$id = c(idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path)
preliminary_experiment_path$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
preliminary_experiment_path$metab = c(6.036, 5.814, 5.676, 5.158, 4.843, 4.727, 4.833)
preliminary_experiment_path$crlb = c(5, 4, 5, 5, 5, 6, 3)
preliminary_experiment_path$metab_kx = c(5.814, 5.814, 5.814, 5.814, 5.814, 5.814, 5.814)
preliminary_experiment_path$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
preliminary_experiment_path$strain = c('path', 'path', 'path', 'path', 'path', 'path', 'path')

nrow_prelimsynth_path = length(preliminary_experiment_path$time)*synthesized_data_n
prelimsynth_path = as.data.frame(matrix(nrow=0,ncol=7))
colnames(prelimsynth_path) = colnames(preliminary_experiment_path)

# Generate synthesized data vectors according to CRLB measured from preliminary data vector 
for (ii in 1:synthesized_data_n) {
  
  if (ii==1){
    prelimsynthnew_path = preliminary_experiment_path
  }
  else{
    prelimsynthnew_path = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(prelimsynthnew_path) = colnames(prelimsynth_path)
    idnum = ii + synthesized_data_n 
    prelimsynthnew_path$id = c(idnum, idnum, idnum, idnum, idnum, idnum, idnum) 
    prelimsynthnew_path$time = preliminary_experiment_path$time
    prelimsynthnew_path$crlb = preliminary_experiment_path$crlb
    prelimsynthnew_path$crlb_kx = preliminary_experiment_path$crlb_kx
    prelimsynthnew_path$strain = preliminary_experiment_path$strain
    
    metabsynth_path = preliminary_experiment_path$metab
    metab_kx_synth_path = preliminary_experiment_path$metab_kx
    
    for (jj in 1:length(preliminary_experiment_path$metab)){
      metab_orig_path = metabsynth_path[jj]; 
      metab_kx_orig_path = metab_kx_synth_path[jj]; 
      
      #Add biological variability 
      sd_bio_variability_path = bio_var*metab_orig_path 
      metab_orig_mean_path = rnorm(1, mean = metab_orig_path, sd = sd_bio_variability_path)
      sd_bio_kx_variability_path = bio_var*metab_kx_orig_path 
      metab_kx_orig_mean_path = rnorm(1, mean = metab_kx_orig_path, sd = sd_bio_kx_variability_path)
      
      # Add measurement variability according to CRLB
      metab_orig_crlb_path = preliminary_experiment_path$crlb[jj]
      metab_orig_sd_path = metab_orig_mean_path*metab_orig_crlb_path/100;
      metabsynth_path[jj] = rnorm(1, mean = metab_orig_mean_path, sd = metab_orig_sd_path)
      
      metab_kx_orig_crlb_path = preliminary_experiment_path$crlb_kx[jj]
      metab_kx_orig_sd_path = metab_kx_orig_mean_path*metab_kx_orig_crlb_path/100;
      metab_kx_synth_path[jj] = rnorm(1, mean = metab_kx_orig_mean_path, sd = metab_kx_orig_sd_path)
      
      # Control for NA metabolite values
      if (is.na(metabsynth_path[jj])) {
        metabsynth_path[jj]=0
      }
      
      # Control for negative metabolite values, which are not output by LCModel
      if (metabsynth_path[jj]<0) {
        metabsynth_path[jj]=0
      }
      
      # Control for NA metabolite values
      if (is.na(metab_kx_synth_path[jj])) {
        metab_kx_synth_path[jj]=0
      }
      
      # Control for negative metabolite values, which are not output by LCModel
      if (metab_kx_synth_path[jj]<0) {
        metab_kx_synth_path[jj]=0
      }
    }
    
    prelimsynthnew_path$metab = metabsynth_path
    prelimsynthnew_path$metab_kx = metab_kx_synth_path
  }
  
  prelimsynth_path = rbind(prelimsynth_path, prelimsynthnew_path)
}

# Convert simulated experimental findings from wide data format to long data format 
data_long <- gather(prelimsynth, condition, conc, metab, metab_kx, factor_key=TRUE)
data_long_path <- gather(prelimsynth_path, condition, conc, metab, metab_kx, factor_key=TRUE)

prelimsynth_btw_group = rbind(prelimsynth, prelimsynth_path)
data_long_btw_group = rbind(data_long, data_long_path)


# Linear model for Aim 1: wt only 
data_long_steady_state =  data_long[which((data_long$condition == 'metab' & data_long$time < 2.9) | (data_long$condition == 'metab_kx' & data_long$time > 6.5)),] # Take only first two isoflurane and last two kx 
data_long_steady_state_means <- aggregate(conc ~ id*condition, data=data_long_steady_state, FUN=mean)

suppressMessages({aim1_lm_sample_synth = lmer(formula = conc ~ 1 + condition + (1|id), data = data_long_steady_state_means)})
aim1_lm_sample_synth_aov = car::Anova(aim1_lm_sample_synth)
aim1_lm_sample_synth_aov_p = aim1_lm_sample_synth_aov[[3]]

p_values[mm, 1] = 0;
if (aim1_lm_sample_synth_aov_p < sig_value){
  p_values[mm, 1] = 1;
}

# Linear model for Aim 2: wt vs. path
data_long_btw_group_steady_state =  data_long_btw_group[which((data_long_btw_group$condition == 'metab' & data_long_btw_group$time < 2.9) | (data_long_btw_group$condition == 'metab_kx' & data_long_btw_group$time > 6.5)),] # Take only first two isoflurane and last two kx 
data_long_btw_group_steady_state_means <- aggregate(conc ~ id*condition*strain, data=data_long_btw_group_steady_state, FUN=mean)

suppressMessages({aim2_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + condition + condition*strain + (1|id), data = data_long_btw_group_steady_state_means)})
aim2_lm_sample_synth_aov = car::Anova(aim2_lm_sample_synth_mixed)
aim2_lm_sample_synth_aov_p = aim2_lm_sample_synth_aov[[3, 3]] 

p_values[mm, 2] = 0;
if (aim2_lm_sample_synth_aov_p < sig_value){
  p_values[mm, 2] = 1;
}

# Linear model for Aim 3: wt vs. path functional 
# Adjust data frame timings so that metab and metab_kx are in series as they will be for the real experiment
data_long_btw_group_time = data_long_btw_group
data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] <- data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] + 11.9

suppressMessages({aim3_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + time + strain*time + (1|id), data = data_long_btw_group_time)})
aim3_lm_sample_synth_mixed_aov = car::Anova(aim3_lm_sample_synth_mixed)
aim3_lm_sample_synth_mixed_aov_p = aim3_lm_sample_synth_mixed_aov[[3, 3]]

p_values[mm, 3] = 0;
if (aim3_lm_sample_synth_mixed_aov_p < sig_value){
  p_values[mm, 3] = 1;
}

}

# Calculate statistical power from model fits to simulated experiments 
power_values[[kk, 1]] = sum(p_values[, 1])/length(p_values[, 1])
power_values[[kk, 2]] = sum(p_values[, 2])/length(p_values[, 2])
power_values[[kk, 3]] = sum(p_values[, 3])/length(p_values[, 3])

}

# Analyze power values matrix to determine appropriate sample sizes
#Aim 1 
number_in_row = 0; 
for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 1]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim1_n = ceiling(kk/margin); 
    break;
  }
}

#Aim 2 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 2]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim2_n = ceiling(kk); 
    break;
  }
}

#Aim 3 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 3]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim3_n = ceiling(kk/margin); 
    break;
  }
}

# Report sample size findings. Only Aim 2 design matters for VR extension.  
print(paste('Aim 1 N for 80% power:', aim1_n))
print(paste('Aim 2 N for 80% power:', aim2_n))
print(paste('Aim 3 N for 80% power:', aim3_n))



#################################### AIM 2 #####################################

# Define basic inputs for experimental simulations 
sig_value = 0.05 # Alpha 
power_value = 0.8 # 1 - beta 
synthesized_data_k = 500 # Number of simulated experiments per N value (higher means more reliable power estimation)
synthesized_data_n_array = 30 # Maximum group N value assessed (should not be too high anyway)
margin = 0.8 # Conservative safeguard for final N, based on experience from participant cancellations
bio_var = 0.1 # Expected between-subject variation based on Swanberg et al., 2021 at 7T

# Initialize array to hold power calculations in 
power_values = data.frame(matrix(nrow = synthesized_data_n_array, ncol = 3)) 

# Calculate statistical power for models fit to experiments simulated at every N value
for (kk in 3:synthesized_data_n_array){
  print(kk);
  synthesized_data_n = kk 
  
  # Initialize array to hold p values for each fit in this round of simulated experiments  
  p_values = data.frame(matrix(nrow = synthesized_data_k, ncol = 3)) 
  
  # It's 1H-MRS o'clock
  for (mm in 1:synthesized_data_k){
    
    ### WILDTYPE 
    # Generate the core data vector based on preliminary empirical findings 
    preliminary_experiment = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    preliminary_experiment$id = c(1, 1, 1, 1, 1, 1, 1)
    preliminary_experiment$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment$metab = c(10, 10, 10, 10, 9, 9, 9) # Based on 10% change in mean seen in Lehmann et al., 2019
    preliminary_experiment$crlb = c(4, 4, 4, 4, 4, 4, 4) # Should not expect much better CRLB than 9.4 T functional data 
    preliminary_experiment$metab_kx = c(9, 9, 9, 9, 9, 9, 9)
    preliminary_experiment$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment$strain = c('wt', 'wt', 'wt', 'wt', 'wt', 'wt', 'wt')
    
    nrow_prelimsynth = length(preliminary_experiment$time)*synthesized_data_n
    prelimsynth = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth) = colnames(preliminary_experiment)
    
    # Generate further synthesized data vectors according to defined biological variability and measured CRLB
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew = preliminary_experiment
      }
      else{
        prelimsynthnew = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew) = colnames(prelimsynth)
        prelimsynthnew$id = c(ii, ii, ii, ii, ii, ii, ii) 
        prelimsynthnew$time = preliminary_experiment$time
        prelimsynthnew$crlb = preliminary_experiment$crlb
        prelimsynthnew$crlb_kx = preliminary_experiment$crlb_kx
        prelimsynthnew$strain = preliminary_experiment$strain
        
        metabsynth = preliminary_experiment$metab
        metab_kx_synth = preliminary_experiment$metab_kx
        
        for (jj in 1:length(preliminary_experiment$metab)){
          metab_orig = metabsynth[jj]; 
          metab_kx_orig = metab_kx_synth[jj]; 
          
          #Add biological variability 
          sd_bio_variability = bio_var*metab_orig 
          metab_orig_mean = rnorm(1, mean = metab_orig, sd = sd_bio_variability)
          sd_bio_kx_variability = bio_var*metab_kx_orig 
          metab_kx_orig_mean = rnorm(1, mean = metab_kx_orig, sd = sd_bio_kx_variability)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb = preliminary_experiment$crlb[jj]
          metab_orig_sd = metab_orig_mean*metab_orig_crlb/100;
          metabsynth[jj] = rnorm(1, mean = metab_orig_mean, sd = metab_orig_sd)
          
          metab_kx_orig_crlb = preliminary_experiment$crlb_kx[jj]
          metab_kx_orig_sd = metab_kx_orig_mean*metab_kx_orig_crlb/100;
          metab_kx_synth[jj] = rnorm(1, mean = metab_kx_orig_mean, sd = metab_kx_orig_sd)
          
          # Control for NA metabolite values
          if (is.na(metabsynth[jj])) {
            metabsynth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth[jj]<0) {
            metabsynth[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth[jj])) {
            metab_kx_synth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth[jj]<0) {
            metab_kx_synth[jj]=0
          }
        }
        
        prelimsynthnew$metab = metabsynth
        prelimsynthnew$metab_kx = metab_kx_synth
      }
      
      prelimsynth = rbind(prelimsynth, prelimsynthnew)
    }
    
    ## AGED
    # Generate the core data vector based on preliminary empirical findings and hypothesis regarding age effects
    preliminary_experiment_path = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment_path) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    idnum_path = synthesized_data_n + 1
    preliminary_experiment_path$id = c(idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path)
    preliminary_experiment_path$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment_path$metab = c(10, 10, 10, 10, 9, 9, 9) # Based on Lehmann et al., 2019
    preliminary_experiment_path$crlb = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment_path$metab_kx = c(9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6) # Intermediate effect for aged
    preliminary_experiment_path$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment_path$strain = c('path', 'path', 'path', 'path', 'path', 'path', 'path')
    
    nrow_prelimsynth_path = length(preliminary_experiment_path$time)*synthesized_data_n
    prelimsynth_path = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth_path) = colnames(preliminary_experiment_path)
    
    # Generate synthesized data vectors according to CRLB measured from preliminary data vector 
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew_path = preliminary_experiment_path
      }
      else{
        prelimsynthnew_path = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew_path) = colnames(prelimsynth_path)
        idnum = ii + synthesized_data_n 
        prelimsynthnew_path$id = c(idnum, idnum, idnum, idnum, idnum, idnum, idnum) 
        prelimsynthnew_path$time = preliminary_experiment_path$time
        prelimsynthnew_path$crlb = preliminary_experiment_path$crlb
        prelimsynthnew_path$crlb_kx = preliminary_experiment_path$crlb_kx
        prelimsynthnew_path$strain = preliminary_experiment_path$strain
        
        metabsynth_path = preliminary_experiment_path$metab
        metab_kx_synth_path = preliminary_experiment_path$metab_kx
        
        for (jj in 1:length(preliminary_experiment_path$metab)){
          metab_orig_path = metabsynth_path[jj]; 
          metab_kx_orig_path = metab_kx_synth_path[jj]; 
          
          #Add biological variability 
          sd_bio_variability_path = bio_var*metab_orig_path 
          metab_orig_mean_path = rnorm(1, mean = metab_orig_path, sd = sd_bio_variability_path)
          sd_bio_kx_variability_path = bio_var*metab_kx_orig_path 
          metab_kx_orig_mean_path = rnorm(1, mean = metab_kx_orig_path, sd = sd_bio_kx_variability_path)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb_path = preliminary_experiment_path$crlb[jj]
          metab_orig_sd_path = metab_orig_mean_path*metab_orig_crlb_path/100;
          metabsynth_path[jj] = rnorm(1, mean = metab_orig_mean_path, sd = metab_orig_sd_path)
          
          metab_kx_orig_crlb_path = preliminary_experiment_path$crlb_kx[jj]
          metab_kx_orig_sd_path = metab_kx_orig_mean_path*metab_kx_orig_crlb_path/100;
          metab_kx_synth_path[jj] = rnorm(1, mean = metab_kx_orig_mean_path, sd = metab_kx_orig_sd_path)
          
          # Control for NA metabolite values
          if (is.na(metabsynth_path[jj])) {
            metabsynth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth_path[jj]<0) {
            metabsynth_path[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth_path[jj])) {
            metab_kx_synth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth_path[jj]<0) {
            metab_kx_synth_path[jj]=0
          }
        }
        
        prelimsynthnew_path$metab = metabsynth_path
        prelimsynthnew_path$metab_kx = metab_kx_synth_path
      }
      
      prelimsynth_path = rbind(prelimsynth_path, prelimsynthnew_path)
    }
    
    # Convert simulated experimental findings from wide data format to long data format 
    data_long <- gather(prelimsynth, condition, conc, metab, metab_kx, factor_key=TRUE)
    data_long_path <- gather(prelimsynth_path, condition, conc, metab, metab_kx, factor_key=TRUE)
    
    prelimsynth_btw_group = rbind(prelimsynth, prelimsynth_path)
    data_long_btw_group = rbind(data_long, data_long_path)
    
    
    # Linear model for Aim 1: wt only 
    data_long_steady_state =  data_long[which((data_long$condition == 'metab' & data_long$time < 2.9) | (data_long$condition == 'metab_kx' & data_long$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_steady_state_means <- aggregate(conc ~ id*condition, data=data_long_steady_state, FUN=mean)
    
    suppressMessages({aim1_lm_sample_synth = lmer(formula = conc ~ 1 + condition + (1|id), data = data_long_steady_state_means)})
    aim1_lm_sample_synth_aov = car::Anova(aim1_lm_sample_synth)
    aim1_lm_sample_synth_aov_p = aim1_lm_sample_synth_aov[[3]]
    
    p_values[mm, 1] = 0;
    if (aim1_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 1] = 1;
    }
    
    # Linear model for Aim 2: wt vs. path
    data_long_btw_group_steady_state =  data_long_btw_group[which((data_long_btw_group$condition == 'metab' & data_long_btw_group$time < 2.9) | (data_long_btw_group$condition == 'metab_kx' & data_long_btw_group$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_btw_group_steady_state_means <- aggregate(conc ~ id*condition*strain, data=data_long_btw_group_steady_state, FUN=mean)
    
    suppressMessages({aim2_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + condition + condition*strain + (1|id), data = data_long_btw_group_steady_state_means)})
    aim2_lm_sample_synth_aov = car::Anova(aim2_lm_sample_synth_mixed)
    aim2_lm_sample_synth_aov_p = aim2_lm_sample_synth_aov[[3, 3]] 
    
    p_values[mm, 2] = 0;
    if (aim2_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 2] = 1;
    }
    
    # Linear model for Aim 3: wt vs. path functional 
    # Adjust data frame timings so that metab and metab_kx are in series as they will be for the real experiment
    data_long_btw_group_time = data_long_btw_group
    data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] <- data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] + 11.9
    
    suppressMessages({aim3_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + time + strain*time + (1|id), data = data_long_btw_group_time)})
    aim3_lm_sample_synth_mixed_aov = car::Anova(aim3_lm_sample_synth_mixed)
    aim3_lm_sample_synth_mixed_aov_p = aim3_lm_sample_synth_mixed_aov[[3, 3]]
    
    p_values[mm, 3] = 0;
    if (aim3_lm_sample_synth_mixed_aov_p < sig_value){
      p_values[mm, 3] = 1;
    }
    
  }
  
  # Calculate statistical power from model fits to simulated experiments 
  power_values[[kk, 1]] = sum(p_values[, 1])/length(p_values[, 1])
  power_values[[kk, 2]] = sum(p_values[, 2])/length(p_values[, 2])
  power_values[[kk, 3]] = sum(p_values[, 3])/length(p_values[, 3])
  
}

# Analyze power values matrix to determine appropriate sample sizes
#Aim 1 
number_in_row = 0; 
for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 1]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim1_n = ceiling(kk/margin); 
    break;
  }
}

#Aim 2 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 2]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim2_n = ceiling(kk); 
    break;
  }
}

#Aim 3 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 3]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim3_n = ceiling(kk/margin); 
    break;
  }
}

# Report sample size findings. Only Aim 2 design matters for VR extension.  
print(paste('Aim 1 N for 80% power:', aim1_n))
print(paste('Aim 2 N for 80% power:', aim2_n))
print(paste('Aim 3 N for 80% power:', aim3_n))




#################################### AIM 3A: Main effect term only #####################################

# Define basic inputs for experimental simulations 
sig_value = 0.05 # Alpha 
power_value = 0.8 # 1 - beta 
synthesized_data_k = 500 # Number of simulated experiments per N value (higher means more reliable power estimation)
synthesized_data_n_array = 15 # Maximum group N value assessed (should not be too high anyway)
margin = 0.8 # Conservative safeguard for final N
bio_var = 0.15 # Expected between-subject variation 

# Initialize array to hold power calculations in 
power_values = data.frame(matrix(nrow = synthesized_data_n_array, ncol = 3)) 

# Calculate statistical power for models fit to experiments simulated at every N value
for (kk in 3:synthesized_data_n_array){
  print(kk);
  synthesized_data_n = kk 
  
  # Initialize array to hold p values for each fit in this round of simulated experiments  
  p_values = data.frame(matrix(nrow = synthesized_data_k, ncol = 3)) 
  
  # It's 1H-MRS o'clock
  for (mm in 1:synthesized_data_k){
    
    ### WILDTYPE 
    # Generate the core data vector based on preliminary empirical findings 
    preliminary_experiment = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    preliminary_experiment$id = c(1, 1, 1, 1, 1, 1, 1)
    preliminary_experiment$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment$metab = c(6.036, 5.814, 5.676, 5.158, 4.843, 4.727, 4.833)
    preliminary_experiment$crlb = c(5, 4, 5, 5, 5, 6, 3)
    preliminary_experiment$metab_kx = c(4.843, 4.843, 5.158, 4.833, 4.833, 4.833, 4.833)
    preliminary_experiment$crlb_kx = c(5, 5, 5, 3, 3, 3, 3)
    preliminary_experiment$strain = c('wt', 'wt', 'wt', 'wt', 'wt', 'wt', 'wt')
    
    nrow_prelimsynth = length(preliminary_experiment$time)*synthesized_data_n
    prelimsynth = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth) = colnames(preliminary_experiment)
    
    # Generate further synthesized data vectors according to defined biological variability and measured CRLB
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew = preliminary_experiment
      }
      else{
        prelimsynthnew = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew) = colnames(prelimsynth)
        prelimsynthnew$id = c(ii, ii, ii, ii, ii, ii, ii) 
        prelimsynthnew$time = preliminary_experiment$time
        prelimsynthnew$crlb = preliminary_experiment$crlb
        prelimsynthnew$crlb_kx = preliminary_experiment$crlb_kx
        prelimsynthnew$strain = preliminary_experiment$strain
        
        metabsynth = preliminary_experiment$metab
        metab_kx_synth = preliminary_experiment$metab_kx
        
        for (jj in 1:length(preliminary_experiment$metab)){
          metab_orig = metabsynth[jj]; 
          metab_kx_orig = metab_kx_synth[jj]; 
          
          #Add biological variability 
          sd_bio_variability = bio_var*metab_orig 
          metab_orig_mean = rnorm(1, mean = metab_orig, sd = sd_bio_variability)
          sd_bio_kx_variability = bio_var*metab_kx_orig 
          metab_kx_orig_mean = rnorm(1, mean = metab_kx_orig, sd = sd_bio_kx_variability)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb = preliminary_experiment$crlb[jj]
          metab_orig_sd = metab_orig_mean*metab_orig_crlb/100;
          metabsynth[jj] = rnorm(1, mean = metab_orig_mean, sd = metab_orig_sd)
          
          metab_kx_orig_crlb = preliminary_experiment$crlb_kx[jj]
          metab_kx_orig_sd = metab_kx_orig_mean*metab_kx_orig_crlb/100;
          metab_kx_synth[jj] = rnorm(1, mean = metab_kx_orig_mean, sd = metab_kx_orig_sd)
          
          # Control for NA metabolite values
          if (is.na(metabsynth[jj])) {
            metabsynth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth[jj]<0) {
            metabsynth[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth[jj])) {
            metab_kx_synth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth[jj]<0) {
            metab_kx_synth[jj]=0
          }
        }
        
        prelimsynthnew$metab = metabsynth
        prelimsynthnew$metab_kx = metab_kx_synth
      }
      
      prelimsynth = rbind(prelimsynth, prelimsynthnew)
    }
    
    ## AGED
    # Generate the core data vector based on preliminary empirical findings and hypothesis regarding age effects
    preliminary_experiment_path = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment_path) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    idnum_path = synthesized_data_n + 1
    preliminary_experiment_path$id = c(idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path)
    preliminary_experiment_path$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment_path$metab = c(6.036, 5.814, 5.676, 5.158, 4.843, 4.727, 4.833)
    preliminary_experiment_path$crlb = c(5, 4, 5, 5, 5, 6, 3)
    preliminary_experiment_path$metab_kx = c(5.814, 5.814, 5.814, 5.814, 5.814, 5.814, 5.814)
    preliminary_experiment_path$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment_path$strain = c('path', 'path', 'path', 'path', 'path', 'path', 'path')
    
    nrow_prelimsynth_path = length(preliminary_experiment_path$time)*synthesized_data_n
    prelimsynth_path = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth_path) = colnames(preliminary_experiment_path)
    
    # Generate synthesized data vectors according to CRLB measured from preliminary data vector 
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew_path = preliminary_experiment_path
      }
      else{
        prelimsynthnew_path = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew_path) = colnames(prelimsynth_path)
        idnum = ii + synthesized_data_n 
        prelimsynthnew_path$id = c(idnum, idnum, idnum, idnum, idnum, idnum, idnum) 
        prelimsynthnew_path$time = preliminary_experiment_path$time
        prelimsynthnew_path$crlb = preliminary_experiment_path$crlb
        prelimsynthnew_path$crlb_kx = preliminary_experiment_path$crlb_kx
        prelimsynthnew_path$strain = preliminary_experiment_path$strain
        
        metabsynth_path = preliminary_experiment_path$metab
        metab_kx_synth_path = preliminary_experiment_path$metab_kx
        
        for (jj in 1:length(preliminary_experiment_path$metab)){
          metab_orig_path = metabsynth_path[jj]; 
          metab_kx_orig_path = metab_kx_synth_path[jj]; 
          
          #Add biological variability 
          sd_bio_variability_path = bio_var*metab_orig_path 
          metab_orig_mean_path = rnorm(1, mean = metab_orig_path, sd = sd_bio_variability_path)
          sd_bio_kx_variability_path = bio_var*metab_kx_orig_path 
          metab_kx_orig_mean_path = rnorm(1, mean = metab_kx_orig_path, sd = sd_bio_kx_variability_path)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb_path = preliminary_experiment_path$crlb[jj]
          metab_orig_sd_path = metab_orig_mean_path*metab_orig_crlb_path/100;
          metabsynth_path[jj] = rnorm(1, mean = metab_orig_mean_path, sd = metab_orig_sd_path)
          
          metab_kx_orig_crlb_path = preliminary_experiment_path$crlb_kx[jj]
          metab_kx_orig_sd_path = metab_kx_orig_mean_path*metab_kx_orig_crlb_path/100;
          metab_kx_synth_path[jj] = rnorm(1, mean = metab_kx_orig_mean_path, sd = metab_kx_orig_sd_path)
          
          # Control for NA metabolite values
          if (is.na(metabsynth_path[jj])) {
            metabsynth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth_path[jj]<0) {
            metabsynth_path[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth_path[jj])) {
            metab_kx_synth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth_path[jj]<0) {
            metab_kx_synth_path[jj]=0
          }
        }
        
        prelimsynthnew_path$metab = metabsynth_path
        prelimsynthnew_path$metab_kx = metab_kx_synth_path
      }
      
      prelimsynth_path = rbind(prelimsynth_path, prelimsynthnew_path)
    }
    
    # Convert simulated experimental findings from wide data format to long data format 
    data_long <- gather(prelimsynth, condition, conc, metab, metab_kx, factor_key=TRUE)
    data_long_path <- gather(prelimsynth_path, condition, conc, metab, metab_kx, factor_key=TRUE)
    
    prelimsynth_btw_group = rbind(prelimsynth, prelimsynth_path)
    data_long_btw_group = rbind(data_long, data_long_path)
    
    
    # Linear model for Aim 1: wt only 
    data_long_steady_state =  data_long[which((data_long$condition == 'metab' & data_long$time < 2.9) | (data_long$condition == 'metab_kx' & data_long$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_steady_state_means <- aggregate(conc ~ id*condition, data=data_long_steady_state, FUN=mean)
    
    suppressMessages({aim1_lm_sample_synth = lmer(formula = conc ~ 1 + condition + (1|id), data = data_long_steady_state_means)})
    aim1_lm_sample_synth_aov = car::Anova(aim1_lm_sample_synth)
    aim1_lm_sample_synth_aov_p = aim1_lm_sample_synth_aov[[3]]
    
    p_values[mm, 1] = 0;
    if (aim1_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 1] = 1;
    }
    
    # Linear model for Aim 2: wt vs. path
    data_long_btw_group_steady_state =  data_long_btw_group[which((data_long_btw_group$condition == 'metab' & data_long_btw_group$time < 2.9) | (data_long_btw_group$condition == 'metab_kx' & data_long_btw_group$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_btw_group_steady_state_means <- aggregate(conc ~ id*condition*strain, data=data_long_btw_group_steady_state, FUN=mean)
    
    suppressMessages({aim2_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + condition + condition*strain + (1|id), data = data_long_btw_group_steady_state_means)})
    aim2_lm_sample_synth_aov = car::Anova(aim2_lm_sample_synth_mixed)
    aim2_lm_sample_synth_aov_p = aim2_lm_sample_synth_aov[[2, 3]] 
    
    p_values[mm, 2] = 0;
    if (aim2_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 2] = 1;
    }
    
    # Linear model for Aim 3: wt vs. path functional 
    # Adjust data frame timings so that metab and metab_kx are in series as they will be for the real experiment
    data_long_btw_group_time = data_long_btw_group
    data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] <- data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] + 11.9
    
    suppressMessages({aim3_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + time + strain*time + (1|id), data = data_long_btw_group_time)})
    aim3_lm_sample_synth_mixed_aov = car::Anova(aim3_lm_sample_synth_mixed)
    aim3_lm_sample_synth_mixed_aov_p = aim3_lm_sample_synth_mixed_aov[[3, 3]]
    
    p_values[mm, 3] = 0;
    if (aim3_lm_sample_synth_mixed_aov_p < sig_value){
      p_values[mm, 3] = 1;
    }
    
  }
  
  # Calculate statistical power from model fits to simulated experiments 
  power_values[[kk, 1]] = sum(p_values[, 1])/length(p_values[, 1])
  power_values[[kk, 2]] = sum(p_values[, 2])/length(p_values[, 2])
  power_values[[kk, 3]] = sum(p_values[, 3])/length(p_values[, 3])
  
}

# Analyze power values matrix to determine appropriate sample sizes
#Aim 1 
number_in_row = 0; 
for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 1]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim1_n = ceiling(kk/margin); 
    break;
  }
}

#Aim 2 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 2]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim2_n = ceiling(kk); 
    break;
  }
}

#Aim 3 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 3]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim3_n = ceiling(kk/margin); 
    break;
  }
}

# Report sample size findings. Only Aim 2 design matters for VR extension.  
print(paste('Aim 1 N for 80% power:', aim1_n))
print(paste('Aim 2 N for 80% power:', aim2_n))
print(paste('Aim 3 N for 80% power:', aim3_n))


#################################### AIM 3B: Main effect term only #####################################

# Define basic inputs for experimental simulations 
sig_value = 0.05 # Alpha 
power_value = 0.8 # 1 - beta 
synthesized_data_k = 500 # Number of simulated experiments per N value (higher means more reliable power estimation)
synthesized_data_n_array = 30 # Maximum group N value assessed (should not be too high anyway)
margin = 0.85 # Conservative safeguard for final N, based on experience from participant cancellations
bio_var = 0.1 # Expected between-subject variation based on Swanberg et al., 2021 at 7T

# Initialize array to hold power calculations in 
power_values = data.frame(matrix(nrow = synthesized_data_n_array, ncol = 3)) 

# Calculate statistical power for models fit to experiments simulated at every N value
for (kk in 3:synthesized_data_n_array){
  print(kk);
  synthesized_data_n = kk 
  
  # Initialize array to hold p values for each fit in this round of simulated experiments  
  p_values = data.frame(matrix(nrow = synthesized_data_k, ncol = 3)) 
  
  # It's 1H-MRS o'clock
  for (mm in 1:synthesized_data_k){
    
    ### WILDTYPE 
    # Generate the core data vector based on preliminary empirical findings 
    preliminary_experiment = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    preliminary_experiment$id = c(1, 1, 1, 1, 1, 1, 1)
    preliminary_experiment$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment$metab = c(10, 10, 10, 10, 9, 9, 9) # Based on 10% change in mean seen in Lehmann et al., 2019
    preliminary_experiment$crlb = c(4, 4, 4, 4, 4, 4, 4) # Should not expect much better CRLB than 9.4 T functional data 
    preliminary_experiment$metab_kx = c(9, 9, 9, 9, 9, 9, 9)
    preliminary_experiment$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment$strain = c('wt', 'wt', 'wt', 'wt', 'wt', 'wt', 'wt')
    
    nrow_prelimsynth = length(preliminary_experiment$time)*synthesized_data_n
    prelimsynth = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth) = colnames(preliminary_experiment)
    
    # Generate further synthesized data vectors according to defined biological variability and measured CRLB
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew = preliminary_experiment
      }
      else{
        prelimsynthnew = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew) = colnames(prelimsynth)
        prelimsynthnew$id = c(ii, ii, ii, ii, ii, ii, ii) 
        prelimsynthnew$time = preliminary_experiment$time
        prelimsynthnew$crlb = preliminary_experiment$crlb
        prelimsynthnew$crlb_kx = preliminary_experiment$crlb_kx
        prelimsynthnew$strain = preliminary_experiment$strain
        
        metabsynth = preliminary_experiment$metab
        metab_kx_synth = preliminary_experiment$metab_kx
        
        for (jj in 1:length(preliminary_experiment$metab)){
          metab_orig = metabsynth[jj]; 
          metab_kx_orig = metab_kx_synth[jj]; 
          
          #Add biological variability 
          sd_bio_variability = bio_var*metab_orig 
          metab_orig_mean = rnorm(1, mean = metab_orig, sd = sd_bio_variability)
          sd_bio_kx_variability = bio_var*metab_kx_orig 
          metab_kx_orig_mean = rnorm(1, mean = metab_kx_orig, sd = sd_bio_kx_variability)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb = preliminary_experiment$crlb[jj]
          metab_orig_sd = metab_orig_mean*metab_orig_crlb/100;
          metabsynth[jj] = rnorm(1, mean = metab_orig_mean, sd = metab_orig_sd)
          
          metab_kx_orig_crlb = preliminary_experiment$crlb_kx[jj]
          metab_kx_orig_sd = metab_kx_orig_mean*metab_kx_orig_crlb/100;
          metab_kx_synth[jj] = rnorm(1, mean = metab_kx_orig_mean, sd = metab_kx_orig_sd)
          
          # Control for NA metabolite values
          if (is.na(metabsynth[jj])) {
            metabsynth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth[jj]<0) {
            metabsynth[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth[jj])) {
            metab_kx_synth[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth[jj]<0) {
            metab_kx_synth[jj]=0
          }
        }
        
        prelimsynthnew$metab = metabsynth
        prelimsynthnew$metab_kx = metab_kx_synth
      }
      
      prelimsynth = rbind(prelimsynth, prelimsynthnew)
    }
    
    ## AGED
    # Generate the core data vector based on preliminary empirical findings and hypothesis regarding age effects
    preliminary_experiment_path = as.data.frame(matrix(nrow=7,ncol=7))
    colnames(preliminary_experiment_path) = c('id', 'time', 'metab', 'crlb', 'metab_kx', 'crlb_kx', 'strain')
    idnum_path = synthesized_data_n + 1
    preliminary_experiment_path$id = c(idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path, idnum_path)
    preliminary_experiment_path$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
    preliminary_experiment_path$metab = c(10, 10, 10, 10, 9, 9, 9) # Based on Lehmann et al., 2019
    preliminary_experiment_path$crlb = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment_path$metab_kx = c(9.6, 9.6, 9.6, 9.6, 9.6, 9.6, 9.6) # Intermediate effect for aged
    preliminary_experiment_path$crlb_kx = c(4, 4, 4, 4, 4, 4, 4)
    preliminary_experiment_path$strain = c('path', 'path', 'path', 'path', 'path', 'path', 'path')
    
    nrow_prelimsynth_path = length(preliminary_experiment_path$time)*synthesized_data_n
    prelimsynth_path = as.data.frame(matrix(nrow=0,ncol=7))
    colnames(prelimsynth_path) = colnames(preliminary_experiment_path)
    
    # Generate synthesized data vectors according to CRLB measured from preliminary data vector 
    for (ii in 1:synthesized_data_n) {
      
      if (ii==1){
        prelimsynthnew_path = preliminary_experiment_path
      }
      else{
        prelimsynthnew_path = as.data.frame(matrix(nrow=7,ncol=7))
        colnames(prelimsynthnew_path) = colnames(prelimsynth_path)
        idnum = ii + synthesized_data_n 
        prelimsynthnew_path$id = c(idnum, idnum, idnum, idnum, idnum, idnum, idnum) 
        prelimsynthnew_path$time = preliminary_experiment_path$time
        prelimsynthnew_path$crlb = preliminary_experiment_path$crlb
        prelimsynthnew_path$crlb_kx = preliminary_experiment_path$crlb_kx
        prelimsynthnew_path$strain = preliminary_experiment_path$strain
        
        metabsynth_path = preliminary_experiment_path$metab
        metab_kx_synth_path = preliminary_experiment_path$metab_kx
        
        for (jj in 1:length(preliminary_experiment_path$metab)){
          metab_orig_path = metabsynth_path[jj]; 
          metab_kx_orig_path = metab_kx_synth_path[jj]; 
          
          #Add biological variability 
          sd_bio_variability_path = bio_var*metab_orig_path 
          metab_orig_mean_path = rnorm(1, mean = metab_orig_path, sd = sd_bio_variability_path)
          sd_bio_kx_variability_path = bio_var*metab_kx_orig_path 
          metab_kx_orig_mean_path = rnorm(1, mean = metab_kx_orig_path, sd = sd_bio_kx_variability_path)
          
          # Add measurement variability according to CRLB
          metab_orig_crlb_path = preliminary_experiment_path$crlb[jj]
          metab_orig_sd_path = metab_orig_mean_path*metab_orig_crlb_path/100;
          metabsynth_path[jj] = rnorm(1, mean = metab_orig_mean_path, sd = metab_orig_sd_path)
          
          metab_kx_orig_crlb_path = preliminary_experiment_path$crlb_kx[jj]
          metab_kx_orig_sd_path = metab_kx_orig_mean_path*metab_kx_orig_crlb_path/100;
          metab_kx_synth_path[jj] = rnorm(1, mean = metab_kx_orig_mean_path, sd = metab_kx_orig_sd_path)
          
          # Control for NA metabolite values
          if (is.na(metabsynth_path[jj])) {
            metabsynth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metabsynth_path[jj]<0) {
            metabsynth_path[jj]=0
          }
          
          # Control for NA metabolite values
          if (is.na(metab_kx_synth_path[jj])) {
            metab_kx_synth_path[jj]=0
          }
          
          # Control for negative metabolite values, which are not output by LCModel
          if (metab_kx_synth_path[jj]<0) {
            metab_kx_synth_path[jj]=0
          }
        }
        
        prelimsynthnew_path$metab = metabsynth_path
        prelimsynthnew_path$metab_kx = metab_kx_synth_path
      }
      
      prelimsynth_path = rbind(prelimsynth_path, prelimsynthnew_path)
    }
    
    # Convert simulated experimental findings from wide data format to long data format 
    data_long <- gather(prelimsynth, condition, conc, metab, metab_kx, factor_key=TRUE)
    data_long_path <- gather(prelimsynth_path, condition, conc, metab, metab_kx, factor_key=TRUE)
    
    prelimsynth_btw_group = rbind(prelimsynth, prelimsynth_path)
    data_long_btw_group = rbind(data_long, data_long_path)
    
    
    # Linear model for Aim 1: wt only 
    data_long_steady_state =  data_long[which((data_long$condition == 'metab' & data_long$time < 2.9) | (data_long$condition == 'metab_kx' & data_long$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_steady_state_means <- aggregate(conc ~ id*condition, data=data_long_steady_state, FUN=mean)
    
    suppressMessages({aim1_lm_sample_synth = lmer(formula = conc ~ 1 + condition + (1|id), data = data_long_steady_state_means)})
    aim1_lm_sample_synth_aov = car::Anova(aim1_lm_sample_synth)
    aim1_lm_sample_synth_aov_p = aim1_lm_sample_synth_aov[[3]]
    
    p_values[mm, 1] = 0;
    if (aim1_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 1] = 1;
    }
    
    # Linear model for Aim 2: wt vs. path
    data_long_btw_group_steady_state =  data_long_btw_group[which((data_long_btw_group$condition == 'metab' & data_long_btw_group$time < 2.9) | (data_long_btw_group$condition == 'metab_kx' & data_long_btw_group$time > 6.5)),] # Take only first two isoflurane and last two kx 
    data_long_btw_group_steady_state_means <- aggregate(conc ~ id*condition*strain, data=data_long_btw_group_steady_state, FUN=mean)
    
    suppressMessages({aim2_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + condition + condition*strain + (1|id), data = data_long_btw_group_steady_state_means)})
    aim2_lm_sample_synth_aov = car::Anova(aim2_lm_sample_synth_mixed)
    aim2_lm_sample_synth_aov_p = aim2_lm_sample_synth_aov[[2, 3]] # Sleep-wake effect only 
    
    p_values[mm, 2] = 0;
    if (aim2_lm_sample_synth_aov_p < sig_value){
      p_values[mm, 2] = 1;
    }
    
    # Linear model for Aim 3: wt vs. path functional 
    # Adjust data frame timings so that metab and metab_kx are in series as they will be for the real experiment
    data_long_btw_group_time = data_long_btw_group
    data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] <- data_long_btw_group_time$time[data_long_btw_group_time$condition == 'metab_kx'] + 11.9
    
    suppressMessages({aim3_lm_sample_synth_mixed <- lmer(formula = conc ~ 1 + strain + time + strain*time + (1|id), data = data_long_btw_group_time)})
    aim3_lm_sample_synth_mixed_aov = car::Anova(aim3_lm_sample_synth_mixed)
    aim3_lm_sample_synth_mixed_aov_p = aim3_lm_sample_synth_mixed_aov[[3, 3]]
    
    p_values[mm, 3] = 0;
    if (aim3_lm_sample_synth_mixed_aov_p < sig_value){
      p_values[mm, 3] = 1;
    }
    
  }
  
  # Calculate statistical power from model fits to simulated experiments 
  power_values[[kk, 1]] = sum(p_values[, 1])/length(p_values[, 1])
  power_values[[kk, 2]] = sum(p_values[, 2])/length(p_values[, 2])
  power_values[[kk, 3]] = sum(p_values[, 3])/length(p_values[, 3])
  
}

# Analyze power values matrix to determine appropriate sample sizes
#Aim 1 
number_in_row = 0; 
for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 1]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim1_n = ceiling(kk/margin); 
    break;
  }
}

#Aim 2 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 2]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim2_n = ceiling(kk/margin); 
    break;
  }
}

#Aim 3 
number_in_row = 0; 

for (kk in 3:synthesized_data_n_array){
  print(kk)
  if (power_values[[kk, 3]] > power_value) {
    number_in_row_old = number_in_row; 
    number_in_row = number_in_row + 1; 
  }
  else {
    number_in_row = 0; # Restart counting 
  }
  if (number_in_row == 2){
    aim3_n = ceiling(kk/margin); 
    break;
  }
}

# Report sample size findings. Only Aim 2 design matters for VR extension.  
print(paste('Aim 1 N for 80% power:', aim1_n))
print(paste('Aim 2 N for 80% power:', aim2_n))
print(paste('Aim 3 N for 80% power:', aim3_n))