library(pwr)

# Define basic inputs for all power analyses 
sig_value = 0.05
power_value = 0.8
v_value = NULL
synthesized_data_n = 10000 
margin = 0.6

# Fit a linear model to the preliminary data vector 
preliminary_experiment = as.data.frame(matrix(nrow=7,ncol=4))
colnames(preliminary_experiment) = c('id', 'time', 'metab', 'crlb')
preliminary_experiment$id = c(1, 1, 1, 1, 1, 1, 1)
preliminary_experiment$time = c(0, 1.1, 2.9, 4.7, 6.5, 8.3, 10.1)
preliminary_experiment$metab = c(0.48305201, 0.188184482, 0.919455138, 0.503200834, 1.926986399, 2.777452416, 4.165521454)
preliminary_experiment$crlb = c(212, 212, 35, 67, 17, 12, 11)

#lm_sample = lm(preliminary_experiment$metab ~ preliminary_experiment$time)
#summary(lm_sample)

nrow_prelimsynth = length(preliminary_experiment$time)*synthesized_data_n
prelimsynth = as.data.frame(matrix(nrow=0,ncol=4))
colnames(prelimsynth) = colnames(preliminary_experiment)

# Generate synthesized data vectors according to CRLB measured from preliminary data vector 
for (ii in 1:synthesized_data_n) {
  print(ii); 
  
  if (ii==1){
    prelimsynthnew = preliminary_experiment
  }
  else{
    prelimsynthnew = as.data.frame(matrix(nrow=7,ncol=4))
    colnames(prelimsynthnew) = colnames(prelimsynth)
    prelimsynthnew$id = c(ii, ii, ii, ii, ii, ii, ii) 
    prelimsynthnew$time = preliminary_experiment$time
    prelimsynthnew$crlb = preliminary_experiment$crlb
  
    metabsynth = preliminary_experiment$metab
    
    for (jj in 1:length(preliminary_experiment$metab)){
      metab_orig = metabsynth[jj]; 
      
      #Add biological variability 
      sd_bio_variability = 0.25*metab_orig #Arbitrary 25% 
      metab_orig_mean = rnorm(1, mean = metab_orig, sd = sd_bio_variability)
      
      # Add measurement variability according to CRLB
      metab_orig_crlb = preliminary_experiment$crlb[jj]
      metab_orig_sd = metab_orig_mean*metab_orig_crlb/100;
      metabsynth[jj] = rnorm(1, mean = metab_orig_mean, sd = metab_orig_sd)
      
      # Control for NA metabolite values
      if (is.na(metabsynth[jj])) {
        metabsynth[jj]=0
      }
      
      # Control for negative metabolite values, which are not output by LCModel
      if (metabsynth[jj]<0) {
        metabsynth[jj]=0
      }
    }
    
    prelimsynthnew$metab = metabsynth
  }
  
  prelimsynth = rbind(prelimsynth, prelimsynthnew)
}

# Fit a linear model to synthesized data vectors
lm_sample_synth = lm(prelimsynth$metab ~ prelimsynth$time)
summary(lm_sample_synth)

lm_sample_synth_r2 = summary(lm_sample_synth)$adj.r.squared
f2_value = lm_sample_synth_r2/(1-lm_sample_synth_r2)
  
# Aim 1 
u_value = 3 - 1 # anesthesia + random subject intercept + offset
v_Aim1 <- pwr.f2.test(u = u_value, v = NULL, f2 = f2_value, power = power_value)
group_sample_size_Aim1 = ceiling((v_Aim1$v + u_value + 1)/(2*margin)) 

# Aim 2
u_value = 5 - 1 # anesthesia + genotype + anesthesia x genotype + random subject intercept + offset
v_Aim2 <- pwr.f2.test(u = u_value, v = NULL, f2 = f2_value, power = power_value)
group_sample_size_Aim2 = ceiling((v_Aim2$v + u_value + 1)/(2*margin))

# Aim 3
u_value = 5 - 1 # time + genotype + time x genotype + random subject intercept + offset
v_Aim3 <- pwr.f2.test(u = u_value, v = NULL, f2 = f2_value, power = power_value)
group_sample_size_Aim3 = ceiling((v_Aim3$v + u_value + 1)/(2*margin))