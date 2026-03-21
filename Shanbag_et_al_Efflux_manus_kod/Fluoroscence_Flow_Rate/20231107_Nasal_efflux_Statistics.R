library('tidyverse')
library('ggplot2')
library('dplyr')
library('ggpubr')
library('see')
library('cowplot')
library('rstatix')
library('R.utils')
library('psych')
library('nlme')
library('purrr')
library('multcomp')

# Create folder for analysis outputs 
directory_name = 'Exp_1' # Should change for each experiment type 
dir.create(directory_name)

# Read in experiment data file 
setwd('C://Users//em6050jo//Documents//Lundgaard_Labb_Swanberg//Efflux_Paper//Data//For_Manuscript//Excel_Files//Fluoroscopy_Flow_Rate//')
data <- read.csv(file = 'C://Users//em6050jo//Documents//Lundgaard_Labb_Swanberg//Efflux_Paper//Data//For_Manuscript//Excel_Files//Fluoroscopy_Flow_Rate//20241107_Static_Fluoroscopy_KMS_Exp_1.csv', header=TRUE)


# PART I ################################ DATA PREPROCESSING ################################# 

# Summarize group statistics for all variables 
data_by_group <- data %>%
                 group_by(Group) %>%
                 summarize(across(where(is.numeric), list(mean=mean, sd=sd, count= ~ n(), se =  ~ sd(.x) / sqrt(n()))),  .groups = 'drop')

data_by_group = arrange(data_by_group, desc(Group))
data$Group <- as.factor(data$Group)
#data$Group <- factor(data$Group, levels=rev(levels(data$Group)))

# PART IIA ################################ DATA PLOTTING ################################# 

num_columns = ncol(data)
data_for_plotting = data 

for (i in 9:num_columns){
   
     data_for_plotting$Group = as.factor(data_for_plotting$Group)
     data_for_plotting$y = data[,i]
     ylabelplot = paste(colnames(data)[i], "\n")
     ylabelplot = gsub("_", " ", ylabelplot)
     ylabel = colnames(data)[i]
     
     # Plot NTe enhancement volume 
     enhancement_volume_hist <- ggdensity(data_for_plotting, x = "y", xlab ="\nGroup", ylab =ylabelplot,
                                add = "mean", rug = TRUE,
                                color = "Group", fill = "Group",
                                palette = c("#00AFBB", "#E7B800", '#A4115c')) + 
                                theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))
  
    filename_png = paste0(directory_name,'//', ylabel, '_hist.png') 
    ggsave(filename_png, enhancement_volume_hist, device = 'png', width = 6, height = 6, dpi = 600)
    filename_eps = paste0(directory_name,'//', ylabel, '_hist.pdf') 
    ggsave(filename_eps, plot = print(enhancement_volume_hist), device = pdf, width = 6, height = 6, dpi = 600)
  
    enhancement_volume_violin <- ggviolin(data_for_plotting, x = "Group", y = "y", fill = "Group", xlab ="\nGroup", ylab =ylabelplot)+ 
                                 geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                                 theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))
  
    filename_png = paste0(directory_name,'//', ylabel, '_violin.png') 
    ggsave(filename_png, enhancement_volume_violin, device = 'png', width = 6, height = 6, dpi = 600)
    filename_eps = paste0(directory_name,'//', ylabel, '_violin.pdf') 
    ggsave(filename_eps, plot = print(enhancement_volume_violin), device = pdf, width = 6, height = 6, dpi = 600)  
}

# Plot select figures for grid 
NT_enhancement_violin <- ggviolin(data, x = "Group", y = "Nasal_40min_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab ="NT mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

OFC_enhancement_violin <- ggviolin(data, x = "Group", y = "Olfactory_Cistern_40min_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab = "OFC mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

Skullcap_Dura_enhancement_violin <- ggviolin(data, x = "Group", y = "Skullcap_Dura_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab = "Skullcap and dura mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

Isolated_Dura_enhancement_violin <- ggviolin(data, x = "Group", y = "Dura_Isolated_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab = "Isolated dura mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

SCLNs_enhancement_violin <- ggviolin(data, x = "Group", y = "SCLNsPostPFA_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab = "SCLN mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

DCLNs_enhancement_violin <- ggviolin(data, x = "Group", y = "DCLNsPostPFA_Mean_Intensity_12b", fill = "Group", xlab ="\nGroup", ylab = "DCLN mean pixel intensity (a.u.)", legend="none")+ 
                         geom_dotplot(binaxis='y', stackdir='center', dotsize=1.5)
                         theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))
                         
                         
## Plot final grids for manuscript 
grid_fluoroscopy_enhancement <- plot_grid(NT_enhancement_violin, OFC_enhancement_violin, Skullcap_Dura_enhancement_violin, Isolated_Dura_enhancement_violin, SCLNs_enhancement_violin, DCLNs_enhancement_violin, ncol=3, align="v")
filename_png = paste0(directory_name,'//Master_Grid_Fluoroscopy_Enhancement.png') 
ggsave(filename_png, grid_fluoroscopy_enhancement, device = 'png', width = 9, height = 9, dpi = 600)
filename_eps = paste0(directory_name,'//Master_Grid_Fluoroscopy_Enhancement.pdf') 
ggsave(filename_eps, plot = print(grid_fluoroscopy_enhancement), device = pdf, width = 9, height = 9, dpi = 600)

# PART IIIA ################################ INFERENTIAL STATISTICS ON PREPROCESSED DATA ################################# 

#Define filename
filename_txt = paste0(directory_name,'//Final_Inferential_Statistics.txt') 
filename_txt = file(filename_txt, 'w')

############################################# BETWEEN-GROUP STATISTICS #############################################

for (i in 9:num_columns){
  data_for_plotting$Group = as.factor(data_for_plotting$Group)
  data_for_plotting$y = data[,i]
  ylabelplot = paste(colnames(data)[i], "\n")
  ylabelplot = gsub("_", " ", ylabelplot)
  ylabel = colnames(data)[i]
  
  # Plot NTe enhancement volume 
  desc_mean_enhancement <- describeBy(y ~ Group, data = data_for_plotting, mat = TRUE)  
  aov_mean_enhancement  = 0
  aov_mean_enhancement <- try(aov(data = data_for_plotting, formula = y ~ Group), silent = TRUE)
  summary(aov_mean_enhancement)
  aov_mean_enhancement_shapiro_residuals = 0 
  aov_mean_enhancement_shapiro_residuals = try(shapiro.test(residuals(aov_mean_enhancement)), silent = TRUE)
  kw_mean_enhancement = kruskal.test(y ~ Group, data_for_plotting)

  
  glht_posthoc <-glht(aov_mean_enhancement, linfct = mcp(Group ="Tukey"))
  glht_posthoc_reported <- summary(glht_posthoc, test = adjusted("none")) 
  
  posthoc_parametric <- pairwise.t.test(data_for_plotting$y, data_for_plotting$Group, paired = FALSE, p.adjust.method = "BH")
  posthoc_nonparametric <- pairwise.wilcox.test(data_for_plotting$y, data_for_plotting$Group, paired = FALSE, p.adjust.method = "BH")

  
  cat("############################################# BETWEEN-GROUP STATISTICS #############################################\n", file = filename_txt, append=TRUE)
  captureOutput(ylabel, file = filename_txt, append = TRUE)
  cat("############################################# ######################## #############################################\n", file = filename_txt, append=TRUE)
  cat("\n\n Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(desc_mean_enhancement, file = filename_txt, append = TRUE)
  cat("\n\n One-Way ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(summary(aov_mean_enhancement), file = filename_txt, append = TRUE)
  cat("\n\n One-Way ANOVA Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(aov_mean_enhancement_shapiro_residuals, file = filename_txt, append = TRUE)
  cat("\n\n Kruskal-Wallis Test\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(kw_mean_enhancement, file = filename_txt, append = TRUE)
  cat("\n\n Parametric post hoc pairwise comparisons (Tukey-corrected p-values)\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(glht_posthoc_reported, file = filename_txt, append = TRUE)
  cat("\n\n Parametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values)\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_parametric, file = filename_txt, append = TRUE)
  cat("\n\n Nonparametric post hoc pairwise comparisons (Benjamini-Hochberg-corrected p-values)\n", file = filename_txt, append=TRUE, sep = "\n")
  captureOutput(posthoc_nonparametric, file = filename_txt, append = TRUE)
}
close(filename_txt)