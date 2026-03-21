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
library("car")
library("lme4")

# Create folder for analysis outputs 
directory_name = 'NTe' # Should change for each experiment type 
dir.create(directory_name)

# Read in experiment data file 
setwd('C://Users//em6050jo//Documents//Lundgaard_Labb_Swanberg//Efflux_Paper//Data//For_Manuscript//Excel_Files//Old_vs_Young')
data <- read.csv(file = 'C://Users//em6050jo//Documents//Lundgaard_Labb_Swanberg//Efflux_Paper//Data//For_Manuscript//Excel_Files//Old_vs_Young//20240530_old_vs_young_NTe_master.csv', header=TRUE, sep=',')
data$Sum_Normalized <- data$Sum / data$Num_Voxels; 


# PART I ################################ DATA PREPROCESSING ################################# 

# Calculate area under NTe enhancement curve as well as slope for each mouse at each time point 
for (row in 1:nrow(data)){
  if (data$Time_Step[row]==75){
    data$Cumulative_Sum[row] = 75*data$Sum[row]
    data$Sum_Slope[row] = data$Sum[row]/75
    data$Cumulative_Sum_Normalized[row] = 75*data$Sum_Normalized[row]
    data$Sum_Normalized_Slope[row] = data$Sum_Normalized[row]/75
  }
  else{
    data$Cumulative_Sum[row] = data$Cumulative_Sum[row-1] + (data$Time_Step[row] - data$Time_Step[row-1]) * data$Sum[row]
    data$Sum_Slope[row] = (data$Sum[row] - data$Sum[row-1])/ (data$Time_Step[row] - data$Time_Step[row-1])
    data$Cumulative_Sum_Normalized[row] = data$Cumulative_Sum_Normalized[row-1] + (data$Time_Step[row] - data$Time_Step[row-1]) * data$Sum_Normalized[row]
    data$Sum_Normalized_Slope[row] = (data$Sum_Normalized[row] - data$Sum_Normalized[row-1])/ (data$Time_Step[row] - data$Time_Step[row-1])
    }
}

filename_csv = paste0(directory_name,'//Data_with_Sums_and_Slopes.csv') 
write.csv(data, filename_csv)

data_by_group <- data %>%
                    group_by(Group, Time_Step) %>%
                    summarize(Num_Voxels_Mean = mean(Num_Voxels), 
                              Mean_Mean = mean(Mean), 
                              Min_Mean = mean(Min), 
                              Max_Mean = mean(Max), 
                              Sum_Mean = mean(Sum), 
                              Standard_Deviation_Mean = mean(Standard_Deviation), 
                              Sum_Normalized_Mean = mean(Sum_Normalized),
                              Cumulative_Sum_Mean = mean(Cumulative_Sum),
                              Cumulative_Sum_Normalized_Mean = mean(Cumulative_Sum_Normalized),
                              Sum_Slope_Mean = mean(Sum_Slope),
                              Sum_Normalized_Slope_Mean = mean(Sum_Normalized_Slope),
                              Num_Voxels_N = n(), 
                              Mean_N = n(), 
                              Min_N = n(), 
                              Max_N = n(), 
                              Sum_N = n(), 
                              Standard_Deviation_N = n(), 
                              Sum_Normalized_N = n(),
                              Cumulative_Sum_N = n(),
                              Cumulative_Sum_Normalized_N = n(),
                              Sum_Slope_N = n(),
                              Sum_Normalized_Slope_N = n(),
                              Num_Voxels_SE = sd(Num_Voxels)/sqrt(Num_Voxels_N), 
                              Mean_SE = sd(Mean)/sqrt(Mean_N), 
                              Min_SE = sd(Min)/sqrt(Min_N), 
                              Max_SE = sd(Max)/sqrt(Max_N), 
                              Sum_SE = sd(Sum)/sqrt(Sum_N), 
                              Standard_Deviation_SE = sd(Standard_Deviation)/sqrt(Sum_Normalized_N), 
                              Sum_Normalized_SE = sd(Sum_Normalized)/sqrt(Sum_Normalized_N), 
                              Cumulative_Sum_SE = sd(Cumulative_Sum)/sqrt(Cumulative_Sum_N), 
                              Cumulative_Sum_Normalized_SE = sd(Cumulative_Sum_Normalized)/sqrt(Cumulative_Sum_Normalized_N), 
                              Sum_Slope_SE = sd(Sum_Slope)/sqrt(Sum_Slope_N), 
                              Sum_Normalized_Slope_SE = sd(Sum_Normalized_Slope)/sqrt(Sum_Normalized_Slope_N)) 

data_by_group = arrange(data_by_group, desc(Group))

data_final_time_step = data[which(data$Time_Step==3838),]
data_final_time_step = arrange(data_final_time_step, desc(Group))

data_max_enhancement <- data %>% 
  group_by(Animal_ID) %>% 
  slice(which.max(Sum)) %>% 
  arrange(desc(Group))

# Calculate slopes before and after max enhancement 
for (row in 1:nrow(data_max_enhancement)){
    analyzed_mouse = data_max_enhancement$Animal_ID[row]
    analyzed_time_step = data_max_enhancement$Time_Step[row]
    analyzed_mouse_data = data[which(data$Animal_ID == analyzed_mouse),]
    analyzed_mouse_final_time_step = max(analyzed_mouse_data$Time_Step)
    analyzed_mouse_final_data = analyzed_mouse_data[which(analyzed_mouse_data$Time_Step == analyzed_mouse_final_time_step),]
    if(data_max_enhancement$Time_Step[row]==0){
      data_max_enhancement$Slope_before_Max_Enhancement[row] =  NaN
    } else{
      data_max_enhancement$Slope_before_Max_Enhancement[row] =  data_max_enhancement$Sum[row] / data_max_enhancement$Time_Step[row]
    }
    if(data_max_enhancement$Time_Step[row]==analyzed_mouse_final_time_step){
      data_max_enhancement$Slope_after_Max_Enhancement[row] =  NaN
    } else {
      data_max_enhancement$Slope_after_Max_Enhancement[row] = (analyzed_mouse_final_data$Sum - data_max_enhancement$Sum[row]) / (analyzed_mouse_final_data$Time_Step - data_max_enhancement$Time_Step[row])
    }
}

data_max_enhancement_normalized <- data %>% 
  group_by(Animal_ID) %>% 
  slice(which.max(Sum_Normalized)) %>% 
  arrange(desc(Group))

# Calculate slopes before and after normalized max enhancement 
for (row in 1:nrow(data_max_enhancement_normalized)){
  analyzed_mouse = data_max_enhancement_normalized$Animal_ID[row]
  analyzed_time_step = data_max_enhancement_normalized$Time_Step[row]
  analyzed_mouse_data = data[which(data$Animal_ID == analyzed_mouse),]
  analyzed_mouse_final_time_step = max(analyzed_mouse_data$Time_Step)
  analyzed_mouse_final_data = analyzed_mouse_data[which(analyzed_mouse_data$Time_Step == analyzed_mouse_final_time_step),]
  if(data_max_enhancement_normalized$Time_Step[row]==0){
    data_max_enhancement_normalized$Slope_before_Max_Enhancement[row] =  NaN
  } else {
    data_max_enhancement_normalized$Slope_before_Max_Enhancement[row] =  data_max_enhancement_normalized$Sum_Normalized[row] / data_max_enhancement_normalized$Time_Step[row]
  }
  if(data_max_enhancement_normalized$Time_Step[row]==analyzed_mouse_final_time_step){
    data_max_enhancement_normalized$Slope_after_Max_Enhancement[row] =  NaN
  } else {
    data_max_enhancement_normalized$Slope_after_Max_Enhancement[row] = (analyzed_mouse_final_data$Sum_Normalized - data_max_enhancement_normalized$Sum_Normalized[row]) / (analyzed_mouse_final_data$Time_Step - data_max_enhancement_normalized$Time_Step[row])
  }
}

# PART IIA ################################ PREPROCESSED DATA PLOTTING ################################# 

# Group average time series: SUM NTe
p_group <- ggplot(data_by_group, aes(x = Time_Step, y = Sum_Mean, color = Group, ymin=Sum_Mean - Sum_SE, ymax=Sum_Mean + Sum_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Total NTe enhancement over time by group",
    x = "\nTime (s)",
    y = "Total NTe enhancement (a.u.)\n",
    color = NULL
  ) +
  #ylim(-15, 400) + 
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Sum_by_Group.png') 
ggsave(filename_png, p_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Sum_by_Group.pdf') 
ggsave(filename_eps, width = 6, height = 6, dpi = 600, plot = print(p_group), device = "pdf")

# Group average time series: Cumulative Sum of SUM NTe
p_sum_group <- ggplot(data_by_group, aes(x = Time_Step, y = Cumulative_Sum_Mean, color = Group, ymin=Cumulative_Sum_Mean - Cumulative_Sum_SE, ymax=Cumulative_Sum_Mean + Cumulative_Sum_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Cumulative NTe enhancement over time by group",
    x = "\nTime (s)",
    y = "NTe enhancement AUC (a.u.)\n",
    color = NULL
  ) +
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Cumulative_Sum_by_Group.png') 
ggsave(filename_png, p_sum_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Cumulative_Sum_by_Group.pdf') 
ggsave(filename_eps, plot = print(p_sum_group), device = pdf, width = 6, height = 6, dpi = 600)

# Group average time series: Slope of SUM NTe
p_slope_group <- ggplot(data_by_group, aes(x = Time_Step, y = Sum_Slope_Mean, color = Group, ymin=Sum_Slope_Mean - Sum_Slope_SE, ymax=Sum_Slope_Mean + Sum_Slope_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Slope of NTe enhancement over time by group",
    x = "\nTime (s)",
    y = expression("NTe enhancement slope (a.u.)"),
    color = NULL
  ) +
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12),  axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Slope_by_Group.png') 
ggsave(filename_png, p_slope_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Slope_by_Group.pdf') 
ggsave(filename_eps, plot = print(p_slope_group), device = pdf, width = 6, height = 6, dpi = 600)

# Plot time to max enhancement
max_enhance_hist <- ggdensity(data_max_enhancement, x = "Time_Step", xlab ="\nGroup", ylab ="Time (s)\n", 
                    add = "mean", rug = TRUE,
                    color = "Group", fill = "Group",
                    palette = c("#00AFBB", "#E7B800"))+ 
                    theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12)) 

filename_png = paste0(directory_name,'//Time_to_Max_Enhancement_Group_Hist.png') 
ggsave(filename_png, max_enhance_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Time_to_Max_Enhancement_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(max_enhance_hist), device = pdf, width = 6, height = 6, dpi = 600)

max_enhance_violin <- ggviolin(data_max_enhancement, x = "Group", y = "Time_Step", fill = "Group", xlab ="\nGroup", ylab ="Time to max NTe enhancement (s)\n", 
                      add = "dotplot", add.params = list(size=3, fill = "black")) + 
                      theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Time_to_Max_Enhancement_Group_Violin.png') 
ggsave(filename_png, max_enhance_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Time_to_Max_Enhancement_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(max_enhance_violin), device = pdf, width = 6, height = 6, dpi = 600)

# Plot slope to max enhancement
max_enhance_slope_hist <- ggdensity(data_max_enhancement, x = "Slope_before_Max_Enhancement", xlab ="\nGroup", ylab ="Slope to max NTe enhancement (a.u.)\n", 
                          add = "mean", rug = TRUE,
                          color = "Group", fill = "Group",
                          palette = c("#00AFBB", "#E7B800"))+ 
                          theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12)) 

filename_png = paste0(directory_name,'//Slope_to_Max_Enhancement_Group_Hist.png') 
ggsave(filename_png, max_enhance_slope_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Slope_to_Max_Enhancement_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(max_enhance_slope_hist), device = pdf, width = 6, height = 6, dpi = 600)

max_enhance_slope_violin <- ggviolin(data_max_enhancement, x = "Group", y = "Slope_before_Max_Enhancement", fill = "Group", xlab ="\nGroup", ylab ="Slope to max NTe enhancement (a.u.)\n", 
                           add = "dotplot", add.params = list(size=3, fill = "black")) + 
                           theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Slope_to_Max_Enhancement_Group_Violin.png') 
ggsave(filename_png, max_enhance_slope_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Slope_to_Max_Enhancement_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(max_enhance_slope_violin), device = pdf, width = 6, height = 6, dpi = 600)


# Plot slope after max enhancement
after_max_enhance_slope_hist <- ggdensity(data_max_enhancement, x = "Slope_after_Max_Enhancement", xlab ="\nGroup", ylab ="Slope after max NTe enhancement (a.u.)\n", 
                                add = "mean", rug = TRUE,
                                color = "Group", fill = "Group",
                                palette = c("#00AFBB", "#E7B800"))+ 
                                theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12)) 

filename_png = paste0(directory_name,'//Slope_after_Max_Enhancement_Group_Hist.png') 
ggsave(filename_png, after_max_enhance_slope_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Slope_after_Max_Enhancement_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(after_max_enhance_slope_hist), device = pdf, width = 6, height = 6, dpi = 600)

after_max_enhance_slope_violin <- ggviolin(data_max_enhancement, x = "Group", y = "Slope_after_Max_Enhancement", fill = "Group", xlab ="\nGroup", ylab ="Slope after max NTe enhancement (a.u.)\n", 
                            add = "dotplot", add.params = list(size=3, fill = "black")) + 
                            theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Slope_after_Max_Enhancement_Group_Violin.png') 
ggsave(filename_png, after_max_enhance_slope_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Slope_after_Max_Enhancement_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(after_max_enhance_slope_violin), device = pdf, width = 6, height = 6, dpi = 600)


# Plot max enhancement values
final_cumulative_enhancement_hist <- ggdensity(data_final_time_step, x = "Cumulative_Sum", xlab ="\nGroup", ylab ="Final NTe enhancement AUC (a.u.)\n",
                                     add = "mean", rug = TRUE,
                                     color = "Group", fill = "Group",
                                     palette = c("#00AFBB", "#E7B800")) + 
                                     theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Final_Cumulative_Enhancement_Group_Hist.png') 
ggsave(filename_png, final_cumulative_enhancement_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Final_Cumulative_Enhancement_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(final_cumulative_enhancement_hist), device = pdf, width = 6, height = 6, dpi = 600)

final_cumulative_enhancement_violin <- ggviolin(data_final_time_step, x = "Group", y = "Cumulative_Sum", fill = "Group", xlab ="\nGroup", ylab ="Final NTe enhancement AUC (a.u.)\n", 
                                       add = "dotplot", add.params = list(size = 3, fill = "black"))+ 
                                       theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Final_Cumulative_Enhancement_Group_Violin.png') 
ggsave(filename_png, final_cumulative_enhancement_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Final_Cumulative_Enhancement_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(final_cumulative_enhancement_violin), device = pdf, width = 6, height = 6, dpi = 600)

# Plot NTe enhancement volume 
enhancement_volume_hist <- ggdensity(data_final_time_step, x = "Num_Voxels", xlab ="\nGroup", ylab ="NTe enhancement volume (voxels)\n",
                           add = "mean", rug = TRUE,
                           color = "Group", fill = "Group",
                           palette = c("#00AFBB", "#E7B800")) + 
                           theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//NTe_Enhancement_Volume_Group_Hist.png') 
ggsave(filename_png, enhancement_volume_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Enhancement_Volume_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(enhancement_volume_hist), device = pdf, width = 6, height = 6, dpi = 600)

enhancement_volume_violin <- ggviolin(data_final_time_step, x = "Group", y = "Num_Voxels", fill = "Group", xlab ="\nGroup", ylab ="NTe enhancement volume (voxels)\n", 
                                       add = "dotplot", add.params = list(size = 3, fill = "black"))+ 
                                       theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//NTe_Enhancement_Volume_Group_Violin.png') 
ggsave(filename_png, enhancement_volume_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Enhancement_Volume_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(enhancement_volume_violin), device = pdf, width = 6, height = 6, dpi = 600)

## Plot final grids for abstract 
grid_time_series <- plot_grid(p_group, p_sum_group, p_slope_group, ncol=1, align="v", labels = c('A', 'B', 'C'), label_size = 20)
filename_png = paste0(directory_name,'//Grid_Time_Series.png') 
ggsave(filename_png, grid_time_series, device = 'png', width = 6, height = 9, dpi = 600)
filename_eps = paste0(directory_name,'//Grid_Time_Series.pdf') 
ggsave(filename_eps, plot = print(grid_time_series), device = pdf, width = 6, height = 9, dpi = 600)

grid_btw_group <- plot_grid(max_enhance_violin, final_cumulative_enhancement_violin, ncol=1, align = "v", labels = c('D', 'E'), label_size = 20) 
filename_png = paste0(directory_name,'//Grid_Btw_Group.png') 
ggsave(filename_png, grid_btw_group, device = 'png', width = 6, height = 9, dpi = 600)
filename_eps = paste0(directory_name,'//Grid_Btw_Group.pdf') 
ggsave(filename_eps, plot = print(grid_btw_group), device = pdf, width = 6, height = 9, dpi = 600)

combined_plot_grid <- plot_grid(grid_time_series, grid_btw_group, ncol = 2, align = "h", rel_widths = c(2, 1))
filename_png = paste0(directory_name,'//Master_Grid__Figure.png') 
ggsave(filename_png, combined_plot_grid, device = 'png', width = 11, height = 10, dpi = 600)
filename_eps = paste0(directory_name,'//Master_Grid__Figure.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid), device = pdf, width = 11, height = 10)

# PART IIB ################################ PREPROCESSED DATA PLOTTING - NORMALIZED TO ENHANCEMENT VOLUME ################################# 

# Group average time series: SUM NTe Normalized
p_group_normalized <- ggplot(data_by_group, aes(x = Time_Step, y = Sum_Normalized_Mean, color = Group, ymin=Sum_Normalized_Mean - Sum_Normalized_SE, ymax=Sum_Normalized_Mean + Sum_Normalized_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Total normalized NTe enhancement over time by group",
    x = "\nTime (s)",
    y = "Total norm. NTe enhancement (a.u.)\n",
    color = NULL
  ) +
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Sum_Normalized_by_Group.png') 
ggsave(filename_png, p_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Sum_Normalized_by_Group.pdf') 
ggsave(filename_eps, width = 6, height = 6, dpi = 600, plot = print(p_group), device = "pdf")

# Group average time series: Cumulative Sum of SUM NTe Normalized
p_sum_group_normalized <- ggplot(data_by_group, aes(x = Time_Step, y = Cumulative_Sum_Normalized_Mean, color = Group, ymin=Cumulative_Sum_Normalized_Mean - Cumulative_Sum_Normalized_SE, ymax=Cumulative_Sum_Normalized_Mean + Cumulative_Sum_Normalized_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Cumulative NTe enhancement normalized over time by group",
    x = "\nTime (s)",
    y = "Norm. NTe enhancement AUC (a.u.)\n",
    color = NULL
  ) +
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Cumulative_Sum_Normalized_by_Group.png') 
ggsave(filename_png, p_sum_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Cumulative_Sum_Normalized_by_Group.pdf') 
ggsave(filename_eps, plot = print(p_sum_group), device = pdf, width = 6, height = 6, dpi = 600)

# Group average time series: Slope of SUM NTe
p_slope_group_normalized <- ggplot(data_by_group, aes(x = Time_Step, y = Sum_Normalized_Slope_Mean, color = Group, ymin=Sum_Normalized_Slope_Mean - Sum_Normalized_Slope_SE, ymax=Sum_Normalized_Slope_Mean + Sum_Normalized_Slope_SE)) +
  geom_ribbon(alpha = 0.15, linetype = 0) + 
  geom_line(size = 0.75) +
  geom_point(size = 2) + 
  labs(
    #title = "Slope of normalized NTe enhancement over time by group",
    x = "\nTime (s)",
    y = expression("Norm. NTe enhancement slope (a.u.)"),
    color = NULL
  ) +
  xlim(0, NA) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12),  axis.text = element_text(size = 12), legend.position = c(0.1, 0.8), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
  scale_color_manual(values = c("#00BFC4", "#F8766D"))

filename_png = paste0(directory_name,'//NTe_Normalized_Slope_by_Group.png') 
ggsave(filename_png, p_slope_group, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//NTe_Normalized_Slope_by_Group.pdf') 
ggsave(filename_eps, plot = print(p_slope_group), device = pdf, width = 6, height = 6, dpi = 600)

# Plot time to max enhancement
max_enhance_hist_normalized <- ggdensity(data_max_enhancement_normalized, x = "Time_Step", xlab ="\nGroup", ylab ="Time (s)\n", 
                               add = "mean", rug = TRUE,
                               color = "Group", fill = "Group",
                               palette = c("#00AFBB", "#E7B800"))+ 
                               theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12)) 

filename_png = paste0(directory_name,'//Time_to_Normalized_Max_Enhancement_Group_Hist.png') 
ggsave(filename_png, max_enhance_hist, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Time_to_Normalized_Max_Enhancement_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(max_enhance_hist), device = pdf, width = 6, height = 6, dpi = 600)

max_enhance_violin_normalized <- ggviolin(data_max_enhancement_normalized, x = "Group", y = "Time_Step", fill = "Group", xlab ="\nGroup", ylab ="Time to max norm. NTe enhancement (s)\n", 
                                 add = "dotplot", add.params = list(size=3, fill = "black")) + 
                                 theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Time_to_Normalized_Max_Enhancement_Group_Violin.png') 
ggsave(filename_png, max_enhance_violin, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Time_to_Normalized_Max_Enhancement_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(max_enhance_violin), device = pdf, width = 6, height = 6, dpi = 600)

# Plot max enhancement values
final_cumulative_enhancement_hist_normalized <- ggdensity(data_final_time_step, x = "Cumulative_Sum_Normalized", xlab ="\nGroup", ylab ="Final normalized NTe enhancement AUC (a.u.)\n",
                                                add = "mean", rug = TRUE,
                                                color = "Group", fill = "Group",
                                                palette = c("#00AFBB", "#E7B800")) + 
                                                theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Final_Cumulative_Enhancement_Normalized_Group_Hist.png') 
ggsave(filename_png, final_cumulative_enhancement_hist_normalized, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Final_Cumulative_Enhancement_Normalized_Group_Hist.pdf') 
ggsave(filename_eps, plot = print(final_cumulative_enhancement_hist_normalized), device = pdf, width = 6, height = 6, dpi = 600)

final_cumulative_enhancement_violin_normalized <- ggviolin(data_final_time_step, x = "Group", y = "Cumulative_Sum_Normalized", fill = "Group", xlab ="\nGroup", ylab ="Final normalized NTe enhancement AUC (a.u.)\n", 
                                                  add = "dotplot", add.params = list(size = 3, fill = "black"))+ 
                                                  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), text = element_text(size = 12), axis.text = element_text(size = 12))

filename_png = paste0(directory_name,'//Final_Cumulative_Enhancement_Normalized_Group_Violin.png') 
ggsave(filename_png, final_cumulative_enhancement_violin_normalized, device = 'png', width = 6, height = 6, dpi = 600)
filename_eps = paste0(directory_name,'//Final_Cumulative_Enhancement_Normalized_Group_Violin.pdf') 
ggsave(filename_eps, plot = print(final_cumulative_enhancement_violin_normalized), device = pdf, width = 6, height = 6, dpi = 600)

## Plot final grids for abstract 
grid_time_series_normalized <- plot_grid(p_group_normalized, p_sum_group_normalized, p_slope_group_normalized, ncol=1, align="v", labels = c('A', 'B', 'C'), label_size = 20)
filename_png = paste0(directory_name,'//Grid_Time_Series_Normalized.png') 
ggsave(filename_png, grid_time_series_normalized, device = 'png', width = 6, height = 9, dpi = 600)
filename_eps = paste0(directory_name,'//Grid_Time_Series_Normalized.pdf') 
ggsave(filename_eps, plot = print(grid_time_series_normalized), device = pdf, width = 6, height = 9, dpi = 600)

grid_btw_group_normalized <- plot_grid(max_enhance_violin_normalized, final_cumulative_enhancement_violin_normalized, ncol=1, align = "v", labels = c('D', 'E'), label_size = 20) 
filename_png = paste0(directory_name,'//Grid_Btw_Group_Normalized.png') 
ggsave(filename_png, grid_btw_group_normalized, device = 'png', width = 6, height = 9, dpi = 600)
filename_eps = paste0(directory_name,'//Grid_Btw_Group_Normalized.pdf') 
ggsave(filename_eps, plot = print(grid_btw_group_normalized), device = pdf, width = 6, height = 9, dpi = 600)

combined_plot_grid_normalized <- plot_grid(grid_time_series_normalized, grid_btw_group_normalized, ncol = 2, align = "h", rel_widths = c(2, 1))
filename_png = paste0(directory_name,'//Master_Grid_Figure_Normalized.png') 
ggsave(filename_png, combined_plot_grid_normalized, device = 'png', width = 11, height = 10, dpi = 600)
filename_eps = paste0(directory_name,'//Master_Grid_Figure_Normalized.pdf') 
ggsave(filename_eps, plot = print(combined_plot_grid_normalized), device = pdf, width = 11, height = 10)

# PART IIIA ################################ INFERENTIAL STATISTICS ON PREPROCESSED DATA ################################# 

############################################# TIME SERIES: SUM #############################################

# Group summary statistics by mean and SD
data_means_sum <- data %>%
                      group_by(Group, Time_Step) %>%
                      get_summary_stats(Sum, type = "mean_sd")

# Check for outliers
data_outliers_sum <- data %>%
                     group_by(Group, Time_Step) %>%
                     identify_outliers(Sum)

# Check for normality 
data_sw_sum <- try(data %>%
               group_by(Group, Time_Step) %>%
               shapiro_test(Sum))

# Run two-way repeated-measures ANOVA 
#rmaov_sum <- data %>% 
#             anova_test(formula = Sum ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Sum, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
#             get_anova_table(rmaov_sum)
modlm_sum <- lmer(Sum ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_sum <- car::Anova(modlm_sum, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_sum_sw <- shapiro_test(residuals(modlm_sum))
          
# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_sum_group <- data %>%
                 group_by(Time_Step) %>%
                 pairwise_t_test(
                 Sum ~ Group, paired = FALSE, 
                 p.adjust.method = "BH"
                 )

# Fit linear model to ranks instead 
modglm_sum <- lmer(rank(Sum) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_sum <- car::Anova(modglm_sum, REML=TRUE, test="F")

# Check rank-transformation linear model residuals for normality 
glm_rmaov_sum_sw <- shapiro_test(residuals(modglm_sum))

# Pairwise comparisons of ranks post hoc with Benjamini-Hochberg correction 
data$Sum_rank <-data$Sum
data$Sum_rank <-rank(data$Sum)
pwc_sum_group_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_rank ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_sum_group_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Sum ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# TIME SERIES: CUMULATIVE INTEGRAL #############################################

# Group summary statistics by mean and SD
data_means_int <- data %>%
                  group_by(Group, Time_Step) %>%
                  get_summary_stats(Cumulative_Sum, type = "mean_sd")

# Check for outliers
data_outliers_int <- data %>%
                     group_by(Group, Time_Step) %>%
                     identify_outliers(Cumulative_Sum)

# Check for normality 
data_sw_int <- try(data %>%
               group_by(Group, Time_Step) %>%
               shapiro_test(Cumulative_Sum))

# Run two-way repeated-measures ANOVA 
#rmaov_int <- data %>% 
#             anova_test(formula = Cumulative_Sum ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Cumulative_Sum, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
#             get_anova_table(rmaov_int)
modlm_int <- lmer(Cumulative_Sum ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_int <- car::Anova(modlm_int, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_int_sw <- shapiro_test(residuals(modlm_int))

# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_int_group <- data %>%
                  group_by(Time_Step) %>%
                  pairwise_t_test(
                  Cumulative_Sum ~ Group, paired = FALSE, 
                  p.adjust.method = "BH"
                  )

# Fit linear model to ranks instead
modglm_int <- lmer(rank(Cumulative_Sum) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_int <- car::Anova(modglm_int, REML=TRUE, test="F")

# Check rank-transformation linear model residuals for normality 
glm_rmaov_int_sw <- shapiro_test(residuals(modglm_int))

# Pairwise comparisons of ranks post hoc with Benjamini-Hochberg correction 
data$Cumulative_Sum_rank <-data$Cumulative_Sum
data$Cumulative_Sum_rank <-rank(data$Cumulative_Sum)
pwc_int_group_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Cumulative_Sum_rank ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )


# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_int_group_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Cumulative_Sum ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# TIME SERIES: SLOPE #############################################

# Group summary statistics by mean and SD
data_means_slope <- data %>%
                  group_by(Group, Time_Step) %>%
                  get_summary_stats(Sum_Slope, type = "mean_sd")

# Check for outliers
data_outliers_slope <- data %>%
                  group_by(Group, Time_Step) %>%
                  identify_outliers(Sum_Slope)

# Check for normality 
data_sw_slope <- try(data %>%
                 group_by(Group, Time_Step) %>%
                 shapiro_test(Sum_Slope))

# Run two-way repeated-measures ANOVA 
# rmaov_slope <- data %>%
#               anova_test(formula = Sum_Slope ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Sum_Slope, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
#               get_anova_table(rmaov_slope)
modlm_slope <- lmer(Sum_Slope ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_slope <- car::Anova(modlm_slope, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_slope_sw <- shapiro_test(residuals(modlm_slope))
               
# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_slope_group <- data %>%
                 group_by(Time_Step) %>%
                 pairwise_t_test(
                 Sum_Slope ~ Group, paired = FALSE, 
                 p.adjust.method = "BH"
                 )

# Fit linear model to rans instead
modglm_slope <- lmer(rank(Sum_Slope) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_slope <- car::Anova(modglm_slope, REML=TRUE, test="F")

# Check rank-transformed linear model residuals for normality 
glm_rmaov_slope_sw <- shapiro_test(residuals(modglm_slope))

# Pairwise comparisons of ranks post hoc with Benjamini-Hochberg correction 
data$Sum_Slope_rank <-data$Sum_Slope
data$Sum_Slope_rank <-rank(data$Sum_Slope)
pwc_slope_group_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_Slope_rank ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_slope_group_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Sum_Slope ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# BETWEEN-GROUP STATISTICS #############################################

# Time to max enhancement
desc_time_to_max_enhancement <- describeBy(data_max_enhancement$Time_Step, data_max_enhancement$Group, mat = TRUE) 
aov_time_to_max_enhancement <- aov(data = data_max_enhancement, formula = Time_Step ~ Group)
summary(aov_time_to_max_enhancement)
aov_time_to_max_enhancement_shapiro_residuals = shapiro_test(residuals(aov_time_to_max_enhancement))

# Slope before max enhancement
desc_slope_to_max_enhancement <- describeBy(data_max_enhancement$Slope_before_Max_Enhancement, data_max_enhancement$Group, mat = TRUE) 
aov_slope_to_max_enhancement <- aov(data = data_max_enhancement, formula = Slope_before_Max_Enhancement ~ Group)
summary(aov_slope_to_max_enhancement)
aov_slope_to_max_enhancement_shapiro_residuals = shapiro_test(residuals(aov_slope_to_max_enhancement))

# Slope after max enhancement
desc_slope_after_max_enhancement <- describeBy(data_max_enhancement$Slope_after_Max_Enhancement, data_max_enhancement$Group, mat = TRUE) 
aov_slope_after_max_enhancement <- aov(data = data_max_enhancement, formula = Slope_after_Max_Enhancement ~ Group)
summary(aov_slope_after_max_enhancement)
aov_slope_after_max_enhancement_shapiro_residuals = shapiro_test(residuals(aov_slope_after_max_enhancement))

# Maximum cumulative enhancement 
desc_final_cumulative_enhancement <- describeBy(data_final_time_step$Cumulative_Sum, data_final_time_step$Group, mat = TRUE) 
aov_final_cumulative_enhancement <- aov(data = data_final_time_step, formula = Cumulative_Sum ~ Group)
summary(aov_final_cumulative_enhancement)
aov_final_cumulative_enhancement_shapiro_residuals = shapiro_test(residuals(aov_final_cumulative_enhancement))

# Enhancement volume
desc_enhancement_volume <- describeBy(data_final_time_step$Num_Voxels, data_final_time_step$Group, mat = TRUE) 
aov_enhancement_volume <- aov(data = data_final_time_step, formula = Num_Voxels ~ Group)
summary(aov_enhancement_volume)
aov_enhancement_volume_shapiro_residuals = shapiro_test(residuals(aov_enhancement_volume))

############################################# PRINT ALL FINDINGS #############################################
#Define filename
filename_txt = paste0(directory_name,'//Final_Inferential_Statistics.txt') 
filename_txt = file(filename_txt, 'w')

cat("############################################# TIME SERIES: SUM #############################################\n", file = filename_txt)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_sum), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_sum), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_sum), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_sum), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_sum_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_sum), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_sum_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# TIME SERIES: CUMULATIVE INTEGRAL #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_int), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_int), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_int), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_int), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_int_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_int), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_int_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# TIME SERIES: SLOPE #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_slope), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_slope), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_slope), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_slope), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_slope_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_slope), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_slope_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# BETWEEN-GROUP STATISTICS #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Time to Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_time_to_max_enhancement, file = filename_txt, append = TRUE)
cat("\n\n Time to Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_time_to_max_enhancement), file = filename_txt, append = TRUE)
cat("\n\n Time to Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_time_to_max_enhancement_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_slope_to_max_enhancement, file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_slope_to_max_enhancement), file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_slope_to_max_enhancement_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_slope_after_max_enhancement, file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_slope_after_max_enhancement), file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_slope_after_max_enhancement_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_final_cumulative_enhancement, file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_final_cumulative_enhancement), file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_final_cumulative_enhancement_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Enhancement Volume Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_enhancement_volume , file = filename_txt, append = TRUE)
cat("\n\n Enhancement Volume\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_enhancement_volume ), file = filename_txt, append = TRUE)
cat("\n\n Enhancement Volume Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_enhancement_volume_shapiro_residuals, file = filename_txt, append = TRUE)

close(filename_txt)

# PART IIIB ################################ INFERENTIAL STATISTICS ON PREPROCESSED DATA - NORMALIZED ################################# 

############################################# TIME SERIES: SUM #############################################

# Group summary statistics by mean and SD
data_means_sum_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  get_summary_stats(Sum_Normalized, type = "mean_sd")

# Check for outliers
data_outliers_sum_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  identify_outliers(Sum_Normalized)

# Check for normality 
data_sw_sum_normalized <- try(data %>%
  group_by(Group, Time_Step) %>%
  shapiro_test(Sum_Normalized))

# Run two-way repeated-measures ANOVA 
# rmaov_sum_normalized <- data %>% 
#  anova_test(formula = Sum_Normalized ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Sum_Normalized, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
# get_anova_table(rmaov_sum_normalized)
modlm_sum_normalized <- lmer(Sum_Normalized ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_sum_normalized <- car::Anova(modlm_sum_normalized, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_sum_normalized_sw <- shapiro_test(residuals(modlm_sum_normalized))

# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_sum_group_normalized <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_Normalized ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Fit rank-transformed linear model in case general linear model residuals are not normal 
modglm_sum_normalized <- lmer(rank(Sum_Normalized) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_sum_normalized <- car::Anova(modglm_sum_normalized, REML=TRUE, test="F")

# Check rank-transformed linear model residuals for normality 
glm_rmaov_sum_normalized_sw <- shapiro_test(residuals(modglm_sum_normalized))

# Rank-transformed pairwise comparisons post hoc with Benjamini-Hochberg correction 
data$Sum_Normalized_rank <-data$Sum_Normalized
data$Sum_Normalized_rank <-rank(data$Sum_Normalized)
pwc_sum_group_normalized_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_Normalized_rank ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_sum_group_normalized_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Sum_Normalized ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# TIME SERIES: CUMULATIVE INTEGRAL #############################################

# Group summary statistics by mean and SD
data_means_int_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  get_summary_stats(Cumulative_Sum_Normalized, type = "mean_sd")

# Check for outliers
data_outliers_int_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  identify_outliers(Cumulative_Sum_Normalized)

# Check for normality 
data_sw_int_normalized <- try(data %>%
  group_by(Group, Time_Step) %>%
  shapiro_test(Cumulative_Sum_Normalized))

# Run two-way repeated-measures ANOVA 
# rmaov_int_normalized <- data %>% 
#  anova_test(formula = Cumulative_Sum_Normalized ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Cumulative_Sum_Normalized, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
# get_anova_table(rmaov_int_normalized)
modlm_int_normalized <- lmer(Cumulative_Sum_Normalized ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_int_normalized <- car::Anova(modlm_int_normalized, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_int_normalized_sw <- shapiro_test(residuals(modlm_int_normalized))

# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_int_group_normalized <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Cumulative_Sum_Normalized ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Fit rank-transformed linear model in case general linear model residuals are not normal 
modglm_int_normalized <- lmer(rank(Cumulative_Sum_Normalized) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_int_normalized <- car::Anova(modglm_int_normalized, REML=TRUE, test="F")

# Check rank-transformed linear model residuals for normality 
glm_rmaov_int_normalized_sw <- shapiro_test(residuals(modglm_int_normalized))

# Rank-transformed pairwise comparisons post hoc with Benjamini-Hochberg correction 
data$Cumulative_Sum_Normalized_rank <-data$Cumulative_Sum_Normalized
data$Cumulative_Sum_Normalized_rank <-rank(data$Cumulative_Sum_Normalized)
pwc_int_group_normalized_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Cumulative_Sum_Normalized_rank ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_int_group_normalized_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Cumulative_Sum_Normalized ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# TIME SERIES: SLOPE #############################################

# Group summary statistics by mean and SD
data_means_slope_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  get_summary_stats(Sum_Normalized_Slope, type = "mean_sd")

# Check for outliers
data_outliers_slope_normalized <- data %>%
  group_by(Group, Time_Step) %>%
  identify_outliers(Sum_Normalized_Slope)

# Check for normality 
data_sw_slope_normalized <- try(data %>%
  group_by(Group, Time_Step) %>%
  shapiro_test(Sum_Normalized_Slope))

# Run two-way repeated-measures ANOVA 
# rmaov_slope_normalized <- data %>%
#  anova_test(formula = Sum_Normalized_Slope ~ Group*Time_Step + Error(Animal_ID/Time_Step), dv = Sum_Normalized_Slope, wid = Animal_ID, within = c(Time_Step), detailed = TRUE)
# get_anova_table(rmaov_slope_normalized)
modlm_slope_normalized <- lmer(Sum_Normalized_Slope ~ Group * Time_Step + (1 | Animal_ID), data = data)
rmaov_slope_normalized <- car::Anova(modlm_slope_normalized, REML=TRUE, test="F")

# Check linear model residuals for normality 
rmaov_slope_normalized_sw <- shapiro_test(residuals(modlm_slope_normalized))

# Pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_slope_group_normalized <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_Normalized_Slope ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Nonparametric pairwise comparisons post hoc with Benjamini-Hochberg correction 
pwc_slope_group_normalized_wilcox <- data %>%
  group_by(Time_Step) %>%
  wilcox_test(
    Sum_Normalized_Slope ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

# Fit log-linked gamma generalized linear model in case general linear model residuals are not normal 
modglm_slope_normalized <- lmer(rank(Sum_Normalized_Slope) ~ Group * Time_Step + (1 | Animal_ID), data = data)
glm_rmaov_slope_normalized <- car::Anova(modglm_slope_normalized, REML=TRUE, test="F")

# Check log-linked gamma generalized linear model residuals for normality 
glm_rmaov_slope_normalized_sw <- shapiro_test(residuals(modglm_slope_normalized))

# Rank-transformed pairwise comparisons post hoc with Benjamini-Hochberg correction 
data$Sum_Normalized_Slope_rank <-data$Sum_Normalized_Slope
data$Sum_Normalized_Slope_rank <-rank(data$Sum_Normalized_Slope)
pwc_slope_group_normalized_rank <- data %>%
  group_by(Time_Step) %>%
  pairwise_t_test(
    Sum_Normalized_Slope ~ Group, paired = FALSE, 
    p.adjust.method = "BH"
  )

############################################# BETWEEN-GROUP STATISTICS #############################################

# Time to max enhancement
desc_time_to_max_enhancement_normalized <- describeBy(data_max_enhancement_normalized$Time_Step, data_max_enhancement_normalized$Group, mat = TRUE) 
aov_time_to_max_enhancement_normalized <- aov(data = data_max_enhancement_normalized, formula = Time_Step ~ Group)
summary(aov_time_to_max_enhancement_normalized)
aov_time_to_max_enhancement_shapiro_residuals_normalized = shapiro_test(residuals(aov_time_to_max_enhancement_normalized))

# Slope before max enhancement
desc_slope_to_max_enhancement_normalized <- describeBy(data_max_enhancement_normalized$Slope_before_Max_Enhancement, data_max_enhancement_normalized$Group, mat = TRUE) 
aov_slope_to_max_enhancement_normalized <- aov(data = data_max_enhancement_normalized, formula = Slope_before_Max_Enhancement ~ Group)
summary(aov_slope_to_max_enhancement_normalized)
aov_slope_to_max_enhancement_normalized_shapiro_residuals = shapiro_test(residuals(aov_slope_to_max_enhancement_normalized))

# Slope after max enhancement
desc_slope_after_max_enhancement_normalized <- describeBy(data_max_enhancement_normalized$Slope_after_Max_Enhancement, data_max_enhancement_normalized$Group, mat = TRUE) 
aov_slope_after_max_enhancement_normalized <- aov(data = data_max_enhancement_normalized, formula = Slope_after_Max_Enhancement ~ Group)
summary(aov_slope_after_max_enhancement_normalized)
aov_slope_after_max_enhancement_normalized_shapiro_residuals = shapiro_test(residuals(aov_slope_after_max_enhancement_normalized))

# Maximum cumulative enhancement 
desc_final_cumulative_enhancement_normalized <- describeBy(data_final_time_step$Cumulative_Sum_Normalized, data_final_time_step$Group, mat = TRUE) 
aov_final_cumulative_enhancement_normalized <- aov(data = data_final_time_step, formula = Cumulative_Sum_Normalized ~ Group)
summary(aov_final_cumulative_enhancement_normalized)
aov_final_cumulative_enhancement_shapiro_residuals_normalized = shapiro_test(residuals(aov_final_cumulative_enhancement_normalized))

#Define filename
filename_txt = paste0(directory_name,'//Final_Inferential_Statistics_Normalized.txt') 
filename_txt = file(filename_txt, 'w')

cat("############################################# TIME SERIES: SUM - NORMALIZED #############################################\n", file = filename_txt)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_sum_normalized), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_sum_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_sum_normalized), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_sum_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_sum_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group_normalized), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_sum_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_sum_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group_normalized_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_sum_group_normalized_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# TIME SERIES: CUMULATIVE INTEGRAL- NORMALIZED #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_int_normalized), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_int_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_int_normalized), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_int_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_int_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group_normalized), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_int_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_int_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group_normalized_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_int_group_normalized_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# TIME SERIES: SLOPE - NORMALIZED #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Means and SD\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_means_slope_normalized), file = filename_txt, append = TRUE)
cat("\n\n Outliers\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_outliers_slope_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(data_sw_slope_normalized), file = filename_txt, append = TRUE)
cat("\n\n Repeated-Measures ANOVA\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(rmaov_slope_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(rmaov_slope_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group_normalized), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")
cat("\n\n Repeated-Measures ANOVA on Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(get_anova_table(glm_rmaov_slope_normalized), file = filename_txt, append = TRUE)
cat("\n\n Shapiro-Wilk Test for Rank-Transformed Linear Model Residuals\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(glm_rmaov_slope_normalized_sw, file = filename_txt, append = TRUE)
cat("\n\n Pairwise Comparisons for Rank-Transformed Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group_normalized_rank), file = filename_txt, append = TRUE)
cat("\n\n Nonparametric Pairwise Comparisons for Linear Model\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(as.data.frame(pwc_slope_group_normalized_wilcox), file = filename_txt, append = TRUE)
cat("\n\n\n", file = filename_txt, append=TRUE, sep = "\n")

cat("############################################# BETWEEN-GROUP STATISTICS - NORMALIZED #############################################\n", file = filename_txt, append=TRUE)
cat("\n\n Time to Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_time_to_max_enhancement_normalized, file = filename_txt, append = TRUE)
cat("\n\n Time to Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_time_to_max_enhancement_normalized), file = filename_txt, append = TRUE)
cat("\n\n Time to Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_time_to_max_enhancement_shapiro_residuals_normalized, file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_slope_to_max_enhancement_normalized, file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_slope_to_max_enhancement_normalized), file = filename_txt, append = TRUE)
cat("\n\n Slope to Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_slope_to_max_enhancement_normalized_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_slope_after_max_enhancement_normalized, file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_slope_after_max_enhancement_normalized), file = filename_txt, append = TRUE)
cat("\n\n Slope after Max Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_slope_after_max_enhancement_normalized_shapiro_residuals, file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement Descriptives\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(desc_final_cumulative_enhancement_normalized, file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(summary(aov_final_cumulative_enhancement_normalized), file = filename_txt, append = TRUE)
cat("\n\n Final Cumulative Enhancement Residual Shapiro-Wilk Test\n", file = filename_txt, append=TRUE, sep = "\n")
captureOutput(aov_final_cumulative_enhancement_shapiro_residuals_normalized, file = filename_txt, append = TRUE)

close(filename_txt)