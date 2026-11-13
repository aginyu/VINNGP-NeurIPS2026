# plot

rm(list = ls())
source('./src/dual_mesh_loc.R')
source('./src/logintensity_vis.R')
load('./data/comparsion_fullcase.RData')

getwd()

logintensity_mean_list = list(logintensity_mean_full_df)
names(logintensity_mean_list) = c("full")
logintensity_sd_list = list(logintensity_sd_full_df)
names(logintensity_sd_list) = c("full")


contour_step_mean = 0.3
contour_step_sd = 0.08

pdf("./LGCP_comparison.pdf", width = 8, height = 9)  # 英寸为单位


layout(matrix(c(1,2,7,3,4,7,5,6,7), nrow=3, ncol=3, byrow=TRUE), widths=c(1, 1, 0.25))
par(mar = c(2,4,2,1), oma = c(3,3,3,2))

val_min_true = min(true_log_intensity_mesh)
val_max_true = max(true_log_intensity_mesh)
zlim = c(val_min_true,val_max_true)

logintensity_vis(mesh = mesh_full, field_value_list = true_log_intensity_mesh, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
mtext("TRUE", side = 2, line = 2.5, cex = 2)

logintensity_vis(mesh = mesh_full, field_value_list = logintensity_mean_list$full$VINNGP, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
mtext("VoGCAM", side = 2, line = 2.5, cex = 2)

logintensity_vis(mesh = mesh_full, field_value_list = logintensity_mean_list$full$INLA, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
mtext("INLA", side = 2, line = 2.5, cex = 2)

logintensity_vis(mesh = mesh_full, field_value_list = logintensity_mean_list$full$VIFRK, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
#mtext("Full case", side = 1, line = 3.5, cex = 2)
mtext("VIFRK", side = 2, line = 2.5, cex = 2)

logintensity_vis(mesh = mesh_full, field_value_list = eta_mincon_n, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
mtext("MINCON", side = 2, line = 2.5, cex = 2)

logintensity_vis(mesh = mesh_full, field_value_list = eta_quasi_n, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = true_log_intensity_mesh)
#mtext("Full case", side = 1, line = 3.5, cex = 2)
mtext("QUASI", side = 2, line = 2.5, cex = 2)

values = true_log_intensity_mesh
val_min = min(values)
val_max = max(values)
quantiles = quantile(values, probs = seq(0, 1, by = 0.01))
color_breaks = seq(val_min,val_max,length.out = 100)
colors = rev(colorRampPalette(brewer.pal(11, "RdBu"))(length(color_breaks)-1))

# 只画图例
plot.new()
image.plot(legend.only = TRUE,
           zlim = range(values),
           col = colors,
           breaks = color_breaks,
           legend.shrink = 1,
           legend.width = 4,
           smallplot = c(0.35, 0.55, 0.022, 0.978),
           axis.args = list(cex.axis = 2, # Increased label size (e.g., 1.2x default)
                            font.axis = 4)
)

dev.off()


#----------------------------------- log-intensity posterior sd --------------------------------


layout(matrix(c(1,2,7,3,4,7,5,6,7), nrow=3, ncol=3, byrow=TRUE), widths=c(1, 1, 0.25))
par(mar=c(2, 2, 2, 1), oma=c(3,3,3,2))

ref_sd = unlist(logintensity_sd_list)

logintensity_vis(mesh = mesh_full, field_value_list = logintensity_sd_list$full$VINNGP, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = ref_sd)
mtext("VoGCAM", side = 2, line = 2.5, cex = 2)


logintensity_vis(mesh = mesh_full, field_value_list = logintensity_sd_list$full$INLA, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = ref_sd)
mtext("INLA", side = 2, line = 2.5, cex = 2)


logintensity_vis(mesh = mesh_full, field_value_list = logintensity_sd_list$full$VIFRK, case = "full", method = "TRUE", contour_step = contour_step_mean,ref_values = ref_sd)
mtext("Full case", side = 1, line = 3.5, cex = 2)
mtext("VIFRK", side = 2, line = 2.5, cex = 1.8)


values_mean = c(
                logintensity_sd_list$full$VINNGP, logintensity_sd_list$full$INLA, logintensity_sd_list$full$VIFRK)
values = sort(values_mean)
val_min = min(values)
val_max = max(values)
quantiles = quantile(values, probs = seq(0, 1, by = 0.01))
color_breaks = c(val_min, log(quantiles + 1 - min(values)))
color_breaks = sort(color_breaks)
colors = colorRampPalette(viridis(200))(length(color_breaks)-1)

# Plot color bar
plot.new()
image.plot(legend.only=TRUE, z = matrix(values, ncol = 1), col=colors, legend.shrink = 1, legend.mar = 10, legend.width = 60,
           smallplot = c(0.35, 0.55, 0.022, 0.978),
           axis.args = list(cex.axis = 2, # Increased label size (e.g., 1.2x default)
                            font.axis = 4))
#dev.off()

table_full

