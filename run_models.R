library('tidyr')
library('tibble')
library('dplyr')
library('ggplot2')
library('mgcv')
library('cowplot')
library('itsadug')
library('grid')

# plotting function
plotter <- function(plot_data, orig_data, fit_xx, fit_yy, orig_xx, orig_yy, fit_grp, orig_grp, lower="lower", upper="upper") {
  ggplot(data=plot_data) +
                 geom_line(aes_string(x = fit_xx, y = fit_yy, group = fit_grp, colour=fit_grp), show.legend = FALSE) +
                 geom_ribbon(aes_string(x = fit_xx, ymin = lower, ymax = upper, fill = fit_grp), alpha=0.075, show.legend=FALSE) +
                 geom_line(data=orig_data, aes_string(x = orig_xx, y = orig_yy, group = orig_grp), alpha=.1, size=.4, linetype=5, show.legend = FALSE) + 
                 expand_limits(y=0) +
                 labs(x = "Age", y = "volume")  +
                 theme_minimal() 
}

# function to create individual data for plotting trajectories
expand_data <- function(data){
  expanded_data = NULL
  for (id in unique(data$ID)){
    age_seq = seq(from=min(data[data$ID==id,][,4]), to=max(data[data$ID==id,][,4]), length.out=50)
    expid = expand.grid(ID=id, Male=data[data$ID==id,][[1,3]], Age=age_seq)
    expanded_data=rbind(expanded_data, expid)
  }
  return(expanded_data)
}

####################################################################################
# Setup and load data
####################################################################################
# default k for GAMs
k=20

# load data
data <- as_tibble(read.csv('./brain_tissue_volumes.csv')[ ,2:14])
data$ID <- as.factor(data$ID)

# convert volume to ml (mm3/1000)
data$ICV <- data$ICV/1000
data$corticalGM <- data$corticalGM/1000
data$subcorticalGM <- data$subcorticalGM/1000
data$WM <- data$WM/1000
data$cerebellum <- data$cerebellum/1000

# Fit models 
##########################################################################
# Intracranial Volume (ICV)
##########################################################################
# Start with a linear mixed effect models with random intercepts, modelling the dependence of tissue volume on sex. 
model0 <- lmerTest::lmer(ICV ~ 1 + (1 + 1|ID) + Male, data=data, REML=F)

# Then add in Age... 
model1 <- lmerTest::lmer(ICV ~ 1 + (1 + 1|ID) + Age + Male, data=data, REML=F)

# Then model age as a smooth function with a GAM...
modelGAM <- gam(ICV ~ 1 + s(ID, bs="re") + s(Age, k=k) + Male, method="ML", data=data)

# compare models with AIC
model_selection = AIC(model0,model1,modelGAM)
print('ICV models')
print(model_selection)

# and finally, compare to the inclusion of random slopes as well
modelGAMslope <- gam(ICV ~ 1 + s(ID, bs="re") + s(ID, Age, bs="re") + s(Age, k=k) + Male, method="ML", data=data)
gam_model_selection <- AIC(modelGAM, modelGAMslope)
compareML(modelGAM, modelGAMslope, suggest.report = T)

####################################################################################
# Plot models
####################################################################################
# choose model to plot: modelGAM, modelGAMslope - check lines 105 and 128 if using modelGAMslope
plot_model <- modelGAM

####################################################################################
# Individual trajectories
expanded_data <- expand_data(data)

# model predictions
predictions <- predict(plot_model, newdata=expanded_data, se.fit = T)
predicted_data <- data.frame(predictions$fit , predictions$se.fit)
colnames(predicted_data) <- c("fit","se")
# CI
predicted_data$upper <- predicted_data$fit  + 1.96*(predicted_data$se)
predicted_data$lower <- predicted_data$fit  - 1.96*(predicted_data$se)

# collect prediction, lower and upper CI together
tempdf <- bind_cols(expanded_data, predicted_data)

# plot individual fits
plot_ICV_1 <- plotter(tempdf, data, "Age", "fit", "Age", "ICV", fit_grp="ID", orig_grp="ID")

####################################################################################
# Average trajectories
age_sequence <- seq(from=min(data$Age),to=max(data$Age), length.out=100)
all_x <- expand.grid(ID=1, #ignored
                     Male=0.5,
                     Age=age_sequence)

# model predictions without random effects
predictions <- predict(plot_model, newdata=all_x, se.fit = T, exclude=c(s(ID)))
#predictions <- predict(plot_model, newdata=all_x, se.fit = T, exclude = c(s(ID), s(ID,Age))) # use for modelGAMslope
predicted_data <- data.frame(fit=predictions$fit , se=predictions$se.fit)
# CI
predicted_data$upper <- predicted_data$fit  + 1.96*(predicted_data$se)
predicted_data$lower <- predicted_data$fit  - 1.96*(predicted_data$se)

# collect prediction, lower and upper CI together
tempdf <- bind_cols(all_x, predicted_data)
tempdf$ID <- factor(tempdf$ID)

# plot group average fit
plot_ICV_2 <- plotter(tempdf, data, "Age", "fit", "Age", "ICV", fit_grp="ID", orig_grp="ID")

####################################################################################
# Male and female trajectories
male_x <- all_x
male_x$Male <- c(1)
female_x <- all_x
female_x$Male <- c(0)
all_x <- rbind(male_x, female_x)

# model predictions without random effects
predictions <- predict(plot_model, newdata=all_x, se.fit = T, exclude = c(s(ID)))
#predictionsm <- predict(plot_model, newdata=all_x, se.fit = T, exclude = c(s(ID), s(ID,Age))) # use for modelGAMslope
predicted_data <- data.frame(fit=predictions$fit , se=predictions$se.fit)
# CI
predicted_data$upper <- predicted_data$fit  + 1.96*(predicted_data$se)
predicted_data$lower <- predicted_data$fit  - 1.96*(predicted_data$se)

# collect prediction, lower and upper CI together
tempdf <- bind_cols(all_x, predicted_data)
tempdf$Male <- factor(tempdf$Male)

# plot trajectories
plot_ICV_3 <-  plotter(tempdf, data, "Age", "fit", "Age", "ICV", fit_grp="Male", orig_grp="ID")

####################################################################################
# save plots
icv_plots <- plot_grid(plot_ICV_1, plot_ICV_2, plot_ICV_3, labels = "auto",nrow=1)
ggplot2::ggsave(file="ICV_plots.pdf",icv_plots, width = 25, height = 10, units = "cm", dpi=200)
####################################################################################



####################################################################################
# Regional tissue volume, correcting for ICV
# in this case using corticalGM, edit as required.
####################################################################################
# run models
model0tissue <- lmerTest::lmer(corticalGM ~ 1 + (1 + 1|ID) + ICV + Male, data=data, REML=F)
model1tissue <- lmerTest::lmer(corticalGM ~ 1 + (1 + 1|ID) + Age + ICV + Male, data=data, REML=F)
modelGAMtissue <- gam(corticalGM ~ 1 + s(ID, bs="re") + s(Age, k=k) + ICV + Male, method="ML", data=data)
modelGAMtissueslope <- gam(corticalGM ~ 1 + s(ID, bs="re") +  s(ID, Age, bs="re") + s(Age, k=k) + ICV + Male, method="ML", data=data)

# compare models with AIC
tissue_model_selection = AIC(model0tissue,model1tissue,modelGAMtissue)
print('tissue models')
print(tissue_model_selection)

####################################################################################
# Plot models 
####################################################################################
# Tissue volume
# choose model to plot: modelGAM, modelGAMslope - check lines 186, 194 and 217 if using modelGAMslope
plot_model <- modelGAMtissue

####################################################################################
# Average trajectories
age_sequence <- seq(from=min(data$Age),to=max(data$Age), length.out=100)
all_x <- expand.grid(ID=1,
                     Male=0.5,
                     ICV=mean(data$ICV),
                     Age=age_sequence)

# refit model on ICV-corrected tissue volume data using same smoothing parameters for Age and ID
data_tmp <- data
data_tmp$corticalGM <- data_tmp$corticalGM/data_tmp$ICV
tmp_model <- gam(formula(plot_model),
                 sp = c(plot_model$sp["s(ID)"],
                        plot_model$sp["s(Age)"]),
                 method="ML", data=data_tmp)
#tmp_model <- gam(formula(plot_model),                       # use this instead for modelGAMslope
#                 sp = c(plot_model$sp["s(ID)"],
#                        plot_model$sp["s(ID,Age)"],
#                        plot_model$sp["s(Age)"]),
#                 method="ML", data=data_tmp)

# model predictions without random effects
predictions <- predict(tmp_model, newdata=all_x, se.fit = T, exclude=c(s(ID)))
#predictions <- predict(tmp_model, newdata=all_x, se.fit = T, exclude=c(s(ID), s(ID,Age)))  # if using modelGAMslope
predicted_data <- data.frame(fit=predictions$fit , se=predictions$se.fit)
# CI
predicted_data$upper <- predicted_data$fit  + 1.96*(predicted_data$se)
predicted_data$lower <- predicted_data$fit  - 1.96*(predicted_data$se)

# collect prediction, lower and upper CI together
tempdf <- bind_cols(all_x, predicted_data)
tempdf$ID <- factor(tempdf$ID)

# plot group average fit
plot_tissue_1 <- plotter(tempdf, data_tmp, "Age", "fit", "Age", "corticalGM", "ID", "ID")
 
####################################################################################
# Male and female trajectories
male_x <- all_x
male_x$Male <- c(1)
female_x <- all_x
female_x$Male <- c(0)
all_x <- rbind(male_x, female_x)

# model predictions without random effects
predictions <- predict(tmp_model, newdata=all_x, se.fit = T, exclude=c(s(ID)))
#predictions <- predict(tmp_model, newdata=all_x, se.fit = T, exclude=c(s(ID), s(ID,Age)))  # if using modelGAMslope
predicted_data <- data.frame(fit=predictions$fit , se=predictions$se.fit)
# CI
predicted_data$upper <- predicted_data$fit  + 1.96*(predicted_data$se)
predicted_data$lower <- predicted_data$fit  - 1.96*(predicted_data$se)

# collect prediction, lower and upper CI together
tempdf <- bind_cols(all_x, predicted_data)
tempdf$Male <- factor(tempdf$Male)

# plot
plot_tissue_2 <- plotter(tempdf, data_tmp, "Age", "fit", "Age", "corticalGM", "Male", "ID")
  
####################################################################################
# save out
tissue_plots <- plot_grid(plot_tissue_1, plot_tissue_2, labels = "auto",nrow=1)
ggplot2::ggsave(file="tissue_volume_plots.pdf", tissue_plots, width = 25, height = 10, units = "cm", dpi=200)
####################################################################################