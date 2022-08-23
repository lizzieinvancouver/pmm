# Code that visualizes speed comparisions between different uberspeed Stan models
# Started by Dinara, August 19 2022

# Housekeeping
rm(list = ls())
setwd("/Users/dindin/Documents/work/pmm/analyses")

# Libraries
install.packages("ggplot2")
install.packages("data.table")
library("ggplot2")
library("data.table")

# Reading in data
filenames <- list.files("output/speedtests", pattern="*.csv", full.names=TRUE)
dfs <- lapply(filenames, function(fn) {
  df <- read.csv(fn)
  df$model <- gsub(".*uberspeed_([^.]+)_compare.csv", "\\1", fn)
  df
})
models <- do.call(rbind, dfs)
models$model <- as.factor(models$model)

# Set global theme for all upcoming plots
theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             strip.background = element_blank(),
             panel.border = element_rect(colour = "black", fill = NA),
             plot.margin = unit(c(1,1,1,1), units = , "cm"),
             axis.text = element_text(size = 10),
             axis.title = element_text(size = 12, face = "plain"),
             plot.title = element_text(size = 16 , hjust = 0.5, face = "bold"))

# Visualize speed vs model type as a box and whisker plot
# for reference
speed_models <- aggregate(runtime ~ model, models, function(x) {
  c(mean = mean(x), quantile(x, probs = c(.25, .75)), max = max(x), min = min(x))
})

plot1 <- ggplot(models, aes(x = model, y = runtime, fill = model)) +
  geom_boxplot(alpha=0.3) +
  labs(y = "Runtime (in minutes)", fill = "Model type") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + 
  ggtitle("Speed vs Model") +
  scale_fill_discrete(labels=c("New function, no cholesky, no priors",
                               "New function, no cholesky, with priors",
                               "New function, with cholesky, no priors",
                               "New function, with cholesky, with priors",
                               "Old function, no cholesky, no priors",
                               "Old function, no cholesky, with priors",
                               "Old function, with cholesky, no priors",
                               "Old function, with cholesky, with priors"))
plot1
ggsave("speedtests_viz/SpeedVsModel_Boxplot.pdf", plot1)


# Visualize true value of params vs their mean value
# True value for parameters
models[models$variable == "b_z", "true_param_val"] = 0.6
models[models$variable == "lam_interceptsb", "true_param_val"] = 0.7
models[models$variable == "sigma_interceptsb", "true_param_val"] = 0.1
models[models$variable == "a_z", "true_param_val"] = 4
models[models$variable == "lam_interceptsa", "true_param_val"] = 0.4
models[models$variable == "sigma_interceptsa", "true_param_val"] = 0.2
models[models$variable == "sigma_y", "true_param_val"] = 0.01

# Create a new variable that would hold value of how different mean param value from the actual
models$param_diff <- models$mean - models$true_param_val


params <- c("b_z", "lam_interceptsb", "sigma_interceptsb", "a_z",
            "lam_interceptsa", "sigma_interceptsa", "sigma_y")

plots2 <- lapply(params, function(pm) {
  plot2 <- ggplot(models[models$variable == pm,], aes(x = model, y = param_diff, fill = model)) +
    geom_boxplot(alpha=0.3) +
    labs(y = "Diff between mean and true values", fill = "Model type") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank()) + 
    ggtitle(paste0("Difference between mean and true ", pm, " values for each model type")) +
    scale_fill_discrete(labels=c("New function, no cholesky, no priors",
                                 "New function, no cholesky, with priors",
                                 "New function, with cholesky, no priors",
                                 "New function, with cholesky, with priors",
                                 "Old function, no cholesky, no priors",
                                 "Old function, no cholesky, with priors",
                                 "Old function, with cholesky, no priors",
                                 "Old function, with cholesky, with priors"))
  ggsave(paste0("speedtests_viz/mean_vs_true_vals/for_", pm, ".pdf"), plot2)
  plot2
})

# Visualize whether the true parameter value is within the 25-75% interval for each model type
for (i in 1:nrow(models)) {
  models[i, "within_interval"] <- if ((models[i, "first_quartile"] < models[i, "true_param_val"])
                                      & (models[i, "third_quartile"] > models[i, "true_param_val"])) 1 else 0
}

models_within_count <- unique(data.table(model_type = models$model, parameter = models$variable, 
                                         in_interval = models$within_interval, seed = models$seed))
models_within_count$count <- nrow(models[models$model == models_within_count$model_type &
                                           models$variable == models_within_count$parameter &
                                           models$within_interval == 1, ])

models_within_count <- aggregate(models$within_interval, by = list(models$model, models$variable), FUN = sum)

plot3 <- ggplot(models_within_count, aes(y = Group.1,  x = x, fill = Group.2)) +
    geom_bar(stat = "identity") + 
    labs(x = "# of model runs", y = "Model type", fill = "Parameter", 
         title = "Number of model runs where true parameter value is within 25-75% interval for each model type",
         subtitle = "Total # of model runs = 70") +
    scale_x_continuous(breaks = seq(0, 70, 10), limits = c(0, 70)) +
    scale_y_discrete(labels=c("New function, no cholesky, no priors",
                              "New function, no cholesky, with priors",
                              "New function, with cholesky, no priors",
                              "New function, with cholesky, with priors",
                              "Old function, no cholesky, no priors",
                              "Old function, no cholesky, with priors",
                              "Old function, with cholesky, no priors",
                              "Old function, with cholesky, with priors"))
plot3
ggsave("speedtests_viz/TrueParamValWithinInterval.pdf", plot3, width = 14)

