##################################################
#                                                #
#  Seminar - Financial Econometrics              #
#                                                #
#  Thomas Theodor Kjølbye                        #
#                                                #
#  The following script produces all output      #
#  used in the paper. On my computer, the        #
#  entire script takes approx. 45 minute to run. #
#                                                #
##################################################

########################## Settings ###############################

# Clear 
rm(list = ls())

# Set working directory to source file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Set language to english
Sys.setenv(LANG = "EN")

# Disable scientific notation (e) 
options(scipen = 999)

# Load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, ggplot2, egg, gridExtra, knitr, kableExtra, glmnet, zeallot, neuralnet, comprehenr, data.table)

# Import own table customization compatible with LaTeX
source("table_theme.R", local = FALSE)

# Initialize warning log container
log <- c()

# Timer start
start_time <- Sys.time()

### 1. Simulate Covariates -----------------------------------------------------

# Simulation Setup 
draws = 50
variables = 40
iterations = 200

simCovariates <- function(d = draws, var = variables, ite = iterations) { # Draws, variables, iterations,
  
  # sim_data() returns the simulated observations (draws) of the covariates (Xs) which, in turn, will be used in the DGP.
  # The covariates are all drawn from the normal distribution. 
  
  # Set seed
  set.seed(007)
  
  # Initialize data holder 
  data <- as_tibble(data.frame(matrix(nrow = 0, ncol = var))) %>% 
    add_column(Iteration = NA, .before = "X1") # Column for each variable + an index
  
  for (i in 1:ite) { 
    
    # Initialize temp data holder
    data_temp = matrix(nrow = d, ncol = var) # Column for each variable 
    
   for (j in 1:var) {
     
     x = rnorm(d, 0, runif(1, min = 1, max = 1.3) ** log(j)) 
     data_temp[, j] = x
     
   }
  
    # Cholesky Decomposition to get independence 
    data_var <- var(data_temp)
    chol <- solve(chol(data_var), tol = 2.30e-18) # Set tolerance to circumvent linear dependence 
    data_temp <- data_temp %*% chol
    
    ### Diagnostics ### 
    if(all(round(var(data_temp)) != diag(var))){
      log <- c(log, "-- Not all simulated variables are independent")
    }

    # Generate random correlation matrix
    R <- matrix(runif(ncol(data_temp)^2,-1,1), ncol=ncol(data_temp)) # square matrix with draws from uniform dist.
    RtR <- R %*% t(R)
    corr <- cov2cor(RtR)
    
    ### Diagnostics ###
    if((sum(eigen(corr)$values > 0)) != ncol(data_temp)){
      log <- c(log, "-- corr is not positive definite")
    }
    
    # Transform data according to correlation matrix
    data_temp <- data_temp %*% chol(corr)
    
    # Transform to tibble, add iteration column, and append 
    data_temp <- as.tibble(data_temp) %>% 
      add_column(Iteration = i, .before = "V1")
    colnames(data_temp)  <- colnames(data)
    
    data <- bind_rows(data, data_temp)
    
  }
  
  return(data)
  
}

# Generate Covariates
covariates <- simCovariates()

# Table 1 - Snip of covariates 
table1_names <- c("Iteration", "X_1", "X_2", "X_3", "X_4", "Dots", "X_{37}", "X_{38}", "X_{39}", "X_{40}")
table1 <- table_theme(rbind(head(covariates, 4), tail(covariates, 3)) %>% select(Iteration, X1:X5, X37:X40), colnames = table1_names, 
                      caption = "Sample of simulated data from simCovariates()", escape = FALSE) %>% 
  kable_styling(latex_options = "scale_down")

# Table 1 - Snip of covariates (Small)
# The table is modified slight in LaTeX
table1_names <- c("Iteration", "X_1", "X_2", "X_3", "X_4", "X_5", "Dots", "X_{36}", "X_{37}", "X_{38}", "X_{39}", "X_{40}")
table1 <- table_theme(rbind(head(covariates, 4), tail(covariates, 3)) %>% select(Iteration, X1:X6, X36:X40), colnames = table1_names, 
                      caption = "Sample of simulated data from simCovariates()", escape = FALSE) %>% 
  kable_styling(latex_options = "scale_down")

# Figure 1 - Distribution of Covariates
norm_dist_df <- as.tibble(to_vec(for (i in 1:variables) runif(1, min = 1, max = 1.3) ** log(i))) %>% 
  `colnames<-`("sd") %>% 
  add_column(mean = 0, .before = "sd")

figure_1 <- ggplot(data = data.frame(x = c(-7, 7)), aes(x)) +            
  map2(.x = t(norm_dist_df)[1, ], .y = t(norm_dist_df)[2, ],
       .f = ~stat_function(fun = dnorm, n = 101, args = list(mean = .x, sd = .y))) +
  labs(y = "f(x)") + 
  theme_article() +
  theme(text = element_text(size = 15))
  #scale_y_continuous(breaks = NULL)

ggsave("figure1.pdf", plot=figure_1, width = 20, height = 18, units= "cm", dpi = 300)

# iterable object where each entry is a ggplot of the gaussian density with varying variance
# map2(.x = t(norm_dist_df)[1, ], .y = t(norm_dist_df)[2, ],
#      .f = ~{
#        
#        ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
#          stat_function(fun = dnorm, n = 101, args = list(mean = .x, sd = .y)) + ylab("") +
#          scale_y_continuous(breaks = NULL) + 
#          theme_article()
#        
#      })


### 2. Create Coefficient Cluster ----------------------------------------------

coefCluster <- function(betas = variables, center = c(1,3,5), radius = 5) {
  
  # coef_cluster returns a vector of betas of which 3 * 2 * radius - 1 are non-zero and split evenly in 3 separate cofficient clusters with varying 
  # center but identical radius. Reducing the radius yields a greater # of non-zero coefficients, though the function adheres to the rule. 
  # Analogous to sim_data(), coef_cluster also carries out some diagnostics. 
  
  # Initialize vector of betas
  beta = matrix(rep(0, times = betas)) #beta = t(matrix(nrow = 1, ncol = betas, dimnames = list("beta_value")))
  
  # Create index of non-zero betas
  index <- sample(seq(betas), size = 3 * (2 * radius - 1), replace = FALSE)
  
  # Clusters
  cluster_temp <- map_df(.x = rep(0, times = 2 * radius -1),
                         function(.x) { 
                           
                           return(data.frame(cluster1 = .x + center[1] + runif(1, -1, 1),
                                           cluster2 = .x + center[2] + runif(1, -1, 1),
                                           cluster3 = .x + center[3] + runif(1, -1, 1))) 
                         }) %>% 
    unlist()
  
  # cluster_temp <- map_df(.x = rep(0, times = 2 * radius -1),
  #                     ~(data.frame(cluster1 = .x + center[1] + runif(1, -1, 1),
  #                                            cluster2 = .x + center[2] + runif(1, -1, 1),
  #                                            cluster3 = .x + centeR[3] + runif(1, -1, 1))))
  
  for (i in enumerate(index)) {

    index_ite = i[[1]]; value = i[[2]]
    beta[value] = cluster_temp[index_ite]

  }
  
  ### Diagnostics ### 
  if(length(unique(beta)) != 3 * (2 * radius - 1) + 1){ # length(beta[, 1] == 0) != betas - 3 * (2 * radius -1) 
    log <- c(log, "-- Incorrect # of non-zero coefficients")
  }
  
  ### Diagnostics ### 
  if(FALSE){
    log <- c(log, "-- Incoerrect # of non-zero coefficients in each cluster")
  }

  return(beta)
  
}

# Betas - 7 set of betas one for each radius
beta1 <- coefCluster(radius = 1, center = c(1,3,5))
beta2 <- coefCluster(radius = 2, center = c(1,3,5))
beta3 <- coefCluster(radius = 3, center = c(1,3,5))
beta4 <- coefCluster(radius = 4, center = c(1,3,5))
beta5 <- coefCluster(radius = 5, center = c(1,3,5))
beta6 <- coefCluster(radius = 6, center = c(1,3,5))
beta7 <- coefCluster(radius = 7, center = c(1,3,5))
beta <- as.tibble(cbind(beta1, beta2, beta3, beta4, beta5, beta6, beta7)) # as.tibble and not matrix to make it compatible with e.g. map_df()

# Figure 2 - 2D plot of the coefficient cluster
figure_2A <- ggplot() + 
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "#00BFC4", alpha = 0.0, color = "white") +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 2, ymax = 4),
            fill = "#00BFC4", alpha = 0.5, color = "black", size = 0.4) + 
  labs(x = "Coefficient Cluster 1", y = "Coefficient Cluster 2") +
  theme_article()

figure_2B <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "#F8766D", alpha = 0.0, color = "white")  +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 4, ymax = 6),
            fill = "#F8766D", alpha = 0.8, color = "black", size = 0.4) +
  labs(x = "Coefficient Cluster 1", y = "Coefficient Cluster 3") +
  theme_article()

figure2 <- gridExtra::grid.arrange(figure_2A, figure_2B, ncol = 2)
ggsave("figure2.pdf", plot=figure2, width = 25, height = 10, units= "cm", dpi = 300)

# Table 2 - non-zero beta coefficients in each coefficient cluster through the radius 
table2_names <- c("Cluster Radius", "# non-zero coefficients", " # zero coefficients", "# of coefficients")
table2 <- table_theme(map_df(.x = beta, 
                             function(.x) {
                               
                               return(data.frame(
                                                 coeff = length(unique(.x)) - 1,
                                                 coeff_zero = 40 - length(unique(.x)) + 1,
                                                 coeff_total = 40)) 
                               
                             }) %>% add_column(seq(7), .before = "coeff"), 
                      colnames = table2_names, caption = "Coeffcient Clusters") %>% 
  kable_styling(latex_options = "scale_down")


### 3. Data Generating Process -------------------------------------------------

# DGP (number corresponds to cluster radius) - there is most likely a neater (more elegant) way to code this up
DGP1 <- as.matrix(covariates[, -1]) %*% beta1 + rnorm(1, mean = 0, sd = 1)
DGP2 <- as.matrix(covariates[, -1]) %*% beta2 + rnorm(1, mean = 0, sd = 1)
DGP3 <- as.matrix(covariates[, -1]) %*% beta3 + rnorm(1, mean = 0, sd = 1)
DGP4 <- as.matrix(covariates[, -1]) %*% beta4 + rnorm(1, mean = 0, sd = 1) 
DGP5 <- as.matrix(covariates[, -1]) %*% beta5 + rnorm(1, mean = 0, sd = 1)
DGP6 <- as.matrix(covariates[, -1]) %*% beta6 + rnorm(1, mean = 0, sd = 1) 
DGP7 <- as.matrix(covariates[, -1]) %*% beta7 + rnorm(1, mean = 0, sd = 1)

# Combine Y and X - there is most likely a neater (more elegant) way to code this up 
data1 <- as.tibble(cbind(DGP1, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP1") %>% 
  rename(DGP = DGP1)

data2 <- as.tibble(cbind(DGP2, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP2") %>% 
  rename(DGP = DGP2)

data3 <- as.tibble(cbind(DGP3, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP3") %>% 
  rename(DGP = DGP3)   

data4 <- as.tibble(cbind(DGP4, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP4") %>% 
  rename(DGP = DGP4)  

data5 <- as.tibble(cbind(DGP5, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP5") %>% 
  rename(DGP = DGP5)

data6 <- as.tibble(cbind(DGP6, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP6") %>% 
  rename(DGP = DGP6)

data7 <- as.tibble(cbind(DGP7, covariates[, -1])) %>% 
  add_column(covariates[, 1], .before = "DGP7") %>% 
  rename(DGP = DGP7)

### 4. Experiment 1 ------------------------------------------------------------

# Experiment #1 
RL_mse <- function(RL = 0, data, folds = 10, output = 1, split = 0.10, radius = 1) { # RL is a dummy indicating ridge or lasso. radius indicates data set. 

  # The function RL_mse returns one MSE the 200 iterations of either the ridge or lasso subset regression methodology. 
  # The chosen regression has been subject to 10 folds cross valiation (CV). 
  
  # RL = 0 --> ridge; RL = 1 --> lasso
  
  # Initialize empty holder for MSEs
  MSEs <- NULL
  
  # # Assign data - old with radius input 
  # ifelse(radius == 1, data <- data1, 
  #        ifelse(radius == 2, data <- data2, 
  #               ifelse(radius == 3, data <- data3, 
  #                      ifelse(radius == 4, data <- data4,
  #                             ifelse(radius == 5, data <- data5,
  #                                    ifelse(radius == 6, data <- data6, 
  #                                           ifelse(radius == 7, data <- data7, data <- NULL)))))))

  for (i in 1:nrow(unique(data[, 1]))) { # 200 iterations 
    
    # Initialize data holder for each iteration
    data_ite <- data %>% 
      filter(Iteration == i)
    
    # Separate data in training and test set 90-10 split 
    test_index <- sample(nrow(data_ite), size = round(split * nrow(data_ite)))
    test_data <- data_ite[test_index, ] %>% select(-Iteration)
    train_data <- data_ite[-test_index, ] %>% select(-Iteration)
    
    # Perform CV on training set
    CV <- cv.glmnet(x = as.matrix(train_data[, 2:41]), y = as.matrix(train_data[, 1]), alpha = RL, nfolds = folds)
    
    # Carry out subset regression
    RLreg <- glmnet::glmnet(y = as.matrix(train_data$DGP), x = as.matrix(train_data[, 2:41]), alpha = RL, lambda = CV$lambda.min)
    
    # Compute MSE
    RL_hat <- predict(RLreg, as.matrix(test_data[, 2:41]))
    RL_MSE <- apply((test_data$DGP - RL_hat) ** 2, MARGIN = 2, FUN = mean)

    # Store MSE
    MSEs <- rbind(MSEs, RL_MSE)
    
  }
  
  # Compute mean MSE of 200 iterations
  MSE <- apply(MSEs, MARGIN = 2, FUN = mean)
  
  if (output == 1) {
    return(MSE)
  }
  
  if (output == 0) {
    
    # Compability with ggplot (figure 3)
    MSEs <- as.tibble(MSEs) %>% 
      mutate(iteration = seq(200), cluster_radius = radius)

    return(MSEs)
    
  }
  
  

  # Det kan være jeg skal udvide funktionen til også at klar forward eller backward regression - jeg kan gøre det med dummy, hvor jeg laver if
  # statement, og så tager hele ridge lasso del af funktionen hvis dummy = 1 og eller så tager jeg back stepwise regression hvis dummy = 0
  
}

# Ridge MSE - Radius 1 through 7
ridge1 <- RL_mse(RL = 0, data = data1, folds = 10)
ridge2 <- RL_mse(RL = 0, data = data2, folds = 10)
ridge3 <- RL_mse(RL = 0, data = data3, folds = 10)
ridge4 <- RL_mse(RL = 0, data = data4, folds = 10)
ridge5 <- RL_mse(RL = 0, data = data5, folds = 10)
ridge6 <- RL_mse(RL = 0, data = data6, folds = 10)
ridge7 <- RL_mse(RL = 0, data = data7, folds = 10)

# Lasso MSE - Radius 1 through 7
lasso1 <- RL_mse(RL = 1, data = data1, folds = 10, split = 0.15)
lasso2 <- RL_mse(RL = 1, data = data2, folds = 10, split = 0.15)
lasso3 <- RL_mse(RL = 1, data = data3, folds = 10, split = 0.15)
lasso4 <- RL_mse(RL = 1, data = data4, folds = 10, split = 0.15)
lasso5 <- RL_mse(RL = 1, data = data5, folds = 10, split = 0.15)
lasso6 <- RL_mse(RL = 1, data = data6, folds = 10, split = 0.15)
lasso7 <- RL_mse(RL = 1, data = data7, folds = 10, split = 0.15)

# Compile data in single tibble 
ridge_data <- rbind(ridge1, ridge2, ridge3, ridge4, ridge5, ridge6, ridge7)
lasso_data <- rbind(lasso1, lasso2, lasso3, lasso4, lasso5, lasso6, lasso7)
RL_MSE <- as.tibble(cbind(ridge_data, lasso_data)) %>% 
  rename(ridge_MSE = s0, lasso_MSE = V2)


# MSEs for figure 3 - Ridge
ridge1_mses <- RL_mse(RL = 0, data = data1, folds = 10, output = 0, radius = 1)
ridge2_mses <- RL_mse(RL = 0, data = data2, folds = 10, output = 0, radius = 2)
ridge3_mses <- RL_mse(RL = 0, data = data3, folds = 10, output = 0, radius = 3)
ridge4_mses <- RL_mse(RL = 0, data = data4, folds = 10, output = 0, radius = 4)
ridge5_mses <- RL_mse(RL = 0, data = data5, folds = 10, output = 0, radius = 5)
ridge6_mses <- RL_mse(RL = 0, data = data6, folds = 10, output = 0, radius = 6)
ridge7_mses <- RL_mse(RL = 0, data = data7, folds = 10, output = 0, radius = 7)

# MSEs for figure 3 - Lasso
lasso1_mses <- RL_mse(RL = 1, data = data1, folds = 10, output = 0, split = 0.15, radius = 1)
lasso2_mses <- RL_mse(RL = 1, data = data2, folds = 10, output = 0, split = 0.15, radius = 2)
lasso3_mses <- RL_mse(RL = 1, data = data3, folds = 10, output = 0, split = 0.15, radius = 3)
lasso4_mses <- RL_mse(RL = 1, data = data4, folds = 10, output = 0, split = 0.15, radius = 4)
lasso5_mses <- RL_mse(RL = 1, data = data5, folds = 10, output = 0, split = 0.15, radius = 5)
lasso6_mses <- RL_mse(RL = 1, data = data6, folds = 10, output = 0, split = 0.15, radius = 6)
lasso7_mses <- RL_mse(RL = 1, data = data7, folds = 10, output = 0, split = 0.15, radius = 7)

# Compile MSEs in single tibble
ridge_data_figure3 <- rbind(ridge1_mses, ridge2_mses, ridge3_mses, ridge4_mses, ridge5_mses, ridge6_mses, ridge7_mses) %>% 
  arrange(iteration, cluster_radius)

lasso_data_figure3 <- rbind(lasso1_mses, lasso2_mses, lasso3_mses, lasso4_mses, lasso5_mses, lasso6_mses, lasso7_mses) %>% 
  arrange(iteration, cluster_radius)

# # Compile data in single tibble for old figure 3
# ridge_data_figure3 <- as.tibble(cbind(ridge1_mses, ridge2_mses, ridge3_mses, ridge4_mses, ridge5_mses, ridge6_mses, ridge7_mses)) %>% 
#   `colnames<-`(c("1", "2", "3", "4", "5", "6", "7"))

# # Figure 3 - MSE through cluster radius 
# figure3_old <- ggplot() +
#   geom_line(data = RL_MSE, aes(y = ridge_MSE, x = seq(7), col = "Ridge"), size = 1) +
#   geom_line(data = RL_MSE, aes(y = lasso_MSE, x = seq(7), col = "Lasso"), size = 1) +
#   geom_point(data = RL_MSE, aes(y = ridge_MSE, x = seq(7))) +
#   geom_point(data = RL_MSE, aes(y = lasso_MSE, x = seq(7))) +
#   theme_article() +
#   labs(x = "Cluster Radius", y = "MSE", col='') + # Legend title
#   theme(legend.position = "bottom") + 
#   theme(text = element_text(size = 15))
# 
# ggsave("figure3_old.pdf", plot=figure3_old, width = 20, height = 18, units= "cm", dpi = 300)

# Figure 3 - Experiment 1 output (MSE through cluster radius)
figure_3A <- ggplot(ridge_data_figure3) +
  geom_line(aes(cluster_radius, s0, color = as.factor(iteration)), alpha = 0.1, size = 0.01) +
  scale_color_discrete(name = "Date") + 
  scale_y_continuous(limits = c(0, 80)) + # Removes extreme observations above 80 MSE
  theme_article() + 
  theme(legend.position = "none") + # Removes legend altogether 
  scale_color_manual(values = rep("#00BFC4", 201)) + 
  labs(x = "Cluster Radius", y = "MSE") +
  geom_line(data = RL_MSE, aes(y = ridge_MSE, x = seq(7), col = "Ridge"), size = 1) # Add average MSE 

figure_3B <- ggplot(lasso_data_figure3) +
  geom_line(aes(cluster_radius, s0, color = as.factor(iteration)), alpha = 0.1, size = 0.01) +
  scale_color_discrete(name = "Date") + 
  scale_y_continuous(limits = c(0, 80)) + # Removes extreme observations above 80 MSE
  theme_article() + 
  theme(legend.position = "none") + # Removes legend altogether 
  scale_color_manual(values = rep("#F8766D", 201)) + 
  labs(x = "Cluster Radius", y = "MSE") +
  geom_line(data = RL_MSE, aes(y = lasso_MSE, x = seq(7), col = "Lasso"), size = 1) # Add average MSE 

figure3 <- gridExtra::grid.arrange(figure_3A, figure_3B, ncol = 2)
ggsave("figure3.pdf", plot=figure3, width = 25, height = 10, units= "cm", dpi = 300)

#Table 3 - Experiment 1 output (MSE through cluster radius) 
# The table is modified slightly in LaTeX to allow for the multi-column
table3_names <- c("1", "2", "3", "4", "5", "6", "7" )
table3 <- table_theme(t(RL_MSE), colnames = table3_names, caption = "Experiment 1: MSE through coefficent cluster radius'", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")
  
### 5. Experiment 2 ------------------------------------------------------------

# Initiate empty list for figure 4
figure4 <- list()
table4_data <- as_tibble(data.frame(matrix(nrow = 0, ncol = 6 )))

# Covariates for experiment 2 
exp2_covariates50 <- simCovariates(d = 50, var = 40, ite = 200) 
exp2_covariates60 <- simCovariates(d = 60, var = 40, ite = 200) 
exp2_covariates70 <- simCovariates(d = 70, var = 40, ite = 200) 
exp2_covariates80 <- simCovariates(d = 80, var = 40, ite = 200) 
exp2_covariates90 <- simCovariates(d = 90, var = 40, ite = 200)
exp2_covariates100 <- simCovariates(d = 100, var = 40, ite = 200)

# Experiment 2 list of MSEs for density plot later 
figure8_data <- list() 

# Experiment 2 loop for Ridge
exp2_radius <- list(beta1, beta3, beta5, beta7)
for (i in enumerate(exp2_radius)) {
  
  index <- i[[1]]; beta <- i[[2]]
  
  # DGP - given cluster size but differnt number of draws 
  exp2_DGP50 <- as.matrix(exp2_covariates50[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP60 <- as.matrix(exp2_covariates60[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP70 <- as.matrix(exp2_covariates70[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP80 <- as.matrix(exp2_covariates80[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP90 <- as.matrix(exp2_covariates90[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP100 <- as.matrix(exp2_covariates100[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  
  # Combine Y and X
  exp2_data50 <- as.tibble(cbind(exp2_DGP50, exp2_covariates50[, -1])) %>% 
    add_column(exp2_covariates50[, 1], .before = "exp2_DGP50") %>% 
    rename(DGP = exp2_DGP50)
  
  exp2_data60 <- as.tibble(cbind(exp2_DGP60, exp2_covariates60[, -1])) %>% 
    add_column(exp2_covariates60[, 1], .before = "exp2_DGP60") %>% 
    rename(DGP = exp2_DGP60)   
  
  exp2_data70 <- as.tibble(cbind(exp2_DGP70, exp2_covariates70[, -1])) %>% 
    add_column(exp2_covariates70[, 1], .before = "exp2_DGP70") %>% 
    rename(DGP = exp2_DGP70)
  
  exp2_data80 <- as.tibble(cbind(exp2_DGP80, exp2_covariates80[, -1])) %>% 
    add_column(exp2_covariates80[, 1], .before = "exp2_DGP80") %>% 
    rename(DGP = exp2_DGP80)
  
  exp2_data90 <- as.tibble(cbind(exp2_DGP90, exp2_covariates90[, -1])) %>% 
    add_column(exp2_covariates90[, 1], .before = "exp2_DGP90") %>% 
    rename(DGP = exp2_DGP90)
  
  exp2_data100 <- as.tibble(cbind(exp2_DGP100, exp2_covariates100[, -1])) %>% 
    add_column(exp2_covariates100[, 1], .before = "exp2_DGP100") %>% 
    rename(DGP = exp2_DGP100)
  
  # Average MSE of Ridge for given cluster radius but with different number of draws 
  exp2_RL50 <- RL_mse(RL = 0, data = exp2_data50, folds = 10)
  exp2_RL60 <- RL_mse(RL = 0, data = exp2_data60, folds = 10)
  exp2_RL70 <- RL_mse(RL = 0, data = exp2_data70, folds = 10)
  exp2_RL80 <- RL_mse(RL = 0, data = exp2_data80, folds = 10)
  exp2_RL90 <- RL_mse(RL = 0, data = exp2_data90, folds = 10)
  exp2_RL100 <- RL_mse(RL = 0, data = exp2_data100, folds = 10)
  
  # Combine average MSE
  exp2_RL_mse <- as.tibble(rbind(exp2_RL50, exp2_RL60, exp2_RL70, exp2_RL80, exp2_RL90, exp2_RL100))
  
  # Bind in tibble for table 4
  table4_data <- rbind(table4_data, t(exp2_RL_mse))
  
  # 200 MSEs of Ridge for given cluster radius but with different number of draws
  exp2_RL50_mses <- RL_mse(RL = 0, data = exp2_data50, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 50)
  exp2_RL60_mses <- RL_mse(RL = 0, data = exp2_data60, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 60)
  exp2_RL70_mses <- RL_mse(RL = 0, data = exp2_data70, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 70)
  exp2_RL80_mses <- RL_mse(RL = 0, data = exp2_data80, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 80)
  exp2_RL90_mses <- RL_mse(RL = 0, data = exp2_data90, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 90)
  exp2_RL100_mses <- RL_mse(RL = 0, data = exp2_data100, folds = 10, output = 0, radius = 1) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 100)
  
  # Append MSEs to list so as to make figure 8 (density plot) later
  list_temp <- list(exp2_RL50_mses, exp2_RL60_mses, exp2_RL70_mses, exp2_RL80_mses, exp2_RL90_mses, exp2_RL100_mses)
  figure8_data[[index]] <- list_temp
  
  # Combine 200 MSEs
  exp2_ridge_data_figure4 <- rbind(exp2_RL50_mses, exp2_RL60_mses, exp2_RL70_mses, exp2_RL80_mses, exp2_RL90_mses, exp2_RL100_mses) %>% 
    arrange(iteration, obs)
  
  # Figure 4 
  figure4_temp <- ggplot(exp2_ridge_data_figure4) +
    geom_line(aes(obs, s0, color = as.factor(iteration)), alpha = 0.1, size = 0.01) +
    scale_color_discrete(name = "Date") + 
    scale_y_continuous(limits = c(0, 80)) + # Removes extreme observations above 80 MSE
    theme_article() + 
    theme(legend.position = "none") + # Removes legend altogether 
    scale_color_manual(values = rep("#00BFC4", 201)) + 
    labs(x = "Observations", y = "MSE") +
    geom_line(data = exp2_RL_mse, aes(y = s0, x = c(50, 60, 70, 80, 90, 100), col = "Ridge"), size = 1) # Add average MSE 
  
  figure4[[index]] <- figure4_temp
  
  print(index)
  
}

# Figure 4
figure4_stacked <- gridExtra::grid.arrange(figure4[[1]], figure4[[2]], figure4[[3]], figure4[[4]], ncol = 2)
ggsave("figure4.pdf", plot=figure4_stacked, width = 20, height = 18, units= "cm", dpi = 300)

# Table 4 
table4_names <- c("50", "60", "70", "80", "90", "100")
table4 <- table_theme(table4_data %>% `rownames<-`(c("1", "3", "5", "7")), colnames = table4_names, 
                            caption = "Experiment 2 output: Ridge MSE", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")

# Holder for table5 data 
table5_data <- as_tibble(data.frame(matrix(nrow = 0, ncol = 6 )))

# Experiment 2 list of MSEs for density plot later 
figure9_data <- list() 

# Experiment 2 loop for Lasso
for (i in enumerate(exp2_radius)) {
  
  index <- i[[1]]; beta <- i[[2]]
  
  # DGP - given cluster size but different number of draws 
  exp2_DGP50 <- as.matrix(exp2_covariates50[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP60 <- as.matrix(exp2_covariates60[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP70 <- as.matrix(exp2_covariates70[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP80 <- as.matrix(exp2_covariates80[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP90 <- as.matrix(exp2_covariates90[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  exp2_DGP100 <- as.matrix(exp2_covariates100[, -1]) %*% beta + rnorm(1, mean = 0, sd = 1)
  
  # Combine Y and X
  exp2_data50 <- as.tibble(cbind(exp2_DGP50, exp2_covariates50[, -1])) %>% 
    add_column(exp2_covariates50[, 1], .before = "exp2_DGP50") %>% 
    rename(DGP = exp2_DGP50)
  
  exp2_data60 <- as.tibble(cbind(exp2_DGP60, exp2_covariates60[, -1])) %>% 
    add_column(exp2_covariates60[, 1], .before = "exp2_DGP60") %>% 
    rename(DGP = exp2_DGP60)   
  
  exp2_data70 <- as.tibble(cbind(exp2_DGP70, exp2_covariates70[, -1])) %>% 
    add_column(exp2_covariates70[, 1], .before = "exp2_DGP70") %>% 
    rename(DGP = exp2_DGP70)
  
  exp2_data80 <- as.tibble(cbind(exp2_DGP80, exp2_covariates80[, -1])) %>% 
    add_column(exp2_covariates80[, 1], .before = "exp2_DGP80") %>% 
    rename(DGP = exp2_DGP80)
  
  exp2_data90 <- as.tibble(cbind(exp2_DGP90, exp2_covariates90[, -1])) %>% 
    add_column(exp2_covariates90[, 1], .before = "exp2_DGP90") %>% 
    rename(DGP = exp2_DGP90)
  
  exp2_data100 <- as.tibble(cbind(exp2_DGP100, exp2_covariates100[, -1])) %>% 
    add_column(exp2_covariates100[, 1], .before = "exp2_DGP100") %>% 
    rename(DGP = exp2_DGP100)
  
  # Average MSE of Lasso for given cluster radius but with different number of draws 
  exp2_RL50 <- RL_mse(RL = 1, data = exp2_data50, folds = 10, split = 0.15)
  exp2_RL60 <- RL_mse(RL = 1, data = exp2_data60, folds = 10, split = 0.15)
  exp2_RL70 <- RL_mse(RL = 1, data = exp2_data70, folds = 10, split = 0.15)
  exp2_RL80 <- RL_mse(RL = 1, data = exp2_data80, folds = 10, split = 0.15)
  exp2_RL90 <- RL_mse(RL = 1, data = exp2_data90, folds = 10, split = 0.15)
  exp2_RL100 <- RL_mse(RL = 1, data = exp2_data100, folds = 10, split = 0.15)
  
  # Combine average MSE
  exp2_RL_mse <- as.tibble(rbind(exp2_RL50, exp2_RL60, exp2_RL70, exp2_RL80, exp2_RL90, exp2_RL100))
  
  # Bind in tibble for table 4
  table5_data <- rbind(table5_data, t(exp2_RL_mse))
  
  # 200 MSEs of Lasso for given cluster radius but with different number of draws
  exp2_RL50_mses <- RL_mse(RL = 1, data = exp2_data50, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 50)
  exp2_RL60_mses <- RL_mse(RL = 1, data = exp2_data60, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 60)
  exp2_RL70_mses <- RL_mse(RL = 1, data = exp2_data70, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 70)
  exp2_RL80_mses <- RL_mse(RL = 1, data = exp2_data80, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 80)
  exp2_RL90_mses <- RL_mse(RL = 1, data = exp2_data90, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 90)
  exp2_RL100_mses <- RL_mse(RL = 1, data = exp2_data100, folds = 10, output = 0, radius = 1, split = 0.15) %>% 
    select(-cluster_radius) %>% 
    mutate(obs = 100)
  
  # Append MSEs to list so as to make figure 8 (density plot) later
  list_temp <- list(exp2_RL50_mses, exp2_RL60_mses, exp2_RL70_mses, exp2_RL80_mses, exp2_RL90_mses, exp2_RL100_mses)
  figure9_data[[index]] <- list_temp
  
  # Combine 200 MSEs
  exp2_ridge_data_figure4 <- rbind(exp2_RL50_mses, exp2_RL60_mses, exp2_RL70_mses, exp2_RL80_mses, exp2_RL90_mses, exp2_RL100_mses) %>% 
    arrange(iteration, obs)
  
  # Figure 4 
  figure4_temp <- ggplot(exp2_ridge_data_figure4) +
    geom_line(aes(obs, s0, color = as.factor(iteration)), alpha = 0.1, size = 0.01) +
    scale_color_discrete(name = "Date") + 
    scale_y_continuous(limits = c(0, 80)) + # Removes extreme observations above 80 MSE
    theme_article() + 
    theme(legend.position = "none") + # Removes legend altogether 
    scale_color_manual(values = rep("#F8766D", 201)) + 
    labs(x = "Observations", y = "MSE") +
    geom_line(data = exp2_RL_mse, aes(y = s0, x = c(50, 60, 70, 80, 90, 100), col = "Ridge"), size = 1) # Add average MSE 
  
  figure4[[index+4]] <- figure4_temp # + 4 so as to make sure they are stored after the 4 ridge plots
  
  print(index)
  
}

# Figure 5
figure5_stacked <- gridExtra::grid.arrange(figure4[[5]], figure4[[6]], figure4[[7]], figure4[[8]], ncol = 2)
ggsave("figure5.pdf", plot=figure5_stacked, width = 20, height = 18, units= "cm", dpi = 300)

# Table 5
table5_names <- c("50", "60", "70", "80", "90", "100")
table5 <- table_theme(table5_data %>% `rownames<-`(c("1", "3", "5", "7")), colnames = table5_names, 
                      caption = "Experiment 2 output: Lasso MSE", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")


### 6. Distribution figures of MSEs of experiment 1 and 2 ----------------------

# Experiment 1 - Density plots for ridge and lasso (figure 6 and 7) 
figure6 <- list()
figure6_data <- list(ridge1_mses, ridge2_mses, ridge3_mses, ridge4_mses, ridge5_mses, ridge6_mses, ridge7_mses)

for (i in enumerate(figure6_data)) {
  
  index <- i[[1]]; data <- i[[2]]
  
  figure6_temp <- ggplot(data, aes(x = s0, fill = variable)) +
    geom_density(alpha = .2, fill = "#00BFC4") + 
    labs(x = "Values", y = "Density") +
    theme_article()
  
  figure6[[index]] <- figure6_temp
  
}

layout <- rbind(c(1,1,2,2), c(3,3,4,4), c(5,5,6,6), c(8,7,7,8))
figure6_stacked <- gridExtra::grid.arrange(figure6[[1]], figure6[[2]], figure6[[3]], figure6[[4]], figure6[[5]], figure6[[6]], figure6[[7]], 
                                           ncol = 2, nrow = 4, layout_matrix = layout)
ggsave("figure6.pdf", plot=figure6_stacked, width = 25, height = 45, units= "cm", dpi = 300)

figure7 <- list()
figure7_data <- list(lasso1_mses, lasso2_mses, lasso3_mses, lasso4_mses, lasso5_mses, lasso6_mses, lasso7_mses)

for (i in enumerate(figure7_data)) {
  
  index <- i[[1]]; data <- i[[2]]
  
  figure7_temp <- ggplot(data, aes(x = s0, fill = variable)) +
    geom_density(alpha = .2, fill = "#FF6666") + 
    labs(x = "Values", y = "Density") +
    theme_article()
  
  figure7[[index]] <- figure7_temp
  
}

layout <- rbind(c(1,1,2,2), c(3,3,4,4), c(5,5,6,6), c(8,7,7,8))
figure7_stacked <- gridExtra::grid.arrange(figure7[[1]], figure7[[2]], figure7[[3]], figure7[[4]], figure7[[5]], figure7[[6]], figure7[[7]], 
                                           ncol = 2, nrow = 4, layout_matrix = layout)
ggsave("figure7.pdf", plot=figure7_stacked, width = 25, height = 45, units= "cm", dpi = 300)

# Table of Standard Deviations for ridge and lasso for experiment 1 
# Table 6 - Compile data in single tibble
colnames <- c("r1", "r2", "r3", "r4", "r5", "r6", "r7")

ridge_data_sd <- cbind(ridge1_mses[, 1], ridge2_mses[, 1], ridge3_mses[, 1], ridge4_mses[, 1], ridge5_mses[, 1], ridge6_mses[, 1], ridge7_mses[, 1]) %>% 
  `colnames<-`(colnames) %>% 
  apply(MARGIN = 2, FUN = sd)
lasso_data_sd <- cbind(lasso1_mses[, 1], lasso2_mses[, 1], lasso3_mses[, 1], lasso4_mses[, 1], lasso5_mses[, 1], lasso6_mses[, 1], lasso7_mses[, 1]) %>% 
  `colnames<-`(colnames) %>% 
  apply(MARGIN = 2, FUN = sd)

table6_volatility_data <- as.tibble(cbind(ridge_data_sd, lasso_data_sd))

table6_names <- c("1", "2", "3", "4", "5", "6", "7" )
table6 <- table_theme(t(table6_volatility_data), colnames = table4_names, caption = "Experiment 1: Std. Dev. versus coefficent cluster radius'", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")

# Experiment 2 - Density plots for ridge and lasso for each complexity 1, 3, 5, and 7 
# Data for density plots of ridge and lasso from experiment 2 loop 
figure8_data
figure9_data 

# First 6 entries are for cluster radius 1 with observations 50, 60, 70, 80, 90, and 100. The next 6 are for cluster radius 2 with observations... (ridge) 
figure8 <- list()

index <- 0 
for (t in 1:length(figure8_data)){
  
  figure_ite <- figure8_data[[t]]
  
  for (i in enumerate(figure_ite)) {
    
    index <- index + 1; data <- i[[2]]
    
    figure8_temp <- ggplot(data, aes(x = s0, fill = variable)) +
      geom_density(alpha = .2, fill = "#00BFC4") + 
      labs(x = "Values", y = "Density") +
      theme_article()
    
    figure8[[index]] <- figure8_temp

  }
  
}

layout_figure8 <- rbind(c(1,2,3), c(4,5,6))
figure8A_stacked <- gridExtra::grid.arrange(figure8[[1]], figure8[[2]], figure8[[3]], figure8[[4]], figure8[[5]], figure8[[6]], ncol = 3, nrow = 2, layout_matrix = layout_figure8)
ggsave("figure8A.pdf", plot=figure8A_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure8B_stacked <- gridExtra::grid.arrange(figure8[[7]], figure8[[8]], figure8[[9]], figure8[[10]], figure8[[11]], figure8[[12]], ncol = 3, nrow = 2, layout_matrix = layout_figure8)
ggsave("figure8B.pdf", plot=figure8B_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure8C_stacked <- gridExtra::grid.arrange(figure8[[13]], figure8[[14]], figure8[[15]], figure8[[16]], figure8[[17]], figure8[[18]], ncol = 3, nrow = 2, layout_matrix = layout_figure8)
ggsave("figure8C.pdf", plot=figure8C_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure8D_stacked <- gridExtra::grid.arrange(figure8[[19]], figure8[[20]], figure8[[21]], figure8[[22]], figure8[[23]], figure8[[24]], ncol = 3, nrow = 2, layout_matrix = layout_figure8)
ggsave("figure8D.pdf", plot=figure8D_stacked, width = 35, height = 21, units= "cm", dpi = 300)

# Lasso Density plots experiment 2 
figure9 <- list()

index <- 0 
for (t in 1:length(figure9_data)){
  
  figure_ite <- figure9_data[[t]]
  
  for (i in enumerate(figure_ite)) {
    
    index <- index + 1; data <- i[[2]]
    
    figure9_temp <- ggplot(data, aes(x = s0, fill = variable)) +
      geom_density(alpha = .2, fill = "#FF6666") + 
      labs(x = "Values", y = "Density") +
      theme_article()
    
    figure9[[index]] <- figure9_temp
    
  }
  
}

layout_figure9 <- rbind(c(1,2,3), c(4,5,6))
figure9A_stacked <- gridExtra::grid.arrange(figure9[[1]], figure9[[2]], figure9[[3]], figure9[[4]], figure9[[5]], figure9[[6]], 
                                            ncol = 3, nrow = 2, layout_matrix = layout_figure8)
ggsave("figure9A.pdf", plot=figure9A_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure9B_stacked <- gridExtra::grid.arrange(figure9[[7]], figure9[[8]], figure9[[9]], figure9[[10]], figure9[[11]], figure9[[12]], ncol = 3, 
                                            nrow = 2, layout_matrix = layout_figure8)
ggsave("figure9B.pdf", plot=figure9B_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure9C_stacked <- gridExtra::grid.arrange(figure9[[13]], figure9[[14]], figure9[[15]], figure9[[16]], figure9[[17]], figure9[[18]], ncol = 3, 
                                            nrow = 2, layout_matrix = layout_figure8)
ggsave("figure9C.pdf", plot=figure9C_stacked, width = 35, height = 21, units= "cm", dpi = 300)

figure9D_stacked <- gridExtra::grid.arrange(figure9[[19]], figure9[[20]], figure9[[21]], figure9[[22]], figure9[[23]], figure9[[24]], ncol = 3, 
                                            nrow = 2, layout_matrix = layout_figure8)
ggsave("figure9D.pdf", plot=figure9D_stacked, width = 35, height = 21, units= "cm", dpi = 300)

# RIDGE: Tables of Standard deviation of ridge and lasso regression in experiment 2 
table7_temp1 <- data.frame(figure8_data[[1]][[1]][, 1], figure8_data[[1]][[2]][, 1], figure8_data[[1]][[3]][, 1], 
                           figure8_data[[1]][[4]][, 1], figure8_data[[1]][[5]][, 1], figure8_data[[1]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table7_temp2 <- data.frame(figure8_data[[2]][[1]][, 1], figure8_data[[2]][[2]][, 1], figure8_data[[2]][[3]][, 1], 
                           figure8_data[[2]][[4]][, 1], figure8_data[[2]][[5]][, 1], figure8_data[[2]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table7_temp3 <- data.frame(figure8_data[[3]][[1]][, 1], figure8_data[[3]][[2]][, 1], figure8_data[[3]][[3]][, 1], 
                           figure8_data[[3]][[4]][, 1], figure8_data[[3]][[5]][, 1], figure8_data[[3]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table7_temp4 <- data.frame(figure8_data[[4]][[1]][, 1], figure8_data[[4]][[2]][, 1], figure8_data[[4]][[3]][, 1], 
                           figure8_data[[4]][[4]][, 1], figure8_data[[4]][[5]][, 1], figure8_data[[4]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)

# Sd of given cluster radius (rows) and observations (columns) for the ridge regression in experiment 2 
table7_data <- rbind(t(table7_temp1), t(table7_temp2), t(table7_temp3), t(table7_temp4))
table7_names <- c("50", "60", "70", "80", "90", "100")
table7 <- table_theme(table7_data %>% `row.names<-`(c("1", "3", "5", "7")), colnames = table7_names,
                      caption = "Experiment 2: Std. Dev. for given radis and number of draws", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")

# LASSO: Tables of Standard deviation of ridge and lasso regression in experiment 2 
table8_temp1 <- data.frame(figure9_data[[1]][[1]][, 1], figure9_data[[1]][[2]][, 1], figure9_data[[1]][[3]][, 1], 
                           figure9_data[[1]][[4]][, 1], figure9_data[[1]][[5]][, 1], figure9_data[[1]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table8_temp2 <- data.frame(figure9_data[[2]][[1]][, 1], figure9_data[[2]][[2]][, 1], figure9_data[[2]][[3]][, 1], 
                           figure9_data[[2]][[4]][, 1], figure9_data[[2]][[5]][, 1], figure9_data[[2]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table8_temp3 <- data.frame(figure9_data[[3]][[1]][, 1], figure9_data[[3]][[2]][, 1], figure9_data[[3]][[3]][, 1], 
                           figure9_data[[3]][[4]][, 1], figure9_data[[3]][[5]][, 1], figure9_data[[3]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)
table8_temp4 <- data.frame(figure9_data[[4]][[1]][, 1], figure9_data[[4]][[2]][, 1], figure9_data[[4]][[3]][, 1], 
                           figure9_data[[4]][[4]][, 1], figure9_data[[4]][[5]][, 1], figure9_data[[4]][[6]][, 1]) %>% 
  map(sd) %>% map(round, digits = 2)

# Sd of given cluster radius (rows) and observations (columns) for the ridge regression in experiment 2 
table8_data <- rbind(t(table8_temp1), t(table8_temp2), t(table8_temp3), t(table8_temp4))
table8_names <- c("50", "60", "70", "80", "90", "100")
table8 <- table_theme(table8_data %>% `row.names<-`(c("1", "3", "5", "7")), colnames = table8_names,
                      caption = "Experiment 2: Std. Dev. for given radis and number of draws", escape = TRUE) %>% 
  kable_styling(latex_options = "scale_down")


### 7. Application: Empirical Asset Pricing ------------------------------------







### 8. Run Diagnostics ---------------------------------------------------------

# Script does not take an input (data) which is why I only run diagnostics throughout as I generate the data myself 

#if(length(unique(table(dat$date))) != 1){
#  log <- c(log, "-- Not all dates have equal number of observations --")
#}






# Timer finished
end_time <- Sys.time()

print(paste("Total time:", end_time - start_time))

bind_rows(data, test)







### 9. OLD CODE  ----------------------

# Ridge & Lasso
ridgeeq <- glmnet::glmnet(y = as.matrix(train_data$DGP), x = as.matrix(train_data[, 3:42]), alpha = 0, lambda = 1)
ridge_hat <- predict(ridgeeq, as.matrix(test_data[, 3:42]))
ridge_mse <- apply((test_data$DGP - ridge_hat) ** 2, MARGIN = 2, FUN = mean)


# Test for at vise matrix regning er korrekt - det var den! 
data_test <- matrix(1:9, nrow  = 3)
covariates_test <- matrix(3:5, nrow = 3)

y <- data_test %*% covariates_test + 1 # +1 er error term holder

# Density Plots - not used 

ggplot(ridge1_mses, aes(x = s0)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white",  bins = 40) +
  geom_density(alpha = .2, fill = "#FF6666") + 
  labs(x = "Values", y = "Density")

ggplot(ridge1_mses, aes(x = s0, fill = variable)) +
  geom_density(alpha = .2, fill = "#00BFC4") + 
  labs(x = "Values", y = "Density") +
  theme_article()

ggplot(lasso6_mses, aes(x = s0, fill = variable)) +
  geom_density(alpha = .2, fill = "#00BFC4") + 
  theme_article()

# Saml data til at lave density plots
data_test <- cbind(ridge6_mses[, 1], lasso6_mses[, 1]) %>% 
  `colnames<-`(c("ridge_mse", "lasso_mse")) %>% 
  stack(c(ridge_mse, lasso_mse))

# Farvekombination 1 
ggplot(data_test, aes(x = values, fill = ind, color = ind)) +
  geom_density(alpha = .4) + 
  scale_color_manual(values = c("#00BFC4", "#FF6666")) + 
  scale_fill_manual(values = c("white", "white")) + 
  theme_article() #+ 
#theme(legend.position = "none") # removes ledger altogether has to come after theme_article()

# Farvekomvinaiton 2 
ggplot(data_test, aes(x = values, fill = ind)) +
  geom_density(alpha = .4) + 
  scale_fill_manual(values = c("#FF6666", "#00BFC4")) + 
  theme_article() #+ 
#theme(legend.position = "none") # removes ledger altogether has to come after theme_article()