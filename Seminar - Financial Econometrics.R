##################################################
#                                                #
#  Seminar - Financial Econometrics              #
#                                                #
#  Thomas Theodor Kjølbye                        #
#                                                #
#  The following script produces all output      #
#  used in the paper. On my computer, the        #
#  entire script takes approx. 15 minute to run. #
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

# Hængepartier / Spørgsmål til del 1:
#   * Hvordan skal variansen defienres for vores covariates, lige nu med log og runif()?
#   * Skal alle covaraites være trukket fra normalfordelingen?
#   * Hvorfor tjekker vi om corr-matrix er positiv definint?
#   * Spørg til mødet om min ide med at generere x'er er god; det med at give dem en korrelation + alle er normalfordelte

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

# Figure 1 - Distribution of Covariates
norm_dist_df <- as.tibble(to_vec(for (i in 1:variables) runif(1, min = 1, max = 1.3) ** log(i))) %>% 
  `colnames<-`("sd") %>% 
  add_column(mean = 0, .before = "sd")

figure_1 <- ggplot(data = data.frame(x = c(-7, 7)), aes(x)) +            
  map2(.x = t(norm_dist_df)[1, ], .y = t(norm_dist_df)[2, ],
       .f = ~stat_function(fun = dnorm, n = 101, args = list(mean = .x, sd = .y))) +
  labs(y = "f(x)") + 
  theme_article()
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

# Betas - 5 set of betas one for each radius
beta1 <- coefCluster(radius = 1)
beta2 <- coefCluster(radius = 2)
beta3 <- coefCluster(radius = 3)
beta4 <- coefCluster(radius = 4)
beta5 <- coefCluster(radius = 5)
beta <- as.tibble(cbind(beta1, beta2, beta3, beta4, beta5)) # as.tibble and not matrix to make it compatible with e.g. map_df()

# Figure 2 - 2D plot of the coefficient cluster
figure_2A <- ggplot() + 
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "blue", alpha = 0.0, color = "white") +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 2, ymax = 4),
            fill = "blue", alpha = 0.5, color = "black", size = 0.4) + 
  labs(x = "Coefficient Cluster 1", y = "Coefficient Cluster 2") +
  theme_article()

figure_2B <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "blue", alpha = 0.0, color = "white")  +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 4, ymax = 6),
            fill = "orange", alpha = 0.8, color = "black", size = 0.4) +
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
                               
                             }) %>% add_column(seq(5), .before = "coeff"), 
                      colnames = table2_names, caption = "Coeffcients in each cluster") %>% 
  kable_styling(latex_options = "scale_down")


### 3. Data Generating Process -------------------------------------------------

# Hængepartier / Spørgsmål til del 3:
#   * jeg bør måske sætte seed her, så error term bliver ens? 
#   * Jeg kan formentlig også skrive det op på en lidt mere elegant måde? 

# DGP 
DGP1 <- as.matrix(covariates[, -1]) %*% beta1 + rnorm(1, mean = 0, sd = 1)
DGP2 <- as.matrix(covariates[, -1]) %*% beta2 + rnorm(1, mean = 0, sd = 1)
DGP3 <- as.matrix(covariates[, -1]) %*% beta3 + rnorm(1, mean = 0, sd = 1)
DGP4 <- as.matrix(covariates[, -1]) %*% beta4 + rnorm(1, mean = 0, sd = 1) 
DGP5 <- as.matrix(covariates[, -1]) %*% beta5 + rnorm(1, mean = 0, sd = 1)

# Combine Y and X
# Jeg tror jeg skal hav Y og X i et samlet datasæt efter at have kigget lidt på gammel kode hvor vi lave ridge og lasso. Jeg får så
# 5 datasæt i alt og jeg skal så lave ridge for hver iteration så jeg skal nok have iteration med som første søjle i hvert datasæt og så
# grouper jeg by iterationen og laver så ridge og lasso

# Test for at vise matrix regning er korrekt - det var den! 
data_test <- matrix(1:9, nrow  = 3)
covariates_test <- matrix(3:5, nrow = 3)

y <- data_test %*% covariates_test + 1 # +1 er error term holder

# Når jeg laver ridge og lasso så, så tror jeg, at jeg skal lave til tibble og så group_by iteration og så for hver group laver jeg det så. 



# Her skal vi have 5 y'er; en for hver radius. Vi får 1 y for hver observation af X, i.e. 50 y'er til hver iteration, vi har så 5 radiuser
# hvor vi gør det med hver af dem 

# 1. Vi simulerer covariates
# 2. Vi laver coefficient cluster - draw med replacement sætter prob op så der er 2 * rc - 1 i hvert cluster, og hvis de trækker du ik nul --> giv den en værdi omkring center + et stochastic shock
# 3. Vi simulerer vores DGP
# 4. Vi har nu data --> vi laver Ridge, lasso, stepwise forward + CV og genererer MSE

# Hvis jeg ikke når at kode hele simulationssetupet inden mødet, så skal jeg have skrrevet ned præcis hvordan jeg vil gøre, så jeg kan spørge! 


### 2. Run Diagnostics -----------------------------------------------------

# Check dataset balanced (all dates have equal number of observations)
if(length(unique(table(dat$date))) != 1){
  log <- c(log, "-- Not all dates have equal number of observations --")
}



# Timer finished
end_time <- Sys.time()

print(paste("Total time:", end_time - start_time))

bind_rows(data, test)
