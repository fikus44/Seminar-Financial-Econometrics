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

sim_data <- function(d = draws, var = variables, ite = iterations) { # Draws, variables, iterations,
  
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
     
     x = rnorm(d, 0, runif(1, min = 0, max = 1) ** log(i)) 
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
covariates <- sim_data()

# Table 1 - Snip of covariates 

# Figure 1 - Distribution of Covariates
norm_dist <- to_vec(for (i in 1:variables) runif(1, min = 0, max = 1) ** log(i))
norm_dist_df <- as.tibble(norm_dist) %>%
  `colnames<-`("sd") %>% 
  add_column(mean = 0, .before = "sd")

# map2(.x = t(norm_dist_df)[1, ], .y = t(norm_dist_df)[2, ],
#      .f = ~{
#        
#        ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
#          stat_function(fun = dnorm, n = 101, args = list(mean = .x, sd = .y)) + ylab("") +
#          scale_y_continuous(breaks = NULL) + 
#          theme_article()
#        
#      })

ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +            
  map2(.x = t(norm_dist_df)[1, ], .y = t(norm_dist_df)[2, ],
       .f = ~stat_function(fun = dnorm, n = 101, args = list(mean = .x, sd = .y))) +
  scale_y_continuous(breaks = NULL) + 
  theme_article()

# Jeg skal have ændret varians i covariates, så den stiger lidt hver gang 



### 2. Create Coefficient Cluster ----------------------------------------------

coef_cluster <- function(betas = variables, center = c(1,3,5), radius = 5) {
  
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

# Generate coefficients
tester = coef_cluster()
tester2 = coef_cluster(radius = 3)

# Figure 2 - 2D plot of the coefficient cluster
figure_1 <- ggplot() + 
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "blue", alpha = 0.0, color = "white") +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 2, ymax = 4),
            fill = "blue", alpha = 0.5, color = "black", size = 0.4) + 
  labs(x = "Coefficient Cluster 1", y = "Coefficient Cluster 2") +
  theme_article()

figure_2 <- ggplot() +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = 0, ymax = 8), 
            fill = "blue", alpha = 0.0, color = "white")  +
  geom_rect(aes(xmin = 0, xmax = 2, ymin = 4, ymax = 6),
            fill = "orange", alpha = 0.8, color = "black", size = 0.4) +
  labs(x = "Coefficient Cluster 1", y = "Coefficient Cluster 3") +
  theme_article()

figure2 <- gridExtra::grid.arrange(figure_1, figure_2, ncol = 2)
ggsave("figure2.pdf", plot=figure2, width = 25, height = 10, units= "cm", dpi = 300)

### 3. Data Generating Process ----------------------------------------------

# Her skal vi have 5 y'er; en for hver radius. 

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
