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
#   * Hvordan skal variansen defienres for vores covariates?
#   * Skal alle covaraites være trukket fra normalfordelingen?
#   * Hvorfor tjekker vi om corr-matrix er positiv definint?
#   * Spørg til mødet om min ide med at generere x'er er god; det med at give dem en korrelation 

draws = 50
variables = 40
iterations = 200

sim_data <- function(draws, variables, iterations) { # Draws, variables, iterations,
  
  # Set seed
  set.seed(007)
  
  # Initialize data holder 
  data <- as_tibble(data.frame(matrix(nrow = 0, ncol = variables))) %>% 
    add_column(Iteration = NA, .before = "X1") # Column for each variable + an index
  
  for (i in 1:iterations) { # Loop f.o.m. 1 t.o.m. # of iterations
    
    # Initialize temp data holder
    data_temp = matrix(nrow = draws, ncol = variables) # Column for each variable 
    
   for (j in 1:variables) {
     
     x = rnorm(draws, 0, runif(1, min = 0, max = 1) ** log(i)) 
     data_temp[, j] = x
     
   }
  
    # Cholesky Decomposition to get independence 
    data_var <- var(data_temp)
    chol <- solve(chol(data_var), tol = 2.30e-18) # Set tolerance to circumvent linear dependence 
    data_temp <- data_temp %*% chol
    
    ### Diagnostics ### 
    if(all(round(var(data_temp)) != diag(variables))){
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
covariates <- sim_data(50, 40, 200)

# Table 1 - Snip of covariates 

# Figure 1 - Distribution of Covariates
norm_dist <- to_df(for (x in 1:variables))

ggplot



### 2. Create Coefficient Cluster ----------------------------------------------

# Ide: Hans coefficient clsuter giver ikke mening, men jeg vil gerne have hans resultater. 

 
# 1. Vi simulerer covariates
# 2. Vi laver coefficient cluster - draw med replacement sætter prob op så der er 2 * rc - 1 i hvert cluster, og hvis de trækker du ik nul --> giv den en værdi omkring center + et stochastic shock
# 3. Vi simulerer vores DGP
# 4. Vi har nu data --> vi laver Ridge, lasso, stepwise forward + CV og genererer MSE


### 2. Run Diagnostics -----------------------------------------------------

# Check dataset balanced (all dates have equal number of observations)
if(length(unique(table(dat$date))) != 1){
  log <- c(log, "-- Not all dates have equal number of observations --")
}



# Timer finished
end_time <- Sys.time()

print(paste("Total time:", end_time - start_time))

bind_rows(data, test)
