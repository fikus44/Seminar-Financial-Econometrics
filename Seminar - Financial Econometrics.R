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

# Table 1 - Snip of covariates (Small)
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
                      colnames = table2_names, caption = "Coeffcient Clusters") %>% 
  kable_styling(latex_options = "scale_down")


### 3. Data Generating Process -------------------------------------------------

# Hængepartier / Spørgsmål til del 3:
#   * jeg bør måske sætte seed her, så error term bliver ens? 
#   * Jeg kan formentlig også skrive det op på en lidt mere elegant måde? 

# DGP (number corresponds to cluster radius) - there is most likely a neater (more elegant) way to code this up
DGP1 <- as.matrix(covariates[, -1]) %*% beta1 + rnorm(1, mean = 0, sd = 1)
DGP2 <- as.matrix(covariates[, -1]) %*% beta2 + rnorm(1, mean = 0, sd = 1)
DGP3 <- as.matrix(covariates[, -1]) %*% beta3 + rnorm(1, mean = 0, sd = 1)
DGP4 <- as.matrix(covariates[, -1]) %*% beta4 + rnorm(1, mean = 0, sd = 1) 
DGP5 <- as.matrix(covariates[, -1]) %*% beta5 + rnorm(1, mean = 0, sd = 1)

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

### 4. Machine Learning (Subsetting) ------------------------------------------

# Experiment #1 
RL_mse <- function(RL = 0, radius, folds = 10) { # RL is a dummy indicating ridge or lasso. radius indicates data set. 

  # The function RL_mse returns one MSE the 200 iterations of either the ridge or lasso subset regression methodology. 
  # The chosen regression has been subject to 10 folds cross valiation (CV). 
  
  # RL = 0 --> ridge; RL = 1 --> lasso
  
  # Initialize empty data holder
  data <- NULL
  MSEs <- NULL
  
  # Assign data
  ifelse(radius == 1, data <- data1, 
         ifelse(radius == 2, data <- data2, 
                ifelse(radius == 3, data <- data3, 
                       ifelse(radius == 4, data <- data4,
                              ifelse(radius == 5, data <- data5, data <- data1)))))
  
  for (i in 1:nrow(unique(data[, 1]))) { # 200 iterations 
    
    # Initialize data holder for each iteration
    data_ite <- data %>% 
      filter(Iteration == i)
    
    # Separate data in training and test set 80-20 split 
    test_index <- sample(nrow(data_ite), size = round(0.1 * nrow(data_ite)))
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
  
  return(MSE)
  

  # Det kan være jeg skal udvide funktionen til også at klar forward eller backward regression - jeg kan gøre det med dummy, hvor jeg laver if
  # statement, og så tager hele ridge lasso del af funktionen hvis dummy = 1 og eller så tager jeg back stepwise regression hvis dummy = 0
  
}

# Ridge - Radius 1 through 5
# * Ridge giver nogle meget store MSES, er det: 1) kodefejl?, 2) prøv at se hvis jeg gør clusters mean større, 3) hvad hvis vi øger obs?
# Jeg ændrede splittet fra 80-20 til 90-10 - hjalp meget! jeg skal måske lige tilføj et par obserevationr! 
Ridge1 <- RL_mse(RL = 0, radius = 1, folds = 10)
Ridge2 <- RL_mse(RL = 0, radius = 2, folds = 10)
Ridge3 <- RL_mse(RL = 0, radius = 3, folds = 10)
Ridge4 <- RL_mse(RL = 0, radius = 4, folds = 10)
Ridge5 <- RL_mse(RL = 0, radius = 5, folds = 10)

# Lasso - Radius 1 through 5
# Lasso har stadig 80-20 splittet 
Lasso1 <- RL_mse(RL = 1, radius = 1, folds = 10)
Lasso2 <- RL_mse(RL = 1, radius = 2, folds = 10)
Lasso3 <- RL_mse(RL = 1, radius = 3, folds = 10)
Lasso4 <- RL_mse(RL = 1, radius = 4, folds = 10)
Lasso5 <- RL_mse(RL = 1, radius = 5, folds = 10)


# Vi ser i det første eksperiment på hvordan ridge og lasso predicter igennem cluster radius. Men vi har ret få observationer; kun lidt flere end
# covariates - typisk vil vi have bare lidt flere og det kunne derfor være interessant at se hvor mange observationer vi egentlig skal bruge før
# vi får en form for konvergens; hvornår bliver ridge og lasso gode; hvor mange observationer før de bliver gode? og hvordan ser konvergensen
# ud? e.g. er de MEGET ringe i starten (som det viste ridge var med 40 obs) og så bliver de hurtigere bedre, eller er konvergensen mere lineær?
# Det kan som sagt tale til hvor mange obserevationer der skal til før vi kan bruge ridge og lasso til noget, i.e. i hvilke casees kan vi 
# bruge den metode til noget? Det kan i en situation sættes i perspektiv til andre ML metoder? f.eks. deep learning, som vist skal brugee
# mere data. 

# Det eksperiment bliver så for en given cluster størrelse - men jeg kan vælge flere, måske 1, 3 og 5 og så se hvordan konvergens er der. Jeg 
# computer så MSE gennem antal observationer

# Ridge & Lasso
ridgeeq <- glmnet::glmnet(y = as.matrix(train_data$DGP), x = as.matrix(train_data[, 3:42]), alpha = 0, lambda = 1)
ridge_hat <- predict(ridgeeq, as.matrix(test_data[, 3:42]))
ridge_mse <- apply((test_data$DGP - ridge_hat) ** 2, MARGIN = 2, FUN = mean)


# Test for at vise matrix regning er korrekt - det var den! 
data_test <- matrix(1:9, nrow  = 3)
covariates_test <- matrix(3:5, nrow = 3)

y <- data_test %*% covariates_test + 1 # +1 er error term holder



# 1. Vi simulerer covariates
# 2. Vi laver coefficient cluster - draw med replacement sætter prob op så der er 2 * rc - 1 i hvert cluster, og hvis de trækker du ik nul --> giv den en værdi omkring center + et stochastic shock
# 3. Vi simulerer vores DGP
# 4. Vi har nu data --> vi laver Ridge, lasso, stepwise forward + CV og genererer MSE



### 2. Run Diagnostics -----------------------------------------------------


# Tilføj noget med at diagnostics er blevet carried out i de forskellige dele af scripted da programmet her ikke taget noget input 
# Check dataset balanced (all dates have equal number of observations)
if(length(unique(table(dat$date))) != 1){
  log <- c(log, "-- Not all dates have equal number of observations --")
}



# Timer finished
end_time <- Sys.time()

print(paste("Total time:", end_time - start_time))

bind_rows(data, test)
