library(tidyverse) # For data wrangling
library(plyr) # For data wrangling
library(readxl) # For reading in data
library(missForest) # For random forest imputation
library(mice) # For imputation with MICE
library(glmnet) # For logistic regression and Cox
library(survminer) # For survival data
library(survival) # For survival data
library(survMisc) # For survival data
library(yardstick) # To calculate metrics
library(probably) # To calculate metrics
library(predRupdate) # For calibration
library(caret) # For model fitting and internal validation
library(DMwR2) # For SMOTE

##### C-index calculation #####
cindex_summary <- function(data, lev = NULL, model = NULL) {
  # Convert factors to numeric (0 and 1) for concordance calculation
  obs_numeric <- ifelse(data$obs == lev[1], 0, 1)
  pred_numeric <- as.numeric(data$pred)
  
  # Calculate Harrell's C-index using predicted probabilities
  cindex <- Hmisc::rcorr.cens(pred_numeric, obs_numeric)[2]
  out <- c(C = cindex)
  names(out) <- "C"
  out
}

# Calibration for binary outcome models
calibration_stats <- function(observed, predicted) {
  mean_observed <- mean(observed)
  mean_predicted <- mean(predicted)
  slope <- cov(predicted, observed) / var(predicted)
  intercept <- mean_observed - slope * mean_predicted
  return(list(intercept = intercept, slope = slope))
}

##### Logistic regression with repeated nested cross validation ##### 
LR_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  dca_clinical <- data.frame()
  dca_general <- data.frame()
  preds <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$chr, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>2) { # If there is more than one predictor, use LASSO
        combined_df[,c(1:16)] <- lapply(combined_df[,c(1:16)], factor)
        nzv <- nearZeroVar(combined_df, saveMetrics = TRUE) # Remove predictors with near zero variance
        combined_df2 <- combined_df[,nzv[,"nzv"] == FALSE] 
        options(na.action='na.pass')
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df2))
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      set.seed(123)
      if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
      imputed_data <- missForest(test)
      test <- imputed_data$ximp
      } else {
      imputed_data <- mice(test)
      test <- complete(imputed_data)
      }
      
      test$chr <- factor(test$chr)
      test$chr <- factor(make.names(levels(test$chr))[test$chr])
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      # If any remaining missing values with the column mean
      test[is.na(test)] <- apply(test, 2, function(col) mean(col, na.rm = TRUE))
      
      make.names(colnames(test), unique=TRUE)
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        set.seed(123)
        if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }

        train$chr <- factor(train$chr)
        train$chr <- factor(make.names(levels(train$chr))[train$chr])
        
        # Set up trainControl with custom summary function
        control <- trainControl(
          method = "cv",
          number = 5,
          summaryFunction = cindex_summary,
          classProbs = TRUE, # Needed if using probabilities
          sampling = "smote"
        )
        
        if (ncol(train)>2) { # If there is more than one predictor, use LASSO
          # Train the model
          tune_grid <- expand.grid(
            alpha = 1,             # Set alpha to 1 for LASSO
            lambda = seq(0.001, 0.1, length = 100) # Define a range for lambda
          )
          
          # Train the model with LASSO and C-index as the metric
          inner_model <- train(
            chr ~ .,
            data = train,
            method = "glmnet",
            trControl = control,
            metric = "C",
            tuneGrid = tune_grid,
            family = "binomial"
          )
        } else {
          # Train the model with logistic regression and C-index as the metric
          inner_model <- train(
            chr ~ .,
            data = train,
            method = "glm",
            trControl = control,
            metric = "C",
            family = "binomial"
          )
        }
        print(inner_model)
        
        if (ncol(train)>2) { 
        if (all(!is.na(inner_model$results$C))){
          # Store the F1 for each inner repeat
          inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                              data.frame(C = inner_model$results$C,
                                                         lambda = inner_model$bestTune$lambda,
                                                         repeat_number = inner_rep))
        }
      }
      }
      
      if (ncol(train)>3) { 
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      
      # Calculate the average C-index for each hyperparameter combination
      avg_results <- all_inner_results %>%
        summarise(avg_C = mean(C, na.rm=TRUE),
                  lambda = mean(lambda, na.rm=TRUE),
                  .groups = "keep") 
      }
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none',
                                    sampling = "smote") 
      
      if (ncol(train)>3) { 
      final_tune <- expand.grid(
        alpha = 1,             # Set alpha to 1 for LASSO
        lambda = avg_results$lambda) # Define best average lambda
      
      final_model <- train(chr ~ .,
                           data = train,
                           metric = "C",
                           method = "glmnet",
                           family = "binomial",
                           tuneGrid = final_tune,
                           trControl = control_final)
      
      } else {
        final_model <- train(chr ~ .,
                             data = train,
                             metric = "C",
                             method = "glm",
                             family = "binomial",
                             trControl = control_final)
      }
      cat(sum(is.na(test))," missing data points \n")
      
      # Predict the linear predictors (PI)
      test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
      test$chr_pred <- predict(final_model, newdata = test, type = "raw")
      
      cm <- confusionMatrix(data = test$chr_pred, reference = test$chr)
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      test$chr_num <- as.numeric(test$chr)-1
      
      # Fit a logistic regression model on the test set using penalized coefficients
      model_test <- glm(chr ~ PI, data = test, family="binomial")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$chr=="X1"),
        n_test = nrow(test),
        events_test = sum(test$chr=="X1"),
        balanced_accuracy = cm$byClass[11],
        sensitivity = cm$byClass[1],
        specificity = cm$byClass[2],
        ppv = cm$byClass[3],
        npv = cm$byClass[4],
        precision = cm$byClass[5],
        recall = cm$byClass[6],
        f1 = cm$byClass[7]
        
      ))
      
      # Calculate calibration slope on the hold-out data (test set)
      cal <- pred_val_probs(binary_outcome=test$chr_num, Prob=test$PI, cal_plot=FALSE)
      calibration_intercept <- cal$CalInt
      calibration_slope <- cal$CalSlope
      brier <- cal$BrierScore
      
      # Calculate ICI
      loess.calibrate <- loess(as.numeric(test$chr)-1 ~ test$PI)
     
      # Estimate loess-based smoothed calibration curve
      P.calibrate <- predict(loess.calibrate, newdata = test)
    
      # This is the point on the loess calibration curve corresponding to a given predicted probability.
      ICI <- mean(abs(P.calibrate - test$PI))
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        intercept = calibration_intercept,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
      
      # Decision curve analysis
      library(dcurves)
      dca <- data.frame(obs=as.numeric(test$chr)-1,
                        pred=as.numeric(test$PI),
                        chr_pred=test$chr_num)
      
      dca_assessment <- dca(obs ~ pred,
                            data = dca,
                            prevalence = 0.192,
                            thresholds = seq(0, 0.5, 0.01)
      ) %>% as_tibble()
      
      dca_clinical <- rbind(dca_clinical, data.frame(label=dca_assessment$label,
                                                     n=dca_assessment$n,
                                                     outer_repeat = rep(outer_rep,nrow(dca_assessment)),
                                                     outer_fold = rep(i,nrow(dca_assessment)),
                                                     threshold=dca_assessment$threshold,
                                                     pos_rate=dca_assessment$pos_rate,
                                                     net_benefit=dca_assessment$net_benefit))
      
      dca_assessment <- dca(obs ~ pred,
                            data = dca,
                            prevalence = 0.017,
                            thresholds = seq(0, 0.5, 0.001)
      ) %>% as_tibble()
      
      dca_general <- rbind(dca_general, data.frame(label=dca_assessment$label,
                                                   n=dca_assessment$n,
                                                   outer_repeat = rep(outer_rep,nrow(dca_assessment)),
                                                   outer_fold = rep(i,nrow(dca_assessment)),
                                                   threshold=dca_assessment$threshold,
                                                   pos_rate=dca_assessment$pos_rate,
                                                   net_benefit=dca_assessment$net_benefit))
      preds <- rbind(preds, dca)
    }
  }
  
  list(
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes,
    dca_clinical = dca_clinical,
    dca_general = dca_general,
    preds = preds
  )
}

##### Cox with repeated nested cross validation ##### 
Cox_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()

  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("OuterRepeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$Transition, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      min_lambda <- data.frame()
      
      if (ncol(combined_df)>6) {
        # If there is more than one predictor, use LASSO
        combined_df[,c(1:16)] <- lapply(combined_df[,c(1:16)], factor)
        options(na.action='na.pass')
        combined_df <- combined_df[, sapply(combined_df, function(x) length(unique(na.omit(x))) > 1)] # Remove columns with zero variance
        
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat[,c(1:21)] <- lapply(PPS_scored_mat[,c(1:21)], factor)
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        # Identify and remove near-zero variance predictors
        near_zero_vars <- nearZeroVar(PPS_scored_mat)
        PPS_scored_mat <- PPS_scored_mat[, -near_zero_vars]
        
        } else if (ncol(combined_df)>3){
        # If there is more than one predictor, use LASSO
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))

        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        } 
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(test)>3) {
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
      if (anyNA(test)) {
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
        imputed_data <- mice(test, m = 1, method = 'pmm', maxit = 50, seed = 123)
        test <- complete(imputed_data)
        }
      } else {
      imputed_data <- mice(test)
      test <- complete(imputed_data)
      }
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        if (ncol(train_base)>3) {
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }
        
        if (ncol(train)>3) { # If there is more than one predictor, tune LASSO
          
          # Convert all variables to numeric for matrix
          train_num <- train %>%
            mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
            mutate(across(where(is.factor), as.numeric))         # Convert factor to numeric
          
          predictors <- as.matrix(train_num[, !(colnames(train_num) %in% c("day_exit", "Transition"))])
          SurvObj <- Surv(time = train$day_exit, event = train$Transition)
          
          # Train the model with LASSO and C-index as the metric
          inner_model <- cv.glmnet(predictors, SurvObj, family = "cox", alpha = 1)
          min_lambda <- rbind(min_lambda, data.frame(lambda=inner_model$lambda.min))
          
          print(inner_model)
          
        } 
      }
        
      cat("Inner model fitted \n")
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Train the final model using the mean minimum lambda
      best_lambda <- mean(min_lambda$lambda)      
      if (ncol(train)>3) {
        
      final_model <- glmnet(predictors, SurvObj, family = "cox", alpha = 1, lambda = best_lambda)
      cat("Final model fitted \n")
      
      # Convert all variables to numeric for matrix
      test_num <- test %>%
        mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
        mutate(across(where(is.factor), as.numeric))         # Convert factor to numeric
      cat("Converted \n")
      
      test_m <- as.matrix(test_num)
      test_m <- test_m[, colnames(predictors)] 
      
      # Predict the linear predictors (PI)
      test$PI <- predict(final_model, newx = test_m, type = "response")
      
      cat("PI predicted \n")
      
      } else {

      final_model <- coxph(Surv(time = day_exit, event = Transition) ~ PRS_resid, data=train)
      cat("Final model fitted \n")
      
      test_m <- test %>% subset(select=c(-day_exit,-Transition))
      test$PI <- predict(final_model, newdata = test_m, type = "risk")
      cat("PI predicted \n")
      }
      
      # Create survival object 
      SurvObj_test <- with(test, Surv(time = day_exit, event = Transition)) 
     
      # Fit a logistic regression model on the test set using penalized coefficients
      model_test <- coxph(SurvObj_test ~ PI, data = test)
      
      cat("Discrimination calculated \n")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train_base$Transition),
        n_test = nrow(test),
        events_test = sum(test$Transition)
      ))
      
      cat("Discrimination saved \n")
      
      # Fit a logistic regression model
      base_surv <- survfit(Surv(day_exit, Transition) ~ 1, data=test)
      test$pred_risk <- 1 - (exp(-base_surv$cumhaz[which.max(base_surv$time)])^exp(test$PI))
      cal <- pred_val_probs(binary_outcome=test$Transition, Prob=test$pred_risk, cal_plot=FALSE)
      calibration_in_large <- cal$CalInt
      calibration_slope <- model_test$coefficients[1]
      brier <- cal$BrierScore
      
      cat("Calibration saved \n")
      
      cat("Outer rep = ", outer_rep, " outer fold = ", i,"\n")

      cat("Brier saved \n")
      
      # Calculate ICI
      if (sd(test$PI)!=0){
      loess.calibrate <- loess(test$Transition ~ test$pred_risk)
      if (sum(!is.na(loess.calibrate$fitted))==length(loess.calibrate$fitted)){
      # Estimate loess-based smoothed calibration curve
      P.calibrate <- predict(loess.calibrate, newdata = test)
      
      # This is the point on the loess calibration curve corresponding to a given predicted probability.
      ICI <- mean(abs(P.calibrate - test$pred_risk))
      } else {
        ICI <- NA
      }
      } else {
        ICI <- NA
        }
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        calibration_in_large = calibration_in_large,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
      cat("Results saved \n")
      
      
      
    }
  }
  
  list(
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

##### RF with repeated nested cross validation #####
## Compute F1-score

F1_score <- function(mat, algoName){
  
  # Remark: left column = prediction // top = real values
  recall <- matrix(1:nrow(mat), ncol = nrow(mat))
  precision <- matrix(1:nrow(mat), ncol = nrow(mat))
  F1_score <- matrix(1:nrow(mat), ncol = nrow(mat))
  
  
  for(i in 1:nrow(mat)){
    recall[i] <- mat[i,i]/rowSums(mat)[i]
    precision[i] <- mat[i,i]/colSums(mat)[i]
  }
  
  for(i in 1:ncol(recall)){
    F1_score[i] <- 2 * ( precision[i] * recall[i] ) / ( precision[i] + recall[i])
  }
  
  # We display the matrix labels
  colnames(F1_score) <- colnames(mat)
  rownames(F1_score) <- algoName
  
  # Display the F1_score for each class
  F1_score
  
  # Display the average F1_score
  mean(F1_score[1,], na.rm=TRUE)
}

# Create a function to perform nested cross-validation with repeats
RF_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  all_inner_results <- list()
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  dca_clinical <- data.frame()
  dca_general <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$chr, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>2) { # If there is more than one predictor, use LASSO
        combined_df[,c(1:15)] <- lapply(combined_df[,c(1:15)], factor)
        nzv <- nearZeroVar(combined_df, saveMetrics = TRUE) # Remove predictors with near zero variance
        combined_df_nzv <- combined_df[,nzv[,"nzv"] == FALSE] 
        options(na.action='na.pass')
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df_nzv))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      test$chr <- factor(test$chr)
      test$chr <- factor(make.names(levels(test$chr))[test$chr])
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        # Tune both alpha and lambda
        for (mtry_value in tuneGrid) {
          library(randomForest)
          library(MLmetrics)
          library(caret)
          
          # Define a custom summary function to calculate F1 score
          customSummary <- function(data, lev = NULL, model = NULL) {
            f1 <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
            out <- c(F1 = f1)
            out
          }
          
          if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
            imputed_data <- missForest(train_base)
            train <- imputed_data$ximp
          } else {
            imputed_data <- mice(train_base)
            train <- complete(imputed_data)
          }
          
          train$chr <- factor(train$chr)
          train$chr <- factor(make.names(levels(train$chr))[train$chr])
          
          tree <- c(50, 100, 250, 500)
          n.tree <- sample(tree,1)
          nodeSize <- seq(1,(nrow(train)/10), by=1)
          node.size <- sample(nodeSize,1)
          
          if (ncol(train)>3){
          tune_grid_temp <- data.frame(mtry=c(NA,NA,NA,NA,NA))
          tune_grid_temp$mtry <- sample(tuneGrid$mtry,5)
          
          # Set up the trainControl with the custom summary function
          control <- trainControl(method = 'cv', 
                                  number = 5, 
                                  classProbs = TRUE,
                                  summaryFunction = customSummary,
                                  search = 'random')
          } else {
          tune_grid_temp <- tuneGrid
          
          # Set up the trainControl with the custom summary function
          control <- trainControl(method = 'cv', 
                                  number = 5, 
                                  classProbs = TRUE,
                                  search = 'random')
          }
          
          if(ncol(train)>2){
            
          # Train the model using the custom F1 metric
            inner_model <- train(chr ~ .,
                                 data = train,
                                 method = "rf",
                                 metric = "F1",
                                 tuneGrid = tune_grid_temp,
                                 tuneLength=10,
                                 ntree = n.tree,
                                 nodesize=node.size,
                                 trControl = control)
          } else {
            inner_model <- train(chr ~ .,
                                 data = train,
                                 method = "rf",
                                 metric = "Accuracy",
                                 tuneGrid = tune_grid_temp,
                                 tuneLength=10,
                                 ntree = n.tree,
                                 nodesize=node.size,
                                 trControl = control)
          }
          print(inner_model)
          cat("Inner model fitted \n")
          
          if (ncol(train)>2){
          if (all(!is.na(inner_model$results$F1))){
            # Store the F1 for each inner repeat
            inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                data.frame(min.node.size = node.size,
                                                           tree = n.tree,
                                                           mtry = inner_model$results$mtry, 
                                                           F1 = inner_model$results$F1, 
                                                           repeat_number = inner_rep))
          }
          } else {
            if (all(!is.na(inner_model$results$Accuracy))){
            # Store the F1 for each inner repeat
            inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                data.frame(min.node.size = node.size,
                                                           tree = n.tree,
                                                           mtry = inner_model$results$mtry, 
                                                           F1 = inner_model$results$Accuracy, 
                                                           repeat_number = inner_rep))
            }
          }
        }
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      print(all_inner_results)
      
      # Calculate the average F1 for each hyperparameter combination
      avg_results <- all_inner_results %>%
        group_by(mtry) %>%
        summarise(avg_F1 = mean(F1, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(tree, na.rm=TRUE),
                  nodesize=mean(min.node.size, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparameter combination based on the highest average F1
      best_inner_result <- avg_results[which.max(avg_results$avg_F1), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_F1 <- best_inner_result$avg_F1
      best_F1_SE <- best_inner_result$SE_F1
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      cat("Inner results saved \n")
      
      # Train the final model using the best hyperparameters on the training data
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      if (ncol(train)>3) {


        final_model <- train(chr ~ .,
                             data = train,
                             method = "rf",
                             metric = "F1",
                             ntree = best_ntree,
                             nodesize=best_min.node.size,
                             trControl=control_final,
                             tuneGrid = repGrid)
      
      } else {
      final_model <- train(chr ~ .,
                           data = train,
                           method = "rf",
                           metric = "Accuracy",
                           ntree = best_ntree,
                           nodesize=best_min.node.size,
                           trControl=control_final,
                           tuneGrid = repGrid)
      }
      
      cat("Final model fitted \n")
      cat("Predicting probabilities and classes...\n")
      
      cat("Outer rep:", outer_rep, " ; Outer fold:", i)
        # Predict the linear predictors (PI) from the model
        if ((i!=2 & outer_rep!=2) & (i!=4 & outer_rep!=3) & (i!=1 & outer_rep!=5) & (i!=5 & outer_rep!=6)){
        test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
        test$chr_pred <- predict(final_model, newdata = test, type = "raw")
        cat("PI generated \n")
        
      cm <- confusionMatrix(data = test$chr_pred, reference = test$chr)
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a logistic regression model on the test set using penalized coefficients
      model_test <- glm(chr ~ PI, data = test, family="binomial")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$chr=="X1"),
        n_test = nrow(test),
        events_test = sum(test$chr=="X1"),
        balanced_accuracy = cm$byClass[11],
        sensitivity = cm$byClass[1],
        specificity = cm$byClass[2],
        ppv = cm$byClass[3],
        npv = cm$byClass[4],
        precision = cm$byClass[5],
        recall = cm$byClass[6],
        f1 = cm$byClass[7]
        
      ))
      
      cat("Discrimination results saved \n")
      
      # Calculate calibration slope on the hold-out data (test set)
      calibration_intercept <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[12])
      calibration_slope <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[13])
      brier <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[11])
      
      # Calculate ICI
      loess.calibrate <- loess(as.numeric(test$chr)-1 ~ test$PI)
      
      # Estimate loess-based smoothed calibration curve
      P.calibrate <- predict(loess.calibrate, newdata = test)
      
      # This is the point on the loess calibration curve corresponding to a given predicted probability.
      ICI <- mean(abs(P.calibrate - test$PI))
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        intercept = calibration_intercept,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
      
      cat("Calibration results saved \n")
      
      # Decision curve analysis
      library(dcurves)
      dca <- data.frame(obs=as.numeric(test$chr)-1,
                        pred=as.numeric(test$PI))
      
      dca_assessment <- dca(obs ~ pred,
                            data = dca,
                            prevalence = 0.192,
                            thresholds = seq(0, 0.5, 0.01)
      ) %>% as_tibble()
      
      dca_clinical <- rbind(dca_clinical, data.frame(label=dca_assessment$label,
                                                     n=dca_assessment$n,
                                                     outer_repeat = rep(outer_rep,nrow(dca_assessment)),
                                                     outer_fold = rep(i,nrow(dca_assessment)),
                                                     threshold=dca_assessment$threshold,
                                                     pos_rate=dca_assessment$pos_rate,
                                                     net_benefit=dca_assessment$net_benefit))
      
      dca_assessment <- dca(obs ~ pred,
                            data = dca,
                            prevalence = 0.017,
                            thresholds = seq(0, 0.5, 0.001)
      ) %>% as_tibble()
      
      dca_general <- rbind(dca_general, data.frame(label=dca_assessment$label,
                                                   n=dca_assessment$n,
                                                   outer_repeat = rep(outer_rep,nrow(dca_assessment)),
                                                   outer_fold = rep(i,nrow(dca_assessment)),
                                                   threshold=dca_assessment$threshold,
                                                   pos_rate=dca_assessment$pos_rate,
                                                   net_benefit=dca_assessment$net_benefit))
      cat("DCA results saved \n")
        }
    }
  }
  
  list(
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes,
    dca_clinical = dca_clinical,
    dca_general = dca_general
  )
}

##### RSF with repeated nested cross validation #####
RSF_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  inner_results <- data.frame()
  all_inner_results <- data.frame()
  c_stat_nested  <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$Transition, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>6) {
        combined_df[,c(1:15)] <- lapply(combined_df[,c(1:15)], factor)
        options(na.action='na.pass')
        combined_df <- combined_df[, sapply(combined_df, function(x) length(unique(na.omit(x))) > 1)] # Remove columns with zero variance
        
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat[,c(1:21)] <- lapply(PPS_scored_mat[,c(1:21)], factor)
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        # Identify and remove near-zero variance predictors
        near_zero_vars <- nearZeroVar(PPS_scored_mat)
        PPS_scored_mat <- PPS_scored_mat[, -near_zero_vars]
        
      } else if (ncol(combined_df)>3){
        # If there is more than one predictor, use LASSO
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        } 
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(test)>3) {
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
        test[,c((ncol(test)-1):ncol(test))] <- round(test[,c((ncol(test)-1):ncol(test))])
        if (anyNA(test)) {
          imputed_data <- missForest(test)
          test <- imputed_data$ximp
          imputed_data <- mice(test, m = 1, method = 'pmm', maxit = 50, seed = 123)
          test <- complete(imputed_data)
        }
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        if (ncol(train_base)>3) {
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
          train[,c((ncol(train)-1):ncol(train))] <- round(train[,c((ncol(train)-1):ncol(train))])
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }
        
        library(randomForestSRC)
        library(survcomp)
        
        tree <- c(50, 100, 250, 500)
        nodeSize <- seq(1,(nrow(train)/10), by=1)

        # Store fold-specific C-index scores
        fold_cindex <- numeric(50)

        for (tune in 1:10) {
        for (fold_idx in 1:5) {
          
          set.seed <- fold_idx
          # Get the training and validation sets for this fold
          folds_cv <- createFolds(train$Transition, k = 5, returnTrain = TRUE)
          
          n.tree <- sample(tree,1)
          node.size <- sample(nodeSize,1)
          mtry <- sample(tuneGrid$mtry,1)
          
          train_indices_cv <- folds_cv[[fold_idx]]
          
          fold_train <- train[train_indices_cv, ]
          fold_valid <- train[-train_indices_cv, ]
          
          predictor_vars <- setdiff(names(fold_train), c("day_exit", "Transition"))
          
          # Dynamically create the formula for the RSF model
          formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))
          
          # Train the RSF model on the training set
          rsf_model <- rfsrc(formula, data = fold_train, ntree = n.tree, mtry = mtry, nodesize = node.size)
          
          # Calculate the C-index on the validation set
          rsf_pred <- predict(rsf_model, fold_valid)
          fold_cindex[tune] <- concordance.index(rsf_pred$survival[,ncol(rsf_pred$survival)], surv.time = fold_valid$day_exit, surv.event = fold_valid$Transition)$c.index
          
        }
        }
        
        # Store the mean C-index for this combination of hyperparameters
        mean_cindex <- mean(fold_cindex)
        inner_results <- rbind(inner_results, data.frame(ntree = n.tree, mtry = mtry, nodesize=node.size, mean_cindex = mean_cindex))
        
        
        # View the results of cross-validation
        print(inner_results)
        
        # Find the best hyperparameter combination based on the highest C-index
        best_params <- inner_results[which.max(inner_results$mean_cindex), ]
        
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      
      # Calculate the average C-index for each hyperparameter combination
      avg_results <- inner_results %>%
        group_by(mtry) %>%
        summarise(avg_C = mean(mean_cindex, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(ntree, na.rm=TRUE),
                  nodesize=mean(nodesize, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparametercombination based on the highest average C-index
      best_inner_result <- avg_results[which.max(avg_results$avg_C), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_C <- best_inner_result$avg_C
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      final_model <- rfsrc(formula, data = train, ntree = best_ntree, mtry = best_mtry, nodesize = best_min.node.size)
      
      # Predict the linear predictors (PI) from the Elastic Net model
      rsf_pred_test <- predict(final_model, test)
      test$PI <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      test$Tx_pred <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a Cox model on the test set using penalized coefficients
      model_test <- coxph(Surv(day_exit, Transition) ~ PI, data = test)
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$Transition==1),
        n_test = nrow(test),
        events_test = sum(test$Transition==1)
        
      ))
      
      # Fit a logistic regression model
      cal <- pred_val_probs(binary_outcome=test$Transition, Prob=test$PI, cal_plot=FALSE)
      calibration_in_large <- cal$CalInt
      calibration_slope <- cal$CalSlope
      brier <- cal$BrierScore
      
      cat("Calibration saved \n")
      
      cat("Outer rep = ", outer_rep, " outer fold = ", i,"\n")
      
      cat("Brier saved \n")
      
      # Calculate ICI
      if (sd(test$PI)!=0){
        loess.calibrate <- loess(test$Transition ~ test$Tx_pred)
        if (sum(!is.na(loess.calibrate$fitted))==length(loess.calibrate$fitted)){
          # Estimate loess-based smoothed calibration curve
          P.calibrate <- predict(loess.calibrate, newdata = test)
          
          # This is the point on the loess calibration curve corresponding to a given predicted probability.
          ICI <- mean(abs(P.calibrate - test$Tx_pred))
        } else {
          ICI <- NA
        }
      } else {
        ICI <- NA
      }
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        calibration_in_large = calibration_in_large,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
      
    }
  }
  
  list(
    best_inner_result_list = best_inner_result_list,
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

##### Harmonised LR with repeated nested cross validation #####
harmonised_lr_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$chr, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>3) { # If there is more than one predictor, use LASSO
        combined_df[,c(1:16)] <- lapply(combined_df[,c(1:16)], factor)
        options(na.action='na.pass')
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      set.seed(123)
      if (ncol(combined_df)>3) { # If there is more than one predictor, impute with random forest, else use mice
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      ### Compute Global Means ###
      global_mean <- colMeans(test[, 1:(ncol(test)-2), drop = FALSE], na.rm = TRUE)
      
      ### Mean Offset Correction ###
      batch_test <- as.factor(test$site)
      test_corrected <- test # Start with the original data
      
      for (b in levels(batch_test)) {
        batch_indices <- which(batch_test == b) # Indices for samples in batch `b`
        
        if (length(batch_indices) == 0) {
          warning(paste("Batch", b, "is empty. Skipping."))
          next
        }
        
        # Extract the subset of test for the current batch
        batch_data <- test[batch_indices, 1:(ncol(test)-2), drop = FALSE] # Exclude site and outcome
        
        # Compute means for each predictor in this batch
        batch_mean <- colMeans(batch_data, na.rm = TRUE)
        
        if (length(batch_mean) == 0) {
          warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
          next
        }
        
        # Compute the offset: batch mean - global mean
        offset <- batch_mean - global_mean
        
        # Subtract the offset to align the batch with the global mean
        test_corrected[batch_indices, 1:(ncol(test)-2)] <- sweep(
          batch_data,
          1,
          offset,
          "-"
        )
      }
      test <- test_corrected # %>% subset(select=c(-site))
      test$chr <- factor(test$chr)
      test$chr <- factor(make.names(levels(test$chr))[test$chr])
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      # If any remaining missing values with the column mean
      test[is.na(test)] <- apply(test, 2, function(col) mean(col, na.rm = TRUE))
      
      make.names(colnames(test), unique=TRUE)
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        set.seed(123)
        if (ncol(combined_df)>3) { # If there is more than one predictor, impute with random forest, else use mice
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }
        
        ### Mean offset correction ###
        batch_train <- as.factor(train$site)
        train_corrected <- train  # Start with the original data
        
        ### Compute Global Means ###
        global_mean <- colMeans(test[, 1:(ncol(test)-2), drop = FALSE], na.rm = TRUE)
        # Iterate over each batch
        for (b in levels(batch_train)) {
          batch_indices <- which(batch_train == b)  # Indices for samples in batch `b`
          
          if (length(batch_indices) == 0) {
            warning(paste("Batch", b, "is empty. Skipping."))
            next
          }
          
          # Extract the subset of test for the current batch
          batch_data <- train[batch_indices, 1:(ncol(train)-2), drop = FALSE]
          
          # Compute row-wise means for this batch
          batch_mean <- colMeans(batch_data, na.rm = TRUE)
          
          # Compute the offset: batch mean - global mean
          offset <- batch_mean - global_mean
          
          if (length(batch_mean) == 0) {
            warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
            next
          }
          
          # Subtract batch means
          train_corrected[batch_indices, 1:(ncol(train)-2)] <- sweep(
            batch_data, 
            1, 
            offset, 
            "-"
          )
        }
        
      train <- train_corrected # %>% subset(select=c(-site))
        
        train$chr <- factor(train$chr)
        train$chr <- factor(make.names(levels(train$chr))[train$chr])
        
        cindex_summary <- function(data, lev = NULL, model = NULL) {
          requireNamespace("survival")
          # Convert factors to numeric (0 and 1) for concordance calculation
          obs_numeric <- ifelse(data$obs == lev[1], 0, 1)
          pred_numeric <- as.numeric(data$pred)
          
          # Calculate Harrell's C-index using predicted probabilities
          cindex <- Hmisc::rcorr.cens(pred_numeric, obs_numeric)[2]
          out <- c(C = cindex)
          names(out) <- "C"
          out
        }
        
        # Set up trainControl with custom summary function
        control <- trainControl(
          method = "cv",
          number = 5,
          summaryFunction = cindex_summary,
          classProbs = TRUE, # Needed if using probabilities
          sampling = "smote"
        )
        
        if (ncol(train)>3) { # If there is more than one predictor, use LASSO
          # Train the model
          tune_grid <- expand.grid(
            alpha = 1,             # Set alpha to 1 for LASSO
            lambda = seq(0.001, 0.1, length = 100) # Define a range for lambda
          )
          
          # Train the model with LASSO and C-index as the metric
          inner_model <- train(
            chr ~ .,
            data = train,
            method = "glmnet",
            trControl = control,
            metric = "C",
            tuneGrid = tune_grid,
            family = "binomial"
          )
        } else {
          # Train the model with logistic regression and C-index as the metric
          inner_model <- train(
            chr ~ .,
            data = train,
            method = "glm",
            trControl = control,
            metric = "C",
            family = "binomial"
          )
        }
        print(inner_model)
        
        if (all(!is.na(inner_model$results$F1))){
          # Store the F1 for each inner repeat
          inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                              data.frame(C = inner_model$results$C, 
                                                         repeat_number = inner_rep))
        }
      }
      
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      
      # Calculate the average F1 for each hyperparameter combination
      avg_results <- all_inner_results %>%
        summarise(avg_C = mean(C, na.rm=TRUE),
                  .groups = "keep") 
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none',
                                    sampling = "smote")
      
      final_model <- train(chr ~ .,
                           data = train,
                           metric = "C",
                           method = "glm",
                           family = "binomial",
                           trControl = control_final)
      
      # Predict the linear predictors (PI) from the Elastic Net model
      test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
      test$chr_pred <- predict(final_model, newdata = test, type = "raw")
      
      cm <- confusionMatrix(data = test$chr_pred, reference = test$chr)
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a logistic regression model on the test set using penalized coefficients
      model_test <- glm(chr ~ PI, data = test, family="binomial")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$chr=="X1"),
        n_test = nrow(test),
        events_test = sum(test$chr=="X1"),
        balanced_accuracy = cm$byClass[11],
        sensitivity = cm$byClass[1],
        specificity = cm$byClass[2],
        ppv = cm$byClass[3],
        npv = cm$byClass[4],
        precision = cm$byClass[5],
        recall = cm$byClass[6],
        f1 = cm$byClass[7]
        
      ))
      
      # Fit a logistic regression model
      test$chr_num <- as.numeric(test$chr)-1
      cal <- pred_val_probs(binary_outcome=test$chr_num, Prob=test$PI, cal_plot=FALSE)
      calibration_intercept <- cal$CalInt
      calibration_slope <- cal$CalSlope
      brier <- cal$BrierScore
      
      # Calculate ICI
      loess.calibrate <- loess(as.numeric(test$chr)-1 ~ test$PI)
      
      # Estimate loess-based smoothed calibration curve
      P.calibrate <- predict(loess.calibrate, newdata = test)
      
      # This is the point on the loess calibration curve corresponding to a given predicted probability.
      ICI <- mean(abs(P.calibrate - test$PI))
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        intercept = calibration_intercept,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
    }
  }
  
  list(
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

##### Harmonised Cox #####
harmonised_Cox_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("OuterRepeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$Transition, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      min_lambda <- data.frame()
      
      if (ncol(combined_df)>7) {
        # If there is more than one predictor, use LASSO
        combined_df[,c(1:16)] <- lapply(combined_df[,c(1:16)], factor)
        options(na.action='na.pass')
        combined_df <- combined_df[, sapply(combined_df, function(x) length(unique(na.omit(x))) > 1)] # Remove columns with zero variance
        
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        # Identify and remove near-zero variance predictors
        near_zero_vars <- nearZeroVar(PPS_scored_mat)
        PPS_scored_mat <- PPS_scored_mat[, -near_zero_vars]
        
      } else if (ncol(combined_df)>4){
        # If there is more than one predictor, use LASSO
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        } 
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(test)>4) {
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
        if (anyNA(test)) {
          imputed_data <- missForest(test)
          test <- imputed_data$ximp
          imputed_data <- mice(test)
          test <- complete(imputed_data)
        }
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      
      ### Mean Offset Correction ###
      batch_test <- as.factor(test$site)
      
      ### Compute Global Means ###
      test_num <- test %>%
        mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
        mutate(across(where(is.factor), as.numeric))   # Convert factor to numeric
      test_num[,c(1:(ncol(test_num)-3))] <- test_num[,c(1:(ncol(test_num)-3))]
      
      global_mean <- colMeans(test_num[, 1:(ncol(test)-3), drop = FALSE], na.rm = TRUE)
      
      test_corrected <- test_num # Start with the original data
      
      for (b in levels(batch_test)) {
        batch_indices <- which(batch_test == b) # Indices for samples in batch `b`
        
        if (length(batch_indices) == 0) {
          warning(paste("Batch", b, "is empty. Skipping."))
          next
        }
        
        # Extract the subset of test for the current batch
        batch_data <- test_num[batch_indices, 1:(ncol(test)-3), drop = FALSE] # Exclude site and outcome
        
        # Compute means for each predictor in this batch
        batch_mean <- colMeans(batch_data, na.rm = TRUE)
        
        if (length(batch_mean) == 0) {
          warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
          next
        }
        
        # Compute the offset: batch mean - global mean
        offset <- batch_mean - global_mean
        
        # Subtract the offset to align the batch with the global mean
        test_corrected[batch_indices, 1:(ncol(test)-3)] <- sweep(
          batch_data,
          1,
          offset,
          "-"
        )
      }
      test <- test_corrected # %>% subset(select=c(-site))
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        if (ncol(train_base)>4) {
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
          if (anyNA(train)) {
            imputed_data <- missForest(train)
            train <- imputed_data$ximp
            imputed_data <- mice(train)
            train <- complete(imputed_data)
          }
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }
        
        if (ncol(train)>4) { # If there is more than one predictor, tune LASSO
          ### Mean offset correction ###
          batch_train <- as.factor(train$site)
          
          # Convert all variables to numeric for matrix
          train_num <- train %>%
            mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
            mutate(across(where(is.factor), as.numeric))  # Convert factor to numeric
            train_num[,c(1:(ncol(train_num)-3))] <- train_num[,c(1:(ncol(train_num)-3))]

            train_corrected <- train_num
            
          ### Compute Global Means ###
          global_mean <- colMeans(train_num[, 1:(ncol(train_num)-3), drop = FALSE], na.rm = TRUE)
          # Iterate over each batch
          for (b in levels(batch_train)) {
            batch_indices <- which(batch_train == b)  # Indices for samples in batch `b`
            
            if (length(batch_indices) == 0) {
              warning(paste("Batch", b, "is empty. Skipping."))
              next
            }
            
            # Extract the subset of test for the current batch
            batch_data <- train_num[batch_indices, 1:(ncol(train_num)-3), drop = FALSE]
            
            # Compute row-wise means for this batch
            batch_mean <- colMeans(batch_data, na.rm = TRUE)
            
            # Compute the offset: batch mean - global mean
            offset <- batch_mean - global_mean
            
            if (length(batch_mean) == 0) {
              warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
              next
            }
            
            # Subtract batch means
            train_corrected[batch_indices, 1:(ncol(train_corrected)-3)] <- sweep(
              batch_data, 
              1, 
              offset, 
              "-"
            )
          }
          
          train <- train_corrected# %>% subset(select=c(-site))
          
          imputed_data <- missForest(train)
          train <- imputed_data$ximp
          
          predictors <- as.matrix(train[, !(colnames(train) %in% c("day_exit", "Transition"))])
          SurvObj <- Surv(time = train$day_exit, event = train$Transition)
          
          # Train the model with LASSO and C-index as the metric
          inner_model <- cv.glmnet(predictors, SurvObj, family = "cox", alpha = 1)
          min_lambda <- rbind(min_lambda, data.frame(lambda=inner_model$lambda.min))
          
          print(inner_model)
          
        } 
      }
      
      cat("Inner model fitted \n")
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Train the final model using the mean minimum lambda
      best_lambda <- mean(min_lambda$lambda)      
      if (ncol(train)>4) {
        
        final_model <- glmnet(predictors, SurvObj, family = "cox", alpha = 1, lambda = best_lambda)
        cat("Final model fitted \n")
        
        # Convert all variables to numeric for matrix
        test_num <- test %>%
          mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
          mutate(across(where(is.factor), as.numeric))         # Convert factor to numeric
        cat("Converted \n")
        
        test_m <- as.matrix(test_num)
        test_m <- test_m[, colnames(predictors)] 
        
        # Predict the linear predictors (PI)
        test$PI <- predict(final_model, newx = test_m, type = "response")
        
        cat("PI predicted \n")
        
      } else {
        
        final_model <- coxph(Surv(time = day_exit, event = Transition) ~ PRS_resid, data=train)
        cat("Final model fitted \n")
        
        test_m <- test %>% subset(select=c(-day_exit,-Transition))
        test$PI <- predict(final_model, newdata = test_m, type = "risk")
        cat("PI predicted \n")
      }
      
      # Create survival object 
      SurvObj_test <- with(test, Surv(time = day_exit, event = Transition)) 
      
      # Fit a logistic regression model on the test set using penalized coefficients
      model_test <- coxph(SurvObj_test ~ PI, data = test)
      
      cat("Discrimination calculated \n")
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train_base$Transition),
        n_test = nrow(test),
        events_test = sum(test$Transition)
      ))
      
      cat("Discrimination saved \n")
      
      # Fit a logistic regression model
      logistic_model <- glm(Transition ~ PI, data = test, family = binomial)
      
      calibration_in_large <- coef(logistic_model)[1]
      calibration_slope <- model_test$coefficients[1]
      
      cat("Calibration saved \n")
      
      cat("Outer rep = ", outer_rep, " outer fold = ", i,"\n")
      
      # Get baseline survival function
      base_surv <- survfit(model_test)
      
      # Predict linear predictor
      test$LP <- predict(model_test, newx = test_m, type = "lp")
      
      # Calculate predicted survival probabilities
      test$surv_prob <-  exp(-base_surv$cumhaz[which.max(base_surv$time)] * exp(test$LP))
      observed <- test$Transition
      preds <- test$surv_prob
      brier <-  mean((preds - observed)^2)
      
      cat("Brier saved \n")
      
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        calibration_in_large = calibration_in_large,
        slope = calibration_slope,
        brier = brier
      ))
      cat("Results saved \n")
      
      
      
    }
  }
  
  list(
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}

##### Harmonised RF with repeated nested cross validation #####
harmonised_RF_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, seed, tuneGrid) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  all_inner_results <- list()
  c_stat_nested <- data.frame()
  calibration_slopes <- data.frame()
  dca_clinical <- data.frame()
  dca_general <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$chr, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>2) { # If there is more than one predictor, use LASSO
        combined_df[,c(1:15)] <- lapply(combined_df[,c(1:15)], factor)
        nzv <- nearZeroVar(combined_df, saveMetrics = TRUE) # Remove predictors with near zero variance
        combined_df_nzv <- combined_df[,nzv[,"nzv"] == FALSE] 
        options(na.action='na.pass')
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df_nzv))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      ### Compute Global Means ###
      global_mean <- colMeans(test[, 1:(ncol(test)-2), drop = FALSE], na.rm = TRUE)
      
      ### Mean Offset Correction ###
      batch_test <- as.factor(test$site)
      test_corrected <- test # Start with the original data
      
      for (b in levels(batch_test)) {
        batch_indices <- which(batch_test == b) # Indices for samples in batch `b`
        
        if (length(batch_indices) == 0) {
          warning(paste("Batch", b, "is empty. Skipping."))
          next
        }
        
        # Extract the subset of test for the current batch
        batch_data <- test[batch_indices, 1:(ncol(test)-2), drop = FALSE] # Exclude site and outcome
        
        # Compute means for each predictor in this batch
        batch_mean <- colMeans(batch_data, na.rm = TRUE)
        
        if (length(batch_mean) == 0) {
          warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
          next
        }
        
        # Compute the offset: batch mean - global mean
        offset <- batch_mean - global_mean
        
        # Subtract the offset to align the batch with the global mean
        test_corrected[batch_indices, 1:(ncol(test)-2)] <- sweep(
          batch_data,
          1,
          offset,
          "-"
        )
      }
      test <- test_corrected # %>% subset(select=c(-site))
      test$chr <- factor(test$chr)
      test$chr <- factor(make.names(levels(test$chr))[test$chr])
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        # Tune both alpha and lambda
        for (mtry_value in tuneGrid) {
          library(randomForest)
          library(MLmetrics)
          library(caret)
          
          # Define a custom summary function to calculate F1 score
          customSummary <- function(data, lev = NULL, model = NULL) {
            f1 <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
            out <- c(F1 = f1)
            out
          }
          
          if (ncol(combined_df)>2) { # If there is more than one predictor, impute with random forest, else use mice
            imputed_data <- missForest(train_base)
            train <- imputed_data$ximp
          } else {
            imputed_data <- mice(train_base)
            train <- complete(imputed_data)
          }
          
          ### Mean offset correction ###
          batch_train <- as.factor(train$site)
          train_corrected <- train  # Start with the original data
          
          ### Compute Global Means ###
          global_mean <- colMeans(test[, 1:(ncol(test)-2), drop = FALSE], na.rm = TRUE)
          # Iterate over each batch
          for (b in levels(batch_train)) {
            batch_indices <- which(batch_train == b)  # Indices for samples in batch `b`
            
            if (length(batch_indices) == 0) {
              warning(paste("Batch", b, "is empty. Skipping."))
              next
            }
            
            # Extract the subset of test for the current batch
            batch_data <- train[batch_indices, 1:(ncol(train)-2), drop = FALSE]
            
            # Compute row-wise means for this batch
            batch_mean <- colMeans(batch_data, na.rm = TRUE)
            
            # Compute the offset: batch mean - global mean
            offset <- batch_mean - global_mean
            
            if (length(batch_mean) == 0) {
              warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
              next
            }
            
            # Subtract batch means
            train_corrected[batch_indices, 1:(ncol(train)-2)] <- sweep(
              batch_data, 
              1, 
              offset, 
              "-"
            )
          }
          
          train <- train_corrected # %>% subset(select=c(-site))
          
          train$chr <- factor(train$chr)
          train$chr <- factor(make.names(levels(train$chr))[train$chr])
          
          tree <- c(50, 100, 250, 500)
          n.tree <- sample(tree,1)
          nodeSize <- seq(1,(nrow(train)/10), by=1)
          node.size <- sample(nodeSize,1)
          
          if (ncol(train)>3){
            tune_grid_temp <- data.frame(mtry=c(NA,NA,NA,NA,NA))
            tune_grid_temp$mtry <- sample(tuneGrid$mtry,5)
            
            # Set up the trainControl with the custom summary function
            control <- trainControl(method = 'cv', 
                                    number = 5, 
                                    classProbs = TRUE,
                                    summaryFunction = customSummary,
                                    search = 'random')
          } else {
            tune_grid_temp <- tuneGrid
            
            # Set up the trainControl with the custom summary function
            control <- trainControl(method = 'cv', 
                                    number = 5, 
                                    classProbs = TRUE,
                                    search = 'random')
          }
          
          if(ncol(train)>2){
            
            # Train the model using the custom F1 metric
            inner_model <- train(chr ~ .,
                                 data = train,
                                 method = "rf",
                                 metric = "F1",
                                 tuneGrid = tune_grid_temp,
                                 tuneLength=10,
                                 ntree = n.tree,
                                 nodesize=node.size,
                                 trControl = control)
          } else {
            inner_model <- train(chr ~ .,
                                 data = train,
                                 method = "rf",
                                 metric = "Accuracy",
                                 tuneGrid = tune_grid_temp,
                                 tuneLength=10,
                                 ntree = n.tree,
                                 nodesize=node.size,
                                 trControl = control)
          }
          print(inner_model)
          cat("Inner model fitted \n")
          
          if (ncol(train)>2){
            if (all(!is.na(inner_model$results$F1))){
              # Store the F1 for each inner repeat
              inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                  data.frame(min.node.size = node.size,
                                                             tree = n.tree,
                                                             mtry = inner_model$results$mtry, 
                                                             F1 = inner_model$results$F1, 
                                                             repeat_number = inner_rep))
            }
          } else {
            if (all(!is.na(inner_model$results$Accuracy))){
              # Store the F1 for each inner repeat
              inner_results[[inner_rep]] <- rbind(inner_results[[inner_rep]], 
                                                  data.frame(min.node.size = node.size,
                                                             tree = n.tree,
                                                             mtry = inner_model$results$mtry, 
                                                             F1 = inner_model$results$Accuracy, 
                                                             repeat_number = inner_rep))
            }
          }
        }
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      print(all_inner_results)
      
      # Calculate the average F1 for each hyperparameter combination
      avg_results <- all_inner_results %>%
        group_by(mtry) %>%
        summarise(avg_F1 = mean(F1, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(tree, na.rm=TRUE),
                  nodesize=mean(min.node.size, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparameter combination based on the highest average F1
      best_inner_result <- avg_results[which.max(avg_results$avg_F1), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_F1 <- best_inner_result$avg_F1
      best_F1_SE <- best_inner_result$SE_F1
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      cat("Inner results saved \n")
      
      # Train the final model using the best hyperparameters on the training data
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      if (ncol(train)>3) {
        
        
        final_model <- train(chr ~ .,
                             data = train,
                             method = "rf",
                             metric = "F1",
                             ntree = best_ntree,
                             nodesize=best_min.node.size,
                             trControl=control_final,
                             tuneGrid = repGrid)
        
      } else {
        final_model <- train(chr ~ .,
                             data = train,
                             method = "rf",
                             metric = "Accuracy",
                             ntree = best_ntree,
                             nodesize=best_min.node.size,
                             trControl=control_final,
                             tuneGrid = repGrid)
      }
      
      cat("Final model fitted \n")
      cat("Predicting probabilities and classes...\n")
      
      cat("Outer rep:", outer_rep, " ; Outer fold:", i)
      # Predict the linear predictors (PI) from the model
      if ((i!=2 & outer_rep!=2) & (i!=4 & outer_rep!=3) & (i!=1 & outer_rep!=5) & (i!=5 & outer_rep!=6)){
        test$PI <- predict(final_model, newdata = test, type = "prob")[,2]
        test$chr_pred <- predict(final_model, newdata = test, type = "raw")
        cat("PI generated \n")
        
        cm <- confusionMatrix(data = test$chr_pred, reference = test$chr)
        test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
        
        # Fit a logistic regression model on the test set using penalized coefficients
        model_test <- glm(chr ~ PI, data = test, family="binomial")
        
        c_stat_nested <- rbind(c_stat_nested, data.frame(
          C_test = concordance(model_test)$concordance,
          SE_test = concordance(model_test)$cvar,
          Fold = i,
          OuterRepeat = outer_rep,
          n_train = nrow(train),
          events_train = sum(train$chr=="X1"),
          n_test = nrow(test),
          events_test = sum(test$chr=="X1"),
          balanced_accuracy = cm$byClass[11],
          sensitivity = cm$byClass[1],
          specificity = cm$byClass[2],
          ppv = cm$byClass[3],
          npv = cm$byClass[4],
          precision = cm$byClass[5],
          recall = cm$byClass[6],
          f1 = cm$byClass[7]
          
        ))
        
        cat("Discrimination results saved \n")
        
        # Calculate calibration slope on the hold-out data (test set)
        calibration_intercept <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[12])
        calibration_slope <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[13])
        brier <- unname(rms::val.prob(p=test$PI, y=as.numeric(test$chr)-1, m=200, pl=F)[11])
        
        # Calculate ICI
        loess.calibrate <- loess(as.numeric(test$chr)-1 ~ test$PI)
        
        # Estimate loess-based smoothed calibration curve
        P.calibrate <- predict(loess.calibrate, newdata = test)
        
        # This is the point on the loess calibration curve corresponding to a given predicted probability.
        ICI <- mean(abs(P.calibrate - test$PI))
        
        # Store calibration results
        calibration_slopes <- rbind(calibration_slopes, data.frame(
          fold = i,
          OuterRepeat = outer_rep,
          intercept = calibration_intercept,
          slope = calibration_slope,
          brier = brier,
          ICI = ICI
        ))
      
      }
    }
  }
  
  list(
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes,
    dca_clinical = dca_clinical,
    dca_general = dca_general
  )
}

##### Harmonised RSF with repeated nested cross validation #####
harmonised_RSF_repeated_nested_cv <- function(combined_df, outerFolds, outerRepeats, innerFolds, innerRepeats, tuneGrid, seed) {
  
  set.seed(seed)
  seeds <- sample(1:10000, outerRepeats) # Generate unique seeds for each outer repeat
  
  best_inner_result_list <- list()
  best_mtry_list <- list()
  best_ntree_list <- list()
  best_nodesize_list <- list()
  inner_results <- data.frame()
  all_inner_results <- data.frame()
  c_stat_nested  <- data.frame()
  calibration_slopes <- data.frame()
  
  # Outer loop for cross-validation
  for (outer_rep in 1:outerRepeats) {
    cat("Outer Repeat:", outer_rep, "\n")
    set.seed(seeds[outer_rep]) # Set a unique seed for each outer repeat
    outer_folds <- createFolds(combined_df$Transition, k = outerFolds, list = TRUE, returnTrain = TRUE)
    
    for (i in seq_along(outer_folds)) {
      cat("  Outer Fold:", i, "\n")
      
      if (ncol(combined_df)>6) {
        combined_df[,c(1:15)] <- lapply(combined_df[,c(1:15)], factor)
        options(na.action='na.pass')
        combined_df <- combined_df[, sapply(combined_df, function(x) length(unique(na.omit(x))) > 1)] # Remove columns with zero variance
        
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        colnames(PPS_scored_mat) <- gsub(" ", "", colnames(PPS_scored_mat))
        PPS_scored_mat[,c(1:21)] <- lapply(PPS_scored_mat[,c(1:21)], factor)
        PPS_scored_mat <- PPS_scored_mat %>% select(-Gender0)
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        # Identify and remove near-zero variance predictors
        near_zero_vars <- nearZeroVar(PPS_scored_mat)
        PPS_scored_mat <- PPS_scored_mat[, -near_zero_vars]
        
      } else if (ncol(combined_df)>3){
        # If there is more than one predictor, use LASSO
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        }        
        
      } else {
        PPS_scored_mat <- as.data.frame(model.matrix(~.-1,combined_df))
        if ("Transition1" %in% colnames(PPS_scored_mat)) {
          PPS_scored_mat <- PPS_scored_mat %>% dplyr::rename(Transition = Transition1)
        } 
      }
      
      train_indices <- outer_folds[[i]]
      train_base <- PPS_scored_mat[train_indices, ] #  This splits the data into train (inner loop)
      test <- PPS_scored_mat[-train_indices, ] #  This splits the data into test (outer loop)
      
      ### Impute missing data in test ###
      if (ncol(test)>3) {
        imputed_data <- missForest(test)
        test <- imputed_data$ximp
        if (anyNA(test)) {
          imputed_data <- missForest(test)
          test <- imputed_data$ximp
          imputed_data <- mice(test, m = 1, method = 'pmm', maxit = 50, seed = 123)
          test <- complete(imputed_data)
        }
      } else {
        imputed_data <- mice(test)
        test <- complete(imputed_data)
      }
      
      if (anyNA(test)) {
        stop("Imputation did not fill all missing values! i=", i, " outer_rep=", outer_rep)
      }
      
      ### Mean Offset Correction ###
      batch_test <- as.factor(test$site)
      
      ### Compute Global Means ###
      if (ncol(test)>3) {
      test_num <- test %>%
        mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
        mutate(across(where(is.factor), as.numeric))   # Convert factor to numeric
      test_num[,c(1:(ncol(test_num)-3))] <- test_num[,c(1:(ncol(test_num)-3))]
      
      global_mean <- colMeans(test_num[, 1:(ncol(test)-3), drop = FALSE], na.rm = TRUE)
      
      test_corrected <- test_num # Start with the original data
      
      for (b in levels(batch_test)) {
        batch_indices <- which(batch_test == b) # Indices for samples in batch `b`
        
        if (length(batch_indices) == 0) {
          warning(paste("Batch", b, "is empty. Skipping."))
          next
        }
        
        # Extract the subset of test for the current batch
        batch_data <- test_num[batch_indices, 1:(ncol(test)-3), drop = FALSE] # Exclude site and outcome
        
        # Compute means for each predictor in this batch
        batch_mean <- colMeans(batch_data, na.rm = TRUE)
        
        if (length(batch_mean) == 0) {
          warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
          next
        }
        
        # Compute the offset: batch mean - global mean
        offset <- batch_mean - global_mean
        
        # Subtract the offset to align the batch with the global mean
        test_corrected[batch_indices, 1:(ncol(test)-3)] <- sweep(
          batch_data,
          1,
          offset,
          "-"
        )
      }
      } else {
        test_num <- test %>%
          mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
          mutate(across(where(is.factor), as.numeric))   # Convert factor to numeric
        test_num[,c(1)] <- test_num[,c(1)]
        
        global_mean <- colMeans(test_num[, 1, drop = FALSE], na.rm = TRUE)
        
        test_corrected <- test_num # Start with the original data
        
        for (b in levels(batch_test)) {
          batch_indices <- which(batch_test == b) # Indices for samples in batch `b`
          
          if (length(batch_indices) == 0) {
            warning(paste("Batch", b, "is empty. Skipping."))
            next
          }
          
          # Extract the subset of test for the current batch
          batch_data <- test_num[batch_indices, 1, drop = FALSE] # Exclude site and outcome
          
          # Compute means for each predictor in this batch
          batch_mean <- colMeans(batch_data, na.rm = TRUE)
          
          if (length(batch_mean) == 0) {
            warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
            next
          }
          
          # Compute the offset: batch mean - global mean
          offset <- batch_mean - global_mean
          
          # Subtract the offset to align the batch with the global mean
          test_corrected[batch_indices, 1] <- sweep(
            batch_data,
            1,
            offset,
            "-"
          )
        }
      }
      test <- test_corrected # %>% subset(select=c(-site))
      
      inner_results <- vector("list", innerRepeats) # Initialize as a list to store results for each inner repeat
      
      # Inner Cross-Validation and Model Training
      for (inner_rep in 1:innerRepeats) {
        cat("    Inner Repeat:", inner_rep, "\n")
        
        # Use a different seed for each inner repeat
        set.seed(seeds[outer_rep] + inner_rep)
        
        ### Impute missing data ###
        if (ncol(train_base)>3) {
          imputed_data <- missForest(train_base)
          train <- imputed_data$ximp
        } else {
          imputed_data <- mice(train_base)
          train <- complete(imputed_data)
        }
        
          ### Mean offset correction ###
          batch_train <- as.factor(train$site)
          
          # Convert all variables to numeric for matrix
          train_num <- train %>%
            mutate(across(where(is.character), as.numeric)) %>%   # Convert character to numeric
            mutate(across(where(is.factor), as.numeric))  # Convert factor to numeric
          
          if (ncol(train)>3) {
          train_num[,c(1:(ncol(train_num)-3))] <- train_num[,c(1:(ncol(train_num)-3))]
          
          train_corrected <- train_num
          
          ### Compute Global Means ###
          global_mean <- colMeans(train_num[, 1:(ncol(train_num)-3), drop = FALSE], na.rm = TRUE)
          # Iterate over each batch
          for (b in levels(batch_train)) {
            batch_indices <- which(batch_train == b)  # Indices for samples in batch `b`
            
            if (length(batch_indices) == 0) {
              warning(paste("Batch", b, "is empty. Skipping."))
              next
            }
            
            # Extract the subset of test for the current batch
            batch_data <- train_num[batch_indices, 1:(ncol(train_num)-3), drop = FALSE]
            
            # Compute row-wise means for this batch
            batch_mean <- colMeans(batch_data, na.rm = TRUE)
            
            # Compute the offset: batch mean - global mean
            offset <- batch_mean - global_mean
            
            if (length(batch_mean) == 0) {
              warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
              next
            }
            
            # Subtract batch means
            train_corrected[batch_indices, 1:(ncol(train_corrected)-3)] <- sweep(
              batch_data, 
              1, 
              offset, 
              "-"
            )
          }
          } else {
            train_num[,c(1)] <- train_num[,c(1)]
            
            train_corrected <- train_num
            
            ### Compute Global Means ###
            global_mean <- colMeans(train_num[, 1, drop = FALSE], na.rm = TRUE)
            # Iterate over each batch
            for (b in levels(batch_train)) {
              batch_indices <- which(batch_train == b)  # Indices for samples in batch `b`
              
              if (length(batch_indices) == 0) {
                warning(paste("Batch", b, "is empty. Skipping."))
                next
              }
              
              # Extract the subset of test for the current batch
              batch_data <- train_num[batch_indices, 1, drop = FALSE]
              
              # Compute row-wise means for this batch
              batch_mean <- colMeans(batch_data, na.rm = TRUE)
              
              # Compute the offset: batch mean - global mean
              offset <- batch_mean - global_mean
              
              if (length(batch_mean) == 0) {
                warning(paste("Batch", b, "has no valid data for mean computation. Skipping."))
                next
              }
              
              # Subtract batch means
              train_corrected[batch_indices, 1] <- sweep(
                batch_data, 
                1, 
                offset, 
                "-"
              )
            }
          }
          train <- train_corrected# %>% subset(select=c(-site))
        
        library(randomForestSRC)
        library(survcomp)
        
        tree <- c(50, 100, 250, 500)
        nodeSize <- seq(1,(nrow(train)/10), by=1)
        
        # Store fold-specific C-index scores
        fold_cindex <- numeric(50)
        
        for (tune in 1:10) {
          for (fold_idx in 1:5) {
            
            set.seed <- fold_idx
            # Get the training and validation sets for this fold
            folds_cv <- createFolds(train$Transition, k = 5, returnTrain = TRUE)
            
            n.tree <- sample(tree,1)
            node.size <- sample(nodeSize,1)
            mtry <- sample(tuneGrid$mtry,1)
            
            train_indices_cv <- folds_cv[[fold_idx]]
            
            fold_train <- train[train_indices_cv, ]
            fold_valid <- train[-train_indices_cv, ]
            
            predictor_vars <- setdiff(names(fold_train), c("day_exit", "Transition"))
            
            # Dynamically create the formula for the RSF model
            formula <- as.formula(paste("Surv(day_exit, Transition) ~", paste(predictor_vars, collapse = " + ")))
            
            # Train the RSF model on the training set
            rsf_model <- rfsrc(formula, data = fold_train, ntree = n.tree, mtry = mtry, nodesize = node.size)
            
            # Calculate the C-index on the validation set
            rsf_pred <- predict(rsf_model, fold_valid)
            fold_cindex[tune] <- concordance.index(rsf_pred$survival[,ncol(rsf_pred$survival)], surv.time = fold_valid$day_exit, surv.event = fold_valid$Transition)$c.index
            
          }
        }
        
        # Store the mean C-index for this combination of hyperparameters
        mean_cindex <- mean(fold_cindex)
        inner_results <- rbind(inner_results, data.frame(ntree = n.tree, mtry = mtry, nodesize=node.size, mean_cindex = mean_cindex))
        
        
        # View the results of cross-validation
        print(inner_results)
        
        # Find the best hyperparameter combination based on the highest C-index
        best_params <- inner_results[which.max(inner_results$mean_cindex), ]
        
      }
      
      # Combine results from all repeats
      all_inner_results <- do.call(rbind, inner_results)
      print(dim(all_inner_results))
      
      # Calculate the average C-index for each hyperparameter combination
      avg_results <- inner_results %>%
        group_by(mtry) %>%
        summarise(avg_C = mean(mean_cindex, na.rm=TRUE),
                  mtry = mean(mtry, na.rm=TRUE),
                  ntree = mean(ntree, na.rm=TRUE),
                  nodesize=mean(nodesize, na.rm=TRUE),
                  .groups = "keep") 
      
      # Determine the best hyperparametercombination based on the highest average C-index
      best_inner_result <- avg_results[which.max(avg_results$avg_C), ]
      best_min.node.size <- best_inner_result$nodesize
      best_mtry <- best_inner_result$mtry
      best_C <- best_inner_result$avg_C
      best_ntree <- best_inner_result$ntree
      
      # Calculate the correct index for the current outer fold and repeat
      index <- (outer_rep - 1) * outerFolds + i
      
      # Store the best hyperparameters for the current outer fold
      best_inner_result_list[[index]] <- best_inner_result
      best_mtry_list[[index]]  <- best_mtry
      best_ntree_list[[index]]  <- best_ntree
      best_nodesize_list[[index]]  <- best_min.node.size
      
      # Train the final model using the best hyperparameters on the inner training set (x_inner_train, y_inner_train)
      control_final <- trainControl(method = 'none')
      repGrid <- data.frame(mtry=best_mtry)  
      
      final_model <- rfsrc(formula, data = train, ntree = best_ntree, mtry = best_mtry, nodesize = best_min.node.size)
      
      # Predict the linear predictors (PI) from the Elastic Net model
      rsf_pred_test <- predict(final_model, test)
      test$PI <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      test$Tx_pred <- 1-rsf_pred_test$survival[,ncol(rsf_pred_test$survival)]
      
      test <- test %>% mutate(PI=case_when(PI==0 ~ 0.001, PI==1 ~ 0.999, TRUE ~ PI))
      
      # Fit a Cox model on the test set using penalized coefficients
      model_test <- coxph(Surv(day_exit, Transition) ~ PI, data = test)
      
      c_stat_nested <- rbind(c_stat_nested, data.frame(
        C_test = concordance(model_test)$concordance,
        SE_test = concordance(model_test)$cvar,
        Fold = i,
        OuterRepeat = outer_rep,
        n_train = nrow(train),
        events_train = sum(train$Transition==1),
        n_test = nrow(test),
        events_test = sum(test$Transition==1)
        
      ))
      
      # Fit a logistic regression model
      cal <- pred_val_probs(binary_outcome=test$Transition, Prob=test$Tx_pred, cal_plot=FALSE)
      calibration_in_large <- cal$CalInt
      calibration_slope <- cal$CalSlope
      brier <- cal$BrierScore
      
      cat("Calibration saved \n")
      
      cat("Outer rep = ", outer_rep, " outer fold = ", i,"\n")
      
      cat("Brier saved \n")
      
      # Calculate ICI
      if (sd(test$PI)!=0){
        loess.calibrate <- loess(test$Transition ~ test$Tx_pred)
        if (sum(!is.na(loess.calibrate$fitted))==length(loess.calibrate$fitted)){
          # Estimate loess-based smoothed calibration curve
          P.calibrate <- predict(loess.calibrate, newdata = test)
          
          # This is the point on the loess calibration curve corresponding to a given predicted probability.
          ICI <- mean(abs(P.calibrate - test$Tx_pred))
        } else {
          ICI <- NA
        }
      } else {
        ICI <- NA
      }
      # Store calibration results
      calibration_slopes <- rbind(calibration_slopes, data.frame(
        fold = i,
        OuterRepeat = outer_rep,
        calibration_in_large = calibration_in_large,
        slope = calibration_slope,
        brier = brier,
        ICI = ICI
      ))
      
    }
  }
  
  list(
    best_inner_result_list = best_inner_result_list,
    best_mtry_list = best_mtry_list,
    best_ntree_list = best_ntree_list,
    best_nodesize_list = best_nodesize_list,
    c_stat_nested = c_stat_nested,
    calibration_slopes = calibration_slopes
  )
}
