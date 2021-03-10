# build base map for fertiliser plot
get_basemap <- function(){
  
  # download full basemap
  base_map <- getMap(resolution = "high")
  
  # convert to correction projection
  proj4string(base_map) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
  
  # return basemap
  return(base_map)
}

### function for extracting covariance matrix and calculating uncertainties
iterate_covar <- function(i, model, prediction_data, factor_no_1, factor_no_2, fixed_effect_no, neg_binom){
  
  # if the family is negative binomial zero inflated, need to index for the fixed effects
  if(neg_binom == TRUE){
    
    # extract fixed effects from covariance matrix 
    coefs <- mvrnorm(n = 1, mu = fixef(object = model)[1]$cond, Sigma = vcov(object = model)[1]$cond)
  }
  
  # if not negative binomial, can just take the fixed effects straight from the model
  else{
    coefs <- mvrnorm(n = 1, mu = fixef(object = model), Sigma = vcov(object = model))
  }
  
  # then extract the terms from the model matrix
  mm <- model.matrix(terms(model), prediction_data)
  
  # set up for loop for adjusting for seperate taxonomic class
  counter <- 1
  subset_predictions <- c()
  
  if(fixed_effect_no == 2){
    
    # remove extra factor for apodiformes
    mm <- mm[,-24]
    y <- mm %*% coefs
    
    # run for loop to take exponential, and adjust as rescaled percentage for each sample
    for(i in 1:(factor_no_1)){
      
      if(i < 6){
        subset_predictions <- c(subset_predictions, (exp(y[counter:(counter+(factor_no_2-1))]) / exp(y)[counter] * 100))
      }
      
      # since apodiformes has one less factor, reduce the last step by 1
      else{
        subset_predictions <- c(subset_predictions, (exp(y[counter:(counter+(factor_no_2-2))]) / exp(y)[counter] * 100))
      }
      
      # step up the counter to move along the vector
      counter <- counter + factor_no_2
    }
  }
  
  # if fixed effect is 1, calculate effect sizes relative to the first value (primary vegetation)
  else{
    y <- mm %*% coefs
    subset_predictions <- (exp(y)/exp(y)[1])*100
  }
  
  # return the vector of adjusted values for that sample
  return(subset_predictions)
}

### run the predicted value function, extract median, upper, and lower confidence, and plot
predict_effects <- function(model, 
                            iterations, 
                            model_data, 
                            response_variable, 
                            factor_number_1,
                            factor_number_2,
                            fixed_number,
                            fixed_column,
                            neg_binom){
  
  # print the fixed effects for the model as a reference for predictions
  print(fixef(model))
  
  # if number of fixed effects = 1, create a dataframe from the model_data with one fixed and response
  if(fixed_number == 1){
    unique_data <- data.frame(factor(as.character(levels(model_data[,fixed_column])),
                                     levels = levels(model_data[,fixed_column])), 0) %>%
      unique()
    
    # assigned the correct column names according to fixed effect and response arguments
    colnames(unique_data) <- c(fixed_column, response_variable)
  }
  
  # if fixed effect is not 1 (2), build dataframe for the 2 fixed effects
  else{
    # set up the dataframe of unique factors and reorder by the model data
    unique_data <- data.frame(factor(as.character(model_data[, fixed_column[1]]),
                                     levels = levels(model_data[, fixed_column[1]])),
                              factor(as.character(model_data[, fixed_column[2]]),
                                     levels = levels(model_data[, fixed_column[2]])), 0) %>%
      unique()
    
    # reorder the data frame by factor order
    colnames(unique_data) <- c(fixed_column[1], fixed_column[2], response_variable)
    col_1 <- unique_data %>% pull(fixed_column[1])
    col_2 <- unique_data %>% pull(fixed_column[2])
    unique_data <- unique_data[with(unique_data, order(col_1, col_2)), ]
  }
  
  # run iterations of extraction of coefficients
  preds.emp <- sapply(X = 1:iterations, iterate_covar, model, prediction_data = unique_data, factor_no_1 = factor_number_1, factor_no_2 = factor_number_2, fixed_effect_no = fixed_number, neg_binom)
  
  # extract the median, upper interval, and lower interval for samples
  preds.emp.summ <- data.frame(Median = apply(X = preds.emp, MARGIN = 1, FUN = median),
                               Upper = apply(X = preds.emp, MARGIN = 1, FUN = quantile, probs = 0.975),
                               Lower = apply(X = preds.emp, MARGIN = 1, FUN = quantile, probs = 0.025))
  
  # bind prediction data back onto the median, upper, and lower intervals, and adjust as percentage change - note factor order
  fin_conf <- cbind(unique_data, preds.emp.summ) %>%
    mutate(Median = Median - 100) %>%
    mutate(Upper = Upper - 100) %>%
    mutate(Lower = Lower - 100)
  
  if(fixed_number == 2){
    # convert the reference factor (primary vegetation) to NA so doesn't plot error bar
    fin_conf$Upper[fin_conf$LUI == "Primary vegetation"] <- NA
    fin_conf$Lower[fin_conf$LUI == "Primary vegetation"] <- NA
    
    # print final dataframe before plotting
    print(fin_conf)
    
    # construct the plot for that biodiversity metric - note the y axis label needs to be specific for that metric
    output_plot <- ggplot(fin_conf) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = fin_conf[,1], y = Median, colour = fin_conf[,2]), position = position_dodge(0.5), size = 2) +
      geom_errorbar(aes(x = fin_conf[,1], ymin = Lower, ymax = Upper, colour = fin_conf[,2]), position = position_dodge(0.5), size = 1.5, width = 1/(factor_number_2)) +
      scale_colour_manual(values = c("black", "yellow2", "orange", "red"), name = "Use intensity") +
      scale_y_continuous(name = gsub(x = paste(response_variable, "difference (%)", sep = " "), pattern = "_", replacement = " ")) +
      theme_bw() +
      theme(panel.grid = element_blank())
  }
  
  # if number of fixed effects =  1 (and single combined LUI) convert LUI into two columns for ggplot 
  else{
    # set up vectors for intensity and land-use types - can be moved as an argument?
    land_use_intensity <- c("MU", "LU", "IU")
    land_use_type <- c("PV", "MSV", "ISV", "YSV", "PF", "P", "C", "U")
    fin_conf$veg_type <- fin_conf[,fixed_column]
    
    # create empty vectors, check for each intensity in veg type, and assign to a new vector
    assign_intensity <- c()
    intensity <- c()
    for(j in 1:length(fin_conf$veg_type)){
      for(i in 1:length(land_use_intensity)){
        assign_intensity <- grep(land_use_intensity[i], fin_conf$veg_type[j])
        if(length(assign_intensity > 0)){
          intensity <-  c(intensity, land_use_intensity[i])
        }
      }
    }
    
    # assign the correctly sorted intensity to a column of dataframe 
    fin_conf$intensity <- intensity 
    
    # amend columns - veg_type
    fin_conf$veg_type <- gsub("MU", "", fin_conf$veg_type)
    fin_conf$veg_type <- gsub("LU", "", fin_conf$veg_type)
    fin_conf$veg_type <- gsub("IU", "", fin_conf$veg_type)
    
    # order factors for intensity and land use type
    fin_conf$intensity <- factor(fin_conf$intensity, levels = land_use_intensity,
                                 labels = c("Minimal use", "Light use", "Intense use"))
    fin_conf$veg_type <- factor(fin_conf$veg_type, levels = land_use_type, 
                                labels = c("Primary", "MSV", "ISV", "YSV", "Plantation", "Pasture", "Cropland", "Urban"))
    
    # print final dataframe before plotting
    print(fin_conf)
    
    # construct the plot for that biodiversity metric - note the y axis label needs to be specific for that metric
    output_plot <- ggplot(fin_conf) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_point(aes(x = veg_type, y = Median, shape = intensity, colour = veg_type), position = position_dodge(0.75), size = 2.5) +
      geom_errorbar(aes(x = veg_type, ymin = Lower, shape = intensity, ymax = Upper, colour = veg_type), position = position_dodge(0.75), size = 1, width = 0.25) +
      scale_colour_manual(values = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")) +
      scale_y_continuous(name = gsub(x = paste(response_variable, "difference (%)", sep = " "), pattern = "_", replacement = " ")) +
      xlab(NULL) +
      labs(shape = NULL) +
      guides(colour = FALSE) +
      theme_bw() +
      theme(legend.position = c(0.15, 0.2), legend.background = element_rect(colour = "black"), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5, size = 13, colour = c("#E69F00", "#009E73", "#F0E442","#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000")))
  }
  
  # return the plot for that biodiversity metric
  return(output_plot)
  
}

# function for predicting continuous variable
predict_continuous <- function(model,
                               model_data, 
                               response_variable,
                               categorical_variable,
                               continuous_variable,
                               continuous_transformation,
                               random_variable,
                               colour_palette){
  
  # set up the prediction dataframe
  prediction_data <- model_data[, c(response_variable, 
                                    random_variable[1], 
                                    random_variable[2], 
                                    random_variable[3], 
                                    categorical_variable[1], 
                                    continuous_variable[1])]
  
  # remove any incomplete rows (NAs) from the prediction data
  prediction_data <- prediction_data[complete.cases(prediction_data),]
  
  # predict the values for the model
  y_value <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[1]]
  y_value_plus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[2]]
  y_value_minus <- c(StatisticalModels::PredictGLMER(model, data = prediction_data, se.fit = TRUE, seMultiplier = 1.96))[[3]]
  
  # bind the predicted values to the prediction data
  bound_values <- data.frame(cbind(prediction_data,
                                   y_value, 
                                   y_value_plus, 
                                   y_value_minus, 
                                   metric = response_variable,
                                   continuous_transformation(prediction_data[, continuous_variable[1]])))
  
  # rename last column after transformation
  colnames(bound_values)[ncol(bound_values)] <- paste(continuous_variable[1], "transform", sep = "_")
  
  # rename the response variable column "response_variable" for later plot function
  bound_values <- bound_values %>%
    rename("response_variable" = all_of(response_variable))
  
  # print the final dataframe of predicted values
  print(bound_values)
  
}

# function for plotting the output from a continuous variable glmer - predict_continuous()
# need to amend to read in categorical variable rather than zone/order
plot_fert_response <- function(data_set, categorical_variable){
  plot_obj <- ggplot(data_set) +
    geom_line(aes_string(x = "fert_transform", y = "y_value", colour = categorical_variable, linetype = "significant"), size = 0.7) +
    geom_ribbon(aes_string(x = "fert_transform", ymin = "y_value_minus", ymax = "y_value_plus", fill = categorical_variable), alpha = 0.4) +
    scale_x_continuous(breaks = c(-1, 0, 1, 2, 2.39794, 2.69897, 3, 3.39794), labels = c(0.1, 1, 10, 100, 250, 500, 1000, 2500)) +
    scale_linetype_manual("", values = c("solid", "dashed"), labels = c("Non-significant", "Significant")) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(plot_obj)
}
