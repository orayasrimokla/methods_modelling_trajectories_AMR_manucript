
# input: 
#   model, 
#   datatable with covariates. Always starts at same time as I0 and is sorted either increasing or decreasing in time per country
# output: datatable with resistance quantiles
run_predictions <- function(model, data,model_name){
  
  if(grepl("spatial", model_name)){
    
    ext_fit <- model
    ext_fit <- model[[1]]$draws(format = "df")
    n_total_draws <- min(1000,length(ext_fit$lp__))
    sample_post <- sample(length(ext_fit$lp__), n_total_draws)
    
    
  }else{
    # extract posterior samples
    ext_fit <- extract(model)
    # get total number of posterior samples
    n_total_draws <- min(1000,nrow(ext_fit$lp__))
    # Sample from posterior (or use all draws)
    sample_post <- sample(nrow(ext_fit$lp__), n_total_draws)
  }
  
  #get dimensions for predictions
  country_list <- unique(data$country)
  
  
  n_pred = nrow(data)
  
  y_pred <- matrix(NA, nrow = length(sample_post), ncol = n_pred)
  
  #check these are the right columns in data for the covariates
  cov_mat <- do.call(cbind, data[,(ncol(data)-length(covar_list) + 1):ncol(data)])
  
  #this is the original loop for the train set 
  for(s in 1:length(sample_post)) {
    
    draw_idx <- sample_post[s]
    
    # Extract parameters for this draw
    phi <- ext_fit$phi[draw_idx, ] #per country
    sigma_phi <- ext_fit$sigma_phi[draw_idx] 
    eta <- ext_fit$eta[draw_idx, ] #per covariate
    
    gamma <- ext_fit$gamma[draw_idx] 
    beta <- ext_fit$beta[draw_idx] 
    
    I0_logit_c <- ext_fit$I0_logit_c[draw_idx,] 
    I0_prob_c <- ext_fit$I0_prob_c[draw_idx,] 
    gamma_c <- ext_fit$gamma_c[draw_idx,] 
    beta_c <- ext_fit$beta_c[draw_idx,] 
    
    #with covariates
    log_beta_c <- ext_fit$log_beta_c[draw_idx,] 
    log_gamma_c <- ext_fit$log_gamma_c[draw_idx,] 
    
    
    #with source
    source_eta_S <- ext_fit$source_eta_S[draw_idx,]
    source_eta <- ext_fit$source_eta[draw_idx]
    
    #with spatial
    sigma_phi_I0 <- ext_fit$sigma_phi_I0[draw_idx] 
    sigma_phi_beta <- ext_fit$sigma_phi_beta[draw_idx] 
    sigma_phi_gamma <- ext_fit$sigma_phi_gamma[draw_idx] 
    
    if(grepl("spatial", model_name)){
      
      I0_logit <- c()
      I0_prob <- c()
      log_beta <- c()
      log_gamma <- c()
      
      k_I0 <- c()
      k_beta <- c()
      k_gamma <- c()
      
      #source and covar
      source_eta_S <- c()
      eta <- c()
      

      regionOrder <- sort(unique(data$owid_region))
      for(region in 1:length(regionOrder)){
        
        ext_fit <- model[[region]]$draws(format = "df")
        
        I0_prob[region] <- ext_fit$I0_prob[draw_idx] 
        
        I0_logit[region] <- ext_fit$I0_logit[draw_idx] 
        log_beta[region] <- ext_fit$log_beta[draw_idx] 
        log_gamma[region] <- ext_fit$log_gamma[draw_idx] 
        
        dataTemp <- data[data$owid_region == regionOrder[region],]
        #need to get length of countries for each region
        num_country_train <- unique(data$country_track)
        num_country_train <- length(num_country_train[num_country_train != -1000])
        
        #spatial with cmdstan - need to map to each region and country that we are extracting - use this to match the actual country ID 
        k_I0[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_I0[",1:num_country_train,"]")])
        k_beta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k[",1:num_country_train,"]")])
        k_gamma[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_gamma[",1:num_country_train,"]")])
        
        #source and covar 
        if(grepl("Source", model_name)){
          source_eta_S[[region]] <- as.numeric(ext_fit[draw_idx,paste0("source_eta_S[",1:length(unique(data$source)),"]")])
        }else{}
        
        if(grepl("Covar", model_name)){
          eta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("eta[",1:ncol(cov_mat),"]")])
        }else{}
        
      }

    }else{
      I0_logit <- ext_fit$I0_logit[draw_idx] 
      log_beta <- ext_fit$log_beta[draw_idx] 
      log_gamma <- ext_fit$log_gamma[draw_idx]
      I0_prob <- ext_fit$I0_prob[draw_idx] 
    }
    
    #initialize tracking variables
    lastCountry <- -1
    lastYear <- -9999
    lastI <- -9999
    
    #predictions for test set - has source level but get single curve for each country - use the other reordered data
    for(i in 1:n_pred) {
      
      country_i <- data$country[i]
      source_i <- data$source[i]
      t_i <- data$t[i]
      
      country_track_i <- data$country_track[i]
      source_track_i <- data$source_track[i]
     
      initial_year_track <- data$init_year_country[i]
      
      neigh_1 <- data$neigh_1[i]
      neigh_2 <- data$neigh_2[i]
      # print(country_i)
      # print(t_i)
      
      #track regions for spatial model
      region_i <- data$owid_region_numeric[i]
      country_spatial_i <- data$countryNumSpatial[i]
      
      
      if(country_i != lastCountry) {
        if(country_track_i != -1000){
          
          if(grepl("spatial", model_name) ){

            y_pred[s, i] <- inv.logit(logit(I0_prob[region_i]) + k_I0[[region_i]][country_spatial_i])
            
           
          }else{
            
            y_pred[s, i] <- I0_prob_c[country_track_i] 
          }
          if(is.na( y_pred[s, i])){
            print(i)
            fail1
          }
          
        }else{
  
          if(grepl("spatial", model_name)){

            y_pred[s, i] <- inv.logit(logit(I0_prob) + (k_I0[neigh_1] +  k_I0[neigh_2])/2)
            
          }else{ 
      
            y_pred[s, i] <- (I0_prob_c[neigh_1] + I0_prob_c[neigh_2])/2
          }
          if(is.na( y_pred[s, i])){
            print(i)
            fail2
          }
        }
        
      } else {
        
        if(abs(t_i - lastYear) < 10^(-4)) {
          y_pred[s, i] <- lastI
        }else{
          
          if(country_track_i != -1000){
            
            if(grepl("spatial", model_name)){
              
              #adding condition for now beta with covariates
              if(grepl("spatialCmdStanSourceCovar", model_name)){
                current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i] + sum(cov_mat[i, ] * eta[[region_i]]))
              }else{
                current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i])
              }
              
              #test with the spatial regional parameters
              current_gamma <- exp(log_gamma[region_i] + k_gamma[[region_i]][country_spatial_i])#uses same neigh struc as beta 
              
            }else if(grepl("covar", model_name)){
              current_beta <- exp(log_beta_c[country_track_i] + sum(cov_mat[i, ] * eta))
              current_gamma <- gamma_c[country_track_i]
            }else{
              current_beta <- beta_c[country_track_i]
              current_gamma <- gamma_c[country_track_i]
            }

            
          }else{
            if(grepl("spatial", model_name)){
              #adding condition for now beta with covariates
              if(grepl("spatialCmdStanSourceCovar", model_name)){
                current_beta <- exp(log_beta +(k_beta[neigh_1]+ k_beta[neigh_2])/2 + sum(cov_mat[i, ] * eta[[region_i]]))
              }else{current_beta <- exp(log_beta + (k_beta[neigh_1]+ k_beta[neigh_2])/2)
              }
              current_gamma <- exp(log_gamma + (k_gamma[neigh_1]+ k_gamma[neigh_2])/2)
            }else if(grepl("covar", model_name)){
              
              current_beta <- exp((log_beta_c[neigh_1] + log_beta_c[neigh_2])/2 + sum(cov_mat[i, ] * eta))
              current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
            }else{
              current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
              current_beta <- (beta_c[neigh_1] + beta_c[neigh_2])/2
            }
            
          }
          
          if(abs(current_beta - current_gamma) < 10^(-4)){
            y_pred[s, i] <- lastI
          }else{
            
            y_pred[s, i] <- ((current_beta - current_gamma) / current_beta) /
              (1 + (((current_beta - current_gamma) / current_beta) / lastI - 1) *
                 exp(-(current_beta - current_gamma)*(t_i - lastYear)))
            
          }}
        
        
        
        if(is.na( y_pred[s, i])){
          print(i)
          fail
        }
      }
      
      
      if(y_pred[s, i] < 0 ){
        
        fail
      }

      # Update tracking
      lastI <- y_pred[s, i]
      lastCountry <- country_i
      lastYear <- t_i
      
      if(grepl("Source", model_name)){
        if(country_track_i != -1000){
          if(grepl("spatial", model_name)){
            current_source <- source_eta_S[[region_i]][source_track_i]
          }else{current_source <- source_eta_S[source_track_i]}
        }else{
          if(grepl("spatial", model_name)){
            current_source <- mean(source_eta_S[[region_i]])
          }else{current_source <- mean(source_eta_S)}
        }
        y_pred[s, i] <- (y_pred[s, i])*current_source
      }else{}
      
      y_pred[s, i] <- max(0.00001, min(0.99999, y_pred[s, i] ))
    }
  }
  

  
  # Calculate summary statistics
  predictions <- data.frame(
    median = apply(y_pred, 2, median),
    q5 = apply(y_pred, 2, quantile, 0.05),
    q95 = apply(y_pred, 2, quantile, 0.95),
    mean = apply(y_pred, 2, mean),
    sd = apply(y_pred, 2, sd))
  
  
  return(list(predictions = predictions))  
}


# input: 
#   model, 
#   datatable with covariates. Always starts at same time as I0 and is sorted either increasing or decreasing in time per country
# output: datatable with resistance quantiles

run_predictions_stacking <- function(model_weights_dt, data, abx_name, sourceData){
  
  total_samples <- 1000 #set to 1000 for nowx
  #set number of samples from the weights and round it 
  model_weights_dt[ , draw_num := round(weight*total_samples)]
  model_weights_dt <- model_weights_dt[ draw_num > 0,] #remove ones that are 0 draws 
  
  list_y_preds <- rbind(llply(1:nrow(model_weights_dt), function(model_i){ 
    print(model_i)
    if(grepl("spatial",model_weights_dt[model_i,1])){
      
      #run for each region
      fit <- readRDS(file = paste0(path_cluster,data_run,paste0(model_weights_dt[model_i,1]),"/",'model: ',paste0(model_weights_dt[model_i,1]),'+2000:2022',abx_name,".rds"))
      model <- fit$draws(format = "df")
      ext_fit <- model
      n_total_draws <- model_weights_dt[model_i,draw_num] 
      sample_post <- sample(length(ext_fit$lp__), n_total_draws)
      
      
    }else{
      #load the model in the dt 
      load(paste0(path_cluster,data_run,paste0(model_weights_dt[model_i,1]),"/",'model: ',paste0(model_weights_dt[model_i,1]),'+2000:2022',abx_name))
      
      # Extract posterior samples from your fitted model
      ext_fit <- extract(model)
      n_total_draws <- model_weights_dt[model_i,draw_num]
      sample_post <- sample(nrow(ext_fit$lp__), n_total_draws)
    }
   
    country_list <- unique(data$country)
    
    n_pred = nrow(data)
    
    y_pred <- matrix(NA, nrow = length(sample_post), ncol = n_pred)
    
    #check these are the right columns in data for covariates
    cov_mat <- do.call(cbind, data[,(ncol(data)-length(covar_list) + 1):ncol(data)])
    
    for(s in 1:length(sample_post)) {
      #print(s)
      draw_idx <- sample_post[s]
      
      # Extract parameters for this draw
      phi <- ext_fit$phi[draw_idx, ] #per country
      sigma_phi <- ext_fit$sigma_phi[draw_idx] 
      eta <- ext_fit$eta[draw_idx, ] #per covariate
     
      gamma <- ext_fit$gamma[draw_idx] 
      beta <- ext_fit$beta[draw_idx] 
      
      I0_logit_c <- ext_fit$I0_logit_c[draw_idx,] 
      I0_prob_c <- ext_fit$I0_prob_c[draw_idx,] 
      gamma_c <- ext_fit$gamma_c[draw_idx,] 
      beta_c <- ext_fit$beta_c[draw_idx,] 
      
      #test with Covar
      log_beta_c <- ext_fit$log_beta_c[draw_idx,] 
      log_gamma_c <- ext_fit$log_gamma_c[draw_idx,] 
      
      
      #test with Source
      source_eta_S <- ext_fit$source_eta_S[draw_idx,]
      source_eta <- ext_fit$source_eta[draw_idx]
      
      #Testing with Spatial
      sigma_phi_I0 <- ext_fit$sigma_phi_I0[draw_idx] 
      sigma_phi_beta <- ext_fit$sigma_phi_beta[draw_idx] 
      sigma_phi_gamma <- ext_fit$sigma_phi_gamma[draw_idx] 
      
      if(grepl("spatial",model_weights_dt[model_i,1])){
        
        I0_logit <- c()
        I0_prob <- c()
        log_beta <- c()
        log_gamma <- c()
        
        k_I0 <- c()
        k_beta <- c()
        k_gamma <- c()
        
        #source and covar
        source_eta_S <- c()
        eta <- c()
        
        regionOrder <- sort(unique(data$owid_region))
        for(region in 1:length(regionOrder)){
          
          ext_fit <- cmdStan_models_region[[region]]$draws(format = "df")
          
          #informative prior
          I0_prob[region] <- ext_fit$I0_prob[draw_idx] 
          
          I0_logit[region] <- ext_fit$I0_logit[draw_idx] 
          log_beta[region] <- ext_fit$log_beta[draw_idx] 
          log_gamma[region] <- ext_fit$log_gamma[draw_idx] 
          
          dataTemp <- data[data$owid_region == regionOrder[region],]
          
          num_country_train <- unique(data$country_track)
          num_country_train <- length(num_country_train[num_country_train != -1000])
          
          #spatial with cmdstan - need to map to each region and country that we are extracting - use this to match the actual country ID 
          k_I0[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_I0[",1:num_country_train,"]")])
          k_beta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k[",1:num_country_train,"]")])
          k_gamma[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_gamma[",1:num_country_train,"]")])
          
          #source and covar 
          if(grepl("Source", model_weights_dt[model_i,1])){
            source_eta_S[[region]] <- as.numeric(ext_fit[draw_idx,paste0("source_eta_S[",1:length(unique(data$source)),"]")])
          }else{}
          
          if(grepl("Covar", model_weights_dt[model_i,1])){
            eta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("eta[",1:ncol(cov_mat),"]")])
          }else{}
          
        }
       
      }else{
        I0_logit <- ext_fit$I0_logit[draw_idx] 
        log_beta <- ext_fit$log_beta[draw_idx] 
        log_gamma <- ext_fit$log_gamma[draw_idx]
        I0_prob <- ext_fit$I0_prob[draw_idx] 
      }
      #initialize tracking variables
      lastCountry <- -1
      lastYear <- -9999
      lastI <- -9999
      
      #predictions for test set - has source level but get single curve for each country - use the other reordered data
      for(i in 1:n_pred) {
        country_i <- data$country[i]
        source_i <- data$source[i]
        t_i <- data$t[i]
        
        country_track_i <- data$country_track[i]
        source_track_i <- data$source_track[i]
      
        initial_year_track <- data$init_year_country[i]
        
        neigh_1 <- data$neigh_1[i]
        neigh_2 <- data$neigh_2[i]
       
        #track regions for spatial model 
        region_i <- data$owid_region_numeric[i]
        country_spatial_i <- data$countryNumSpatial[i]
        
        if(country_i != lastCountry) {
          if(country_track_i != -1000){
            
            if(grepl("spatial", model_weights_dt[model_i,1])){
              
             
              y_pred[s, i] <- inv.logit(logit(I0_prob[region_i]) + k_I0[[region_i]][country_spatial_i])
              
            }else{
              
              y_pred[s, i] <- I0_prob_c[country_track_i] 
            }
            #just a safety check 
            if(is.na( y_pred[s, i])){
              print(i)
              fail1
            }
          }else{
            
            if(grepl("spatial", model_weights_dt[model_i,1])){
               y_pred[s, i] <- inv.logit(logit(I0_prob) + (k_I0[neigh_1] +  k_I0[neigh_2])/2)
            }else{ 
              
              y_pred[s, i] <- (I0_prob_c[neigh_1] + I0_prob_c[neigh_2])/2
            }
            #just a safety check 
            if(is.na( y_pred[s, i])){
              print(i)
              fail2
            }
          }} else {
            
            if(abs(t_i - lastYear) < 10^(-4)) {
              y_pred[s, i] <- lastI
            }else{
              
              if(country_track_i != -1000){
                if(grepl("spatial", model_weights_dt[model_i,1])){
                  #test with the spatial regional parameters
                  
                  #adding condition for now beta with covariates
                  if(grepl("spatialCmdStanSourceCovar", model_weights_dt[model_i,1])){
                    current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i] + sum(cov_mat[i, ] * eta[[region_i]]))
                  }else{
                    current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i])
                  }
                  
                  
                  current_gamma <- exp(log_gamma[region_i] + k_gamma[[region_i]][country_spatial_i])#uses same neigh struc as beta 
                  
                  
                }else if(grepl("covar", model_weights_dt[model_i,1])){
                  current_beta <- exp(log_beta_c[country_track_i] + sum(cov_mat[i, ] * eta))
                  current_gamma <- gamma_c[country_track_i]
                }else{
                  current_beta <- beta_c[country_track_i]
                  current_gamma <- gamma_c[country_track_i]
                }
                
                
                
              }else{
                if(grepl("spatial", model_weights_dt[model_i,1])){
                  
                  #adding condition for now beta with covariates
                  if(grepl("spatialCmdStanSourceCovar",model_weights_dt[model_i,1] )){
                    current_beta <- exp(log_beta +(k_beta[neigh_1]+ k_beta[neigh_2])/2 + sum(cov_mat[i, ] * eta[[region_i]]))
                  }else{current_beta <- exp(log_beta + (k_beta[neigh_1]+ k_beta[neigh_2])/2)
                  }
                  current_gamma <- exp(log_gamma + (k_gamma[neigh_1]+ k_gamma[neigh_2])/2)
                }else if(grepl("covar", model_weights_dt[model_i,1])){
                  
                  current_beta <- exp((log_beta_c[neigh_1] + log_beta_c[neigh_2])/2 + sum(cov_mat[i, ] * eta))
                  current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
                }else{
                  current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
                  current_beta <- (beta_c[neigh_1] + beta_c[neigh_2])/2
                }
                }
              
              
              if(abs(current_beta - current_gamma) < 10^(-4)){
                y_pred[s, i] <- lastI
              }else{
                
                y_pred[s, i] <- ((current_beta - current_gamma) / current_beta) /
                  (1 + (((current_beta - current_gamma) / current_beta) / lastI - 1) *
                     exp(-(current_beta - current_gamma)*(t_i - lastYear)))
                
                
              }}
            
            
            
            if(is.na( y_pred[s, i])){
              print(i)
              fail
            }
          }
       
        
        if(y_pred[s, i] < 0 ){
          
          fail
        }
        
        # Update tracking
        lastI <- y_pred[s, i]
        lastCountry <- country_i
        lastYear <- t_i
        
        if(grepl("Source", model_weights_dt[model_i,1])){
          if(country_track_i != -1000){
            if(grepl("spatial", model_weights_dt[model_i,1])){
              current_source <- source_eta_S[[region_i]][source_track_i]
            }else{current_source <- source_eta_S[source_track_i]}
          }else{
            if(grepl("spatial" , model_weights_dt[model_i,1])){
              current_source <- mean(source_eta_S[[region_i]])
            }else{current_source <- mean(source_eta_S)}
          }
          y_pred[s, i] <- (y_pred[s, i])*current_source
        }else{}
        
        y_pred[s, i] <- max(0.00001, min(0.99999, y_pred[s, i] ))
      }
    }
    
    
    return(y_pred)
    
  }))
  
  y_pred_allModels <-  do.call(rbind,list_y_preds)
  
  # Calculate summary statistics
  predictions <- data.frame(
    median = apply(y_pred_allModels, 2, median),
    q5 = apply(y_pred_allModels, 2, quantile, 0.05),
    q95 = apply(y_pred_allModels, 2, quantile, 0.95),
    mean = apply(y_pred_allModels, 2, mean),
    sd = apply(y_pred_allModels, 2, sd))
  
  return(list(predictions = predictions))  
  
  
}

# output: datatable with resistance quantiles
run_predictions_test_GPR <- function(model, data,model_name, GPR_list,covar_list){
  
  if(grepl("spatial", model_name)){
    
    ext_fit <- model$draws(format = "df")
    n_total_draws <- min(1000,length(ext_fit$lp__))
    sample_post <- sample(length(ext_fit$lp__), n_total_draws)
    
  }else{
    # Extract posterior samples from your fitted model
    ext_fit <- extract(model)
    # Get total number of posterior samples
    n_total_draws <- min(1000,nrow(ext_fit$lp__))
    # Sample from posterior (or use all draws)
    sample_post <- sample(nrow(ext_fit$lp__), n_total_draws)
  }
  
  # Get dimensions for predictions
  country_list <- unique(data$country)
  
  n_pred = nrow(data)
  
  y_pred <- matrix(NA, nrow = length(sample_post), ncol = n_pred)
  
  cov_mat <- do.call(cbind, data[,(ncol(data)-length(covar_list) + 1):ncol(data)])
  
  for(s in 1:length(sample_post)) {
    print(s)
    
    draw_idx <- sample_post[s]
    
    # Extract parameters for this draw
    phi <- ext_fit$phi[draw_idx, ] #per country
    sigma_phi <- ext_fit$sigma_phi[draw_idx] 
    eta <- ext_fit$eta[draw_idx, ] #per covariate
    
    gamma <- ext_fit$gamma[draw_idx] 
    beta <- ext_fit$beta[draw_idx] 
    
    I0_logit_c <- ext_fit$I0_logit_c[draw_idx,] 
    I0_prob_c <- ext_fit$I0_prob_c[draw_idx,] 
    gamma_c <- ext_fit$gamma_c[draw_idx,] 
    beta_c <- ext_fit$beta_c[draw_idx,] 
    
    
    log_beta_c <- ext_fit$log_beta_c[draw_idx,] 
    
    log_gamma_c <- ext_fit$log_gamma_c[draw_idx,] 
    
    
   
    source_eta_S <- ext_fit$source_eta_S[draw_idx,]
    source_eta <- ext_fit$source_eta[draw_idx]
    
    
    I0_prob <- ext_fit$I0_prob[draw_idx] 
   
    sigma_phi_I0 <- ext_fit$sigma_phi_I0[draw_idx] 
    sigma_phi_beta <- ext_fit$sigma_phi_beta[draw_idx] 
    sigma_phi_gamma <- ext_fit$sigma_phi_gamma[draw_idx] 
    
    if(grepl("spatial", model_name)){
      
      #from preds function below
      I0_logit <- c()
      I0_prob <- c()
      log_beta <- c()
      log_gamma <- c()
      
      k_I0 <- c()
      k_beta <- c()
      k_gamma <- c()
      
      k_I0_new <- c()
      k_beta_new <- c()
      k_gamma_new <- c()
      
      #source and covar
      source_eta_S <- c()
      eta <- c()
      
      regionOrder <- sort(unique(data$owid_region))
      for(region in 1:length(regionOrder)){
        
        ext_fit <- model[[region]]$draws(format = "df")
        
        rhosq <- ext_fit$rhosq[draw_idx]
        etasq <- ext_fit$etasq[draw_idx]
        
        
        I0_prob[region] <- ext_fit$I0_prob[draw_idx] 
        
        
        etasq_I0 <- ext_fit$etasq_I0[draw_idx]
        
      
        etasq_gamma <- ext_fit$etasq_gamma[draw_idx]
        
        I0_logit[region] <- ext_fit$I0_logit[draw_idx] 
        log_beta[region] <- ext_fit$log_beta[draw_idx] 
        log_gamma[region] <- ext_fit$log_gamma[draw_idx] 
        
        dataTemp <- data[data$owid_region == regionOrder[region],]
       
        num_country_train <- unique(dataTemp$country_track)
        num_country_train <- length(num_country_train[num_country_train != -1000])
        
        num_source_train <- unique(dataTemp$source_track)
        num_source_train <- length(num_source_train[num_source_train != -1000])
        
        #spatial with cmdstan - need to map to each region and country that we are extracting - use this to match the actual country ID 
        k_I0[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_I0[",1:num_country_train,"]")])
        k_beta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k[",1:num_country_train,"]")])
        k_gamma[[region]] <- as.numeric(ext_fit[draw_idx,paste0("k_gamma[",1:num_country_train,"]")])
        
        print(paste0("s",s))
        print(region)
        
        
        k_beta_new[[region]] <- calc_cov_mat(GPR_list$dist_mat,etasq, rhosq, 0.05, GPR_list$old_country[[region]],GPR_list$new_country[[region]],k_beta[[region]])
        k_I0_new[[region]] <- calc_cov_mat(GPR_list$dist_mat,etasq_I0, rhosq, 0.05, GPR_list$old_country[[region]],GPR_list$new_country[[region]],k_I0[[region]])
        k_gamma_new[[region]] <- calc_cov_mat(GPR_list$dist_mat,etasq_gamma, rhosq, 0.05, GPR_list$old_country[[region]],GPR_list$new_country[[region]],k_gamma[[region]])
        
        #source and covar 
        if(grepl("Source", model_name)){
          source_eta_S[[region]] <- as.numeric(ext_fit[draw_idx,paste0("source_eta_S[",1:num_source_train,"]")])
        }else{}
        
        if(grepl("Covar", model_name)){
          eta[[region]] <- as.numeric(ext_fit[draw_idx,paste0("eta[",1:ncol(cov_mat),"]")])
        }else{}
        
      }
      
    }else{
      I0_logit <- ext_fit$I0_logit[draw_idx] 
      I0_prob <- ext_fit$I0_prob[draw_idx] 
      log_beta <- ext_fit$log_beta[draw_idx] 
      log_gamma <- ext_fit$log_gamma[draw_idx] 
    }
    
    # Initialize tracking variables
    lastCountry <- -1
    lastYear <- -9999
    lastI <- -9999
    
    #predictions for test set - has source level but get single curve for each country - use the other reordered data
    for(i in 1:n_pred) {
      
      
      country_i <- data$country[i]
      source_i <- data$source[i]
      t_i <- data$t[i]

      country_track_i <- data$country_track[i]
      source_track_i <- data$source_track[i]
      
      initial_year_track <- data$init_year_country[i]
      
      neigh_1 <- data$neigh_1[i]
      neigh_2 <- data$neigh_2[i]
      
      #track regions for spatial model 
      region_i <- data$owid_region_numeric[i]
      country_spatial_i <- data$countryNumSpatial[i]
      
      
      if(country_i != lastCountry) {
        if(country_track_i != -1000){
          
          if(grepl("spatial", model_name) ){
            
            y_pred[s, i] <- inv.logit(logit(I0_prob[region_i]) + k_I0[[region_i]][country_spatial_i])
            
          }else{
            
            y_pred[s, i] <- I0_prob_c[country_track_i] #informative prior
          }
          if(is.na( y_pred[s, i])){
            print(i)
            fail1
          }
          
        }else{
          
          if(grepl("spatial", model_name)){
            #test with GPR - Spatial region -informative
            y_pred[s, i] <- inv.logit(logit(I0_prob[region_i]) + k_I0_new[[region_i]][country_ID == country_i]$value)
            
            
            
          }else{ 
            
           #placeholder to test GPR not called 
            y_pred[s, i] <- (I0_prob_c[neigh_1] + I0_prob_c[neigh_2])/2
          }
          if(is.na( y_pred[s, i])){
            print(i)
            fail2
          }
        }
        
      } else {
        
        if(abs(t_i - lastYear) < 10^(-4)) {
          y_pred[s, i] <- lastI
        }else{
          
          if(country_track_i != -1000){
            
            if(grepl("spatial", model_name)){
              
              #adding condition for now beta with covariates
              if(grepl("spatialCmdStanSourceCovar", model_name)){
                current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i] + sum(cov_mat[i, ] * eta[[region_i]]))
              }else{
                current_beta <- exp(log_beta[region_i] + k_beta[[region_i]][country_spatial_i])
              }
              
              #test with the spatial regional parameters
              current_gamma <- exp(log_gamma[region_i] + k_gamma[[region_i]][country_spatial_i])#uses same neigh struc as beta 
              
             
            }else if(grepl("covar", model_name)){
              current_beta <- exp(log_beta_c[country_track_i] + sum(cov_mat[i, ] * eta))
              current_gamma <- gamma_c[country_track_i]
            }else{
              current_beta <- beta_c[country_track_i]
              current_gamma <- gamma_c[country_track_i]
            }
           
            
            
          }else{
            if(grepl("spatial", model_name)){

              #adding condition for now beta with covariates
              if(grepl("spatialCmdStanSourceCovar", model_name)){
                current_beta <- exp(log_beta[region_i] +  k_beta_new[[region_i]][country_ID == country_i]$value + sum(cov_mat[i, ] * eta[[region_i]]))
              }else{
                current_beta <- exp(log_beta[region_i] +  k_beta_new[[region_i]][country_ID == country_i]$value)
              }
              
              #test with the spatial regional parameters
              current_gamma <- exp(log_gamma[region_i] + k_gamma_new[[region_i]][country_ID == country_i]$value)#uses same neigh struc as beta 
              
              
            }else if(grepl("covar", model_name)){
              
              current_beta <- exp((log_beta_c[neigh_1] + log_beta_c[neigh_2])/2 + sum(cov_mat[i, ] * eta))
              current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
            }else{
              current_gamma <- (gamma_c[neigh_1] + gamma_c[neigh_2])/2
              current_beta <- (beta_c[neigh_1] + beta_c[neigh_2])/2
            }
           
          }
          
          
          if(abs(current_beta - current_gamma) < 10^(-4)){
            y_pred[s, i] <- lastI
          }else{
            
           
            y_pred[s, i] <- ((current_beta - current_gamma) / current_beta) /(1 + (((current_beta - current_gamma) / current_beta) / lastI - 1) *
                 exp(-(current_beta - current_gamma)*(t_i - lastYear)))
            
           
          }}
        
        
        
        if(is.na( y_pred[s, i])){
          print(i)
          fail
        }
      }
     
      
      if(y_pred[s, i] < 0 ){
        
        fail
      }
      
      lastI <- y_pred[s, i]
      lastCountry <- country_i
      lastYear <- t_i
      
     
      if(grepl("Source", model_name)){
        if(country_track_i != -1000){
          if(grepl("spatial", model_name)){
            current_source <- source_eta_S[[region_i]][source_track_i]
          }else{current_source <- source_eta_S[source_track_i]}
        }else{
          if(grepl("spatial", model_name)){ #must check if new country has same or new source
            if(source_track_i != -1000)
              current_source <- source_eta_S[[region_i]][source_track_i]
            else{
              current_source <- mean(source_eta_S[[region_i]])
            }
            
          }else{current_source <- mean(source_eta_S)}
        }
        y_pred[s, i] <- (y_pred[s, i])*current_source
      }else{}
      
      y_pred[s, i] <- max(0.00001, min(0.99999, y_pred[s, i] ))
    }
  }

  N_vec <- data$cases
  country_median_N <- tapply(data$cases, data$country_name, median, na.rm = TRUE)
  N_vec[is.na(N_vec) | N_vec == 0] <- country_median_N[data$country_name[is.na(N_vec) | N_vec == 0]]
  N_vec[is.na(N_vec)] <- median(data$cases, na.rm = TRUE)
  
 
  y_pred_obs <- matrix(NA, nrow = nrow(y_pred), ncol = ncol(y_pred))
  for(s in 1:nrow(y_pred)){
    y_pred_obs[s, ] <- rbinom(n = ncol(y_pred),  size = round(N_vec),  prob = y_pred[s, ]) / round(N_vec)
  }
  
  predictions <- data.frame(
    median = apply(y_pred, 2, median),
    q5 = apply(y_pred, 2, quantile, 0.05),
    q95 = apply(y_pred, 2, quantile, 0.95),
    mean = apply(y_pred, 2, mean),
    sd = apply(y_pred, 2, sd),
    obs_mean = apply(y_pred_obs, 2, mean),
    obs_q5 = apply(y_pred_obs, 2, quantile, 0.05),
    obs_q95 = apply(y_pred_obs, 2, quantile, 0.95))
  
  
  return(list(predictions = predictions))  
}


#function to take in the old country, new country, rhosq, etasq, and old_country_GP_values, and distance matrix
#to get new country estimates

calc_cov_mat <- function(dist_mat, etasq, rhosq, delta, old_country, new_country,old_country_GP_values){

 #get number of old and new countries
  N_old <- length(old_country)
  N_new <- length(new_country)
  
  K_old <- matrix(NA,nrow = N_old, ncol = N_old)
  K_new <- matrix(NA,nrow = N_new, ncol =  N_new)
  K_old_new <- matrix(NA,nrow = N_old, ncol =  N_new)
  
  #K normal with old countries or countries found in the train set same in Stan model file 
  for (i in 1:(N_old-1)) {
    K_old[i, i] <- etasq + delta
    for (j in (i + 1):N_old) {
      K_old[i, j] <- etasq * exp(-rhosq * (dist_mat[old_country[i],old_country[j]])^2)
      K_old[j, i] <- K_old[i, j]
    }
  }
  K_old[N_old, N_old] = etasq + delta;
  
  #K** with only new countries in the test set   
  if(N_new == 1){
    K_new[N_new, N_new] <- etasq + delta;
  }else{
    for (i in 1:(N_new-1)) {
      K_new[i, i] <- etasq + delta
      for (j in (i + 1):N_new) {
        K_new[i, j] <- etasq * exp(-rhosq * (dist_mat[new_country[i],new_country[j]])^2)
        K_new[j, i] <- K_new[i, j]
      }
    }
    K_new[N_new, N_new] <- etasq + delta;
  }
  #K* with only new countries and old countries -not mirroring so don't need previous type of indexing
  for (i in 1:N_old) {
    for (j in 1:N_new) {
      #print((dist_mat[old_country[i],new_country[j]])^2)
      K_old_new[i, j] <- etasq * exp(-rhosq * (dist_mat[old_country[i],new_country[j]])^2)
    }
  }
  
  #now we calculate the conditional distribution
  mu_star <- t(K_old_new) %*% solve(K_old) %*% old_country_GP_values
  sigma_star <- K_new - t(K_old_new) %*% solve(K_old) %*% K_old_new
  sigma_star <- sigma_star + diag(0.02, N_new)
  
  #new values for one parameter at a time for each sample
  new_country_GP_values <- mvrnorm(1, mu = mu_star, Sigma = sigma_star)
  new_GP <- data.table(value = new_country_GP_values, country_ID = new_country)
 
  return(new_GP)
 
}


