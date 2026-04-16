#generate estimates for countries with whole missing data using Gaussian processes from the spatial models

#compare for whole countries
path_cluster = "yourPath"
data_run <- "modelRunExamples/" 
folder_name_save <- "all_models/"
dir.create(paste0(path_cluster,data_run))
dir.create(paste0(path_cluster,data_run,folder_name_save))

#load the functions and packages
source(paste0(path_cluster,"amr_testing/load_packages.R"))
source(paste0(path_cluster,"amr_testing/load_functions.R"))
source(paste0(path_cluster,"amr_testing/generatePredictions_functions.R"))
options(mc.cores = parallel::detectCores())

data_grouped_covar_source <-data.frame(read.csv(file=paste0(path_cluster,"data_grouped_covar_source.csv")))

#this should be the same for both data aggregates 
country_id <- sort(unique(data_grouped_covar_source$location_id))

#all countries 
data_test_source <- data_grouped_covar_source

data_test_source <- data_grouped_covar_source[data_grouped_covar_source$year_id >= 2000,] %>% filter(owid_region %in% c(paste0("Asia")))

data_test_source_train <- data_test_source[!is.na(data_test_source$cases_sum),]

list_models_run <- c("spatial_cmdStan","spatialcmdStanSource","spatialcmdStanSourceCovar")

#remove the [1] to run all abx - right now just set to first one 
list_abx <- sort(unique(data_test_source$abx_class))[1]

#test with just the abx consumption for that class - for A baum - select the covariate you want to test
covar_list <- c("ddd_aminoglycoside","ddd_per_1000_fitted_ihme","ddd_betaLactam_penicilins","ddd_quinolone","ddd_othrBetaLactam")


model_label = c(
  "Spatial Effects",
  "Spatial Effects + Source",
  "Spatial Effects + Source + Covar") #just spatial

#now separate loop to run for countries with missing data for GPR
llply(1:length(list_abx), function(i){
  
  #get full dt across all years 
  data_country_source <- prep_data_abxFilt(data_test_source, source = 1, abx_name = list_abx[i], calc_lik = 1, model_name = "all_country_sorting")
  data_country_source$cov_mat <- do.call(cbind, data_country_source[paste0(covar_list)])
  data_country_source$K <- length(covar_list)
  
  data_country_source_train <- prep_data_abxFilt(data_test_source_train, source = 1, abx_name = list_abx[i], calc_lik = 1, model_name = "all_country_sorting")
  data_country_source_train$cov_mat <- do.call(cbind, data_country_source_train[paste0(covar_list)])
  data_country_source_train$K <- length(covar_list)
  
  data_model <- data_country_source_train
  
  #add the source tracking and  initial country tracking to match the train set indexing 
  data_country_source$source_track <- data_country_source_train$source[match( data_country_source$source_name,data_country_source_train$source_name)]
  data_country_source$source_track[is.na(data_country_source$source_track)] <- -1000
  data_country_source$country_track <- data_country_source_train$country[match( data_country_source$country_name,data_country_source_train$country_name)]
  data_country_source$country_track[is.na(data_country_source$country_track)] <- -1000
  
  #table of unique country IDs and names 
  df_country_all <- unique(data.frame(data_country_source$country, data_country_source$country_name))
  df_country_train <- unique(data.frame(data_country_source_train$country, data_country_source_train$country_name))
  
  #get neighbours 
  neigh_1 <- c()
  neigh_2 <- c()
  #make list for for neigh 1 and 2
  for(m in 1:nrow(df_country_all)){
    
    rep_neigh <- as.numeric(data.table(table(data_country_source$country))[m,2])
    
    #now what is country ID of x_name in the all country for example 
    all_country_ID <- df_country_all$data_country_source.country[m]
    
    distance_mat_country_all_ID <- data_country_source$D[all_country_ID,]
    
    list_neigh_distance <- c()
    for(k in 2:length(distance_mat_country_all_ID)){
      #get the country ID for the ones with the closest distances to the one we are looking at ( get three of them)
      index_low_2 <- order(distance_mat_country_all_ID)[k]
      
      #let's check what the names of these countries are that are close to our country of interest
      name_low_2 <- df_country_all$data_country_source.country_name[df_country_all$data_country_source.country== index_low_2]
      #check if that country is in the train set. If it is then append
      if(sum(df_country_train$data_country_source_train.country_name == paste0(name_low_2)) != 0){
        list_neigh_distance <- append(list_neigh_distance,name_low_2)
      }else{}}
    
    #now get the country ID of the neighbours from the train set
    train_country_ID_connect1 <- df_country_train$data_country_source_train.country[df_country_train$data_country_source_train.country_name == paste0(list_neigh_distance[1])]
    train_country_ID_connect2 <- df_country_train$data_country_source_train.country[df_country_train$data_country_source_train.country_name ==  paste0(list_neigh_distance[2])]
    
    neigh_1 <- append(neigh_1,rep(train_country_ID_connect1,rep_neigh))
    neigh_2 <- append(neigh_2,rep(train_country_ID_connect2, rep_neigh))
    
  }
  data_country_source$neigh_1 <- neigh_1
  data_country_source$neigh_2 <- neigh_2
  
 
  #redundant code from before when separating by region - keeping for now since predictions use these values 
  data_country_source$countryNumSpatial <- data_country_source$country_track
  data_model$countryNumSpatial <- data_model$country
  
  #prepare the data for plotting
  train <- data.table(country = data_model$country,country_track = data_model$country, t = data_model$t, 
                      source_track = data_model$source,
                      country_name = data_model$country_name,resis_prop = data_model$resis_prop,year_id = data_model$year_id,source_name = data_model$source_name,
                      owid_region = data_model$owid_region, owid_region_numeric = data_model$owid_region_numeric, countryNumSpatial = data_model$countryNumSpatial,
                      cases = data_model$cases, neigh_1 = NA,neigh_2 = NA,initial_year_track = NA, init_year_country = -9999, # todo: maybe fix or maybe it's fine
                      data_model$cov_mat)
  train[, init_year_country := min(t), by = list(country)]
  
  #test has all country and combo
  test<- data.table(t = data_country_source$t,country = data_country_source$country, resis_prop = data_country_source$resis_prop,
                    country_name = data_country_source$country_name,year_id = data_country_source$year_id,
                    source_track = data_country_source$source_track,
                    source_name = data_country_source$source_name, country_track = data_country_source$country_track,
                    owid_region = data_country_source$owid_region, owid_region_numeric = data_country_source$owid_region_numeric, countryNumSpatial = data_country_source$countryNumSpatial,
                    cases = data_country_source$cases,neigh_1 = data_country_source$neigh_1,neigh_2 = data_country_source$neigh_2,
                    initial_year_track = data_country_source$initial_year_track,init_year_country = -9999,data_country_source$cov_mat)
  
  test <- unique(test, by = c("country","t"))
  test[, init_year_country := train$init_year_country[match( test$country_name,train$country_name)]]
  
  
  preds <- rbindlist(llply(1:length(list_models_run), function(j){
    print(j)
  
    if(grepl("spatial", list_models_run[j])){
    
      fit <- readRDS(file = paste0(path_cluster,data_run,list_models_run[j],"/",'model: ',list_models_run[j],'+2000:2022',list_abx[i],".rds"))
      model <- fit$draws(format = "df")
      cmdStan_models_region <- fit
    }else{
      load(paste0(path_cluster,data_run,list_models_run[j],"/",'model: ',list_models_run[j],'+2000:2022',list_abx[i]))
    }
    model
    
    #merge in test data or could just use test
    test_temp <- test[test$country_track == -1000,] #test_data
    test_temp$t <- test_temp$year_id - 2000
    colnames(test_temp)[colnames(test_temp) =="COUNTRY"] <- "country_name"
    colnames(test_temp)[colnames(test_temp) =="source"] <- "source_name"
    colnames(test_temp)[colnames(test_temp) =="cases_sum"] <- "cases"
    colnames(test_temp)[colnames(test_temp) =="resistance_prop"] <- "resis_prop"
    test_temp$country_track <- -1000
    
    combDat <- rbind(train, test_temp, fill = TRUE)
    
    if(grepl("Source", list_models_run[j])){

      new_countries <- test[test$country_name %in% unique(test_data$COUNTRY)]#test[country_track ==-1000]
      sourceID <- unique(data.frame(train$source_name,train$source_track))
      neigh_ID <- unique(test[,c("country_name","neigh_1","neigh_2")])
      colnames(sourceID) <- c("source_name","source_track")
      
      sourceTestData <- combDat %>%
        dplyr::group_by(country_name) %>%
        expand(year_id = 2000:2022, source_name)
      
      #need to match the country ID from old  and new based on the large test set with all countries since the distanace is organised that way TO DO 
      
      sourceTestData$source_track <- sourceID$source_track[match( sourceTestData$source_name,sourceID$source_name)]
      sourceTest <- merge(sourceTestData,combDat[, c("country_name","year_id","source_name","resis_prop","cases")], by = c("country_name","year_id","source_name"),all.x = TRUE)
      
      sourceTest$neigh_1 <- neigh_ID$neigh_1[match( sourceTest$country_name,neigh_ID$country_name)]
      sourceTest$neigh_2 <- neigh_ID$neigh_2[match(sourceTest$country_name,neigh_ID$country_name)]
      
      sourceTest$country_track <- train$country_track[match(sourceTest$country_name,train$country_name)]
      sourceTest$country_track[is.na(sourceTest$country_track)] <- -1000
      sourceTest$owid_region <- test$owid_region[match(sourceTest$country_name,test$country_name)]
      sourceTest$countryNumSpatial <- train$countryNumSpatial[match(sourceTest$country_name,train$country_name)]
      sourceTest$owid_region_numeric <- test$owid_region_numeric[match(sourceTest$country_name,test$country_name)]
      
      sourceTest$t <- sourceTest$year_id -2000 
      sourceTest$country <- test$country[match(sourceTest$country_name,test$country_name)]    #sourceTest$country_track #need to change to match test 
      sourceTest$source <- sourceTest$source_track
      sourceTest$source_track[is.na(sourceTest$source_track)] <- -1000
      
      #covar
      col_select <- c("country_name","year_id", covar_list)
      sourceTest <- merge(sourceTest,test[, ..col_select], by = c("country_name","year_id"),all.x = TRUE)
      
      
      sourceTest <- sourceTest[order(sourceTest$country_name,sourceTest$t), ] 
      
      #get region information
      regionOrder <- sort(unique(sourceTest$owid_region))
      
      #separate countries by new and old per region in list based on region( nested list)
      oldTest <- c()
      newTest <- c()
      for(region in 1:length(regionOrder)){
        #old countries in train 
        oldTest[[region]] <- unique(sourceTest[sourceTest$country_track != -1000 & sourceTest$owid_region == regionOrder[region],"country"])
        #newcountries 
        newTest[[region]] <- unique(sourceTest[sourceTest$country_track == -1000 & sourceTest$owid_region == regionOrder[region],"country"])
      }
      
      GPR_list <- list(oldTest,newTest,data_country_source$D/(6378100*pi))
     
      
      names(GPR_list) <- c("old_country","new_country","dist_mat")
      
      #data with old and new countries
      if(grepl("CmdStanSource", list_models_run[j])){
        modelSource <- cmdStan_models_region
      }else{modelSource <- model}
      
      #gives the predictions for each source in that country
      predictions_forward <- run_predictions_test_GPR(modelSource,sourceTest, list_models_run[j], GPR_list)
      predictions <-  rbind(cbind(predictions_forward$predictions,sourceTest[,c("t","country_name","resis_prop","year_id","source_name","cases")])
      )
      #averaging the predictions across the sources for plotting
      predDT <- data.table(predictions)
      predDT <- predDT[ , lapply(.SD, mean) , by=c("t", "country_name"),  .SDcols = sapply(predDT, is.numeric) ][,-1]
      predictions <- merge(predDT,combDat, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
      colnames(predictions)[colnames(predictions) == "source_name"] <- "source_name.y"
    }else{
    
      data_forward <- test[order(test$country_name,test$t), ]
      
      #get region information
      regionOrder <- sort(unique(data_forward$owid_region))
      
      #TO DO: need to separate countries by new and old per region in list based on region( nested list)
      oldTest <- c()
      newTest <- c()
      for(region in 1:length(regionOrder)){
        #old countries in train 
        oldTest[[region]] <- unique(data_forward[data_forward$country_track != -1000 & data_forward$owid_region == regionOrder[region],country])
        #newcountries 
        newTest[[region]] <- unique(data_forward[data_forward$country_track == -1000 & data_forward$owid_region == regionOrder[region],country])
      }
      
      GPR_list <- list(oldTest,newTest,data_country_source$D/(6378100*pi))
      #GPR_list <- list(trainCountryID_inTest,newCountryID,data_country_source$D)
      
      names(GPR_list) <- c("old_country","new_country","dist_mat")
      predictions_forward <- run_predictions_test_GPR(cmdStan_models_region,data_forward, list_models_run[j], GPR_list)

      predictions <-  rbind(cbind(predictions_forward$predictions,data_forward[,c("t","country_name","resis_prop","year_id","source_name","cases")])
      )
      
      
      predictions <- merge(predictions,combDat, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
      
    }
    predictions$ModelType = model_label[j]
    
    # #save some diagnostic plots
    # a <- mcmc_pairs(model, pars = c("beta","gamma"))
    # b <- mcmc_trace(model, pars = c("beta","gamma"))
    return(data.table(predictions))
    
  }), fill = TRUE)
  
  preds
  write.csv(preds,paste0(path_cluster,data_run,folder_name_save,'preds: ',list_abx[i],'GPRPred.csv'), row.names = FALSE)
  
  
  model_colors <- viridis(length(unique(preds$ModelType)),begin = 0.15, end = 0.95,option = "D")
  source_colors <- mako(length(unique(na.omit(preds$source_name.y))),begin = 0.15, end = 0.85)
  names(source_colors) <- unique(na.omit(preds$source_name.y))
  names(model_colors) <- unique(preds$ModelType)
  
  p <- ggplot(preds, mapping = aes(x = t + 2000)) + geom_ribbon(aes(ymin = q5 , ymax = q95, fill = ModelType), alpha = 0.20) +
    geom_point(mapping = aes(y = resis_prop.y, size = cases.y, fill = factor(source_name.y),shape=factor(source_name.y)),colour = "grey",stroke = 0.2)  +
    geom_line(mapping = aes(x = t + 2000, y = mean, colour = ModelType), size = 0.7) +
    scale_shape_manual(name = "Data Sources",values = rep(c(21:25), length.out = length(source_colors)), breaks = c(names(source_colors))) +
    scale_fill_manual(name = "Model CrI",values = c(model_colors, source_colors), breaks = c(names(model_colors))) +
    scale_colour_manual(name = "Model Type", values=model_colors) +
    scale_size_continuous(range = c(1.5,8),
                          breaks = c(100, 500 , 5000, 10000,20000),
                          #breaks = c(1, 30, 50, 80, 100, 150),
                          labels = ~paste0(scales::comma(x = .x)),
                          name = "Total Observations (R+S)") +
    labs(x = "Year", y = "Proportion Resistant",fill="Source + Model CrI", colour="Model") + ylim(0,1) + xlim(2000, 2022) + theme_minimal() +
    theme(legend.title = element_text(size=8, face = "bold"),legend.text = element_text(size = 7),legend.key.size = unit(0.6, "cm"), legend.position = "right",
          strip.text = element_text(size = 8,face = "bold"), axis.title = element_text(size = 13),axis.text = element_text(size = 10),axis.text.x = element_text(angle = -45) ,strip.background = element_rect(colour = "gray96", fill = "gray96"))  +
    facet_wrap_paginate(~ country_name , ncol = 4, nrow=4) +  guides(fill = guide_legend(ncol = 2),shape = guide_legend(ncol = 2, override.aes = list(fill = source_colors, size = 3)))
  
  p
  
  for(k in 1:n_pages(p)){
    
    p_save <-  p +
      facet_wrap_paginate(~country_name, ncol = 4, nrow = 4, page = k, scales = "free")
    ggsave(paste0(path_cluster,data_run,folder_name_save, k, 'abx: ',list_abx[i],'_individual.pdf'), p_save, width = 15, height = 8, device = "pdf") 

  }
  
}) 
