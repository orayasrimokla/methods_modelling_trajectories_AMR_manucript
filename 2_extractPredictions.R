#after fitting the models we need to extract the predictions 

#for local 
path_cluster = "yourPath"
data_run <- "modelRunExamples/" 
folder_name_save <- "all_models/"
#create new folder to 
dir.create(paste0(path_cluster,data_run,folder_name_save))

#load the functions and packages
source(paste0(path_cluster,"load_packages.R"))
source(paste0(path_cluster,"load_functions.R"))
source(paste0(path_cluster,"generatePredictions_functions.R")) #script with functions for generating predictions 
options(mc.cores = parallel::detectCores()) #stan model will still run the number of cores equal to chain 
#so for our code we have 3 chains so it will run at most 3 cores 

data_grouped_covar_source <- data.frame(read.csv(file=paste0("pathTo/dataFile.csv")))

#this should be the same for both data aggregates 
country_id <- sort(unique(data_grouped_covar_source$location_id))

#uncomment this line and comment the line below to just test with this for location id less than 200 for example
#data_test_source <- data_grouped_covar_source[data_grouped_covar_source$location_id < 200,]

#all countries 
data_test_source <- data_grouped_covar_source

#selecting baseline year to be 2000 -  filter per region -  lots of ways to define what countries are in what region
data_test_source <- data_grouped_covar_source[data_grouped_covar_source$year_id >= 2000 & data_grouped_covar_source$owid_region == "Asia" ,]

#get the data that is not NA
data_test_source_train <- data_test_source[!is.na(data_test_source$cases_sum),]

list_models_run <- c("betaGammaI0_country","betaGammaI0_covar",
                     "betaGammaI0_AddSource","spatial_cmdStan","spatialcmdStanSource",
                     "spatialcmdStanSourceCovar" )


#remove the [1] to run all abx - right now just set to first one 
list_abx <- sort(unique(data_test_source$abx_class))[1]


#test with just the abx consumption for that class - for A baum - select the covariate you want to test
covar_list <- c("ddd_aminoglycoside","ddd_per_1000_fitted_ihme","ddd_betaLactam_penicilins","ddd_quinolone","ddd_othrBetaLactam")

#can index of model we want
modelNum <- c(1,6)
model_label = c("Country RE",
                "Country RE+Covar",
                "Country RE+Source",
                "Spatial Effects",
                "Spatial Effects+Source",
                "Spatial Effects+Source+Covar"
)
model_label <- model_label[modelNum[1]:modelNum[2]]


#this is an example loop - the logic is to load the individual models after fitting and then extract predictions 
#this is just adding a few lines to keep track of the neighbours for the predictions 
#in this prediction extraction - simplify extractions for countries with existing resistance data used in train -> another script will go 
#over GPR for countries without resistance data used to train model
llply(1:length(list_abx), function(i){
    
    data_country_source_train <- prep_data_abxFilt(data_test_source_train, source = 1, abx_name = list_abx[i], calc_lik = 1, list_models_run[j])
    data_country_source_train$cov_mat <- do.call(cbind, data_country_source_train[paste0(covar_list)])
    data_country_source_train$K <- length(covar_list)
    
    data_model <- data_country_source_train
    
    #adding tests for running models 
    #get full dt across all years 
    data_country_source <- prep_data_abxFilt(data_test_source, source = 1, abx_name = list_abx[i], calc_lik = 1, list_models_run[j])
    data_country_source$cov_mat <- do.call(cbind, data_country_source[paste0(covar_list)])
    data_country_source$K <- length(covar_list)
    
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
                        countryNumSpatial = data_model$countryNumSpatial,
                        cases = data_model$cases, neigh_1 = NA,neigh_2 = NA,initial_year_track = NA, init_year_country = -9999, # todo: maybe fix or maybe it's fine
                        data_model$cov_mat)
    train[, init_year_country := min(t), by = list(country)]
    
    
    #we just want rows with no prop resis_reads? so make data frame with all the countrya nd source combinations for all time
    #test has all country and combo
    test<- data.table(t = data_country_source$t,country = data_country_source$country, resis_prop = data_country_source$resis_prop,
                      country_name = data_country_source$country_name,year_id = data_country_source$year_id,
                      source_track = data_country_source$source,
                      source_name = data_country_source$source_name, country_track = data_country_source$country_track,
                      cases = data_country_source$cases,neigh_1 = data_country_source$neigh_1,neigh_2 = data_country_source$neigh_2,
                      countryNumSpatial = data_country_source$countryNumSpatial,
                      initial_year_track = data_country_source$initial_year_track,init_year_country = -9999,data_country_source$cov_mat)
    
    test <- unique(test, by = c("country","t"))
    test[, init_year_country := train$init_year_country[match( test$country_name,train$country_name)]]
    
    loos <- list()
    #this runs all models but you can specify j = x to set it to get predictions for a single model for this abx
    #the current version of this loop generates predictions for countries that have observational data in the model fitting and 
    #is filtered by the country_track != -1000 line.
    preds <- rbindlist(llply(1:length(list_models_run), function(j){
      
      #need to load the spatial models different since it was fit with cmdstanr
      if(grepl("spatial", list_models_run[j])){
        fit <- readRDS(file = paste0(path_cluster,data_run,list_models_run[j],"/",'model: ',list_models_run[j],'+2000:2022',list_abx[i],".rds"))
        model <- fit$draws(format = "df")
        loos[[j]] <<- fit$loo()
      }else{
        load(paste0(path_cluster,data_run,list_models_run[j],"/",'model: ',list_models_run[j],'+2000:2022',list_abx[i]))
        loos[[j]] <<- extract_log_lik(model) %>% loo()
          }
    model
    
    #For source models we need to generate predictions depending on the source information and then rejoin them later
    #this includes the country RE + source, spatial + source, and spatial + source + covar
    if(grepl("Source", list_models_run[j])){
      #testing - these ARE ALL HACKS NEED TO DELETE LATER
      
      #TESTING COUNTRIES WITH EXISTING RESISTANCE DATA FOR NOW ANOTHER SCRIPT WILL GO OVER PREDICTIONS FOR COUNTRIES
      #WITH WHOLE COUNTRY MISSING DATA
      sourceID <- unique(data.frame(train$source_name,train$source_track))
      neigh_ID <- unique(test[,c("country_name","neigh_1","neigh_2")])
      colnames(sourceID) <- c("source_name","source_track")

      #get sources specific for each country and then expand to have the years we want for our analysis
      sourceTestData <- train %>%
        dplyr::group_by(country_name) %>%
        expand(year_id = 2000:2022, source_name)
      
      sourceTestData$source_track <- sourceID$source_track[match( sourceTestData$source_name,sourceID$source_name)]
      sourceTest <- merge(sourceTestData,train[, c("country_name","year_id","source_name","resis_prop","cases")], by = c("country_name","year_id","source_name"),all.x = TRUE)
      
      sourceTest$neigh_1 <- neigh_ID$neigh_1[match( sourceTest$country_name,neigh_ID$country_name)]
      sourceTest$neigh_2 <- neigh_ID$neigh_2[match(sourceTest$country_name,neigh_ID$country_name)]
      sourceTest$country_track <- train$country_track[match(sourceTest$country_name,train$country_name)]
      sourceTest$countryNumSpatial <- train$countryNumSpatial[match(sourceTest$country_name,train$country_name)]
      
      sourceTest$t <- sourceTest$year_id -2000 #1990
      sourceTest$country <- sourceTest$country_track
      sourceTest$source <- sourceTest$source_track
      
      #covar
      col_select <- c("country_name","year_id", covar_list)
      sourceTest <- merge(sourceTest,test[, ..col_select], by = c("country_name","year_id"),all.x = TRUE)
      
      #sourceTest <- sourceTest[order(sourceTest$country_name,sourceTest$source,sourceTest$t), ]
      sourceTest <- sourceTest[order(sourceTest$country_name,sourceTest$t), ] #need to order by year annd country
      
      #gives the predictions for each source in that country
      predictions_forward <- run_predictions(model,sourceTest, list_models_run[j])
      predictions <-  rbind(cbind(predictions_forward$predictions,sourceTest[,c("t","country_name","resis_prop","year_id","source_name","cases")])
      )
      #averaging the predictions across the sources
      predDT <- data.table(predictions)
      predDT <- predDT[ , lapply(.SD, mean) , by=c("t", "country_name"),  .SDcols = sapply(predDT, is.numeric) ][,-1]
      predictions <- merge(predDT,train, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
      colnames(predictions)[colnames(predictions) == "source_name"] <- "source_name.y"
    }else if(grepl("cmdStan", list_models_run[j])){
      #another way to get the model predictions but predictions generated in R instead of Stan forward predictions
      #TESTING COUNTRIES WITH EXISTING RESISTANCE DATA FOR NOW ANOTHER SCRIPT WILL GO OVER PREDICTIONS FOR COUNTRIES
      #WITH NO RESISTANCE DATA
      data_forward <- test[test$country_track != -1000, ]#test[order(test$country_name,test$t), ]
      
      data_forward <- data_forward[order(data_forward$country_name,data_forward$t), ]
      predictions_forward <- run_predictions(model,data_forward, list_models_run[j])
      predictions <-  rbind(cbind(predictions_forward$predictions,data_forward[,c("t","country_name","resis_prop","year_id","source_name","cases")])
      )
      predictions <- merge(predictions,train, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
    }else{
      
      ##TESTING COUNTRIES WITH EXISTING RESISTANCE DATA FOR NOW ANOTHER SCRIPT WILL GO OVER PREDICTIONS FOR COUNTRIES
      #WITH NO RESISTANCE DATA
      data_forward <- test[test$country_track != -1000, ]#test[order(test$country_name,test$t), ]
      
      #another way to get the model predictions but predictions generated in R instead of Stan forward predictions
      data_forward <- data_forward[order(data_forward$country_name,data_forward$t), ]
      predictions_forward <- run_predictions(model,data_forward, list_models_run[j])
      predictions <-  rbind(cbind(predictions_forward$predictions,data_forward[,c("t","country_name","resis_prop","year_id","source_name","cases")])
      )
      predictions <- merge(predictions,train, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
    }
    
    predictions$ModelType = model_label[j]
    return(data.table(predictions))
    
    }), fill = TRUE)
    
    names(loos) = list_models_run
    #save some things/ get model weights
    if(length(list_models_run) > 1 ){
      weights_abx_model <- loo_model_weights(loos)
      write.csv(weights_abx_model,paste0(path_cluster,data_run,folder_name_save,'model_comparison_WEIGHTS_abx: ',list_abx[i],'.csv'), row.names = TRUE)
    }else{}
    
    preds
    #save the predictions
    write.csv(preds,paste0(path_cluster,data_run,folder_name_save,'preds: ',list_abx[i],'.csv'), row.names = FALSE)
    
    #can add custom plotting code to save plots based on SOP
    p <- ggplot(preds, mapping = aes(x = t + 2000)) + geom_ribbon(aes(ymin = q5 , ymax = q95, fill = country_name), alpha = 0.35) + 
      geom_line(mapping = aes(x = t + 2000, y = mean), colour = "#808080") +  geom_point(mapping = aes(y = resis_prop.y, size = cases.y, colour = factor(source_name.y)))  + 
      #scale_fill_manual(name = "Source",values = colorRampPalette(brewer.pal(8, "BuPu"))(length(unique(data_model$source_name)))) + 
      scale_colour_manual(name = "Model Type", values=colorRampPalette(brewer.pal(8, "PRGn"))(length(unique(predictions$source_name.y)))) +
      scale_size_area(max_size = 5, 
                      #breaks = c(50, 100, 500 , 1000, 2000),
                      #breaks = c(1, 30, 50, 80, 100, 150),
                      labels = ~paste0(scales::comma(x = .x)), 
                      name = "Total Observations (R+S)") + 
      labs(title = paste0("Model Fits for  ",list_abx[i]),x = "Year", y = "Proportion Resistant",fill="Source + Model CI", colour="Model") + ylim(0,1) + xlim(2000, 2022) + theme_minimal() + 
      theme(legend.title = element_text(size=3),legend.text = element_text(size = 3),legend.key.size = unit(0.2, "cm"), legend.position = "right") + 
      facet_wrap_paginate(~country_name, ncol = 4, nrow=4) 
    
    p
    for(k in 1:n_pages(p)){
      p_save <-  p +
        facet_wrap_paginate(~country_name, ncol = 4, nrow = 4, page = k, scales = "free")
      ggsave(paste0(path_cluster,data_run,folder_name_save, k, 'abx: ',list_abx[i],'_individual.pdf'), p_save, width = 15, height = 8, device = "pdf") #width = 10 , height = 6
    }
  
})



