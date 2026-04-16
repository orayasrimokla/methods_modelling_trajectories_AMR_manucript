
#after extracting the individual predictions, we can stack the models

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

data_grouped_covar_source <- data.frame(read.csv(file=paste0("/Users/orayasrim/Documents/AMR_AB/data_grouped_covar_source.csv")))

#this should be the same for both data aggregates 
country_id <- sort(unique(data_grouped_covar_source$location_id))

#uncomment this line and comment the line below to just test with this for location id less than 200 for example
#data_test_source <- data_grouped_covar_source[data_grouped_covar_source$location_id < 200,]

#all countries 
data_test_source <- data_grouped_covar_source

#selecting baseline year to be 2000 -  filter per region -  lots of ways to define what countries are in what region
data_test_source <- data_grouped_covar_source[data_grouped_covar_source$year_id >= 2000 & data_grouped_covar_source$owid_region == "Asia" ,]

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
  data_country_source_train <- prep_data_abxFilt(data_test_source_train, source = 1, abx_name = list_abx[i], calc_lik = 1, "all_country_sorting")
  data_country_source_train$cov_mat <- do.call(cbind, data_country_source_train[paste0(covar_list)])
  data_country_source_train$K <- length(covar_list)
  
  data_model <- data_country_source_train
  
  #adding tests for running models 
  #get full dt across all years 
  data_country_source <- prep_data_abxFilt(data_test_source, source = 1, abx_name = list_abx[i], calc_lik = 1, "all_country_sorting")
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
  
  
 #test has all country and combo
  test<- data.table(t = data_country_source$t,country = data_country_source$country, resis_prop = data_country_source$resis_prop,
                    country_name = data_country_source$country_name,year_id = data_country_source$year_id,
                    source_track = data_country_source$source_track,
                    source_name = data_country_source$source_name, country_track = data_country_source$country_track,
                    cases = data_country_source$cases,neigh_1 = data_country_source$neigh_1,neigh_2 = data_country_source$neigh_2,
                    countryNumSpatial = data_country_source$countryNumSpatial,
                    initial_year_track = data_country_source$initial_year_track,init_year_country = -9999,data_country_source$cov_mat)
  
  test <- unique(test, by = c("country","t"))
  test[, init_year_country := train$init_year_country[match( test$country_name,train$country_name)]]
  
  model_weights_dt <- data.table(read.csv(file=paste0(path_cluster,data_run,folder_name_save,'model_comparison_WEIGHTS_abx: ',list_abx[i],'.csv') ))
  colnames(model_weights_dt) <- c("model_name","weight")
  
  #call function 
  #with source level data 
  new_countries <- test[country_track ==-1000] 
  sourceID <- unique(data.frame(train$source_name,train$source_track))
  neigh_ID <- unique(test[,c("country_name","neigh_1","neigh_2")])
  colnames(sourceID) <- c("source_name","source_track")
  
  #informative prior
  sourceTestData <- train %>%
    dplyr::group_by(country_name) %>%
    expand(year_id = 2000:2022, source_name)
  
  sourceTestData$source_track <- sourceID$source_track[match( sourceTestData$source_name,sourceID$source_name)]
  col_select <- c("country_name","year_id","source_name","resis_prop","cases")
  sourceTest <- merge(sourceTestData,train[, ..col_select], by = c("country_name","year_id","source_name"),all.x = TRUE)
  
  sourceTest$neigh_1 <- neigh_ID$neigh_1[match( sourceTest$country_name,neigh_ID$country_name)]
  sourceTest$neigh_2 <- neigh_ID$neigh_2[match(sourceTest$country_name,neigh_ID$country_name)]
  sourceTest$country_track <- train$country_track[match(sourceTest$country_name,train$country_name)]
  sourceTest$owid_region <- train$owid_region[match(sourceTest$country_name,train$country_name)]
  sourceTest$countryNumSpatial <- train$countryNumSpatial[match(sourceTest$country_name,train$country_name)]
  sourceTest$owid_region_numeric <- train$owid_region_numeric[match(sourceTest$country_name,train$country_name)]
  
  sourceTest$t <- sourceTest$year_id -2000 #1990
  sourceTest$country <- sourceTest$country_track
  sourceTest$source <- sourceTest$source_track
  
  #covar
  col_select <- c("country_name","year_id", covar_list)
  sourceTest <- merge(sourceTest,test[, ..col_select], by = c("country_name","year_id"),all.x = TRUE)
  
  
  #sourceTest <- sourceTest[order(sourceTest$country_name,sourceTest$source,sourceTest$t), ]
  sourceTest <- sourceTest[order(sourceTest$country_name,sourceTest$t), ] #only sort by country and time since we need it ordered for the other models 
  
  #run per region? 
  predictions_forward <- run_predictions_stacking(model_weights_dt,sourceTest, list_abx[i])
  predictions <-  rbind(cbind(predictions_forward$predictions,sourceTest[,c("t","country_name","resis_prop","year_id","source_name","cases")])
  )
  
  write.csv(predictions,paste0(path_cluster,data_run,folder_name_save,'preds_stacked: ',list_abx[i],'.csv'), row.names = FALSE)
  
  #averaging the predictions across the sources 
  predDT <- data.table(predictions)
  predDT <- predDT[ , lapply(.SD, mean) , by=c("t", "country_name"),  .SDcols = sapply(predDT, is.numeric) ][,-1]
  predictions <- merge(predDT,train, by = c("t", "country_name"), all.y = TRUE,all.x = TRUE)
  colnames(predictions)[colnames(predictions) == "source_name"] <- "source_name.y"
  
  predictions$ModelType <- "Stacked Model"
  model_colors <- "deepskyblue3" #mako(length(unique(predictions$ModelType)),begin = 0.15, end = 0.50,option = "D")
  source_colors <- mako(length(unique(na.omit(predictions$source_name.y))),begin = 0.15, end = 0.85)
  names(source_colors) <- unique(na.omit(predictions$source_name.y))
  names(model_colors) <- unique(predictions$ModelType)
  
  #plots and save as pdf 
  p <- ggplot(predictions, mapping = aes(x = t + 2000)) + geom_ribbon(aes(ymin = q5 , ymax = q95, fill = ModelType), alpha = 0.35) + 
    geom_line(mapping = aes(x = t + 2000, y = mean, colour = ModelType)) +  
    geom_point(mapping = aes(y = resis_prop.y, size = cases.y, fill = factor(source_name.y),shape=factor(source_name.y)),colour = "grey",stroke = 0.2)  +
    #scale_fill_manual(name = "Source",values = colorRampPalette(brewer.pal(8, "BuPu"))(length(unique(data_model$source_name)))) + 
    scale_shape_manual(name = "Data Sources",values = rep(c(21:25), length.out = length(source_colors)), breaks = c(names(source_colors))) +
    scale_colour_manual(name = "Model Type", values=model_colors)+
    scale_fill_manual(name = "Model CrI",values = c(model_colors, source_colors), breaks = c(names(model_colors))) + 
    scale_size_continuous(range =  c(1.5,5),
                          breaks = c(100, 500 , 5000, 10000,20000),
                          #breaks = c(1, 30, 50, 80, 100, 150),
                          labels = ~paste0(scales::comma(x = .x)), 
                          name = "Total Observations (R+S)") + 
    labs(x = "Year", y = "Proportion Resistant",fill="Model CrI/ Source", colour="Model Type") + ylim(0,1) + xlim(2000, 2022) + theme_minimal() +
    theme(legend.title = element_text(size=3),legend.text = element_text(size = 3),legend.key.size = unit(0.2, "cm"), legend.position = "right") + 
    facet_wrap_paginate(~country_name, ncol = 4, nrow=4) 
  p
  
  for(k in 1:n_pages(p)){
    
    p_save <-  p +
      facet_wrap_paginate(~country_name, ncol = 4, nrow = 4, page = k, scales = "free")
    ggsave(paste0(path_cluster,data_run,folder_name_save, k, 'abx: ',list_abx[i],"_modelStacked",'.pdf'), p_save, width = 15, height = 8, device = "pdf") #width = 10 , height = 6

    
  }
  
})
  
  
  
  
  
  
  
  
  
  
  
  
  