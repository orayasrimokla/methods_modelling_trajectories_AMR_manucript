#file for fitting the model in R 

#for local 
path_cluster = "yourPath"
data_run <- "modelRunExamples/" 
dir.create(paste0(path_cluster,data_run))

#load the functions and packages
source(paste0(path_cluster,"load_packages.R"))
source(paste0(path_cluster,"load_functions.R"))
options(mc.cores = parallel::detectCores()) #stan model will still run the number of cores equal to chain 
#so for our code we have 3 chains so it will run at most 3 cores 

#run and save the prepare_data_covar -> read from CSV (faster this way)
#just going to read the file for now since it was taking too long to rerun MICE 
#data_grouped_covar_source <- prepare_data_covar(source_filter = 1)
#write.csv(data_grouped_covar_source, paste0(path_cluster,"data_grouped_covar_source.csv"), row.names = FALSE)
data_grouped_covar_source <-data.frame(read.csv(file=paste0("yourPath/data_grouped_covar_source.csv")))

#this should be the same for both data aggregates 
country_id <- sort(unique(data_grouped_covar_source$location_id))



#uncomment this line and comment the line below to just test with this for location id less than 200 for example
#data_test_source <- data_grouped_covar_source[data_grouped_covar_source$location_id < 200,]

#all countries 
data_test_source <- data_grouped_covar_source

#selecting baseline year to be 2000 - filter per region -  lots of ways to define what countries are in what region
data_test_source <- data_grouped_covar_source[data_grouped_covar_source$year_id >= 2000 & data_grouped_covar_source$owid_region == "Asia" ,]

#get the data that is not NA
data_test_source_train <- data_test_source[!is.na(data_test_source$cases_sum),]

list_models_run <- c("betaGammaI0_country","betaGammaI0_covar",
                     "betaGammaI0_AddSource","spatial_cmdStan","spatialcmdStanSource",
                     "spatialcmdStanSourceCovar" )

list_models_run <- c("betaGammaI0_country" )

#remove the [1] to run all abx - right now just set to first one 
list_abx <- sort(unique(data_test_source$abx_class))[1]


#test with just the abx consumption for that class - for A baum - select the covariate you want to test
covar_list <- c("ddd_aminoglycoside","ddd_per_1000_fitted_ihme","ddd_betaLactam_penicilins","ddd_quinolone","ddd_othrBetaLactam")

#this is just normal fitting the models but added the a few lines to keep track of the neighbours for the predictions 
llply(1:length(list_models_run), function(i){
  
  #create folder to save based on model we want 
  dir.create(paste0(path_cluster,data_run,list_models_run[i],"/"))
  
  #get model run times
  #times <- matrix(nrow= length(list_abx), ncol = 2)
  
  llply(1:length(list_abx), function(j){
    
    data_country_source_train <- prep_data_abxFilt(data_test_source_train, source = 1, abx_name = list_abx[j], calc_lik = 1, list_models_run[i])
    data_country_source_train$cov_mat <- do.call(cbind, data_country_source_train[paste0(covar_list)])
    data_country_source_train$K <- length(covar_list)
    
    data_model <- data_country_source_train
    
    niter <- 5000 #can test with others
    
    #create model
    model_name_path <- paste0(path_cluster,"Models/",list_models_run[i],".stan")
    
    #spatial models run with regional data filtered since global one needs different priors/scaling
    if(grepl("spatial", list_models_run[i])){
      
        mod <- cmdstan_model(model_name_path)
        #simplifying for now without regional data - redundant
        data_model <- prep_data_abxFilt(data_test_source_train, source = 1, abx_name = list_abx[j], calc_lik = 1, list_models_run[i])
        
        data_model$D <- data_model$D/(6378100*pi) #furthest point on globe 
        
        #add in covariate data information
        data_model$cov_mat <- do.call(cbind, data_model[paste0(covar_list)])
        data_model$K <- length(covar_list)
        
        data_model$country_numTrain <- data_country_source_train$country[match(data_model$country_name,data_country_source_train$country_name)]
        data_cmdstan <- data_model[sapply(data_model, is.numeric)]
        data_cmdstan <- data_cmdstan[!sapply(data_cmdstan, function(x) any(is.na(x)))]
        
        fit <-mod$sample(
          iter_warmup =niter/2,
          iter_sampling =niter/2,
          data=data_cmdstan,
          chains=4,
          parallel_chains = 4,
          seed=100, 
          refresh=100,
          adapt_delta = 0.95,
          max_treedepth = 30) 
        
       
        fit$save_object(file = paste0(path_cluster,data_run,list_models_run[i],"/",'model: ',list_models_run[i],"+2000:2022",list_abx[j],".rds"))
      
    }else{
    model <-stan(model_name_path,
                 iter=niter,
                 data=data_model,
                 chains=4,
                 seed=100, 
                 refresh=100,
                 
                 control = list(
                   adapt_delta = 0.95, 
                   max_treedepth = 30)) 
    
    save(model, ascii=FALSE, file=paste0(path_cluster,data_run,list_models_run[i],"/",'model: ',list_models_run[i],'+2000:2022',list_abx[j]))
    
    }
    
  })

})





