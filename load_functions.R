# Purpose: Functions for Antimicrobial Prevalence Modelling using Logistic Growth Modelling
# @author:Oraya S
# @Date: Developed: Jan 31 2025, Semi-Cleaned and Semi-Commented: May 13th 2025
# 

initfun <- function(...) {
 
  list(beta=runif(1,0.2,1), gamma=runif(1,0,1), I0_logit=runif(1,-1,1),eta=runif(length(covar_list),-0.05,0.05))
}

inv.logit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

logit <- function(x) {
  return(log(x/(1-x)))
}

#Test function  for scaling and 0in data when merging to the larger covar data frame
scale_format <- function(total_data,list_var_scale_filter){
  for(i in 1:length(list_var_scale_filter)){
    total_data[,paste0(list_var_scale_filter[i])] <- scale(total_data[,paste0(list_var_scale_filter)])[,1]
    total_data[,paste0(list_var_scale_filter[i])][is.na(total_data[,paste0(list_var_scale_filter[i])])] <- 0 
  } 
  return(total_data)}

#just pasting the prepare data covar function 
prepare_data_covar <- function(source_filter){
  
  #explore the data from A. bau - replace with your pathogen data
  a_bau <- data.frame(read.csv(file=paste0(path_cluster,"ab.csv")))
  a_bau$resistance_prop <- a_bau$resistant_unadj/a_bau$cases
  
  #load the shape file to merge your data 
  world_shape <- st_read(paste0(path_cluster,"World Country Boundaries_20250116/geo_export_39720b2d-42d1-46a5-8d79-f5c21ec46ccf.shp"))
  colnames(world_shape)[colnames(world_shape) =="country"] <- "COUNTRY"
  world_shape <- change_name(world_shape)
  
  unique_map_country <- unique(world_shape$COUNTRY)
  df_all <- a_bau %>% expand(loc_name = unique_map_country,year_id = 1990:2022, source, abx_class)

  #summarise your data , skip this if your data is already grouped
  if(source_filter == 1){
    data_grouped <- a_bau %>% group_by(loc_name,location_id, year_id, abx_class, source) %>% dplyr::summarise(cases_sum = sum(cases), resistant_sum = sum(resistant_unadj), susceptible_sum = sum(susceptible_unadj))
  }else{data_grouped <- a_bau %>% group_by(loc_name,location_id, year_id, abx_class) %>% dplyr::summarise(cases_sum = sum(cases), resistant_sum = sum(resistant_unadj), susceptible_sum = sum(susceptible_unadj))}
  
  data_grouped$resistance_prop <- data_grouped$resistant_sum/data_grouped$cases_sum
  data_grouped$susceptible_prop <- data_grouped$susceptible_sum/data_grouped$cases_sum
  

  #merge your covariates to the data based on country and year. An example is shown. This process is repeated for all covariates data that 
  #we wanted merged:
  x <- data.frame(read.csv(file=paste0(path_cluster,"covar_joined.csv")))
  colnames(data_grouped)[colnames(data_grouped) =="loc_name"] <- "COUNTRY"
  data_grouped_covar <- merge(x,data_grouped, by =c("year_id", "location_id", "COUNTRY"), all.y = TRUE) # merge to match the available pathogen data
  

  #get region classification from OWID
  owid <- data.frame(read.csv(file=paste0(path_cluster,"continents-according-to-our-world-in-data/continents-according-to-our-world-in-data.csv")))
  colnames(owid)[colnames(owid) == "Entity"] <- "COUNTRY"
  owid <- change_name(owid)
  data_grouped_covar$owid_region <- owid$World.regions.according.to.OWID[match( data_grouped_covar$COUNTRY,owid$COUNTRY)]
  data_grouped_covar$owid_region[data_grouped_covar$COUNTRY %in%  c("Bonaire","Saba","Saint Eustatius","Virgin Islands")] <- "North America"
  
  return(data_grouped_covar)
}

#reformat country names to match shape file - this will be different depending on naming: example
change_name <- function(df){
 
  df$COUNTRY[df$COUNTRY == "Vietnam"] <- "Viet Nam"
  #....
  
  
  return(df)
}

#this function prepares the data for the models - it creates a huge list depending on the model type.
#the model can then extract specific parameters in the list to use for the runs 
#the function takes in the model and antibiotic name so it can be used in a loop to run all models and antibiotics 
prep_data_abxFilt <- function(df, source, abx_name, calc_lik, model_name){
  
  world_shape <- st_read(paste0(path_cluster,"World Country Boundaries_20250116/geo_export_39720b2d-42d1-46a5-8d79-f5c21ec46ccf.shp"))
  colnames(world_shape)[colnames(world_shape) =="country"] <- "COUNTRY"
  world_shape <- change_name(world_shape)
  
  #filters data based on antibiotic name
  data_with_country <- df %>% filter(abx_class %in% c(paste0(abx_name))) 
  data_with_country = data.table(data_with_country)
  
  #since our code now takes a look at either order country or the source+country combination,
  #we need to order our data to match the order that we use for the random effects
  
  #first we look at if source is mentioned in the model/file name -> if it is, we then we need to order the table by the source and 
  #country combination over time
  
  #if "source" is not present in the file/model name, it just gets ordered by country and year
  #current version of the model does not need this sorting
  if(grepl("source", model_name)){
    ### TEST FOR SORTING 
    country_numeric <- as.numeric(factor(as.vector(data_with_country$COUNTRY)))
    data_with_country$source_numeric <- as.numeric(factor(as.vector(data_with_country$source)))
    data_with_country$new_source <- as.integer(as.factor(paste(data_with_country$source_numeric, country_numeric, sep = "_")))
    
    setkeyv(data_with_country, c("COUNTRY","source_numeric","year_id"))
    
    }else{
  
    # needs to be in order of year for each country
    setkeyv(data_with_country, c("COUNTRY","year_id"))
    }
  
  #get the number of unique countries
  N_country <- length(unique(data_with_country$COUNTRY))
  #get the number of observations or rows within our subset data
  N_data <- nrow(data_with_country)
  #format time to start from 1 for model runs 
  t_data <- data_with_country$year_id - 2000 #TO DO CHANGE LATER SO NO HARDCODING # min(data_with_country$year_id)
  #get the resistance proportion we calculated earlier 
  y_data <- data_with_country$resistance_prop
  
  #code the numbers in numeric form 
  country_numeric <- as.numeric(factor(as.vector(data_with_country$COUNTRY)))
  #get the country name
  country_name <- factor(as.vector(data_with_country$COUNTRY))
  
  #get the total number of cases tested 
  cases_sum <- round(data_with_country$cases_sum)
  #the the number of resistant cases 
  resis_cases <- round(data_with_country$resistant_sum)
  
  #now get your covariate data example: 
  frac_oop_hexp_ihme <- data_with_country$frac_oop_hexp
  ddd_per_1000_fitted_ihme <- data_with_country$ddd_per_1000_fitted
  cov_temp <- data_with_country$Temperature

  
  ###FOR SPATIAL CORRELATION EDGELIST
  filt_loc <- unique(data.table(country_numeric,country_name))
  colnames(filt_loc)[colnames(filt_loc) =="country_name"] <- "COUNTRY"
  #filt_loc <- unique(data_with_country[,c("COUNTRY","location_id")])
  filt_loc$row_id <- 1:nrow(filt_loc)
  map_filtered <- left_join(filt_loc,world_shape, by = "COUNTRY")
  map_filtered <- st_as_sf(map_filtered)
  
  #nb <- poly2nb(map_filtered, row.names = map_filtered$location_id)
  nb <- poly2nb(map_filtered)
  adj <- nb2mat(nb,style="B", zero.policy=TRUE)
  
  #we can to change the node values from the edge list to the actual country values used for analysis
  edge_list <- data.frame(adjmat_to_edgelist(adj))
  edge_list$node1_locId <- filt_loc$country_numeric[match(edge_list[,"ego"], filt_loc$row_id)]
  edge_list$node2_locId <- filt_loc$country_numeric[match(edge_list[,"alter"], filt_loc$row_id, nomatch = NA)]
  
  #test creating distance matrix in m
  #shp_centroid <- st_point_on_surface(map_filtered)
  # calculate the centroids of each country
  map_filtered_centroids <- st_centroid(map_filtered)
  
  #calculate the distance matrix between these centroids
  distance_mat <- cbind(st_distance(map_filtered_centroids))
  
  
  #get owid regions 
  N_regions_owid <- sum(!is.na(unique(data_with_country$owid_region)))
  owid_region <- data_with_country$owid_region
  owid_region_numeric <-  as.numeric(factor(as.vector(data_with_country$owid_region)))
  
  # #test distance-based neighbours -knearneigh ensures all areas have k neighbours
  dist_nb <- knn2nb(knearneigh(st_centroid(map_filtered), k=1), row.names=map_filtered$location_id)
  dist_adj <- nb2mat(dist_nb,style="B", zero.policy=TRUE)
  edge_list_dist <- data.frame(adjmat_to_edgelist(dist_adj))
  edge_list_dist$node1_locId <- filt_loc$country_numeric[match(edge_list_dist[,"ego"], filt_loc$row_id)]
  edge_list_dist$node2_locId <- filt_loc$country_numeric[match(edge_list_dist[,"alter"], filt_loc$row_id, nomatch = NA)]

   N_edges <- nrow(edge_list_dist)
   node1 <- edge_list_dist[,"node1_locId"]
   node2 <- edge_list_dist[,"node2_locId"]

  
  #add the scaling factor - should have the same number as the number of countries
  scaling_mat <- -Matrix(adj, sparse = TRUE)
  diag(scaling_mat) = 0
  diag(scaling_mat) = -rowSums(scaling_mat)
  n = dim(scaling_mat)[1]
  scaled <- inla.scale.model(scaling_mat, constr=list(A=matrix(1,1,n), e=0))
  
  #Extract the scaling factor - there should be a scaling factor for each country
  scaling_factor <- diag(scaled)
  ###END SPATIAL CORRELATION EDGELIST
  
  ##Format the SOURCE data 
  N_source <- length(unique(data_with_country$source))
  source_numeric <- as.numeric(factor(as.vector(data_with_country$source)))
  new_source = as.integer(as.factor(paste(source_numeric, country_numeric, sep = "_")))
  new_source_country = data.table(Country = country_numeric, New_Source = new_source)
  new_source_country = unique(new_source_country)
  setkey(new_source_country, "New_Source")
  new_source_country = new_source_country$Country
  N_source_country_combin <- length(new_source_country)
  
  #if the original data is separated by source
  if(source == 1){
    #source data -add in covariate data later 
    data_country <- list(N = N_data, t = t_data, y = resis_cases, calc_likelihood = calc_lik, C = N_country, country = country_numeric, cases = cases_sum, S = N_source, 
                         source = source_numeric, N_SC_comb = N_source_country_combin, new_source = new_source,new_source_country = new_source_country, temp = cov_temp,
                          N_edges = N_edges, node1 = node1, node2 = node2, 
                         scaling_factor = scaling_factor, country_name = country_name, year_id = data_with_country$year_id,resis_prop = y_data, source_name = data_with_country$source, D = distance_mat,
                         N_regions_owid = N_regions_owid, owid_region = owid_region, owid_region_numeric = owid_region_numeric,
                         mean_temp_ihme =  mean_temp_ihme, haqi_ihme = haqi_ihme, hospital_beds_per1000_ihme = hospital_beds_per1000_ihme,
                         frac_oop_hexp_ihme = frac_oop_hexp_ihme, ddd_per_1000_fitted_ihme = ddd_per_1000_fitted_ihme
    )
  }else{
    
    #just country 
    data_country <- list(N = N_data, t = t_data, y = resis_cases, calc_likelihood = calc_lik, C = N_country, country = country_numeric, cases = cases_sum, temp = cov_temp,
                         precip = cov_precip, WHOhandwash = cov_WHOhandwash, ddd_per_1000 = cov_abxConsumption, N_edges = N_edges, node1 = node1, node2 = node2, 
                         scaling_factor = scaling_factor, country_name = country_name,year_id = data_with_country$year_id,resis_prop = y_data,D = distance_mat,
                         N_regions_owid = N_regions_owid, owid_region = owid_region, owid_region_numeric = owid_region_numeric,
                         frac_oop_hexp_ihme = frac_oop_hexp_ihme, ddd_per_1000_fitted_ihme = ddd_per_1000_fitted_ihme
                         )#,source_name = data_with_country$source)
  }
  return(data_country)
}
