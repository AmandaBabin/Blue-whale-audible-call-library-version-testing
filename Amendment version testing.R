## Blue whale audible call library amendement version testing

# Authors: G. Macklin, J. Stanistreet, and A. Babin 2025

# Download and install packages if not already installed: pacman
library(pacman)

# Then open the packages
p_load(tidyverse, purrr, ggpattern, here)

#BmA vs. other species detectors----

# Set the path to your main folder
main_folder <- r"(R:\Science\CetaceanOPPNoise\CetaceanOPPNoise_5\BaleenWhale_AcousticAnalysis\Deployments\MAR)" # direct to all deployments

# List all subfolders
subfolders <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE) #list all subfolders (aka Validation, Results)
subfolders <-subfolders[!str_detect(subfolders, "TRAIN")]

# Get all results

# Initialize an empty list to store dataframes
dfs <- list()

# Iterate over each subfolder, read the CSV files, and combine into a single dataframe
for (subfolder in subfolders) {
  deployment <- basename(subfolder) #extract deployment name
  
  # Read the CSV file
  csv_file <- list.files(paste0(subfolder,"\\Results\\"), pattern = "_FINAL.csv", #find the presence csv
                         full.names = TRUE)
  
  df <- read_csv(csv_file) #read presence csv
  
  df$deployment <- deployment
  
  df<- df %>% 
    mutate(deployment = deployment,
           station = str_extract(deployment, "^[^_]+")) %>%
    relocate(station,deployment)
  
  # Store the dataframe in the list
  dfs[[length(dfs) + 1]] <- df
  
  dfs <- keep(dfs, ~ nrow(.) > 0)
  
}

dfs <- lapply(dfs, function(x) {
  select(x,
         station,
         deployment,
         filename,
         species = Species,
         call_type = `Call type`,
         call_cat = V25)
})

data.input <- bind_rows(dfs)

sp_list <- c("BW","FW","SW","HB","RW")

all.results <- data.input %>% 
  filter(species %in% sp_list) %>% 
  mutate(sp_code = case_when((species == "BW" & call_cat =="A")~ "BmA",
                             (species == "BW" & call_cat =="IF") ~ "BmT",
                             species == "FW"~ 'Bp',
                             species == "SW"~ "Bb",
                             species == "HB"~ 'Mn',
                             species == "RW" ~ 'Eg'))

## Get all detections

# Initialize an empty list to store dataframes
dfs <- list()

# Iterate over each subfolder, read the CSV files, and combine into a single dataframe
for (subfolder in subfolders) {
  deployment <- basename(subfolder) #extract deployment name
  
  # Read the CSV file
  csv_file <- list.files(paste0(subfolders,"\\Validation\\Arklite_Inputs\\"), pattern = paste0(deployment,"_MBDetections.csv"), #find the detection csv
                         full.names = TRUE)
  
  df <- read_csv(csv_file) #read presence csv
  
  df$deployment <- deployment
  
  df<- df %>% 
    mutate(deployment = deployment,
           station = str_extract(deployment, "^[^_]+")) %>%
    relocate(station,deployment)
  
  # Store the dataframe in the list
  dfs[[length(dfs) + 1]] <- df
  
  dfs <- keep(dfs, ~ nrow(.) > 0)
  
}

detection.input <- bind_rows(dfs)

tier_list <- c(BmT="MD4",BmA="MD4",Bp="MD3",Bb="MD4",Mn="MD3", Eg="MD4")

single_sp = "BmA" #Blue whale audible

## Isolate detections of interest

detections.sing <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == single_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(station,deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection_sing = 1) %>%
  unique() 

## Isolate results of interest

results <- all.results %>% 
  filter(sp_code == single_sp) %>% 
  mutate(presence = 1) %>% 
  select(station,deployment,filename, sp_code,presence) %>% 
  unique() 

#for each species/calltype(s) of interest, create the detections dataframe

second_sp = "Bp" # Blue whale infrasonic: BmT, Sei whale": Bb, Humpback whale: Mn, Fin whale: Bp

detections.BmT <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == second_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(station,deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection_BmT = 1) %>% 
  unique() %>% 
  select(-sp_code)

detections.Bb <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == second_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(station,deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection_Bb = 1) %>% 
  unique() %>% 
  select(-sp_code)

detections.Mn <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == second_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(station,deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection_Mn = 1) %>% 
  unique() %>% 
  select(-sp_code)

detections.Bp <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == second_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(station,deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection_Bp = 1) %>% 
  unique() %>% 
  select(-sp_code)

full_join <- Reduce(function (...) { merge(..., all = TRUE) },
               list(results,detections.sing,detections.BmT,detections.Bb,detections.Mn,detections.Bp))

day.results <- full_join %>%
  mutate(presence = replace_na(presence,0),
         detection_sing = replace_na(detection_sing,0),
         detection_BmT = replace_na(detection_BmT,0),
         detection_Bb = replace_na(detection_Bb,0),
         detection_Mn = replace_na(detection_Mn,0),
         detection_Bp = replace_na(detection_Bp,0)) %>%
    arrange(station,deployment,filename) %>% 
  
  mutate(datestring = str_extract(filename, "\\d{8}\\w\\d{6}\\w")) %>% 
  mutate(filedate = as_date(datestring, format="%Y%m%dT%H%M%SZ")) %>%
  
  group_by(station,deployment,sp_code,filedate) %>% 
  summarise(presence = sum(presence),
            detection_sing = sum(detection_sing),
            detection_BmT = sum(detection_BmT),
            detection_Bb = sum(detection_Bb),
            detection_Mn = sum(detection_Mn),
            detection_Bp = sum(detection_Bp)) %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  ungroup() %>% 
  mutate(sp_code = case_when(detection_sing >0~"BmA",
                            (detection_sing ==0 & detection_BmT>0)~"BmT",
                            (detection_sing ==0 & detection_Bb>0)~"Bb",
                            (detection_sing ==0 & detection_Mn>0)~"Mn",
                            (detection_sing ==0 & detection_Bp>0)~"Bp"))

#this table shows the percent of days that contained files with blue whale audible call detections divided by the number of days with blue whale audible call presence,
#followed by the percent of the days that contained files with blue whale infratonic call detections divided by the number of REMAINING days with BmA call presence,
#THEN sei whale call detections, humpback whale call detections, and fin whale call detections, in the order that they are generally validated
percent.results <- day.results %>% 
  group_by(station,deployment) %>% 
  summarise(days_no_BmA_detection = sum(presence>0 & detection_sing ==0),
            days_BmA_detection = sum(presence >0 & sp_code=="BmA" & detection_sing > 0),
            days_BmT_detection = sum(presence>0 & sp_code=="BmT" & detection_BmT >0),
            days_Bb_detection = sum(presence>0 & sp_code=="Bb" & detection_Bb >0),
            days_Mn_detection = sum(presence>0 & sp_code=="Mn" & detection_Mn >0),
            days_Bp_detection = sum(presence>0 & sp_code=="Bp" & detection_Bp >0)) %>% 
  ungroup() %>%
  
  group_by(deployment) %>%
  mutate(total_presence_days = sum(days_no_BmA_detection + days_BmA_detection),
         percent_BmA = round((days_BmA_detection / total_presence_days) * 100,0),
         percent_BmT = round((days_BmT_detection / (total_presence_days-days_BmA_detection)) * 100,0),
         percent_Bb = round((days_Bb_detection / (total_presence_days-days_BmA_detection)) * 100,0),
         percent_Mn = round((days_Mn_detection / (total_presence_days-days_BmA_detection)) * 100,0),
         percent_Bp = round((days_Bp_detection / (total_presence_days-days_BmA_detection)) * 100,0)) %>%
  select(deployment,percent_BmA,percent_BmT,percent_Bb,percent_Mn,percent_Bp) %>%
  ungroup()

#number of files----

# Set the path to your main folder
main_folder <- r"(E:\)" # direct to working harddrive same setup of files as the Deployment folder (Results, Validation, etc.)
                        # organized in two folders of the 'original version' and 'new versions' (v2-v7) of the blue whale audible call library

version <- "new versions\\v3_pooled"  # original version
                                                    # new versions\\v2_additional exemplars
                                                    # new versions\\v3_pooled
                                                    # new versions\\v4_reduced
                                                    # new versions\\v5_12clusters
                                                    # new versions\\v6_12clusters_MD4
                                                    # new versions\\v7_17clusters_MD4

# List all subfolders
subfolders <- list.dirs(paste0(main_folder,version), full.names = TRUE, recursive = FALSE) #list all subfolders
#subfolders <-subfolders[!str_detect(subfolders, "TRAIN")] #remove any deployments, if wanted

# Get all results

# Initialize an empty list to store dataframes
dfs <- list()

# Iterate over each subfolder, read the CSV files, and combine into a single dataframe
for (subfolder in subfolders) {
  deployment <- basename(subfolder) #extract deployment name
  
  # Read the CSV file
  csv_file <- list.files(paste0(subfolder,"\\Validation\\ArkLite_Inputs\\"), pattern = "_BmA4.csv",
                         full.names = TRUE)
  
  df <- read_csv(csv_file,col_names=FALSE)
  
  colnames(df) <- "filename"
  
  df$deployment <- deployment
  
  df<- df %>% 
    mutate(deployment = deployment)

  # Store the dataframe in the list
  dfs[[length(dfs) + 1]] <- df
  
  dfs <- keep(dfs, ~ nrow(.) > 0)
  
}

dfs <- lapply(dfs, function(x) {
  select(x,
         deployment,
         filename)
})

deployment_filename <- bind_rows(dfs)

#run the number_file and colnames corresponding to the version below
number_files_v1 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v1) <- c("deployment","v1_n_BmA4_files")

number_files_v2 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v2) <- c("deployment","v2_n_BmA4_files")

number_files_v3 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v3) <- c("deployment","v3_n_BmA4_files")

number_files_v4 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v4) <- c("deployment","v4_n_BmA4_files")

number_files_v5 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v5) <- c("deployment","v5_n_BmA4_files")

number_files_v6 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v6) <- c("deployment","v6_n_BmA4_files")

number_files_v7 <- deployment_filename %>%
  group_by(deployment) %>%
  count()
colnames(number_files_v7) <- c("deployment","v7_n_BmA4_files")

#then merge them all together
number_files <- Reduce(function (...) { merge(..., all = TRUE) },
                    list(number_files_v1,number_files_v2,number_files_v3,number_files_v4,number_files_v5,number_files_v6,number_files_v7))

#the number of BmA4 files higher or lower than v1 for each version 
number_files_vs_v1 <- number_files %>%
  mutate(
    v2_n_higher_or_lower_than_v1 = v2_n_BmA4_files-v1_n_BmA4_files,
    v3_n_higher_or_lower_than_v1 = v3_n_BmA4_files-v1_n_BmA4_files,
    v4_n_higher_or_lower_than_v1 = v4_n_BmA4_files-v1_n_BmA4_files,
    v5_n_higher_or_lower_than_v1 = v5_n_BmA4_files-v1_n_BmA4_files,
    v6_n_higher_or_lower_than_v1 = v6_n_BmA4_files-v1_n_BmA4_files,
    v7_n_higher_or_lower_than_v1 = v7_n_BmA4_files-v1_n_BmA4_files)

new_versions_files_vs_v1 <- number_files_vs_v1 %>%
  mutate(
    v2_percent_files_higher_or_lower_than_v1 = round(((v2_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0),
    v3_percent_files_higher_or_lower_than_v1 = round(((v3_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0),
    v4_percent_files_higher_or_lower_than_v1 = round(((v4_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0),
    v5_percent_files_higher_or_lower_than_v1 = round(((v5_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0),
    v6_percent_files_higher_or_lower_than_v1 = round(((v6_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0),
    v7_percent_files_higher_or_lower_than_v1 = round(((v7_n_BmA4_files-v1_n_BmA4_files)/v1_n_BmA4_files)*100,0))

#relative recall----

single_sp = "BmA"

main_folder <- r"(E:\)" # direct to working harddrive same setup of files as the Deployment folder (Results, Validation, etc.)
# organized in two folders of the 'original version' and 'new versions' (v2-v7) of the blue whale audible call library

version <- "new versions\\v7_17clusters_MD4"  # original version
                                                    # new versions\\v2_additional exemplars
                                                    # new versions\\v3_pooled
                                                    # new versions\\v4_reduced
                                                    # new versions\\v5_12clusters
                                                    # new versions\\v6_12clusters_MD4
                                                    # new versions\\v7_17clusters_MD4

subfolders <- list.dirs(paste0(main_folder,version), full.names = TRUE, recursive = FALSE) #list all subfolders
#subfolders <-subfolders[!str_detect(subfolders, "TRAIN")] #remove any deployments, if wanted

## Get all results

# Initialize an empty list to store dataframes
dfs <- list()

# Iterate over each subfolder, read the CSV files, and combine into a single dataframe
for (subfolder in subfolders) {
  deployment <- basename(subfolder) #extract deployment name
  
  # Read the CSV file
  csv_file <- list.files(paste0(subfolder,"\\Results\\"), pattern = "_FINAL.csv", #find the presence csv
                         full.names = TRUE)
  
  df <- read_csv(csv_file) #read presence csv
  
  df$deployment <- deployment
  
  df<- df %>% 
    mutate(deployment = deployment,
           station = str_extract(deployment, "^[^_]+")) %>%
    relocate(station,deployment)
  
  # Store the dataframe in the list
  dfs[[length(dfs) + 1]] <- df
  
  dfs <- keep(dfs, ~ nrow(.) > 0)
  
}

dfs <- lapply(dfs, function(x) {
  select(x,
         deployment,
         filename,
         species = Species,
         call_type = `Call type`,
         call_cat = V25)
})

data.input <- bind_rows(dfs)

sp_list <- c("BW","FW","SW","HB","RW")

all.results <- data.input %>% 
  filter(species %in% sp_list) %>% 
  mutate(sp_code = case_when((species == "BW" & call_cat =="A")~ "BmA",
                             (species == "BW" & call_cat =="IF") ~ "BmT",
                             species == "FW" & call_cat =="IS"~ 'Bp',
                             species == "SW" & call_cat =="FF"~ "Bb",
                             species == "HB"~ 'Mn',
                             species == "RW" ~ 'Eg')) %>% 
  filter(!is.na(sp_code))


## Get all detections

# Initialize an empty list to store dataframes
dfs <- list()

# Iterate over each subfolder, read the CSV files, and combine into a single dataframe
for (subfolder in subfolders) {
  deployment <- basename(subfolder) #extract deployment name
  
  # Read the CSV file
  csv_file <- list.files(paste0(subfolders,"\\Validation\\Arklite_Inputs\\"), pattern = paste0(deployment,"_MBDetections.csv"), #find the detection csv
                         full.names = TRUE)
  
  df <- read_csv(csv_file) #read presence csv
  
  df$deployment <- deployment
  
  df<- df %>% 
    mutate(deployment = deployment,
           station = str_extract(deployment, "^[^_]+")) %>%
    relocate(station,deployment)
  
  # Store the dataframe in the list
  dfs[[length(dfs) + 1]] <- df
  
  dfs <- keep(dfs, ~ nrow(.) > 0)
  
}

detection.input <- bind_rows(dfs)

## Isolate detections of interest

tier_list <- c(BmT="MD4",BmA="MD4",Bp="MD3",Bb="MD4",Mn="MD3", Eg="MD4")

detections <- detection.input %>%
  filter(deployment %in% unique(all.results$deployment)) %>% 
  filter(Species == single_sp) %>% 
  filter(.data[[tier_list[[single_sp]]]] != 0) %>% 
  select(deployment,'filename'=Filename,'sp_code'=Species) %>% 
  mutate(detection = 1) %>% 
  unique() 

## Isolate results of interest

results <- all.results %>% 
  filter(sp_code == single_sp) %>% 
  mutate(presence = 1) %>% 
  select(deployment,filename,presence) %>% 
  unique() 

day.results <- full_join(results,detections) %>% 
  mutate(presence = replace_na(presence,0),
         detection = replace_na(detection,0)) %>% 
  
  mutate(datestring = str_extract(filename, "\\d{8}\\w\\d{6}\\w")) %>% 
  mutate(filedate = as_date(datestring, format="%Y%m%dT%H%M%SZ")) %>%
  
  group_by(deployment,filedate) %>% 
  summarise(presence = sum(presence),
            detection = sum(detection)) %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) %>%
  ungroup()

#run the relative.recall corresponding to the version below
relative.recall_v1 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v1_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v1_relative_recall)

relative.recall_v2 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v2_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v2_relative_recall)

relative.recall_v3 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v3_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v3_relative_recall)

relative.recall_v4 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v4_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v4_relative_recall)

relative.recall_v5 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v5_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v5_relative_recall)

relative.recall_v6 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v6_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v6_relative_recall)

relative.recall_v7 <- day.results %>% 
  group_by(deployment) %>% 
  summarise(days_detect = sum(presence >0 & detection > 0),
            days_nodetect = sum(presence>0 & detection ==0)) %>% 
  mutate(v7_relative_recall = round((days_detect/(days_detect+days_nodetect))*100,0)) %>%
  select(deployment,v7_relative_recall)

#then merge them all together
relative.recall <- Reduce(function (...) { merge(..., all = TRUE) },
                       list(relative.recall_v1,relative.recall_v2,relative.recall_v3,relative.recall_v4,relative.recall_v5,relative.recall_v6,relative.recall_v7))

#the relative recall higher or lower than v1 for each version 
relative_recall_vs_v1 <- relative.recall %>%
  mutate(
    v2_recall_higher_or_lower_than_v1 = v2_relative_recall-v1_relative_recall,
    v3_recall_higher_or_lower_than_v1 = v3_relative_recall-v1_relative_recall,
    v4_recall_higher_or_lower_than_v1 = v4_relative_recall-v1_relative_recall,
    v5_recall_higher_or_lower_than_v1 = v5_relative_recall-v1_relative_recall,
    v6_recall_higher_or_lower_than_v1 = v6_relative_recall-v1_relative_recall,
    v7_recall_higher_or_lower_than_v1 = v7_relative_recall-v1_relative_recall)

new_versions_recall_vs_v1 <- relative_recall_vs_v1 %>%
  mutate(
    v2_percent_recall_higher_or_lower_than_v1 = round(((v2_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0),
    v3_percent_recall_higher_or_lower_than_v1 = round(((v3_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0),
    v4_percent_recall_higher_or_lower_than_v1 = round(((v4_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0),
    v5_percent_recall_higher_or_lower_than_v1 = round(((v5_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0),
    v6_percent_recall_higher_or_lower_than_v1 = round(((v6_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0),
    v7_percent_recall_higher_or_lower_than_v1 = round(((v7_relative_recall-v1_relative_recall)/v1_relative_recall)*100,0))

#trade-off between the number of files and relative recall----

#this is an index of files/recall where lower values represent a more efficient trade-off between the number of files you need to include
#in the analysis for the difference in recall

files_recall <- Reduce(function (...) { merge(..., all = TRUE) },
                          list(new_versions_files_vs_v1,new_versions_recall_vs_v1))

new_versions_vs_v1 <- files_recall %>%
  mutate(
    v2_tradeoff = ifelse(v2_percent_recall_higher_or_lower_than_v1==0,0,round(v2_percent_files_higher_or_lower_than_v1/v2_percent_recall_higher_or_lower_than_v1,1)),
    v3_tradeoff = ifelse(v3_percent_recall_higher_or_lower_than_v1==0,0,round(v3_percent_files_higher_or_lower_than_v1/v3_percent_recall_higher_or_lower_than_v1,1)),
    v4_tradeoff = ifelse(v4_percent_recall_higher_or_lower_than_v1==0,0,round(v4_percent_files_higher_or_lower_than_v1/v4_percent_recall_higher_or_lower_than_v1,1)),
    v5_tradeoff = ifelse(v5_percent_recall_higher_or_lower_than_v1==0,0,round(v5_percent_files_higher_or_lower_than_v1/v5_percent_recall_higher_or_lower_than_v1,1)),
    v6_tradeoff = ifelse(v6_percent_recall_higher_or_lower_than_v1==0,0,round(v6_percent_files_higher_or_lower_than_v1/v6_percent_recall_higher_or_lower_than_v1,1)),
    v7_tradeoff = ifelse(v7_percent_recall_higher_or_lower_than_v1==0,0,round(v7_percent_files_higher_or_lower_than_v1/v7_percent_recall_higher_or_lower_than_v1,1)))

#mean and SD of selected deployments and columns

new_versions_vs_v1_selected_deployments <- new_versions_vs_v1 %>%
  filter(deployment %in% c("FCH_2019_10",
                      "GMB_2019_04",
                      "JOBW_2019_04",
                      "MGL_2019_10",
                      "ROB_2019_10",
                      "WSS_2019_10",
                      "COC_2020_09",
                      "EMBS_2020_09",
                      "FCH_2020_09",
                      "GBK_2020_09",
                      "WSS_2020_09",
                      "COC_2021_08",
                      "CSE_2021_05",
                      "EFC_2021_08",
                      "FCM_2021_08",
                      "MBK_2021_09",
                      "MGL_2021_08",
                      "SAB_2021_09",
                      "SBVC3_2021_09",
                      "SFD_2021_08",
                      "CCU_2022_09",
                      "CS1_2022_10",
                      "CS2_2022_10",
                      "CS3_2022_10",
                      "CSE_2022_10",
                      "ELC_2022_10",
                      "FCD_2022_09",
                      "GDSE_2022_10",
                      "GLNE_2022_10",
                      "GLSW_2022_10",
                      "MGE_2022_10",
                      "ROBV_2022_09",
                      "SABW_2022_10",
                      "SBVC1_2022_10",
                      "CBN_2023_08",
                      "CS1_2023_08",
                      "CS2_2023_08",
                      "CS3_2023_08",
                      "FCM_2023_08",
                      "MBK_2023_08",
                      "MGL_2023_08",
                      "MLC_2023_08"))

new_versions_vs_v1_selected_columns<-new_versions_vs_v1_selected_deployments[,c("v2_percent_files_higher_or_lower_than_v1",
                                        "v3_percent_files_higher_or_lower_than_v1",
                                        "v4_percent_files_higher_or_lower_than_v1",
                                        "v5_percent_files_higher_or_lower_than_v1",
                                        "v6_percent_files_higher_or_lower_than_v1",
                                        "v7_percent_files_higher_or_lower_than_v1",
                                        "v1_relative_recall",
                                        "v2_relative_recall",
                                        "v3_relative_recall",
                                        "v4_relative_recall",
                                        "v5_relative_recall",
                                        "v6_relative_recall",
                                        "v7_relative_recall",
                                        "v2_recall_higher_or_lower_than_v1",
                                        "v3_recall_higher_or_lower_than_v1",
                                        "v4_recall_higher_or_lower_than_v1",
                                        "v5_recall_higher_or_lower_than_v1",
                                        "v6_recall_higher_or_lower_than_v1",
                                        "v7_recall_higher_or_lower_than_v1",
                                        "v2_percent_recall_higher_or_lower_than_v1",
                                        "v3_percent_recall_higher_or_lower_than_v1",
                                        "v4_percent_recall_higher_or_lower_than_v1",
                                        "v5_percent_recall_higher_or_lower_than_v1",
                                        "v6_percent_recall_higher_or_lower_than_v1",
                                        "v7_percent_recall_higher_or_lower_than_v1",
                                        "v2_tradeoff",
                                        "v3_tradeoff",
                                        "v4_tradeoff",
                                        "v5_tradeoff",
                                        "v6_tradeoff",
                                        "v7_tradeoff"
                                        )]

new_versions_vs_v1_selected_columns_mean <- sapply(new_versions_vs_v1_selected_columns,mean,na.rm=TRUE)
new_versions_vs_v1_selected_columns_SD <- sapply(new_versions_vs_v1_selected_columns,sd,na.rm=TRUE)
new_versions_vs_v1_selected_columns_mean_SD <-bind_rows(new_versions_vs_v1_selected_columns_mean,new_versions_vs_v1_selected_columns_SD)
