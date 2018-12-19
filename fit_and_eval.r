fit_and_eval <- function(list_of_modalities, outcome, fold_to_evaluate, fold_range = NULL, ...) {
  
  if (length(fold_range) == 0) {up_to_fold <- 1:max(fold_to_evaluate)} else {up_to_fold <- fold_range}
  
  SMO_classifier <- make_Weka_classifier("weka/classifiers/functions/SMOreg")
  
  out <- foreach(fold_index = up_to_fold, .inorder = FALSE, 
                 .packages = c("tidyverse","dplyr", "CORElearn", "spatstat", "numDeriv", "quantmod", "Biocomb", "RWeka"),
                 .export = c("sd_thresholding_for_categorical_outcome_variables_vec", "select_features_relieff_derivatives_threshold_CORElearn",
                             "extract_weights_from_SMO", "SMO_classifier", "list_of_modalities", "outcome", "fold_to_evaluate")) %do% {
                               
                               all_mods_train <- list()
                               all_relieff_features <- list()
                               all_coordinates <- list()
                                 
                               
                               print(paste("working on TRAINING SET fold",fold_index, sep = " "))
                               for (mod in 1:length(list_of_modalities)) {
                                
                                
                                train <- list_of_modalities[[mod]]$matrix[fold != fold_index,]
                                outcome_train <- outcome[fold != fold_index]
                                img_dim <- list_of_modalities[[mod]]$img_dim
                                name_of_mod <- names(list_of_modalities)[[mod]]
                               
                                #variance thresholding
                                print(paste("variance thresholding, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                 var_thr <- sd_thresholding_for_categorical_outcome_variables_vec(train, .25)
                               
                                var_thr$outcome <- outcome_train
                               
                                #relieff
                                print(paste("relieff, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                relieff <- select_features_relieff_derivatives_threshold_CORElearn(var_thr, "outcome", 
                                                                                                     estimator = "RReliefFequalK")
                                rm(var_thr)
                               
                                #coordinates finding 
                                print(paste("coordinates finding, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                coordinates_from_features_colnames <- relieff[,-1] %>%
                                 colnames(.) %>%
                                 str_split(., pattern = "X") %>%
                                 map_chr(~`[`(.,2)) %>%
                                 as.numeric(.) %>%
                                 arrayInd(., img_dim) %>%
                                 as_data_frame(.)
                               
                                coordinates_from_features_colnames$index <- relieff[,-1] %>%
                                 colnames(.) %>%
                                 str_split(., pattern = "X") %>%
                                 map_chr(~`[`(.,2)) %>%
                                 as.numeric(.)
                               
                                print(paste("clustering, modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                coordinates_from_features_colnames <- cluster_voxels(coordinates_from_features_colnames, ...)
                               
                               
                                #averaging of clusters
                                df_clusters <- relieff %>%
                                 select(., -outcome) %>%
                                 mutate(subject = row_number()) %>%
                                 gather(., feature, value, -subject) %>%
                                 mutate(., index = str_split(feature, pattern = "X") %>% map(~`[`(.,2)) %>% unlist() %>% as.numeric()) %>%
                                 inner_join(., select(coordinates_from_features_colnames, index, cluster_id), by = "index") %>%
                                 group_by(subject, cluster_id) %>%
                                 summarise(mean_of_cluster = mean(value)) %>%
                                 spread(cluster_id, mean_of_cluster) %>%
                                 ungroup(.) %>%
                                 select(-subject)
                                
                                relieff_features <- colnames(relieff)[-1]
                                rm(relieff)
                               
                                #change names of features for easiness of recognition later on
                                colnames(df_clusters) <- paste(name_of_mod, colnames(df_clusters),sep ="_")
                                all_mods_train[[mod]] <- df_clusters
                                names(all_mods_train)[mod] <- name_of_mod 
                                all_relieff_features[[mod]] <- relieff_features
                                names(all_relieff_features)[mod] <- name_of_mod 
                                all_coordinates[[mod]] <- coordinates_from_features_colnames
                                names(all_coordinates)[mod] <- name_of_mod }
                               
                                merged_modalities_df <- Reduce(bind_cols, all_mods_train)
                               
                                merged_modalities_df$outcome <- outcome_train
                               
                                merged_modalities_df_selected <- merged_modalities_df %>%
                                 select(., select.cfs(merged_modalities_df)$Index, outcome)
                               
                                rm(merged_modalities_df)
                               
                               #cluster and select test set
                               
                              all_mods_test <- list()
                              print("working on TEST SET")
                              for (mod in 1:length(list_of_modalities)) {
                                print(paste("modality is", name_of_mod, "modality", mod, "of", length(list_of_modalities), sep = " "))
                                test <- list_of_modalities[[mod]]$matrix[fold == fold_index,]
                                img_dim <- list_of_modalities[[mod]]$img_dim
                                name_of_mod <- names(list_of_modalities)[mod]
                                outcome_test <- outcome[fold == fold_index]
                                
                               
                               test_selected <- test %>%
                                 select(., all_relieff_features[[mod]]) %>%
                                 mutate(subject = row_number()) %>%
                                 gather(., feature, value, -subject) %>%
                                 mutate(., index = str_split(feature, pattern = "X") %>% map(~`[`(.,2)) %>% unlist() %>% as.numeric()) %>%
                                 inner_join(., select(all_coordinates[[mod]], index, cluster_id), by = "index") %>%
                                 group_by(subject, cluster_id) %>%
                                 summarise(mean_of_cluster = mean(value)) %>%
                                 spread(cluster_id, mean_of_cluster) %>%
                                 ungroup(.) %>%
                                 select(-subject)
                                 
                               colnames(test_selected) <- paste(name_of_mod, colnames(test_selected),sep ="_")
                               all_mods_test[[mod]] <- test_selected
                               names(all_mods_test)[mod] <- name_of_mod}
                              
                              
                               merged_modalities_df_test <- Reduce(bind_cols, all_mods_test) %>%
                                 select(., head(colnames(merged_modalities_df_selected),-1))
                               
                               
                               
                               model_SMO <- SMO_classifier(outcome ~ ., data = merged_modalities_df_selected)
                               
                               SMO_weights <- extract_weights_from_SMOreg(model_SMO)
                               
                               classification <- predict(model_SMO, merged_modalities_df_test)
                               
                               
                               accuracy <- data_frame(classification = classification, ground = outcome_test)
                               out <- list(all_coordinates, accuracy = accuracy, weights = SMO_weights)
                               
                             }
  
}
  
  