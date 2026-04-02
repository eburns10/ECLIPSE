
ECLIPSE <- function(folders = NULL, file_type = "all", file_paths = NULL, seurat_object, batch, donor, group, original_barcode, write_folder, mode = "standard", format = "Blank") {
                
  #### THIS FIRST PART IS READING THE DATA IN AND DOING QC
                ### Checking the input data is usable
                if (write_folder %in% folders & file_type != "manual") {
                  message("Folders cannot be the same as the write folder. This would result in the writing of a new .csv file over the raw data, removing the raw data from the computer.")
                  return()
                }
  
                if (is.null(folders) & file_type != "manual") {
                  message("Please provide folders that have the contigs annotations in them.")
                  return()
                }
                
                if (is.null(file_paths) & file_type == "manual") {
                  message("No file paths provided. Please provide a vector of file path names to file_paths")
                  return()
                }
  
                if (!batch %in% colnames(seurat_object@meta.data)) {
                  message("Batch column name provided is not in the Seurat object. Please try again.")
                  return()
                }
  
                if (!donor %in% colnames(seurat_object@meta.data)) {
                  message("Donor column name provided is not in the Seurat object. Please try again.")
                  return()
                }
  
                if (!group %in% colnames(seurat_object@meta.data)) {
                  message("Group column name provided is not in the Seurat object. Please try again.")
                  return()
                }
                
                if (!original_barcode %in% colnames(seurat_object@meta.data)) {
                  message("Original barcode column name provided is not in the Seurat object. Please try again.")
                  return()
                }
                
                if (!file_type %in% c("ALL", "All", "all", "Filtered", "FILTERED", "filtered", "all_contig_annotations.csv", "filtered_contig_annotations.csv", "manual")) {
                  print("Incorrect file type. Should be all or filtered. Returning nothing")
                  return()
                }
                
                ### Reading the data in
                
    if (file_type != "manual") {
                
                contigs_list <- vector("list", length(folders))
      
        for (i in 1:length(folders)) {
          
                if (file_type %in% c("ALL", "All", "all", "all_contig_annotations.csv")) {
                  files <- list.files(folders[i], pattern = "all_contig_annotations.csv", full.names = TRUE)
                }
                else if (file_type %in% c("Filtered", "filtered", "FILTERED", "filtered_contig_annotations.csv")) {
                  files <- list.files(folders[i], pattern = "filtered_contig_annotations.csv", full.names = TRUE)
                }
                
          
                if (length(files) == 0) {
                  print("No contig files found, returning nothing")
                  return()
                }
                else if (length(files) > 1) {
                  print("Multiple contig files found, returning nothing")
                  return()
                }
                else if (length(files) == 1) {
                  contigs_list[[i]] <- read.csv(files[1])
                  contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
                }
  
        }
      
    } else if (file_type == "manual") {
      
               contigs_list <- vector("list", length(file_paths))
        
        for (i in 1:length(file_paths)) {
                
                contigs_list[[i]] <- read.csv(file_paths[i])

                if (nrow(contigs_list[[i]]) == 0) {
                  
                  message(paste("File path ", i, " did not lead to any contigs being loaded in. Check the path name and that it is a .csv file", sep = ""))
                  return()
                  
                }
                
                contigs_list[[i]]$barcode <- paste(i, "_", contigs_list[[i]]$barcode, sep = "")
          
          }
          
    }
                
                ### Joins together all batches and adds data to them about which sample/donor/batch each contig was from
                contigs_full <- bind_rows(contigs_list)
                
                if("sample" %in% colnames(contigs_full)) {
                  
                  contigs_full <- contigs_full %>% select(-sample)
                  
                }
                
                if (format %in% c("None", "none")) {
                  contigs_full[contigs_full == "None"] <- ""
                }
                
                seurat_object$tcrEclipse_barcode <- paste(seurat_object@meta.data[[batch]], "_", seurat_object@meta.data[[original_barcode]], sep = "")
                seurat_object$seurat_object_barcode <- rownames(seurat_object@meta.data)
               
                cell_data <- seurat_object@meta.data %>% select(tcrEclipse_barcode, unique(c(batch, donor, group)), seurat_object_barcode)
                cell_data <- cell_data %>% filter(!is.na(index))
                
                if (nrow(cell_data) != n_distinct(cell_data$tcrEclipse_barcode)) {
                  message("Warning: there are likely multiple sequencing runs labeled as one batch. Please relabel the batch column such that each batch is a unique sequncing run. Returning nothing")
                  return()
                }
                
                all <- contigs_full %>% left_join(cell_data, by = c("barcode" = "tcrEclipse_barcode"), keep = FALSE)
                all <- all %>% mutate(barcode = seurat_object_barcode) %>% select(-seurat_object_barcode)
                if("barcode" %in% colnames(seurat_object@meta.data)) {
                  seurat_object@meta.data <- seurat_object@meta.data %>% select(-barcode)
                }
                
                seurat_object@meta.data <- seurat_object@meta.data %>% dplyr::rename("barcode" = seurat_object_barcode)
                message("QC Filtering:")
                message(paste(nrow(all), " intial contigs found before filtering", sep = ""))
                
                ### Renaming barcodes, filtering to only cells in seurat object, removing cells without a CDR3, and renaming QC metrics from cellranger so these contigs aren't thrown out
                all <- all %>% filter(barcode %in% seurat_object$barcode)
                message(paste(nrow(all), " contigs were from true cells", sep = ""))
                if(nrow(all) == 0){
                  message("Check that your barcodes are the original (i.e. likely end with -1 and match the contig files)")
                  return()
                }
                all <- all %>% filter(!cdr3 %in% c("", "None") & !cdr3_nt %in% c("", "None"))
                all <- all %>% filter(chain %in% c("TRA", "TRB", "TRG", "TRD"))
  
                message(paste(sum(all$high_confidence %in% c("true", "True", "TRUE") == FALSE), " contigs were removed for not being high confidence", sep = ""))
                all <- all %>% filter(high_confidence %in% c("true", "True", "TRUE"))
                all <- all %>% mutate(productive = "true",
                                      is_cell = "true",
                                      full_length = "true")
          
                ### Removing cells with abnormally long CDR3s
                all <- all %>% mutate(cdr3_length = str_length(cdr3))
                plot <- all %>% ggplot() + geom_histogram(aes(x = cdr3_length, fill = chain), bins = 30) + geom_vline(aes(xintercept = 30)) + xlab("CDR3 Length (aa)") + ylab("Count")
                print(plot)
                message(paste(sum(all$cdr3_length > 30), " contigs had a CDR3 longer than 30 amino acids and were removed", sep = ""))
                all <- all %>% filter(cdr3_length <= 30)
                
                ### There are TCRs with TRDV and then TRAJ and TRAC. Cellranger calls them as TRD, but they are actually TRA since they pair with a beta chain and have no J segment.
                all <- all %>% mutate(ad_hybrid_tcr = case_when(
                                                                chain %in% c("TRD", "Multi") & str_detect(v_gene, "TRDV") & (str_detect(j_gene, "TRAJ") | str_detect(c_gene, "TRAC")) ~ TRUE,
                                                                TRUE ~ FALSE))
                message(paste(sum(all$ad_hybrid_tcr), " contigs were found that were hybrid TRD/TRA chains. These were converted to TRA chains", sep = ""))
                all <- all %>% mutate(chain = case_when(
                                                        chain %in% c("TRD", "Multi") & str_detect(v_gene, "TRDV") & (str_detect(j_gene, "TRAJ") | str_detect(c_gene, "TRAC")) ~ "TRA",
                                                        TRUE ~ chain))
                message(paste(sum(all$chain %in% c("TRG", "TRD", "Multi")), " true TRG/TRD contigs were found and removed", sep = ""))
                all <- all %>% filter(chain %in% c("TRA", "TRB"))
                
                ### Removing cells with a * at the end of their cdr3. These indicate stop codons, so removing these contigs
                if (sum(str_detect(all$cdr3, "\\*")) > 0) {
                  message(paste("Warning, ",  sum(str_detect(all$cdr3, "\\*")), " contigs have a stop codon at the end of their CDR3 sequence, so removing them", sep = ""))
                  all <- all %>% filter(!str_detect(cdr3, "\\*"))
                }
                
                ### Removing duplicated chains. Some datasets (particularly contig files made with older version of Cellranger) have lots of duplicated rows that have identical CDR3 and gene segments for the same barcode. This messes up future analysis
                duplicates <- all %>% arrange(desc(umis), desc(reads)) %>% group_by(barcode, chain, cdr3) %>% filter(row_number() > 1) %>% ungroup()
                message(paste(nrow(duplicates), " duplicated contigs were removed", sep = ""))
                all <- all %>% arrange(desc(umis), desc(reads)) %>% group_by(barcode, chain, cdr3) %>% filter(row_number() == 1) %>% ungroup()
                
                ### Noting cells with only 1 type of chain (alpha or beta) and that chain isn't seen in any other cell. We don't want these cells to be predicted.
                all2 <- all %>% group_by(barcode) %>% mutate(chainTypes = n_distinct(chain)) %>% ungroup()
                all2 <- all2 %>% group_by(cdr3, chain, .data[[donor]]) %>% mutate(maxChains = max(chainTypes)) %>% ungroup()
                message(paste(sum(all2$maxChains == 1), " orphan contigs were found with no observed pairing.", sep = ""))
                orphanContigs <- all2 %>% filter(maxChains == 1) %>% select(barcode, chain, cdr3)
                
                message(paste("After all filtering, ", nrow(all2), " contigs from ", n_distinct(all2$barcode), " true cells remained", sep = ""))
                
                rm(contigs_full, contigs_list)
                   
                ### Gives better spacing
                message()
                message()
                message()
                
                message("Finding clones with extra chain...")
                message()
                
                donors <- unique(as.vector(all[, donor])[[1]])
                new_contigs <- vector("list", length(donors))
         for (i in 1:length(donors)) {
                  
                ### THIS PART IS FINDING THE COMBO CLONES AND REWRITING THE CONTIG FILES
                message(paste("Running donor ", i, "/", length(donors), " (", donors[i], ")", sep = ""))
                slim <- all2 %>% filter(.data[[donor]] == donors[i])
                
                ### Calling clones with multiple alpha chains and then rewriting the contig files to address this
                fused_a <- findThreeChainClones(contigs = slim, chain = "a")

                ### Calling clones with multiple beta chains and then rewriting the contig files to address this
                new_contigs[[i]] <- findThreeChainClones(contigs = fused_a, chain = "b")
                
         }
                
                contigs_post <- bind_rows(new_contigs)
                
   ### THIS PART IS FEEDING THE EDITED CONTIG FILES INTO VDJdive WHICH CALLS AMBIGUOUS CELLS
                
                message("Running VDJdive clonoStats() for clone EM assignment...")
                ### Writing a .csv file as filtered_contig_annotations.csv
                ### Reading that into VDJdive
                ### Using clonoStats for EM assignment
                ### Extracting output and changing format
                
                if (!write_folder %in% folders) {
                  write.csv(contigs_post, paste(write_folder, "/filtered_contig_annotations.csv", sep = ""), row.names = FALSE)
                }
                contigs <- readVDJcontigs(write_folder)
                vdj <- suppressMessages(clonoStats(contigs, method = "EM", type = "TCR", assignment = TRUE, group = donor))
                
                message("Extracting results from VDJdive...")
                df <- vdj@assignment
                df <- as(df, "RsparseMatrix")
                colnames(df) <- clonoNames(vdj)
                
                rm(vdj)
                
                #### Making a loop that that finds the most likely clone for each cell as well as the second most likely and the ratio of likelihoods between the two
                ### Putting all that data into a df called em
                barcode <- rownames(df)
                clone <- vector(length = nrow(df))
                maxes <- vector(length = nrow(df))
                second_max <- vector(length = nrow(df))
                
                for (i in 1:nrow(df)) {
                  
                  maxes[i] <- max(df[i, ])
                  clone[i] <- names(which.max(df[i, ]))
                  second_max[i] <- max(df[i, colnames(df) != clone[i]])
                  
                }
                
                rm(df)
                
                ### Summarizing the VDJdive EM data and then filtering cells that don't meet QC thresholds
                em <- data.frame(barcode, maxes, clone, second_max)
                em$ratio <- em$maxes / em$second_max
                em2 <- em %>% arrange(maxes) %>% mutate(nm = row_number(),
                                                        assigned_by_VDJdive = maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))
                plot2 <- em2 %>% ggplot() + geom_point(aes(x = nm, y = maxes, color = assigned_by_VDJdive)) +
                                            geom_hline(yintercept = 0.8) +
                                            xlab("Ranked Position") + ylab("Proportion Assigned to Most Likely Clone") + ggtitle("VDJdive EM Assignment QC")
                print(plot2)
                final <- em %>% filter(maxes >= 0.8 | (maxes >= 0.7 & ratio > 5) | (maxes > 0.5 & ratio > 10))  
                
                ### Removing cells where an orphan TCR was assigned. There's no way for the algorithm to accurately predict these so we don't want them.
                orphanContigs <- orphanContigs %>% left_join(final, by = "barcode")
                orphanContigs <- orphanContigs %>% mutate(remove = case_when(
                                                                             is.na(clone) ~ FALSE,
                                                                             chain == "TRA" & str_detect(clone, paste("^", cdr3, " ", sep = "")) ~ TRUE,
                                                                             chain == "TRB" & str_detect(clone, paste(" ", cdr3, "$", sep = "")) ~ TRUE,
                                                                             TRUE ~ FALSE))
                orphanContigs <- orphanContigs %>% filter(remove) %>% select(barcode, remove)
                orphanContigs <- unique(orphanContigs)
                
                final <- final %>% left_join(orphanContigs, by = "barcode")
                message(paste(sum(final$remove, na.rm = T), " cells were removed for being assigned an orphan chain", sep = ""))
                final <- final %>% filter(is.na(remove))
                
                ### Adding the VDJdive and scRepertoire calls to the Seurat object
                seurat_object$tcrEclipseGroup <- seurat_object@meta.data[[group]]
                contig_list <- createHTOContigList(all, 
                                                   seurat_object,
                                                   group.by = "tcrEclipseGroup")
                no_vdjd <- combineTCR(contig_list)
                obj <- combineExpression(no_vdjd, seurat_object, cloneCall = "aa", proportion = TRUE)
                rm(seurat_object)
                
                obj@meta.data <- obj@meta.data %>% rename(scR_CTgene = CTgene,
                                                          scR_CTnt = CTnt,
                                                          scR_CTaa = CTaa,
                                                          scR_CTstrict = CTstrict,
                                                          scR_clonalProportion = clonalProportion,
                                                          scR_clonalFrequency = clonalFrequency,
                                                          scR_cloneSize = cloneSize)
                final <- final %>% mutate(clone2 = str_replace(clone, " ", "_"))
                final <- final %>% mutate(clone2 = str_replace_all(clone2, ";", ":::"))
                final <- final %>% select(barcode,
                                          "VDJdive_clone" = clone2,
                                          "EM_max_prop" = maxes,
                                          "EM_2ndmax_prop" = second_max,
                                          "VDJ_assign_ratio" = ratio)
                
                obj@meta.data <- obj@meta.data %>% left_join(final, by = "barcode")

                
                ### Making CTaa our final clone call which uses VDJdive if available but if not the scRepertoire call
                obj@meta.data <- obj@meta.data %>% mutate(cloneCallSource = case_when(
                                                                                      !is.na(VDJdive_clone) ~ "VDJdive",
                                                                                      is.na(VDJdive_clone) & !is.na(scR_CTaa) ~ "scRepertoire",
                                                                                      TRUE ~ NA))

                obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                           !is.na(VDJdive_clone) ~ VDJdive_clone,
                                                                            is.na(VDJdive_clone) & !is.na(scR_CTaa) ~ scR_CTaa,
                                                                            TRUE ~ NA))

                
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = str_replace_all(CTaa, ";", ":::"))
                
                ### Combining clones. Cells that were ambiguous may be assigned to a a combo clone by VDJdive but only be listed with one chain
                ### Here those cells have their CDR3 renamed as the combo CDR3
                obj@meta.data <- obj@meta.data %>% ungroup() %>% mutate(CTaa_donor = case_when(
                                                                                               !is.na(CTaa) & !is.na(.data[[donor]]) ~ paste(CTaa, .data[[donor]], sep = " "),
                                                                                                TRUE ~ NA))
 
                custom_a <- obj@meta.data %>% filter(!is.na(CTaa) & !is.na(.data[[donor]])) %>% group_by(CTaa, .data[[donor]], cloneCallSource) %>% summarise(count = n()) %>% arrange(desc(count))
                custom_a <- custom_a %>% separate(CTaa, into = c("alpha", "beta"), remove = FALSE, sep = "_")
                custom_a <- custom_a %>% mutate(special_a = case_when(
                                                                      str_detect(alpha, ":::") & cloneCallSource == "VDJdive" ~ alpha,
                                                                      !str_detect(alpha, ":::") | cloneCallSource != "VDJdive" ~ NA))
                              
                custom_a <- custom_a %>% group_by(beta, .data[[donor]]) %>% arrange(desc(count)) %>% mutate(special_aChain = first(special_a)) %>% ungroup()
                custom_a <- custom_a %>% mutate(new_clone = case_when(
                                                                      str_detect(special_aChain, paste("^", alpha, ":::", "|", ":::", alpha, "$", sep = "")) & alpha != special_aChain ~ paste(special_aChain, beta, sep = "_"),
                                                                      TRUE ~ NA))
                
                custom_a <- custom_a %>% unite(col = "CTaa_donor", c(CTaa, .data[[donor]]), sep = " ", remove = FALSE)
                custom_a <- custom_a %>% filter(!is.na(new_clone)) %>% select(CTaa_donor, new_clone)
                obj@meta.data <- obj@meta.data %>% left_join(custom_a, by = "CTaa_donor")
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                           !is.na(new_clone) ~ new_clone,
                                                                            is.na(new_clone) ~ CTaa)) %>% select(-new_clone)
                
                custom_b <- obj@meta.data %>% filter(!is.na(CTaa) & !is.na(.data[[donor]])) %>% group_by(CTaa, .data[[donor]], cloneCallSource) %>% summarise(count = n()) %>% arrange(desc(count))
                custom_b <- custom_b %>% separate(CTaa, into = c("alpha", "beta"), remove = FALSE, sep = "_")
                custom_b <- custom_b %>% mutate(special_b = case_when(
                                                                      str_detect(beta, ":::") & cloneCallSource == "VDJdive" ~ beta,
                                                                      !str_detect(beta, ":::") | cloneCallSource != "VDJdive" ~ NA))
                custom_b <- custom_b %>% group_by(alpha, .data[[donor]]) %>% arrange(desc(count)) %>% mutate(special_bChain = first(special_b)) %>% ungroup()
                custom_b <- custom_b %>% mutate(new_clone = case_when(
                                                                      str_detect(special_bChain, paste("^", beta, ":::", "|", ":::", beta, "$", sep = "")) & beta != special_bChain ~ paste(alpha, special_bChain, sep = "_"),
                                                                      TRUE ~ NA))
                
                custom_b <- custom_b %>% unite(col = "CTaa_donor", c(CTaa, .data[[donor]]), sep = " ", remove = FALSE)
                custom_b <- custom_b %>% filter(!is.na(new_clone)) %>% select(CTaa_donor, new_clone)
                obj@meta.data <- obj@meta.data %>% left_join(custom_b, by = "CTaa_donor")
                obj@meta.data <- obj@meta.data %>% mutate(CTaa = case_when(
                                                                           !is.na(new_clone) ~ new_clone,
                                                                            is.na(new_clone) ~ CTaa)) %>% select(-new_clone, -CTaa_donor)
                
                ### Making columns so that scRepertoire visualizations can be used. Grouping by the group specified in the arguments, not the donor
                obj@meta.data <- obj@meta.data %>% group_by(.data[[group]]) %>% mutate(eclipseTcrGroupSize = sum(!is.na(CTaa))) %>% ungroup()
                obj@meta.data <- obj@meta.data %>% group_by(CTaa, .data[[group]]) %>% mutate(clonalProportion = case_when(
                                                                                                                          !is.na(CTaa) ~ n() / first(eclipseTcrGroupSize),
                                                                                                                           is.na(CTaa) ~ NA),
                                                                                             clonalFrequency = case_when(
                                                                                                                         !is.na(CTaa) ~ n(),
                                                                                                                          is.na(CTaa) ~ NA)) %>% ungroup() %>% select(-eclipseTcrGroupSize)
                obj@meta.data <- obj@meta.data %>% mutate(cloneSize = factor(case_when(
                                                                                       clonalProportion <= 1 & clonalProportion > 0.1 ~ "Hyperexpanded (0.1 < X <= 1)",
                                                                                       clonalProportion <= 0.1 & clonalProportion > 0.01 ~ "Large (0.01 < X <= 0.1)",
                                                                                       clonalProportion <= 0.01 & clonalProportion > 0.001 ~ "Medium (0.001 < X <= 0.01)",
                                                                                       clonalProportion <= 0.001 & clonalProportion > 1e-04 ~ "Small (1e-04 < X <= 0.001)",
                                                                                       clonalProportion <= 1e-04 & clonalProportion > 0 ~ "Rare (0 < X <= 1e-04)",
                                                                                       clonalProportion <= 0 ~ "None ( < X <= 0)",
                                                                                       is.na(clonalProportion) ~ NA)))
                
                
                ### Making a column that notes whether the cells in the clone appear to be MAIT or iNKT based on TCR segment usage
                obj@meta.data <- obj@meta.data %>% mutate(scR_mait_conv = case_when(
                                                                                    str_detect(scR_CTgene, "TRAV1-2\\.TRAJ33") ~ "MAIT Conventional (TRAV1-2 TRAJ33)",
                                                                                    TRUE ~ "Normal"),
                                                          scR_mait_unconv = case_when(
                                                                                      str_detect(scR_CTgene, "TRAV1-2\\.TRAJ12|TRAV1-2\\.TRAJ20") ~ "MAIT Unconventional (TRAV1-2 TRAJ12 or TRAJ20)",
                                                                                      TRUE ~ "Normal"),
                                                          scR_inkt = case_when(
                                                                               str_detect(scR_CTgene, "TRAV10\\.TRAJ18") & str_detect(scR_CTgene, "TRBV25-1") ~ "iNKT",
                                                                               TRUE ~ "Normal"))
                obj@meta.data <- obj@meta.data %>% mutate(scR_unconventional_count = 3 - (scR_mait_conv == "Normal") - (scR_mait_unconv == "Normal") - (scR_inkt == "Normal"))
                obj@meta.data <- obj@meta.data %>% mutate(cell_unconventional_subset = case_when(
                                                                                                 scR_unconventional_count == 0 ~ "Conventional",
                                                                                                 scR_unconventional_count == 1 & scR_mait_conv != "Normal" ~ scR_mait_conv,
                                                                                                 scR_unconventional_count == 1 & scR_mait_unconv != "Normal" ~ scR_mait_unconv,
                                                                                                 scR_unconventional_count == 1 & scR_inkt != "Normal" ~ scR_inkt,
                                                                                                 scR_unconventional_count > 1 ~ "Multiple Unconventional Types"))
                                                                          
                df <- obj@meta.data %>% group_by(CTaa, cell_unconventional_subset, .data[[group]]) %>% summarise(count = n()) %>% ungroup() %>% filter(cell_unconventional_subset != "Conventional") %>% arrange(desc(count))
                
                if (nrow(df) > 0) {
                     obj@meta.data <- obj@meta.data %>% ungroup() %>% mutate(CTaa_group = case_when(
                                                                                                    !is.na(CTaa) & !is.na(.data[[group]]) ~ paste(CTaa, .data[[group]], sep = " "),
                                                                                                     TRUE ~ NA))
                     df <- df %>% group_by(CTaa, .data[[group]]) %>% summarise(clone_unconventional_subset = first(cell_unconventional_subset))
                     df <- df %>% ungroup() %>% mutate(CTaa_group = case_when(
                                                                              !is.na(CTaa) & !is.na(.data[[group]]) ~ paste(CTaa, .data[[group]], sep = " "),
                                                                               TRUE ~ NA))
                     df <- df %>% ungroup() %>% select(CTaa_group, clone_unconventional_subset)
                     obj@meta.data <- obj@meta.data %>% left_join(df, by = "CTaa_group")
                     obj@meta.data <- obj@meta.data %>% mutate(clone_unconventional_subset = case_when(
                                                                                                        is.na(clone_unconventional_subset) ~ "Conventional",
                                                                                                       !is.na(clone_unconventional_subset) ~ clone_unconventional_subset)) %>% select(-CTaa_group)
                     message(paste(sum(obj$clone_unconventional_subset != "Conventional"), " cells were detected in clones that appear to be MAIT or iNKT cells", sep = ""))
                }
                
                obj@meta.data <- obj@meta.data %>% select(-scR_mait_conv, -scR_mait_unconv, -scR_inkt, -scR_unconventional_count, -tcrEclipseGroup, -tcrEclipse_barcode)
                
                obj@meta.data <- as.data.frame(obj@meta.data)
                rownames(obj@meta.data) <- obj$barcode
                
                message(paste(round(sum(obj$cloneCallSource == "VDJdive", na.rm = T) / nrow(obj@meta.data) * 100, 1), "% of cells in the Seurat object were assigned with high confidence by VDJdive clonoStats()", sep = ""))
                message(paste(round(sum(!is.na(obj$CTaa)) / nrow(obj@meta.data) * 100, 1), "% of cells in the Seurat object were assigned in the end (VDJdive or scRepertoire assignment)", sep = ""))
                message("These final clonal annotations are stored in CTaa")
                message("Summary of the top 5 largest clones:")

                big_clones <- obj@meta.data %>% group_by(.data[[donor]]) %>% mutate(donor_count = sum(!is.na(CTaa))) %>% ungroup() %>% 
                                                group_by(CTaa, .data[[donor]]) %>% summarise(Clone_Size = n(), Percent_of_Repertoire = 100 * n() / first(donor_count)) %>%
                                                select(CTaa, Clone_Size, Percent_of_Repertoire) %>% arrange(desc(Clone_Size)) %>% filter(!is.na(CTaa)) %>% head(n = 5)
                print(big_clones)
                
                
                
                return(obj)
          
}
















findThreeChainClones <- function(contigs, chain, mode = "standard") {
  
  if (chain %in% c("alpha", "Alpha", "a", "A", "ALPHA")) {
    
      ### THE FIRST PART HERE IS FINDING CLONES WITH 2 ALPHA CHAINS AND A SHARED BETA CHAIN     
    
              message("findThreeChainClones for alpha chain:")       
    
              ### Uses scRepertoire to call clones at baseline
              tcrs <- combineTCR(contigs)[[1]]
              tcrs <- tcrs %>% separate(col = "CTaa", into = c("alpha", "beta"), sep = "_", remove = FALSE)
              tcrs <- tcrs %>% filter(beta != "NA" & alpha != "NA")
              
              ### Makes data frame called combos that has all cells with 3+ chains and lists the different alphas with each beta
              combos_orig <- tcrs %>% filter(str_count(alpha, ";") < 2  & !str_detect(beta, ";"))
              combos_full <- combos_orig %>% group_by(beta, alpha) %>% summarise(count = n()) %>% ungroup() %>% 
                                             group_by(beta) %>% filter(!((count / sum(count)) < 0.1 & str_detect(alpha, ";"))) %>% ungroup() 
              combos <- combos_full %>% group_by(beta) %>% arrange(desc(count)) %>% 
                                        summarise(alphas = paste(alpha, collapse = " "),
                                                  counts = paste(count, collapse = " "),
                                                  total = sum(count),
                                                  chains = n())
              if (nrow(combos) > 0) {
                
                    combos <- suppressWarnings(combos %>% filter(chains >= 3) %>% select(chains, counts, total, beta, alphas) %>% 
                                         separate(alphas, into = paste("A", seq(1, max(combos$chains)), sep = ""), sep = " ", remove = FALSE))
                    
              } else if (nrow(combos) == 0) {
                return(contigs)
              }
              
              combos <- as.data.frame(combos)
              combos <- combos %>% mutate(combo = NA)
              
        if (nrow(combos) > 0) {
              
              
              ### Scans if there are any clones that have multiple chains and assigns all those cells to 1 clone
              ### If there is 1+ each of A1A2-B, A1-B, and A2-B then it is called as a double. Or also 3+ of A1A2-B.
              for(i in 1:nrow(combos)){
                
                end <- combos[i, 1] + 5
                for(j in 6:end - 1){
                  
                  if (str_detect(combos[i, j], ";") == FALSE){
                    c1 <- combos[i, j]
                    c2 <- NA
                    for(k in (j + 1):end){
                      
                      if (str_detect(combos[i, k], ";") == FALSE){
                        c2 <- combos[i, k]
                        chains <- as.character(combos[i, 6:end])
                        
                        combo1 <- paste(c1, c2, sep = ";")
                        combo2 <- paste(c2, c1, sep = ";")
                        
                        if (combo1 %in% chains & is.na(combos[i, "combo"])){
                          combos[i, "combo"] <- combo1
                        } 
                        else if (combo2 %in% chains & is.na(combos[i, "combo"])){
                          combos[i, "combo"] <- combo2
                        }
                        
                      }
                    }
                    
                  }
                  
                }
              }
              
        }     
              if (mode == "showSpecial") {
                combo_table <- combos
              }
              combos <- combos %>% select(beta, combo)
              combos <- combos %>% filter(!is.na(combo))
              
              ### Uses combos_full to find more clones that may not have all 3 conditions met, but still have cells with 2 alpha that are >= 10% and also 3+ in count
              combos2 <- combos_full %>% filter(str_detect(alpha, ";") & count >= 3)
              combos2 <- combos2 %>% arrange(desc(count)) %>% group_by(beta) %>% filter(row_number() == 1) %>% ungroup()
              combos2 <- combos2 %>% select(beta, "combo" = alpha) %>% filter(!beta %in% combos$beta)
                                                              
              combos <- rbind(combos, combos2)
              
              
              #### This finds the beta chain that is most common for each alpha chain. Note: this is not just 3 TCR clones, this is for all clones
              #### We want to be able to filter so that for each 2A 1B clone, each A must be most commonly seen with that B.
              pairings <- combos_orig %>% group_by(alpha, beta) %>% summarise(count = n())
              suppressWarnings(pairings <- pairings %>% separate(alpha, into = c("a1", "a2"), sep = ";"))
              pairings <- pairings %>% pivot_longer(cols = c(a1, a2), names_to = "original", values_to = "alpha") %>% 
                                       select(-original) %>% filter(!is.na(alpha))
              pairings <- pairings %>% group_by(alpha, beta) %>% summarise(Count = sum(count))
              pairings <- pairings %>% arrange(desc(Count)) %>% group_by(alpha) %>% mutate(rank = min_rank(desc(Count))) %>% 
                                                                group_by(alpha, rank) %>% filter(sum(rank) == 1) %>% ungroup() %>% select(-rank)
              
              ### This part compares the most common beta for each alpha and sees if for every 3 TCR clone both alphas are most commonly with the beta
              ### If this isn't the case, they are filtered out here
              pairings2 <- pairings %>% left_join(combos, by = "beta")
              pairings2 <- pairings2 %>% filter(str_detect(combo, paste("^", alpha, ";|;", alpha, "$", sep = "")))
              pairings2 <- pairings2 %>% group_by(beta) %>% mutate(combos = n_distinct(combo)) %>% 
                                         group_by(beta, combo) %>% summarise(count = n(), combos = first(combos)) %>% 
                                         filter(count == 2)
              
              if (nrow(pairings2) > 0) {
                if (max(pairings2$combos) > 1) {
                  message("Warning! Muliple combos with a beta remain")
                }
              }
              combos <- pairings2 %>% select(beta, combo) %>% ungroup()
              if (mode == "showSpecial") {
                combos3 <- combos
              }
              
              ### Makes a vector special that has all of the alpha CDR3s of chains that have 2 alphas
              special <- as.character((combos %>% filter(!is.na(combo)) %>% select(combo))[[1]])
              message(paste("There were ", length(special), " special combo clones detected", sep = ""))
              
              ### Makes a data frame tcr_mod where combos (beta CDR3 + combo alpha CDR3s) is joined to to tcrs based off of shared beta CDR3
              ### This allows us to know which barcodes need to change.
              ### Then, if the current alpha CDR3 is detected within the combo alpha CDR3 that has the same beta CDR3 then it will be assigned as then new alpha CDR3
              ### Lastly, only the barcode and combo alpha CDR3 are retained and cells that aren't getting assigned with a combo alpha CDR3 are removed
              tcr_mod <- tcrs %>% left_join(combos, by = "beta")
  
              tcr_mod <- tcr_mod %>% mutate(alpha_mod = case_when(
                                                                  str_detect(combo, paste("^", alpha, ";", "|", ";", alpha, "$", sep = "")) | combo == alpha ~ combo,
                                                                  TRUE ~ NA))
              
              tcr_mod <- tcr_mod %>% select(barcode, alpha_mod) %>% filter(!is.na(alpha_mod))
              
              ### A df called new_contigs is made that takes the original contigs data and left joins in tcr_mod by shared barcode
              ### Then, if the row is an alpha chain and has a new alpha combo it is called as being a combo
              new_contigs <- contigs %>% left_join(tcr_mod, by = "barcode")
              new_contigs <- new_contigs %>% mutate(new_cdr3 = case_when(
                                                                         !is.na(alpha_mod) & chain == "TRA" ~ alpha_mod))
              
              ### new_contigs is then filtered by reads/UMI and then, for newly assigned clones only, only the chain with the most amount of UMIs/reads for a given chain/barcode pair is retained.
              ### This filtering is necessary because essentially each row = 1 chain, so if you assign 1 row to a combo (2 chains), then you would be duplicating the data.
              ### Removing duplicates rows ensures that downstream the clone will not get duplicated.
              new_contigs <- new_contigs %>% group_by(barcode, chain) %>% arrange(desc(umis), desc(reads)) %>% mutate(row_num = row_number()) %>% ungroup()
              new_contigs <- new_contigs %>% filter(!(!is.na(new_cdr3) & row_num > 1))
              new_contigs <- new_contigs %>% mutate(cdr3 = case_when(
                                                                     !is.na(new_cdr3) ~ new_cdr3,
                                                                     is.na(new_cdr3) ~ cdr3))
              
              ### Makes a new df ab with the modified contigs.
              ### Feeds the first 32 columns of ab into scRepertoire so show find new clone calls. This makes a new df called post
              ab <- new_contigs
              message(paste("This led to ", sum(str_detect(new_contigs$cdr3, ";") & new_contigs$chain == "TRA"), " cells being called as in a combo alpha clone", sep = ""))
              post <- combineTCR(ab)[[1]]
              
              if (mode == "showSpecial"){
                message("Post-filtering combo alpha clones:")
                print(combos3)
                return(combo_table)
              }
              
              
              
              
              
     ### THIS SECOND PART IS FINDING ALL THE CELLS WITH AN ALPHA AND A MISSING BETA THAT ARE IN CLONES WITH COMBO ALPHA CHAINS
     ### THEN THESE CELLS ARE ASSIGNED TO THEIR COMBO CLONE SO THEY ARE CALLED DOWNSTREAM AS 1 CLONE
          
          if (length(special) > 0) {
              ### Makes a df vdj which summarizes the new clone sizes
              ### Then, if the alpha isn't missing, calls the number of combo CDR3s it could be from. 
              vdj <- post %>% group_by(CTaa) %>% summarise(count = n()) %>% arrange(desc(count)) %>% ungroup()
              vdj <- vdj %>% separate(CTaa, remove = FALSE, into = c("alpha", "beta"), sep = "_")
              vdj <- vdj %>% rowwise() %>% mutate(special_count = case_when(
                                                                            alpha != "NA" ~ sum(str_detect(special, paste("^", alpha, ";", "|", ";", alpha, "$", sep = "")) | special == alpha),
                                                                            alpha == "NA" ~ 0))
              
              
              ### Makes a new column called special_tcr that has the combo alpha that the alpha chain belongs to.
              vdj$special_tcr <- NA
              for (i in 1:nrow(vdj)) {
                
                if (vdj$special_count[i] == 1) {
                   vdj$special_tcr[i] <- special[which(str_detect(special, paste("^", vdj$alpha[i], ";", "|", ";", vdj$alpha[i], "$", sep = "")) | special == vdj$alpha[i])]
                }
                
              }
              
              vdj_bad <- vdj %>% filter(special_count > 1 & beta == "NA")
              if (nrow(vdj_bad) > 0) {
                message(paste("There are some alpha chains that are in multiple combo clones. Because of this, ", sum(vdj_bad$count), " cells with these alphas and no beta won't be able to be assinged to a combo clone", sep = ""))
              }
              
              ### Makes a df call special_cells that calls the cells with beta NA and that belong to combo alpha clones
              ### Begins by only retaining cells that could be called and then calculates the percentage of cells with each beta and the same combo alpha.
              ### You can think of these percentages as the percentage likelihood of each beta for cells with a NA beta and the same alpha
              special_cells <- vdj %>% filter(!is.na(special_tcr)) %>% mutate(mod_count = case_when(
                                                                                                    beta == "NA" ~ 0,
                                                                                                    beta != "NA" ~ count))
              special_cells <- special_cells %>% group_by(special_tcr) %>% mutate(group_count = sum(mod_count)) %>% ungroup()
              special_cells <- special_cells %>% mutate(group_pct = 100 * mod_count / group_count)
              
              ### Next we try to assign for NA beta cells. This starts with making a column new_alpha which is only there for rows where the beta is 80%+ likely.
              ### We then add that same alpha (that pairs with the likely beta) to a new column new_alpha2 for all cells in the same alpha group.
              ### Lastly, makes a new column new_alpha3 that, for cells with no beta and a paired alpha, has the combo alpha chain CDR3
              special_cells <- special_cells %>% mutate(new_alpha = case_when(
                                                                              group_pct >= 80 ~ alpha,
                                                                              group_pct < 80 ~ NA))
              special_cells <- special_cells %>% arrange(desc(group_pct)) %>% group_by(special_tcr) %>% mutate(new_alpha2 = dplyr::first(new_alpha)) %>% ungroup()
              special_cells <- special_cells %>% rowwise() %>% mutate(new_alpha3 = case_when(
                                                                                             is.na(new_alpha2) ~ NA,
                                                                                             beta == "NA" & (str_detect(new_alpha2, paste("^", alpha, ";", "|", ";", alpha, "$", sep = "")) | new_alpha2 == alpha) ~ new_alpha2,
                                                                                             TRUE ~ NA)) %>% ungroup()
              message(paste(sum(!is.na(special_cells$new_alpha)), " of the ", n_distinct(special_cells$special_tcr), " clones had cells with NA alphas assigned to them", sep = ""))
              
              ### Makes a df special_cells2 that selects only the CDR3s and the new alpha combo CDR3 and then filters rows with no new assignment
              special_cells2 <- special_cells %>% dplyr::select(CTaa, new_alpha3) %>% filter(!is.na(new_alpha3))
              
              ### Joins special_cells with post (the scRepertoire output with combo alpha chains) by CDR3s to make vdj2
              ### This allow us to know which barcodes need to change.
              ### Barcodes that don't need to change are then removed
              vdj2 <- post %>% left_join(special_cells2, by = "CTaa")
              vdj2 <- vdj2 %>% select(barcode, new_alpha3) %>% filter(!is.na(new_alpha3))
              message(paste("This led to ", nrow(vdj2), " cells being reassigned", sep = ""))
              
              ### ab (contigs after joining alpha chains) is joined to vdj2 to make ab2
              ### ab2 is filtered to retain only the row of that barcode/chain combo with the most UMIs/reads. This is only for alpha chains that are changing
              ### This is necessary because, again, not doing so would result in some duplications down stream
              ### Lastly, the alpha cdr3 for these combo alpha - missing beta cells are changed
              ab2 <- ab %>% left_join(vdj2, by = "barcode")
              ab2 <- ab2 %>% group_by(barcode, chain) %>% arrange(desc(umis), desc(reads)) %>% mutate(row_num = row_number()) %>% ungroup()
              ab2 <- ab2 %>% filter(!(!is.na(new_alpha3) & row_num > 1 & chain == "TRA"))
              ab2 <- ab2 %>% mutate(cdr3 = case_when(
                                                     chain == "TRA" & !is.na(new_alpha3) ~ new_alpha3,
                                                     TRUE ~ cdr3))
          }
          else if (length(special) == 0) {
            ab2 <- ab
          }
              message(paste("In total, this means ", sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRA"), " cells were altered", sep = ""))
              message()
              message()
              message()
              
              return(ab2)
  }
  else if (chain %in% c("beta", "Beta", "b", "B", "BETA")) {
    
    
    ### THE FIRST PART HERE IS FINDING CLONES WITH 2 BETA CHAINS AND A SHARED ALPHA CHAIN     
    
              message("findThreeChainClones for beta chain:")          
    
              ### Uses scRepertoire to call clones at baseline
              tcrs <- combineTCR(contigs)[[1]]
              tcrs <- tcrs %>% separate(col = "CTaa", into = c("alpha", "beta"), sep = "_", remove = FALSE)
              tcrs <- tcrs %>% filter(beta != "NA" & alpha != "NA")
              
              ### Makes data frame called combos that has all cells with 3+ chains and lists the different betas with each alpha
              combos_orig <- tcrs %>% filter(str_count(beta, ";") < 2 & !str_detect(alpha, ";"))
              combos_full <- combos_orig %>% group_by(beta, alpha) %>% summarise(count = n()) %>% ungroup() %>% 
                                             group_by(alpha) %>% filter(!((count / sum(count)) < 0.1 & str_detect(beta, ";"))) %>% ungroup()
              combos <- combos_full %>% group_by(alpha) %>% arrange(desc(count)) %>% 
                                        summarise(betas = paste(beta, collapse = " "),
                                                  counts = paste(count, collapse = " "),
                                                  total = sum(count),
                                                  chains = n())
              if (nrow(combos) > 0) {
                
                  combos <- suppressWarnings(combos %>% filter(chains >= 3) %>% select(chains, counts, total, alpha, betas) %>% 
                                             separate(betas, into = paste("A", seq(1, max(combos$chains)), sep = ""), sep = " ", remove = FALSE))
              
              } else if (nrow(combos) == 0) {
                    return(contigs)
              }
              
              combos <- as.data.frame(combos)
              combos <- combos %>% mutate(combo = NA)
              
              
        if (nrow(combos) > 0) {
              
              ### Scans if there are any clones that have multiple chains and assigns all those cells to 1 clone
              ### If there is 1+ each of A-B1B2, A-B1, and A-B2 then it is called as a double. Or also 3+ of A-B1B2.
              for(i in 1:nrow(combos)){
                
                end <- combos[i, 1] + 5
                for(j in 6:end - 1){
                  
                  if (str_detect(combos[i, j], ";") == FALSE){
                    c1 <- combos[i, j]
                    c2 <- NA
                    for(k in (j + 1):end){
                      
                      if (str_detect(combos[i, k], ";") == FALSE){
                        c2 <- combos[i, k]
                        chains <- as.character(combos[i, 6:end])
                        
                        combo1 <- paste(c1, c2, sep = ";")
                        combo2 <- paste(c2, c1, sep = ";")
                        
                        if(combo1 %in% chains & is.na(combos[i, "combo"])){
                          combos[i, "combo"] <- combo1
                        } 
                        else if (combo2 %in% chains & is.na(combos[i, "combo"])){
                          combos[i, "combo"] <- combo2
                        }
                        
                      }
                    }
                    
                  }
                  
                }
              }
        
        }
              if (mode == "showSpecial") {
                combo_table <- combos
              }
              combos <- combos %>% select(alpha, combo)
              combos <- combos %>% filter(!is.na(combo))
              
              ### Uses combos_full to find more clones that may not have all 3 conditions met, but still have cells with 2 beta that are >= 10% and also 3+ in count
              combos2 <- combos_full %>% filter(str_detect(beta, ";") & count >= 3)
              combos2 <- combos2 %>% arrange(desc(count)) %>% group_by(alpha) %>% filter(row_number() == 1) %>% ungroup()
              combos2 <- combos2 %>% select(alpha, "combo" = beta) %>% filter(!alpha %in% combos$alpha)
                                                              
              combos <- rbind(combos, combos2)
              
              #### This finds the alpha chain that is most common for each beta chain. Note: this is not just 3 TCR clones, this is for all clones
              #### We want to be able to filter so that for each 2B 1A clone, each B must be most commonly seen with that A.
              pairings <- combos_orig %>% group_by(alpha, beta) %>% summarise(count = n())
              suppressWarnings(pairings <- pairings %>% separate(beta, into = c("b1", "b2"), sep = ";"))
              pairings <- pairings %>% pivot_longer(cols = c(b1, b2), names_to = "original", values_to = "beta") %>% 
                                       select(-original) %>% filter(!is.na(beta))
              pairings <- pairings %>% group_by(alpha, beta) %>% summarise(Count = sum(count))
              pairings <- pairings %>% arrange(desc(Count)) %>% group_by(beta) %>% mutate(rank = min_rank(desc(Count))) %>% 
                                                                group_by(beta, rank) %>% filter(sum(rank) == 1) %>% ungroup() %>% select(-rank)
              
              ### This part compares the most common alpha for each beta and sees if for every 3 TCR clone both betas are most commonly with the alpha
              ### If this isn't the case, they are filtered out here
              pairings2 <- pairings %>% left_join(combos, by = "alpha")
              pairings2 <- pairings2 %>% filter(str_detect(combo, paste("^", beta, ";|;", beta, "$", sep = "")))
              pairings2 <- pairings2 %>% group_by(alpha) %>% mutate(combos = n_distinct(combo)) %>% 
                                         group_by(alpha, combo) %>% summarise(count = n(), combos = first(combos)) %>% 
                                         filter(count == 2)
              
              if (nrow(pairings2) > 0) {
                if (max(pairings2$combos) > 1) {
                  message("Warning! Muliple combos with an alpha remain")
                }
              }
              combos <- pairings2 %>% select(alpha, combo) %>% ungroup()
              if (mode == "showSpecial") {
                combos3 <- combos
              }
              
              ### Makes a vector special that has all of the beta CDR3s of chains that have 2 betas
              special <- as.character((combos %>% filter(!is.na(combo)) %>% select(combo))[[1]])
              message(paste("There were ", length(special), " special combo clones detected", sep = ""))
              
              ### Makes a data frame tcr_mod where combos (alpha CDR3 + combo beta CDR3s) is joined to to tcrs based off of shared alpha CDR3
              ### This allows us to know which barcodes need to change.
              ### Then, if the current beta CDR3 is detected within the combo beta CDR3 that has the same alpha CDR3 then it will be assigned as then new beta CDR3
              ### Lastly, only the barcode and combo beta CDR3 are retained and cells that aren't getting assigned with a combo beta CDR3 are removed
              tcr_mod <- tcrs %>% left_join(combos, by = "alpha")
              tcr_mod <- tcr_mod %>% mutate(beta_mod = case_when(
                                                                  str_detect(combo, paste("^", beta, ";", "|", ";", beta, "$", sep = "")) | combo == beta ~ combo,  
                                                                  TRUE ~ NA))
              
              tcr_mod <- tcr_mod %>% select(barcode, beta_mod) %>% filter(!is.na(beta_mod))
              
              ### A df called new_contigs is made that takes the original contigs data and left joins in tcr_mod by shared barcode
              ### Then, if the row is a beta chain and has a new beta combo it is called as being a combo
              new_contigs <- contigs %>% left_join(tcr_mod, by = "barcode")
              new_contigs <- new_contigs %>% mutate(new_cdr3 = case_when(
                                                                         !is.na(beta_mod) & chain == "TRB" ~ beta_mod))
              
              ### new_contigs is then filtered by reads/UMI and then, for newly assigned clones only, only the chain with the most amount of UMIs/reads for a given chain/barcode pair is retained.
              ### This filtering is necessary because essentially each row = 1 chain, so if you assign 1 row to a combo (2 chains), then you would be duplicating the data.
              ### Removing duplicates rows ensures that downstream the clone will not get duplicated.
              new_contigs <- new_contigs %>% group_by(barcode, chain) %>% arrange(desc(umis), desc(reads)) %>% mutate(row_num = row_number()) %>% ungroup()
              new_contigs <- new_contigs %>% filter(!(!is.na(new_cdr3) & row_num > 1))
              new_contigs <- new_contigs %>% mutate(cdr3 = case_when(
                                                                     !is.na(new_cdr3) ~ new_cdr3,
                                                                     is.na(new_cdr3) ~ cdr3))
              
              ### Makes a new df ab with the modified contigs.
              ### Feeds the first 32 columns of ab into scRepertoire so show find new clone calls. This makes a new df called post
              ab <- new_contigs
              message(paste("This led to ", sum(str_detect(new_contigs$cdr3, ";") & new_contigs$chain == "TRB"), " cells being called as in a combo beta clone", sep = ""))
              post <- combineTCR(ab)[[1]]
              
              if (mode == "showSpecial"){
                message("Post-filtering combo beta clones:")
                print(combos3)
                return(combo_table)
              }
              
              
              
     ### THIS SECOND PART IS FINDING ALL THE CELLS WITH AN BETA AND A MISSING ALPHA THAT ARE IN CLONES WITH COMBO BETA CHAINS
     ### THEN THESE CELLS ARE ASSIGNED TO THEIR COMBO CLONE SO THEY ARE CALLED DOWNSTREAM AS 1 CLONE
            if (length(special) > 0) {  
              ### Makes a df vdj which summarizes the new clone sizes
              ### Then, if the beta isn't missing, calls the number of combo CDR3s it could be from.
              vdj <- post %>% group_by(CTaa) %>% summarise(count = n()) %>% arrange(desc(count)) %>% ungroup()
              vdj <- vdj %>% separate(CTaa, remove = FALSE, into = c("alpha", "beta"), sep = "_")
              vdj <- vdj %>% rowwise() %>% mutate(special_count = case_when(
                                                                            beta != "NA" ~ sum(str_detect(special, paste("^", beta, ";", "|", ";", beta, "$", sep = "")) | special == beta),
                                                                            beta == "NA" ~ 0))
              
              ### Makes a new column called special_tcr that has the combo beta that the beta chain belongs to.
              vdj$special_tcr <- NA
              for (i in 1:nrow(vdj)) {
                
                if (vdj$special_count[i] == 1) {
                  vdj$special_tcr[i] <- special[which(str_detect(special, paste("^", vdj$beta[i], ";", "|", ";", vdj$beta[i], "$", sep = "")) | special == vdj$beta[i])]
                }
                
              }
              
              vdj_bad <- vdj %>% filter(special_count > 1 & alpha == "NA")
              if (nrow(vdj_bad) > 0) {
                message(paste("There are some beta chains that are in multiple combo clones. Because of this, ", sum(vdj_bad$count), " cells with these betas and no alpha won't be able to be assinged to a combo clone", sep = ""))
              }
              
              ### Makes a df call special_cells that calls the cells with alpha NA and that belong to combo beta clones
              ### Begins by only retaining cells that could be called and then calculates the percentage of cells with each alpha and the same combo beta
              ### You can think of these percentages as the percentage likelihood of each alpha for cells with a NA alpha and the same beta
              special_cells <- vdj %>% filter(!is.na(special_tcr)) %>% mutate(mod_count = case_when(
                                                                                                    alpha == "NA" ~ 0,
                                                                                                    alpha != "NA" ~ count))
              special_cells <- special_cells %>% group_by(special_tcr) %>% mutate(group_count = sum(mod_count)) %>% ungroup()
              special_cells <- special_cells %>% mutate(group_pct = 100 * mod_count / group_count)
              
              ### Next we try to assign for NA alpha cells. This starts with making a column new_beta which is only there for rows where the alpha is 80%+ likely.
              ### We then add that same beta (that pairs with the likely alpha) to a new column new_beta2 for all cells in the same beta group. 
              ### Lastly, makes a new column new_beta3 that, for cells with no alpha and a paired beta, has the combo beta chain CDR3
              special_cells <- special_cells %>% mutate(new_beta = case_when(
                                                                             group_pct >= 80 ~ beta,
                                                                             group_pct < 80 ~ NA))
              special_cells <- special_cells %>% arrange(desc(group_pct)) %>% group_by(special_tcr) %>% mutate(new_beta2 = dplyr::first(new_beta)) %>% ungroup()
              special_cells <- special_cells %>% rowwise() %>% mutate(new_beta3 = case_when(
                                                                                            is.na(new_beta2) ~ NA,
                                                                                            alpha == "NA" & (str_detect(new_beta2, paste("^", beta, ";", "|", ";", beta, "$", sep = "")) | new_beta2 == beta) ~ new_beta2,
                                                                                            TRUE ~ NA)) %>% ungroup()
              message(paste(sum(!is.na(special_cells$new_beta)), " of the ", n_distinct(special_cells$special_tcr), " clones had cells with NA alphas assigned to them", sep = ""))
              
              ### Makes a df special_cells2 that selects only the CDR3s and the new beta combo CDR3 and then filters rows with no new assignment
              special_cells2 <- special_cells %>% dplyr::select(CTaa, new_beta3) %>% filter(!is.na(new_beta3))
              
              ### Joins special_cells with post (the scRepertoire output with combo beta chains) by CDR3s to make vdj2
              ### This allow us to know which barcodes need to change.
              ### Barcodes that don't need to change are then removed
              vdj2 <- post %>% left_join(special_cells2, by = "CTaa")
              vdj2 <- vdj2 %>% select(barcode, new_beta3) %>% filter(!is.na(new_beta3))
              message(paste("This led to ", nrow(vdj2), " cells being reassigned", sep = ""))
              
              ### ab (contigs after joining beta chains) is joined to vdj2 to make ab2
              ### ab2 is filtered to retain only the row of that barcode/chain combo with the most UMIs/reads. This is only for beta chains that are changing
              ### This is necessary because, again, not doing so would result in some duplications down stream
              ### Lastly, the beta cdr3 for these combo beta - missing alpha cells are changed
              ab2 <- ab %>% left_join(vdj2, by = "barcode")
              ab2 <- ab2 %>% group_by(barcode, chain) %>% arrange(desc(umis), desc(reads)) %>% mutate(row_num = row_number()) %>% ungroup()
              ab2 <- ab2 %>% filter(!(!is.na(new_beta3) & row_num > 1 & chain == "TRB"))
              ab2 <- ab2 %>% mutate(cdr3 = case_when(
                                                     chain == "TRB" & !is.na(new_beta3) ~ new_beta3,
                                                     TRUE ~ cdr3))
            }
            else if (length(special) == 0) {
              ab2 <- ab
            }
              
              message(paste("In total, this means ", sum(str_detect(ab2$cdr3, ";") & ab2$chain == "TRB"), " cells were altered", sep = ""))
              message()
              message()
              message()
        
              return(ab2)
              
  }
  else {
    print("Unknown chain, returning original contigs")
    return(contigs)
  }
  
  
}
  
  








tcrDoubletDetect <- function(object, filter = TRUE, singleChainLimit = 2, totalChainLimit = -1) {
      
      ### Makes column in meta data that has the number of alpha/beta chains both originally and after ECLIPSE.
      object@meta.data <- object@meta.data %>% separate(scR_CTaa, remove = FALSE, sep = "_", into = c("scr_alpha", "scr_beta"))
      object@meta.data <- object@meta.data %>% separate(CTaa, remove = FALSE, sep = "_", into = c("ctaa_alpha", "ctaa_beta"))
      
      object@meta.data <- object@meta.data %>% mutate(scr_alphaCount = case_when(
                                                                                 scr_alpha == "NA" | is.na(scr_alpha) ~ 0,
                                                                                 scr_alpha != "NA" ~ str_count(scr_alpha, ";") + 1),
                                                      scr_betaCount = case_when(
                                                                                scr_beta == "NA" | is.na(scr_beta) ~ 0,
                                                                                scr_beta != "NA" ~ str_count(scr_beta, ";") + 1))
      
      plot <- VlnPlot(object, group.by = "scr_alphaCount", features = "nCount_RNA") + xlab("Number of Alpha Chains")
      suppressWarnings(print(plot))
      plot2 <- VlnPlot(object, group.by = "scr_betaCount", features = "nCount_RNA") + xlab("Number of Beta Chains")
      suppressWarnings(print(plot2))
      
      ### This adjusts for if the clones are called by my pipeline as supposed to be having 2 of a chain. Still if something has 2 that doesn't mean it's 100% a doublet because not everything can be called by my pipeline
      object@meta.data <- object@meta.data %>% mutate(adj_scr_alphaCount = case_when(
                                                                                     scr_alpha == str_replace(ctaa_alpha, ":::", ";") & str_detect(ctaa_alpha, ":::") & cloneCallSource == "VDJdive" ~ scr_alphaCount - 1,
                                                                                     TRUE ~ scr_alphaCount),
                                                      adj_scr_betaCount = case_when(
                                                                                    scr_beta == str_replace(ctaa_beta, ":::", ";") & str_detect(ctaa_beta, ":::") & cloneCallSource == "VDJdive" ~ scr_betaCount - 1,
                                                                                    TRUE ~ scr_betaCount))

      ### Anything with more than the desired number of chains is removed
      cell_count <- nrow(object@meta.data)
      
      if (filter == TRUE) {
        
           message(paste(sum(object$scr_alphaCount > singleChainLimit | object$scr_betaCount > singleChainLimit), " cells (", 
                         round(sum(object$scr_alphaCount > singleChainLimit | object$scr_betaCount > singleChainLimit) / cell_count * 100, 1), 
                         "% of total) were removed due to having >", singleChainLimit, " alpha or beta chains originally", sep = ""))
           object <- subset(object, scr_alphaCount <= singleChainLimit & scr_betaCount <= singleChainLimit)
           
           if (singleChainLimit < 2) {
             message("")
             message(paste("Warning:\n", "Filtering cells with more than 1 alpha or beta chain is not recommended as some clones express 2 alpha and/or beta chains", sep = ""))
             message("")
           }
      }
      
      
      ### Finds the total amount of chains, both raw and adjusted for clones called as being 2 of a chain in my pipeline
      object@meta.data <- object@meta.data %>% mutate(scr_count = scr_alphaCount + scr_betaCount,
                                                      adj_scr_count = adj_scr_alphaCount + adj_scr_betaCount)
      
      if (filter == TRUE) {
        
          if (totalChainLimit >= 0) {
              message(paste(sum(object$scr_count > totalChainLimit), " cells (",
                            round(sum(object$scr_count > totalChainLimit) / cell_count * 100, 1),
                            "% of total) were removed due to having >", totalChainLimit, " total chains originally", sep = ""))
              object <- subset(object, scr_count <= totalChainLimit)
          }
          
          if (totalChainLimit <= 2 & totalChainLimit > -1) {
              message("")
              message(paste("Warning:\n", "Filtering all cells with more than 2 total chains is not recommended as some clones express 2 alpha and/or beta chains", sep = ""))
          }
            
      }
      
      plot3 <- VlnPlot(object, group.by = "scr_count", features = "nCount_RNA") + xlab(paste("Number of TCR Chains", "\n", "originally detected", sep = ""))
      suppressWarnings(print(plot3))
      plot4 <- VlnPlot(object, group.by = "adj_scr_count", features = "nCount_RNA") + xlab(paste("Number of TCR Chains", "\n", "(adjusted to account for clones", "\n" ,"that express 2 of one chain type)", sep = ""))
      suppressWarnings(print(plot4))
      
      object@meta.data <- object@meta.data %>% select(-scr_alpha, -scr_beta, -ctaa_alpha, -ctaa_beta, -scr_alphaCount, -scr_betaCount, -adj_scr_alphaCount, -adj_scr_betaCount)
      return(object)
      
}

