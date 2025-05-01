
# For debugging
# data <- Data
# row_start <- RowStartForFeatureTableData - 1
# adduct_col <- end_col - 5
# annotation_col <- end_col - 8
# window <- RT_Window
# RT_col <- RTColumn
# AppendFrags_From_LibHits

remove_duplicates<-function(data,row_start,adduct_col,annotation_col,window,RT_col){
  #remove any commas in annotations
  data[,annotation_col]<-gsub(",","",data[,annotation_col])
  # create a vector with just the ranking score (1_, 2_, 3_, or 4_)
  ranking_class<-substr(as.vector(data[,annotation_col]),start=1,stop=2)
  # Identify and remove rows with 4_ or 5_
  not_MSMS<-sapply(ranking_class,function(x){
    if(x=="4_"||x=="5_"){
      return(FALSE)
    }else{
      return(TRUE)
    }
  })
  if (length(which(not_MSMS)) == 1) {
    message("only one compound ranking 1_, 2_, or 3_, stopping remove_duplicates function --> No CombinedIDed output") # PJS 1/18/2023 - having only one causes lots of errors because single-row matrices become vectors
    return(data[0,])
  }
  data<-data[which(not_MSMS),]
  # Calculate median peak area or peak height
  ## Potential Place of error since text strings can change for samples
  peak_cols<-which(grepl("Peak.area|Peak.height|.mzXML|Neg|Pos|neg|pos",colnames(data)))
  
  if (length(peak_cols)>1) {
    m <- as.matrix(apply(data[,peak_cols], 2, as.numeric))
    median_intensities<-apply(m,1,median, na.rm = TRUE)
  }
  
  # A specific case where there is only one sample
  if (length(peak_cols)==1){
    median_intensities<-as.numeric(data[,peak_cols])
  }
  
  # Take top ranked lipid for each row
  top_ranked_lipid<-sapply(as.vector(data[,annotation_col]), function(x){
    strsplit(x,split = " | ")[[1]][1]
  })
  top_ranked_lipid<-as.vector(top_ranked_lipid)
  
  # Remove first 2 digits from lipid names
  top_ranked_lipid<-substring(top_ranked_lipid,3)
  
  # Remove adducts from names
  top_ranked_lipid<-gsub("\\+.*| .*|\\-H.*|\\-2H.*","",top_ranked_lipid)
  # top_ranked_lipid <- top_ranked_lipid[-1]
  names(median_intensities)<-top_ranked_lipid
  # Identify sodiated rows
  not_sodiated<-sapply(as.vector(data[,adduct_col]),function(x){
    # Yang: if x is missing return TRUE, need to verify with Jeremy
    if (is.na(x)){
      return(TRUE)
    }
    else{
      if(x=="[M+Na]+"){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
  })
  
  # Identify row # of most abundant duplicate for each group. Make a string of adducts within the window.
  keepers<-c()
  adducts<-c()
  for(i in unique(top_ranked_lipid)){
    # print(i)
    rows = which(top_ranked_lipid==i)
    if(length(rows)==1){
      keepers<-c(keepers,rows)
      adducts<-c(adducts,as.vector(data[,adduct_col])[rows])
    }else{
      candidates<-median_intensities[rows]
      ranks<-order(candidates, decreasing = TRUE)
      
      # Prefer negative adducts if they're there
      all_adducts<-as.vector(data[,adduct_col])[rows]
      ### EDITED JPK (added new adducts for negative mode)
      if(grepl("+", paste(all_adducts,sep = "",collapse = "|")) && ("[M-H]-" %in% all_adducts || "[M-2H]-" %in% all_adducts || "[M+C2H3O2]-" %in% all_adducts || "[M+HCO2]-" %in% all_adducts || "[M+Cl]-" %in% all_adducts)){
        negative_adduct_indices<-c(match("[M-H]-",all_adducts),match("[M-2H]-",all_adducts),match("[M-2H]-",all_adducts),match("[M+C2H3O2]-",all_adducts),match("[M+HCO2]-",all_adducts),match("[M+Cl]-",all_adducts))
        negative_adduct_indices<-negative_adduct_indices[!is.na(negative_adduct_indices)]
        candidates_neg<-candidates[negative_adduct_indices]
        best_abundance<-max(candidates_neg)
        best_index<-intersect(which(median_intensities==best_abundance),rows)[1]
      } else {
        best_abundance<-max(candidates)
        best_index<-intersect(which(median_intensities==best_abundance),rows)[1]
      }
      
      # Use the runner-up if the top ranked adduct is sodium
      if(not_sodiated[best_index]==FALSE){
        best_index<-intersect(which(median_intensities==candidates[match(2,ranks)]),rows)[1]
      }
      # Find rows within RT window
      rows<-rows[ranks]
      
      in_window<-intersect(which(as.numeric(as.vector(data[,RT_col])[rows]) <= as.numeric(as.vector(data[,RT_col])[best_index]) + window/2),
                           which(as.numeric(as.vector(data[,RT_col])[rows]) >= as.numeric(as.vector(data[,RT_col])[best_index]) - window/2))
      
      in_window_candidates<-candidates[in_window]
      
      # Get adducts within that window
      adduct_indices<-rows[in_window]
      adduct_names<-unique(as.vector(data[,adduct_col])[adduct_indices])
      
      # Move negative adducts to the front of the list and sort by intensity
      if(grepl("+", paste(all_adducts,sep = "",collapse = "|")) && ("[M-H]-" %in% adduct_names || "[M-2H]-" %in% adduct_names || "[M+C2H3O2]-" %in% all_adducts || "[M+HCO2]-" %in% all_adducts || "[M+Cl]-" %in% all_adducts)){
        negative_adduct_indices<-c(match("[M-H]-",adduct_names),match("[M-2H]-",adduct_names),match("[M+C2H3O2]-",adduct_names),match("[M+HCO2]-",adduct_names),match("[M+Cl]-",adduct_names))
        negative_adduct_indices<-negative_adduct_indices[!is.na(negative_adduct_indices)]
        negative_adduct_indices<-negative_adduct_indices[order(in_window_candidates[negative_adduct_indices],decreasing = TRUE)]
        adduct_names<-c(adduct_names[negative_adduct_indices],adduct_names[-negative_adduct_indices])
      }
      adducts<-c(adducts,paste(adduct_names,sep = "",collapse = " || "))
      keepers<-c(keepers,best_index)
      
      ### EDITED PJS
      ncol_data <- ncol(data)
      rows <- adduct_indices[which(adduct_indices != best_index)]
      if (length(rows) == 0) {
        next
      }
      data[best_index,  ncol_data - 8] <- paste(data[best_index,  ncol_data - 8], paste(data[rows, ncol_data - 8],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data - 7] <- paste(data[best_index,  ncol_data - 7], paste(data[rows, ncol_data - 7],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data - 4] <- paste(data[best_index,  ncol_data - 4], paste(data[rows, ncol_data - 4],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data - 3] <- paste(data[best_index,  ncol_data - 3], paste(data[rows, ncol_data - 3],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data - 2] <- paste(data[best_index,  ncol_data - 2], paste(data[rows, ncol_data - 2],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data - 1] <- paste(data[best_index,  ncol_data - 1], paste(data[rows, ncol_data - 1],sep = "",collapse = " || "), sep = " || ")
      data[best_index,  ncol_data] <- paste(data[best_index,  ncol_data], paste(data[rows, ncol_data],sep = "",collapse = " || "), sep = " || ")
    }
  }
  
  # Return a frame with no duplicates
  data<-data.frame(data[keepers,],top_ranked_lipid[keepers],adducts)
  data<-data[which(as.vector(data[,ncol(data)])!="[M+Na]+"),]
  colnames(data)[ncol(data)-1] = "Molecular"
  colnames(data)[ncol(data)] = "Adducts Confirmed by MS/MS"
  directory <- "C:/Users/pstel/Documents/PFAS Lab/remove_duplicates_debugging"
  # write.table(data, file.path(directory, "CombinedIDed_Fragments.csv"), sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")  
  return(data)
}