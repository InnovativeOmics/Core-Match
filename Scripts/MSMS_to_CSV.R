

final_MSMS_export <- function(NegPosIDed_dir, FeatureList_in_dir, NegPos) {
  #Mandatory Parameters
  #NEW PAUL, ALWAYS Frag_Annotations.csv, place in library folder (should always be there)
  if (NegPos == "Neg") {
    Frag_lib_dataFrame_file_directory <- file.path(InputLibrary, "Frag_Annotations_Neg.csv")
  }
  if (NegPos == "Pos") {
    Frag_lib_dataFrame_file_directory <- file.path(InputLibrary, "Frag_Annotations_Pos.csv")
  }
  #PAUL: next 3 lines (including commented) MSMSFile, should be fixed / always the same?
  Frag_MZ_col <- 6
  MSMS_scan_RT_col <- 4
  Data_Start_Row_rawMSMS <- 1
  #PAUL: These will be fixed, next 4 lines (ToBeAppended can be changed throughout to Frag_lib)
  ID_Column_Frag_lib <- 2
  MZ_Column_Frag_lib <- 1
  # RT_Column_Frag_lib <- as.numeric(ginput(message="What is the column containing retention times in the feature table \n(the table containing information to append) \nInput should be numeric", title="Retention Time Column",icon="question"))
  Data_Start_Row_Frag_lib <- 1
  #PAUL: user inputs from FluoroMatch Modular / Flow: link up
  ppm_Window <- ppm_Window
  # RT_Window <- as.numeric(ginput(message="What is the retention time window for matching features from the two files? \n(e.g. 0.3 => +/- 0.15 minutes) \nInput should be numeric", title="retention time window",icon="question"))
  #PAUL: This "ID_Method" can just be set to "Fragments" always
  ID_Method <- "Fragments"
  #PAUL: The output should be Neg_rawMSMS.csv in the same directory as NegIDed_FIN.csv, name wont change
  output_file <- paste(ID_Method,"_Appended.csv",sep="")
  
  #Optional inputs specific for aligning features mass using feature numbers for Visualizer platform, turned off for normal use
  
  #NEW PAUL, ALWAYS NL_Annotations.csv, found in library folder (place there)
  if (NegPos == "Neg") {
    NL_dir <- file.path(InputLibrary, "NL_Annotations_Neg.csv")
  }
  if (NegPos == "Pos") {
    NL_dir <- file.path(InputLibrary, "NL_Annotations_Pos.csv")
  }
  #PAUL: Link from user inputs in FluoroMatch Modular / Flow for next3 lines
  MZcol_NegPosIDed <- MZColumn
  FeatureCol_NegPosIDed <- CommentColumn
  Intensity_Threshold <- intensityCutOff
  #PAUL: MSMS columns should always be the same
  FeatureCol_MSMS <- 2
  PrecursorMZCol_MSMS <- 5
  IntensityCol_MSMS <- 7
  #PAUL: Next two lines these files should already be readin and existing in R, or you can re-read in
  NegPosIDed <- read.csv(NegPosIDed_dir, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE, check.names = FALSE)
  NegPosIDed <- as.matrix(NegPosIDed)
  #PAUL: New file to read in
  NL_lib <- read.csv(NL_dir, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE, check.names = FALSE)
  NL_lib <- as.matrix(NL_lib)
  
  
  ############Code##############
  #PAUL: MSMS file should already exist as a data frame no nead to read in
  if (NegPos == "Neg") {
    rawMSMS_df <- Neg_rawMSMS
  } else if (NegPos == "Pos") {
    rawMSMS_df <- Pos_rawMSMS
  }
  #PAUL: New file to read in "Frag_Annotations.csv"
  Frag_lib_df <- read.csv(Frag_lib_dataFrame_file_directory, sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=FALSE, check.names = FALSE)
  rawMSMS_df <- as.matrix(rawMSMS_df)
  
  #Filters for the MSMS Data (only keep fragments below precursor masses, intensity threshold)
  #keep only those below precursor mass + threshold
  Frags_Below_Precursor_Index<-c(1,which(as.numeric(rawMSMS_df[,Frag_MZ_col])<(as.numeric(rawMSMS_df[,PrecursorMZCol_MSMS])+FilterAbovePrecursor))) ## Remove that first 1??
  rawMSMS_df<-rawMSMS_df[Frags_Below_Precursor_Index,]
  #keep only those above intensity threshold
  Intensity_Threshold_Index<-c(1,which(as.numeric(rawMSMS_df[,IntensityCol_MSMS])>Intensity_Threshold))
  rawMSMS_df<-rawMSMS_df[Intensity_Threshold_Index,]
  
  ## Align features to MSMS Table
  
  FeatureList_in <- as.matrix(read.csv(FeatureList_in_dir, sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=TRUE, check.names = FALSE))
  nrow_FeatureList_in <- nrow(FeatureList_in)
  MZ_Append <- as.numeric(FeatureList_in[, 6]) # Column should be constant (future error?)
  RT_Append <- as.numeric(FeatureList_in[, 7]) # Column should be constant
  # IDs_Append <- as.numeric(FeatureList_in[, 12]) # Column should be constant
  rowID_col_new<-which(colnames(FeatureList_in) == "row.ID")
  IDs_Append <- as.numeric(FeatureList_in[, rowID_col_new]) # Column should be constant
  
  Precursor_rawMSMS <- as.numeric(rawMSMS_df[, PrecursorMZCol_MSMS])
  RT_rawMSMS <- as.numeric(rawMSMS_df[, MSMS_scan_RT_col])
  #reverse because the last one is going to be the one that has an actual match and replace the others
  for (j in nrow_FeatureList_in:(RowStartForFeatureTableData-1)) { # Row should be constant (future error?)
    ## Old, more nuanced approach
    
    # message(paste(j, "/", nrow_FeatureList_in, sep = ""))
    # MZ_feature <- as.numeric(as.character(FeatureList_in[j, 1]))
    # MZ_MS2 <- as.numeric(Neg_rawMSMS[, 5])
    # MZConditional <- (MZ_feature - SelectionAccuracy/2) <= MZ_MS2 & MZ_MS2 <= (MZ_feature + SelectionAccuracy/2)
    # MZConditional <- which(MZConditional == TRUE)
    # RT_feature <- as.numeric(as.character(FeatureList_in[j, 2]))
    # RT_MS2 <- as.numeric(Neg_rawMSMS[, 4])
    # RTConditional <- (RT_feature - RT_Window/2) <= RT_MS2 & RT_MS2 <= (RT_feature + RT_Window/2)
    # RTConditional <- which(RTConditional == TRUE)
    # # Add features to Neg_rawMSMS and duplicate rows when already features
    # feat <- as.character(FeatureList_in[j, 3])
    # MZ_RT_int <- intersect(MZConditional, RTConditional)
    # rows_to_dup <- MZ_RT_int[which(Neg_rawMSMS[MZ_RT_int, 2] != "")]
    # if (length(rows_to_dup) > 0) {
    #   Neg_rawMSMS <- rbind(Neg_rawMSMS, Neg_rawMSMS[rows_to_dup,])
    #   nrow_Neg_rawMSMS <- nrow(Neg_rawMSMS)
    #   Neg_rawMSMS[(nrow_Neg_rawMSMS - length(rows_to_dup) + 1):nrow_Neg_rawMSMS, 2] <- feat
    # }
    # rows_to_update <- MZ_RT_int[which(Neg_rawMSMS[MZ_RT_int, 2] == "")]
    # if (length(rows_to_update) > 0) {
    #   Neg_rawMSMS[rows_to_update, 2] <- feat
    # }
    
    ## Simplified approach
    
    MZ_Append_current <- MZ_Append[j]
    RT_Append_current <- RT_Append[j]
    IDs_Append_current <- IDs_Append[j]
    if (alignMS2_ppm==TRUE) {
      #calculate conditionals of whether or not a feature will have an associated MS/MS spectrum
      ppm_error <- ((Precursor_rawMSMS - MZ_Append_current)*10^6)/Precursor_rawMSMS
      MZ_Conditional <- abs(ppm_error)<(ppm_Window/2)
    } else {
      #calculate conditionals of whether or not a feature will have an associated MS/MS spectrum
      Da_error <- Precursor_rawMSMS - MZ_Append_current
      MZ_Conditional <- abs(Da_error)<(SelectionAccuracy/4)
    }
    MZ_Index <- which(MZ_Conditional==TRUE)
    RT_Conditional <- ((RT_Append_current - RT_Window/2) < RT_rawMSMS) & (RT_rawMSMS < (RT_Append_current + RT_Window/2))
    RT_Index <- which(RT_Conditional==TRUE)
    MZ_RT_int <- intersect(MZ_Index, RT_Index)
    rawMSMS_df[MZ_RT_int,FeatureCol_MSMS] <- IDs_Append_current
  }
  
  #Filters remove all MSMS scans without features
  #which rows have features
  Features_MSMSfile<-as.numeric(rawMSMS_df[,FeatureCol_MSMS])
  AllFeature_Index<-c(which(!is.na(Features_MSMSfile)))
  rawMSMS_df<-rawMSMS_df[AllFeature_Index,]
  
  ##consolidate the Fragment Screening list first before searching, otherwise values will just be overwritten, still values will be over written if within ppm but not exactly the same
  Frag_lib_df<-aggregate(Frag_lib_df, list(Frag_lib_df[,MZ_Column_Frag_lib]), FUN=function(x) paste(x,collapse=";"))
  Frag_lib_df <- as.matrix(Frag_lib_df)
  #check if works for other formats COULD ERROR
  Frag_lib_df<-Frag_lib_df[,-2]
  
  ##(Same as above for NL) consolidate the Fragment Screening list first before searching, otherwise values will just be overwritten, still values will be over written if within ppm but not exactly the same
  NL_lib<-aggregate(NL_lib, list(NL_lib[,1]), FUN=function(x) paste(x,collapse=";"))
  NL_lib <- as.matrix(NL_lib)
  #check if works for other formats COULD ERROR
  NL_lib<-NL_lib[,-2]
  
  #Add 6 columns for data to be append on to: Feature m/z, Name, and ppm error
  nrowrawMSMS <- nrow(rawMSMS_df)
  ID_and_ppm_cols<-matrix("",nrowrawMSMS,6)
  rawMSMS_df <- cbind(rawMSMS_df,ID_and_ppm_cols)
  colnames(rawMSMS_df)[(ncol(rawMSMS_df)-5):ncol(rawMSMS_df)] <- c("Feature m/z","NL mz",ID_Method,"ppm_Error","NL","ppm_Error_NL")
  
  MSMS_PrecursorMZ_Col<-ncol(rawMSMS_df)-5
  NL_Col<-ncol(rawMSMS_df)-4
  Name_col<-ncol(rawMSMS_df)-3
  ppm_col<-ncol(rawMSMS_df)-2
  NL_Name_col<-ncol(rawMSMS_df)-1
  NL_ppm_col<-ncol(rawMSMS_df)
  
  nrowFrag_lib <- nrow(Frag_lib_df)
  ncolrawMSMS <- ncol(rawMSMS_df)
  
  ##Appending Precursor Masses to MSMS Table for NL Searching (inputs for ease below)
  # MZcol_NegPosIDed
  # FeatureCol_NegPosIDed
  # FeatureCol_MSMS
  # Matching numbers in vectors is a lot faster then indexing a matrix and matching characters!
  Features_MSMSfile<-as.numeric(rawMSMS_df[,FeatureCol_MSMS])
  Features_NegPosIDed<-as.numeric(NegPosIDed[, FeatureCol_NegPosIDed])
  for(i in 2:nrow(NegPosIDed)){
    Feature_Conditional <- Features_NegPosIDed[i]==Features_MSMSfile
    Feature_Index <- which(Feature_Conditional==TRUE)
    rawMSMS_df[Feature_Index,MSMS_PrecursorMZ_Col] <- NegPosIDed[i, MZcol_NegPosIDed]
  }
  #Calculate Neutral loss
  #which rows have features
  AllFeature_Index<-which(!is.na(Features_MSMSfile))
  #Precursor masses and index
  PrecursorMZs<-as.numeric(rawMSMS_df[AllFeature_Index,MSMS_PrecursorMZ_Col])
  FragmentMZs<-as.numeric(rawMSMS_df[AllFeature_Index,Frag_MZ_col])
  rawMSMS_df[AllFeature_Index,NL_Col] <- PrecursorMZs-FragmentMZs
  
  ## Identify neutral losses (NLs) (could make this a function since used twice, also will be used in Pos and Neg... and for multiple files... but then the whole things needs to be a function)
  MZ_NL_lib <- as.numeric(NL_lib[, 1])
  MZ_NL_MSMS <- as.numeric(rawMSMS_df[, NL_Col])
  
  for(o in 1:nrow(NL_lib)){
    MZ_Append_current <- MZ_NL_lib[o]
    IDs_Append <- NL_lib[o, 2]
    ppm_error <- ((MZ_NL_MSMS - MZ_Append_current)*10^6)/MZ_NL_MSMS
    MZ_Conditional <- abs(ppm_error)<(ppm_Window/2)
    MZ_Index <- which(MZ_Conditional==TRUE)
    rawMSMS_df[MZ_Index,NL_Name_col] <- IDs_Append
    rawMSMS_df[MZ_Index,NL_ppm_col] <- ppm_error[MZ_Index]
    rawMSMS_df[MZ_Index,Name_col] <- IDs_Append
    rawMSMS_df[MZ_Index,ppm_col] <- ppm_error[MZ_Index]
  }
  
  MZ_Append <- as.numeric(Frag_lib_df[, MZ_Column_Frag_lib])
  # Frag_lib_df[, RT_Column_Frag_lib] <- as.numeric(as.character(Frag_lib_df[, RT_Column_Frag_lib]))
  MZ_rawMSMS <- as.numeric(rawMSMS_df[, Frag_MZ_col])
  # rawMSMS_df[, MSMS_scan_RT_col] <- as.numeric(as.character(rawMSMS_df[, MSMS_scan_RT_col]))
  
  ## Identify Fragments
  for(o in Data_Start_Row_Frag_lib:nrowFrag_lib){
    MZ_Append_current <- MZ_Append[o]
    # RT_rawMSMS <- rawMSMS_df[o, MSMS_scan_RT_col]
    IDs_Append <- Frag_lib_df[o, ID_Column_Frag_lib]
    # RT_Frag_lib <- Frag_lib_df[tba, RT_Column_Frag_lib]
    ppm_error <- ((MZ_rawMSMS - MZ_Append_current)*10^6)/MZ_rawMSMS
    MZ_Conditional <- abs(ppm_error)<(ppm_Window/2)
    MZ_Index <- which(MZ_Conditional==TRUE)
    # RT_Conditional <- ((RT_Frag_lib - RT_Window/2) < RT_rawMSMS) && (RT_rawMSMS < (RT_Frag_lib + RT_Window/2))
    rawMSMS_df[MZ_Index,Name_col] <- IDs_Append
    rawMSMS_df[MZ_Index,ppm_col] <- ppm_error[MZ_Index]
  }
  
  # Remove obsolete columns
  rawMSMS_df <- cbind(rawMSMS_df[, 1:7], rawMSMS_df[, 10:15])
  
  if (NegPos == "Neg") {
    write.csv(rawMSMS_df, file.path(InputDirectory, "Output/Neg_rawMSMS.csv"), row.names = FALSE, na = "")
  } else if (NegPos == "Pos") {
    write.csv(rawMSMS_df, file.path(InputDirectory, "Output/Pos_rawMSMS.csv"), row.names = FALSE, na = "")
  }
}