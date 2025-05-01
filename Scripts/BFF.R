#InputDirectory<-"C:/Users/Jeremy Koelmel/Desktop/Desktop/Innovative_Omics/CLIENTS/2024_USGS/JBB_QCs/Annotated/"; CommentColumn<-1; RowStartForFeatureTableData<-2; NegPos = "Neg"; OutDir="/Output/NegIDed_FIN.csv"; OutDirOnlyFrags="/Output/Neg_OnlyIDs_PredFrags.csv"; InputDir_Append="Output/ddMS/NegByClass/Additional_Files"; ID_name="PredictedFrag_IDs"; FragName="Frags"; nFrag="Num_Frags"; fileNames="Files"; ImportTable="/Output/NegIDed_Fragments.csv"
#testing parameters
# InputDirectory<-"D:/FluoroMatch_Data/2024_06_26_Testing_OUTPUTS/2024_12_31_FLUOROMATCH55/McGill_Rapid_MZ3_BFF/Annotated/"
# runPosddMS<-FALSE
# runNegddMS<-TRUE
# Polarity<-"neg"

#####################################Blank Feature Filtering (BFF), currently only for MZMine 3 and Flow#########################################
read_extra_parameters <- function(InputDirectory) {
  
  ExtraParDir<-paste0(dirname(InputDirectory),"/Run_Details/ExtraParameters.txt")
  
  # Read the text file
  text <- readLines(ExtraParDir)
  
  # Extract MZMine version
  mzmine_version_line <- grep("MZmine Version:", text, value = TRUE)
  
  if (length(mzmine_version_line) > 0) {
    mzmine_version <- as.integer(gsub("MZmine Version: ", "", mzmine_version_line))
    print(paste("MZMine Version =", mzmine_version))
  } else {
    print("MZMine Version not found.")
    mzmine_version <- NA
  }
  
  # Extract blank filtering parameters
  blank_filtering_line <- grep("Blank Filtering Values:", text, value = TRUE)
  
  if (length(blank_filtering_line) > 0) {
    blank_filtering_params <- strsplit(gsub("Blank Filtering Values: ", "", blank_filtering_line), ", ")[[1]]
    a_BFF <- as.integer(strsplit(blank_filtering_params[1], " = ")[[1]][2])
    b_BFF <- as.integer(strsplit(blank_filtering_params[2], " = ")[[1]][2])
    c_BFF <- as.integer(strsplit(blank_filtering_params[3], " = ")[[1]][2])
    print("parsing blank filtering parameters")
    print(paste("a =", a_BFF))
    print(paste("b =", b_BFF))
    print(paste("c =", c_BFF))
  } else {
    print("Blank Filtering Values not found.")
    a_BFF <- NA
    b_BFF <- NA
    c_BFF <- NA
  }
  # Return extracted parameters
  return(list(
    mzmine_version = mzmine_version,
    a_BFF = a_BFF,
    b_BFF = b_BFF,
    c_BFF = c_BFF
  ))
}

BFF_Filter <- function(InputDirectory, Polarity, a_BFF, b_BFF, c_BFF) {
  
  # List files ending with "pos.csv" (case-insensitive)
  FileDir <- dir(path = InputDirectory, pattern = paste0(Polarity,"\\.csv$"), ignore.case = TRUE, full.names = TRUE)[1]
  
  # Read the CSV file
  FeatureTable<-read.csv(FileDir)
  
  # Find columns with "blank" and desired keywords (case-insensitive)
  blank_cols <- grep("blank", names(FeatureTable), ignore.case = TRUE)
  sample_cols <- grep(paste0("(",Polarity,"|mzxml|mzml|peak\\.area|peak\\.height)"), names(FeatureTable), ignore.case = TRUE)
  
  # Replace non-numeric values with 0 in specified columns
  FeatureTable[, sample_cols] <- lapply(FeatureTable[, sample_cols], function(x) {
    as.numeric(as.character(x))
  })
  
  # Replace NA (resulting from conversion) with 0
  FeatureTable[, sample_cols][is.na(FeatureTable[, sample_cols])] <- 0
  
  blank_cols <- intersect(blank_cols, sample_cols)
  
  BFF <- TRUE
  if(length(blank_cols)<3) {
    print(paste0("only ",length(blank_cols)," blanks found, skipping BFF, a minimum of 3 blanks needed"))
    BFF<-FALSE
  }
  
  if(BFF==TRUE){
    sample_cols <- setdiff(sample_cols, blank_cols)
    
    # Calculate percentiles for each row in other_cols
    if (length(sample_cols) == 1) {
      percentiles <- as.numeric(FeatureTable[, sample_cols])
    } else {
      percentiles <- apply(FeatureTable[, sample_cols], 1, function(row) {
        quantile(row, probs = a_BFF, na.rm = TRUE)
      })
    }
    
    
    # Calculate average + 3 * standard deviation for blank_cols
    BFFs <- apply(FeatureTable[, blank_cols], 1, function(row) {
      b_BFF * (mean(row, na.rm = TRUE) + c_BFF * sd(row, na.rm = TRUE))
    })
    
    Treshold<-percentiles>BFFs
    
    print(paste("Of ",length(Treshold), "Features, The following number of features were retained after BFF: ",length(Treshold[Treshold])," (",round((length(Treshold[Treshold])/length(Treshold))*100)," percent )"))
    
    # Filter rows where percentile is greater than threshold
    filtered_df <- FeatureTable[Treshold,]
    
    write.csv(write.csv(filtered_df,FileDir,row.names = FALSE))
  }
}

