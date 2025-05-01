##Jeremy Koelmel  jeremykoelmel@gmail.com
##Michael Kummer  mk.code.na@gmail.com
##Paul Stelben    paul.stelben@yale.edu
##Nick Kroeger    nkroeger.cs@gmail.com
##David Schiessel david@schiessel.us 
##Jonathan Sparks jonathansparks0758@gmail.com 

rm(list = ls()) #remove all R objects from memory, this program can be memory intensive as it is dealing with huge datasets
closeAllConnections()

#### Mandatory Parameter to Change ####
#If you want to manually input your variables...
# 1. Set ManuallyInputVariables <- TRUE (all caps)
# 2. Assign variables under the next if statement, "if (ManuallyInputVariables==TRUE)"
FLOW <- FALSE
csvInput <- TRUE
ManuallyInputVariables <- FALSE
RT_flagging <- TRUE #JPK: for PFAS analysis
ParallelComputing <- TRUE
Lipid <- FALSE
Tween_pos <- FALSE #PJS: for PolyMatch
FilterAbovePrecursor <- 1 #how far from the precursor should fragment masses be kept (e.g. if precursor is 700, should 702 be considered?)
TargetEIC_Only <- TRUE
alignMS2_ppm <- FALSE #selection accuracy will be divided by 4, eg. 1 Da Window = +/- 0.25 difference, if FALSE, Otherwise the mass accuracy in ppm will be used to align MSMS scans to the feature table
EICcpp <- FALSE #Whether to use the C++ version of the EIC / MS1 generation code, seems to give the same results, faster but less well tested
IMfirst<-FALSE #For running Modular after Running FluoroMatch IM using PNNL processed DIA to DDA workflow or similar
#### END Mandatory Parameter to Change ####
#### END "Read Me" Section ####

if (ManuallyInputVariables){
  #######################MANUAL INPUTS#######################
  # Input Directory for Libraries ...must have forward slashes but nothing at the end of the directory
  InputLibrary<-"C:/Users/Jeremy Koelmel/Downloads/LipidMatch_Flow_3.5/LipidMatch_Flow_3.5/LipidMatch_Modular/LipidMatch_Libraries_Acetate"
  source(paste(InputLibrary,"/Scripts/Input_Manual.R",sep=""))
} else if (csvInput){
  #######################CSV INPUTS (Recommended)#######################
  #parametersDirectory
parametersDir <- "C:/SOFTWARE/FluoroMatch_5.7/Flow/Distribution/"
parametersFile <- paste(parametersDir, "PARAMETERS.csv", sep="")
  parametersInput_csv <- read.csv(parametersFile, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  parametersInput_csv <- as.matrix(parametersInput_csv)
  # Input Directory for Libraries ...must have \\ (double backslash) or / at the end of the directory
  InputLibrary<-as.character(parametersInput_csv[15,2])
  InputLibrary <- substring(InputLibrary,1, nchar(InputLibrary)-1)
  source(paste(InputLibrary,"/Scripts/Input_CSV.R",sep=""))
}else{
  #######################POPUPs for Modular#######################
  if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets", repos = "http://cran.us.r-project.org")}
  if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk", repos = "http://cran.us.r-project.org")}
  # if (FLOW == FALSE && csvInput == FALSE && ManuallyInputVariables == FALSE) {
  require(gWidgets)
  require(gWidgetstcltk)
  options(guiToolkit="tcltk")
  
  ## Input Directory for Libraries
  InputLibrary<-tk_choose.dir(caption="Input Directory of libraries (came with LipidMatch)")
  if(is.na(InputLibrary)){
    stop()
  }
  source(paste(InputLibrary,"/Scripts/Input_Popups.R",sep=""))
} 
####################### END VARIABLE DECLARATIONS SECTION ####################

#############################################Libraries and Dependencies########################################################################
source(paste(InputLibrary,"/Scripts/Dependencies_Libraries_And_Parallelization.R",sep=""))

#Error handling for manual input section... specifically the Input Directory. Checks for Output folder, and stops.
source(paste(InputLibrary,"/Scripts/CheckOutputFolder.R",sep=""))

#Library Information
# NAME AND DIRECTORY for exact mass library

ImportLibPOS<-file.path(InputLibrary, "Precursor_Library_POS.csv")
ImportLibNEG<-file.path(InputLibrary, "Precursor_Library_NEG.csv")

LibCriteria<- file.path(InputLibrary, "ID_CRITERIA.csv")

LibColInfo<-1  #Look at your Library. What column are your IDs? (Columns to retrieve ID information and align with data)
ParentMZcol_in<-2 #Look at your Library. What column are your precursor masses in? (Columns to retrieve ID information and align with data)
ddMS2_Code<-"1"
AIF_Code<-"2"
Class_Code<-"3"
ExactMass_Code<-"4"
NoID_Code<-"5"

PPM_CONST <- (10^6 + ppm_Window/2) / 10^6

RepeatingUnits_dir <- file.path(InputLibrary, "REPEATING_UNITS_INPUT.csv")

################### FUNCTIONS ###################

source(paste(InputLibrary,"/Scripts/BFF.R",sep=""))
source(paste(InputLibrary,"/Scripts/AIF_OG_Nick.R",sep=""))
source(paste(InputLibrary,"/Scripts/MS2_RuleBased_Match.R",sep=""))
source(paste(InputLibrary,"/Scripts/Gen_MZMine_Target_List_From_Feature_Table.R",sep=""))
source(paste(InputLibrary,"/Scripts/Load_MS2_as_RObject.R",sep=""))
source(paste(InputLibrary,"/Scripts/Link_IDs_to_FeatureTable.R",sep=""))
source(paste(InputLibrary,"/Scripts/AppendFrags_From_LibHits.R",sep=""))
source(paste(InputLibrary,"/Scripts/Remove_Duplicates.R",sep=""))
source(paste(InputLibrary,"/Scripts/HomologousSeries_Detection_And_FeatureScoring.R",sep=""))
source(paste(InputLibrary,"/Scripts/MSMS_to_CSV.R",sep=""))
source(paste(InputLibrary,"/Scripts/IM_And_MSMS_Combine.R",sep=""))
source(paste(InputLibrary,"/Scripts/Column_Name_Manipulation.R",sep=""))
source(paste(InputLibrary,"/Scripts/IsotopePercentages.R",sep=""))
source(paste(InputLibrary,"/Scripts/Kaufmann.R",sep=""))
source(paste(InputLibrary,"/Scripts/kaufmann_kde_percentiles.R",sep=""))
source(paste(InputLibrary,"/Scripts/Manual_Review.R",sep=""))
source(paste(InputLibrary,"/Scripts/Stats.R",sep=""))
source(paste(InputLibrary,"/Scripts/Formula_Prediction.R",sep=""))
source(paste(InputLibrary,"/Scripts/homologousVoting.R",sep=""))
if (EICcpp==TRUE) {
  source(paste(InputLibrary,"/Scripts/EIC_MS1_fns_Cpp.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genEIC_Cpp.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genIsoTable_Cpp.R",sep=""))
  source(paste(InputLibrary,"/Scripts/MS1Spectragen_Cpp.R",sep=""))
} else {
  source(paste(InputLibrary,"/Scripts/EIC_MS1_fns.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genEIC.R",sep=""))
  source(paste(InputLibrary,"/Scripts/genIsoTable.R",sep=""))
  source(paste(InputLibrary,"/Scripts/MS1Spectragen.R",sep=""))
}

####################################End functions##########################################

########################### Run BFF ##############################

FeatureTable_NEG <- list.files(path = InputDirectory, pattern = "Neg\\.csv$", ignore.case = TRUE)
if(length(FeatureTable_NEG) > 1){
  stop(paste("ERROR: You should only have 1 Negative mode Feature Table... we detected", length(FeatureTable_NEG)," Feature Tables in the folder:", InputDirectory))
}else if(length(FeatureTable_NEG) == 0){
  print(paste("CAUTION: Could not find any negative mode Feature Tables... we detected", length(FeatureTable_NEG)," Feature Tables in the folder: ", InputDirectory," ...If this incorrect, check that you have 'Neg.csv' at the end of the file name."))
}

FeatureTable_POS <- list.files(path = InputDirectory, pattern = "Pos\\.csv$", ignore.case = TRUE)
if(length(FeatureTable_POS) > 1){
  stop(paste("ERROR: You should only have 1 Positive mode Feature Table... we detected", length(FeatureTable_POS)," Feature Tables in the folder:", InputDirectory))
}else if(length(FeatureTable_POS) == 0){
  print(paste("CAUTION: Could not find any Positive mode Feature Tables... we detected", length(FeatureTable_POS)," Feature Tables in the folder: ", InputDirectory," ...If this incorrect, check that you have 'Neg.csv' at the end of the file name."))
}

runPosddMS <- FALSE
runNegddMS <- FALSE
# runPosAIF <- FALSE
# runNegAIF <- FALSE
#Run Negative mode analysis if there are negative .csv files
if(length(FeatureTable_NEG) == 1){  
  runNegddMS <- TRUE  
}
#Run Positive mode analysis if there are positive .csv files
if(length(FeatureTable_POS) == 1) {
  runPosddMS <- TRUE
}

if (FLOW == TRUE) {
  BFF_params <- read_extra_parameters(InputDirectory)
  if (BFF_params$mzmine_version == 3){
    if(runPosddMS) { 
      BFF_Filter(InputDirectory, Polarity = "pos", a_BFF = BFF_params$a_BFF, b_BFF = BFF_params$b_BFF, c_BFF = BFF_params$c_BFF)
    }
    if(runNegddMS) { 
      BFF_Filter(InputDirectory, Polarity = "neg", a_BFF = BFF_params$a_BFF, b_BFF = BFF_params$b_BFF, c_BFF = BFF_params$c_BFF)
    }
  }
}

####################### IonDecon Integration ####################
source(paste(InputLibrary,"/Scripts/AIF.R",sep=""))
source(paste(InputLibrary,"/Scripts/EIC_MS1_AIF_QRAI_fns.R",sep=""))
source(paste(InputLibrary,"/Scripts/run_all_ions.R",sep=""))
if (runNegddMS==TRUE) {
  MZMLs <- list.files(path=InputDirectory, pattern="All_Ions.*[nNgG]\\.mzX*ML", ignore.case=TRUE) #All_Ions _Neg _Pos.mzML or _P _N case insensitive
  IonDecon_AIF(
    InputDirectory, FeatureTable_NEG, MZMLs,
    MZColumn, RTColumn, CommentColumn,
    RT_Window/2, PrecursorMassAccuracy,
    corrMin, minNumberOfAIFScans, intensityCutOff)
}
if (runPosddMS==TRUE) {
  MZMLs <- list.files(path=InputDirectory, pattern="All_Ions.*[pPsS]\\.mzX*ML", ignore.case=TRUE) #All_Ions _Neg _Pos.mzML or _P _N case insensitive
  IonDecon_AIF(
    InputDirectory, FeatureTable_POS, MZMLs,
    MZColumn, RTColumn, CommentColumn,
    RT_Window/2, PrecursorMassAccuracy,
    corrMin, minNumberOfAIFScans, intensityCutOff)
}

####Read in files, create folder structure, error handle if outputs exist and run functions####

##sets up which files (negative or positive) to run from which directories, etc.
fileName <- basename(InputDirectory)

##For testing Parsing
# ddMS2NEG_in <- list.files("C:/Users/User/Downloads/TEST/", pattern="[nNgG]\\.(mzML|ms2)", ignore.case=FALSE)
# ddMS2NEG_in[grep("[dD][dD][mM][sS]", ddMS2NEG_in)]
ddMS2NEG_in <- list.files(path=InputDirectory, pattern="[nNgG]\\.(mzML|ms2)", ignore.case=FALSE) #_Neg _Pos.mzML or _P _N case insensitive
ddMS2NEG_in <- ddMS2NEG_in[grep("[dD][dD][mM][sS]", ddMS2NEG_in)]
ddMS2NEG_in_ms2 <- ddMS2NEG_in[grep("\\.ms2", ddMS2NEG_in)]
ddMS2NEG_in_mzML <- ddMS2NEG_in[grep("\\.mzML", ddMS2NEG_in)]
ddMS2POS_in <- list.files(path=InputDirectory, pattern="[pPsS]\\.(mzML|ms2)", ignore.case=FALSE)
ddMS2POS_in <- ddMS2POS_in[grep("[dD][dD][mM][sS]", ddMS2POS_in)]
ddMS2POS_in_ms2 <- ddMS2POS_in[grep("\\.ms2", ddMS2POS_in)]
ddMS2POS_in_mzML <- ddMS2POS_in[grep("\\.mzML", ddMS2POS_in)]

#user info outputted for error handling
if(length(ddMS2POS_in) == 0){
  print(paste("CAUTION: We detected", length(ddMS2POS_in),"positive ddMS .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have 'p', 'P', 'pos', or 'POS' at the end of the file name and you must have a 'dd' within the name. OR Remove the folder: ", fileName))
}
if(length(ddMS2NEG_in) == 0){
  print(paste("CAUTION: We detected", length(ddMS2NEG_in),"negative ddMS .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have 'n', 'N', 'neg', or 'NEG' at the end of the file name and you must have a 'dd' within the name. OR Remove the folder: ", fileName))
}

#Negative/Positive mode sample names (took .ms2 files and dropped the ".ms2")
ExtraSampleNameddMSNEG_in <- vector()
ExtraSampleNameddMSPOS_in <- vector()
# ExtraSampleNameAIFNEG_in <- vector()
# ExtraSampleNameAIFPOS_in <- vector()
for(j in seq_len(length(ddMS2NEG_in))){   ExtraSampleNameddMSNEG_in[j] <- sub("\\.\\w+", "", ddMS2NEG_in[j])    }
for(j in seq_len(length(ddMS2POS_in))){   ExtraSampleNameddMSPOS_in[j] <- sub("\\.\\w+", "", ddMS2POS_in[j])    }
# for(j in seq_len(length(AIFMS2NEG_in))){   ExtraSampleNameAIFNEG_in[j] <- sub("\\.\\w+", "", AIFMS2NEG_in[j])    }
# for(j in seq_len(length(AIFMS2POS_in))){   ExtraSampleNameAIFPOS_in[j] <- sub("\\.\\w+", "", AIFMS2POS_in[j])    }

runPosddMS <- FALSE
runNegddMS <- FALSE
# runPosAIF <- FALSE
# runNegAIF <- FALSE
#Run Negative mode analysis if there are negative .ms2 and .csv files
if(length(ddMS2NEG_in) != 0 && length(FeatureTable_NEG) == 1){  
  runNegddMS <- TRUE  
}
#Run Positive mode analysis if there are positive .ms2 and .csv files
if(length(ddMS2POS_in) != 0 && length(FeatureTable_POS) == 1) {
  runPosddMS <- TRUE
}
# if(length(AIFMS2NEG_in) != 0 && length(AIFMS1NEG_in) != 0 && length(FeatureTable_NEG) == 1){  runNegAIF <- TRUE   }
# if(length(AIFMS2POS_in) != 0 && length(AIFMS1POS_in) != 0 && length(FeatureTable_POS) == 1){  runPosAIF <- TRUE   }

OutputDirectory<-file.path(InputDirectory, "Output",sep="")
if(!dir.exists(OutputDirectory)){ dir.create(OutputDirectory) }

##Generate output directories and whether to run Positive or Negative MS/MS matching or Both
if(runPosddMS || runNegddMS){#ddMS
  # OutputDirectoryddMS <- file.path(InputDirectory, "Output", "ddMS") # PJS 11/13/2022
  OutputDirectoryddMS <- file.path(InputDirectory, "Output", "ddMS") # PJS 11/13/2022
  if(!dir.exists(OutputDirectoryddMS)){ dir.create(OutputDirectoryddMS) }
  if(runPosddMS){#pos ddMS
    OutputDirectoryddMSPos_in <- file.path(OutputDirectoryddMS,"Pos")
    if(!dir.exists(OutputDirectoryddMSPos_in)){ dir.create(OutputDirectoryddMSPos_in) }
    OutputDirectoryddMSPosByClass_in <- file.path(OutputDirectoryddMS,"PosByClass")
    if(!dir.exists(OutputDirectoryddMSPosByClass_in)){ dir.create(OutputDirectoryddMSPosByClass_in) }
  }
  if(runNegddMS){#negByClass ddMS
    OutputDirectoryddMSNegByClass_in <- file.path(OutputDirectoryddMS,"NegByClass")
    if(!dir.exists(OutputDirectoryddMSNegByClass_in)){ dir.create(OutputDirectoryddMSNegByClass_in) }
    OutputDirectoryddMSNeg_in <- file.path(OutputDirectoryddMS,"Neg")
    if(!dir.exists(OutputDirectoryddMSNeg_in)){ dir.create(OutputDirectoryddMSNeg_in) }
  }
}

####################### Run the libraries and input data ######
NegClassDDLib <- FALSE
PosDDLib <- FALSE
NegDDLib <- FALSE
PosClassDDLib <- FALSE
NegAIFLib <- FALSE
PosAIFLib <- FALSE

#NEG
LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "NEG" & toupper(LibraryCriteria[,6]) == "FALSE",] #subset LibraryCriteria to find negative libraries
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,4]) == "TRUE",] #subset LibraryCriteria to find ddMS libraries to run

if(runNegddMS && nrow(LibraryCriteria)>0){
  NegDDLib <- TRUE
  FeatureTable_dir_in<-file.path(InputDirectory, FeatureTable_NEG)
  cat(paste0("Reading in file:\t", FeatureTable_NEG,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
  FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
  # Create list of dataframes with all ms2 data
  
  length_ddMS2NEG_in <- length(ddMS2NEG_in)
  MS2_df_list <- vector("list", length_ddMS2NEG_in)
  nrow_all_scans <- 0
  #decipher .MS2 files OR (else) mzML files
  
  if (length(ddMS2NEG_in_ms2)>0){
    for (c in 1:length(ddMS2NEG_in_ms2)){
      MS2_dir_in <- file.path(InputDirectory, ddMS2NEG_in_ms2[c])
      MS2_df_in <- createDataFrame(MS2_dir_in)
      MS2_df_list[[c]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
      }
    }
  }
  
  if (length(ddMS2NEG_in_mzML)>0){
    for (c in 1:length(ddMS2NEG_in_mzML)){
      MS2_dir_in <- file.path(InputDirectory, ddMS2NEG_in_mzML[c])
      MS2_df_in_peaks <- openMSfile(MS2_dir_in)
      MS2_df_in_h <- header(MS2_df_in_peaks)
      MS2_df_in_h <- MS2_df_in_h[MS2_df_in_h$msLevel==2,c("precursorMZ","retentionTime")]
      MS2_df_in_h[,2] <- MS2_df_in_h[,2]/60
      PEAKS = lapply(as.numeric(row.names(MS2_df_in_h)), function(i) peaks(MS2_df_in_peaks, i) )
      MS2_df_in<-cbind(MS2_df_in_h,as.matrix(PEAKS))
      colnames(MS2_df_in) <- c("precursor", "rt","mz_intensity")
      NonEmptyMS2spectra<-as.logical(lapply(1:nrow(MS2_df_in), function(i) length(unlist(MS2_df_in[i,3]))>2)) #throw away things with only 1 row might not want this (:
      MS2_df_in<-MS2_df_in[NonEmptyMS2spectra,]
      MS2_df_list[[c + length(ddMS2NEG_in_ms2)]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
      }
    }
  }
  
  
  # Add scans to master MSMS export matrix
  Neg_rawMSMS <- matrix("", nrow_all_scans, 9)
  colnames(Neg_rawMSMS)[1:9] <- c("Arbitrary_Identifier", "Feature", "File", "Selected_RT", "Selected_mz", "mz", "Intensity", "Fragments", "LibraryFile")
  start <- 1
  for (c in 1:length_ddMS2NEG_in){
    # message(paste(c, "/", length_ddMS2NEG_in), sep = "")
    nrow_MS2_df_list_c <- nrow(MS2_df_list[[c]])
    MS2_RT <- as.numeric(MS2_df_list[[c]][,2])
    MS2_MZ <- as.numeric(MS2_df_list[[c]][,1])
    for(j in 1:nrow_MS2_df_list_c){
      scan_rows <- nrow(MS2_df_list[[c]][j,][[3]][[1]])
      # start <- which(Neg_rawMSMS[, 1] == "")[1]
      stop <- start + scan_rows - 1
      Neg_rawMSMS[start:stop, 1] <- j
      Neg_rawMSMS[start:stop, 2] <- ""
      Neg_rawMSMS[start:stop, 3] <- ddMS2NEG_in[c]
      Neg_rawMSMS[start:stop, 4] <- MS2_RT[j]
      Neg_rawMSMS[start:stop, 5] <- MS2_MZ[j]
      Neg_rawMSMS[start:stop, 6:7] <- MS2_df_list[[c]][j,][[3]][[1]]
      Neg_rawMSMS[start:stop, 8] <- ""
      Neg_rawMSMS[start:stop, 9] <- ""
      start <- stop + 1
    }
  }
  write.csv(Neg_rawMSMS, file.path(InputDirectory, "Output/Neg_rawMSMS_NoAnnotations.csv"), row.names = FALSE, col.names = TRUE, na = "")
  # Run RunTargeted
  for (c in 1:length_ddMS2NEG_in){
    # message(paste("c = ", c, sep = ""))
    cat(paste0("Reading in file:\t", ddMS2NEG_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
    ExtraSample<-ExtraSampleNameddMSNEG_in[c]
    OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in)
    Neg_rawMSMS_List <- c()
    nrow_LibraryCriteria <- nrow(LibraryCriteria)
    
    for(i in seq_len(nrow_LibraryCriteria)){
      # message(paste(i, "/", nrow_LibraryCriteria, sep = ""))
      LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
      OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
      ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
      ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
      if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
      if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
      RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSNeg_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2NEG_in[c], "Neg")
    }
  }
  print(paste("Finished Negative ddMS analysis", timestamp(), sep = " --> "))
  # Export MSMS combined file
  # write.csv(Neg_rawMSMS, file.path(InputDirectory, "Output/Neg_rawMSMS_Example2.csv"), row.names = FALSE, col.names = TRUE, na = "")
}

#NEG BY CLASS
LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "NEG" & toupper(LibraryCriteria[,6]) == "TRUE",] #subset LibraryCriteria to find negative libraries
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,4]) == "TRUE",] #subset LibraryCriteria to find ddMS libraries to run
if(runNegddMS && nrow(LibraryCriteria)>0){
  NegClassDDLib <- TRUE
  FeatureTable_dir_in<-file.path(InputDirectory, FeatureTable_NEG)
  cat(paste0("Reading in file:\t", FeatureTable_NEG,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
  FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
  length_ddMS2NEG_in <- length(ddMS2NEG_in)
  # Run RunTargeted
  for (c in 1:length_ddMS2NEG_in){
    # message(paste("c = ", c, sep = ""))
    cat(paste0("Reading in file:\t", ddMS2NEG_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
    ExtraSample<-ExtraSampleNameddMSNEG_in[c]
    OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in)
    nrow_LibraryCriteria <- nrow(LibraryCriteria)
    if (ParallelComputing == TRUE) {
      foreach (i = seq_len(nrow(LibraryCriteria)), .packages = c("sqldf")) %dopar% {
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSNegByClass_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2NEG_in[c], "NegByClass")
      }
    } else {
      for(i in seq_len(nrow_LibraryCriteria)){
        # message(paste(i, "/", nrow_LibraryCriteria, sep = ""))
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSNegByClass_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2NEG_in[c], "NegByClass")
      }
    }
  }
  print(paste("Finished Negative by class ddMS analysis", timestamp(), sep = " --> "))
  # Export MSMS combined file
  # write.csv(Neg_rawMSMS, file.path(InputDirectory, "Output/Neg_rawMSMS_Example3.csv"), row.names = FALSE, col.names = TRUE, na = "")
}

#POS ddMS
LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "POS" & toupper(LibraryCriteria[,6]) == "FALSE",] #subset LibraryCriteria to find positive libraries(not pos-by-class)
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,4]) == "TRUE",] #subset LibraryCriteria to find ddMS libraries to run
if(runPosddMS && nrow(LibraryCriteria)>0){
  PosDDLib <- TRUE
  FeatureTable_dir_in<-file.path(InputDirectory, FeatureTable_POS)
  cat(paste0("Reading in file:\t", FeatureTable_POS,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
  FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
  # Create list of dataframes with all ms2 data
  
  length_ddMS2POS_in <- length(ddMS2POS_in)
  MS2_df_list <- vector("list", length_ddMS2POS_in)
  nrow_all_scans <- 0
  #decipher .MS2 files OR (else) mzML files
  
  if (length(ddMS2POS_in_ms2)>0){
    for (c in 1:length(ddMS2POS_in_ms2)){
      MS2_dir_in <- file.path(InputDirectory, ddMS2POS_in_ms2[c])
      MS2_df_in <- createDataFrame(MS2_dir_in)
      MS2_df_list[[c]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
      }
    }
  }
  
  if (length(ddMS2POS_in_mzML)>0){
    for (c in 1:length(ddMS2POS_in_mzML)){
      MS2_dir_in <- file.path(InputDirectory, ddMS2POS_in_mzML[c])
      MS2_df_in_peaks <- openMSfile(MS2_dir_in)
      MS2_df_in_h <- header(MS2_df_in_peaks)
      MS2_df_in_h <- MS2_df_in_h[MS2_df_in_h$msLevel==2,c("precursorMZ","retentionTime")]
      MS2_df_in_h[,2] <- MS2_df_in_h[,2]/60
      PEAKS = lapply(as.numeric(row.names(MS2_df_in_h)), function(i) peaks(MS2_df_in_peaks, i) )
      MS2_df_in<-cbind(MS2_df_in_h,as.matrix(PEAKS))
      colnames(MS2_df_in) <- c("precursor", "rt","mz_intensity")
      NonEmptyMS2spectra<-as.logical(lapply(1:nrow(MS2_df_in), function(i) length(unlist(MS2_df_in[i,3]))>2)) #throw away things with only 1 row might not want this (:
      MS2_df_in<-MS2_df_in[NonEmptyMS2spectra,]
      MS2_df_list[[c + length(ddMS2POS_in_ms2)]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
      }
    }
  }
  
  # Add scans to master MSMS export matrix
  Pos_rawMSMS <- matrix("", nrow_all_scans, 9)
  colnames(Pos_rawMSMS)[1:9] <- c("Arbitrary_Identifier", "Feature", "File", "Selected_RT", "Selected_mz", "mz", "Intensity", "Fragments", "LibraryFile")
  start <- 1
  for (c in 1:length_ddMS2POS_in){
    # message(paste(c, "/", length_ddMS2NEG_in), sep = "")
    nrow_MS2_df_list_c <- nrow(MS2_df_list[[c]])
    MS2_RT <- as.numeric(MS2_df_list[[c]][,2])
    MS2_MZ <- as.numeric(MS2_df_list[[c]][,1])
    for(j in 1:nrow_MS2_df_list_c){
      scan_rows <- nrow(MS2_df_list[[c]][j,][[3]][[1]])
      # start <- which(Pos_rawMSMS[, 1] == "")[1]
      stop <- start + scan_rows - 1
      Pos_rawMSMS[start:stop, 1] <- j
      Pos_rawMSMS[start:stop, 2] <- ""
      Pos_rawMSMS[start:stop, 3] <- ddMS2NEG_in[c]
      Pos_rawMSMS[start:stop, 4] <- MS2_RT[j]
      Pos_rawMSMS[start:stop, 5] <- MS2_MZ[j]
      Pos_rawMSMS[start:stop, 6:7] <- MS2_df_list[[c]][j,][[3]][[1]]
      Pos_rawMSMS[start:stop, 8] <- ""
      Pos_rawMSMS[start:stop, 9] <- ""
      start <- stop + 1
    }
  }
  write.csv(Pos_rawMSMS, file.path(InputDirectory, "Output/Pos_rawMSMS_NoAnnotations.csv"), row.names = FALSE, col.names = TRUE, na = "")
  # Run RunTargeted
  for (c in 1:length_ddMS2POS_in){
    # message(paste("c = ", c, sep = ""))
    cat(paste0("Reading in file:\t", ddMS2POS_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
    ExtraSample<-ExtraSampleNameddMSPOS_in[c]
    OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in)
    Pos_rawMSMS_List <- c()
    nrow_LibraryCriteria <- nrow(LibraryCriteria)
    if (ParallelComputing == TRUE) {
      foreach (i = seq_len(nrow(LibraryCriteria)), .packages = c("sqldf")) %dopar% {
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSPos_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2POS_in[c], "Pos")
      }
    } else {
      for(i in seq_len(nrow_LibraryCriteria)){
        # message(paste(i, "/", nrow_LibraryCriteria, sep = ""))
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSPos_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2POS_in[c], "Pos")
      }
    }
  }
  print(paste("Finished Posative ddMS analysis", timestamp(), sep = " --> "))
  # Export MSMS combined file
  # write.csv(Pos_rawMSMS, file.path(InputDirectory, "Output/Pos_rawMSMS_Example2.csv"), row.names = FALSE, col.names = TRUE, na = "")
}


#POS BY CLASS
LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "POS" & toupper(LibraryCriteria[,6]) == "TRUE",] #subset LibraryCriteria to find positive class libraries
LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,4]) == "TRUE",] #subset LibraryCriteria to find ddMS libraries to run
if(runPosddMS && nrow(LibraryCriteria)>0){
  PosClassDDLib <- TRUE
  FeatureTable_dir_in<-file.path(InputDirectory, FeatureTable_POS)
  cat(paste0("Reading in file:\t", FeatureTable_POS,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
  FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
  length_ddMS2POS_in <- length(ddMS2POS_in)
  # Run RunTargeted
  for (c in 1:length_ddMS2POS_in){
    # message(paste("c = ", c, sep = ""))
    cat(paste0("Reading in file:\t", ddMS2POS_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
    ExtraSample<-ExtraSampleNameddMSPOS_in[c]
    OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in)
    nrow_LibraryCriteria <- nrow(LibraryCriteria)
    if (ParallelComputing == TRUE) {
      foreach (i = seq_len(nrow(LibraryCriteria)), .packages = c("sqldf")) %dopar% {
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSPosByClass_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2POS_in[c], "PosByClass")
      }
    } else {
      for(i in seq_len(nrow_LibraryCriteria)){
        # message(paste(i, "/", nrow_LibraryCriteria, sep = ""))
        LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
        OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
        ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
        ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
        if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
        if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
        RunTargeted(MS2_df_list[[c]], FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryddMSPosByClass_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo, ddMS2POS_in[c], "PosByClass")
      }
    }
  }
  print(paste("Finished Positive by class ddMS analysis", timestamp(), sep = " --> "))
  # Export MSMS combined file
  # write.csv(Pos_rawMSMS, file.path(InputDirectory, "Output/Pos_rawMSMS_Example3.csv"), row.names = FALSE, col.names = TRUE, na = "")
}

if(runPosddMS){
  print(paste("Creating Identifications for Positive Mode", timestamp(), sep = " --> "))
}
if(runPosddMS) { #& !runPosAIF){
  ddMS2directory<-file.path(OutputDirectoryddMSPos_in,"Confirmed_Compounds")
  Classdirectory<-file.path(OutputDirectoryddMSPosByClass_in,"Confirmed_Compounds")
  AIFdirectory<-"Nothing"
  ForceColumnNames(file.path(InputDirectory,FeatureTable_POS))
  CreateIDs(file.path(InputDirectory,FeatureTable_POS), ddMS2directory, Classdirectory, AIFdirectory, ImportLibPOS, OutputDirectory, PosDDLib, PosClassDDLib, PosAIFLib, "Pos")
  AppendFrag(CommentColumn, RowStartForFeatureTableData, InputDirectory, NegPos = "Pos", OutDir="Output/PosIDed_Fragments.csv", OutDirOnlyFrags="Output/Pos_OnlyIDs_Fragments.csv", InputDir_Append="Output/ddMS/Pos/Additional_Files", ID_name="Potential_IDs", FragName="Frags", nFrag="Num_Frags", fileNames="Files", ImportTable="Output/PosIDed.csv",firstFragAppender=TRUE)
  ##append by class (EPA MASTER in the case of PFAS as well) - in new columns
  if (Lipid == FALSE) {
    AppendFrag(CommentColumn, RowStartForFeatureTableData, InputDirectory, NegPos = "Pos", OutDir="Output/PosIDed_insilico.csv", OutDirOnlyFrags="Output/Pos_OnlyIDs_PredFrags.csv", InputDir_Append="Output/ddMS/PosByClass/Additional_Files", ID_name="PredictedFrag_IDs", FragName="Frags", nFrag="Num_Frags", fileNames="Files", ImportTable="Output/PosIDed_Fragments.csv",firstFragAppender=FALSE)
    ##Scores all the features, finds homologous series, and sorts data
    Scoring(file.path(InputDirectory,"Output/PosIDed_insilico.csv"), OutDir = file.path(InputDirectory,"Output"), Mass_col = MZColumn, Retention_col = RTColumn, RowStartForFeatureTableData-1, RT_flagging, KMD_Bin_Window = (PrecursorMassAccuracy*2), upper, lower, RepeatingUnits_dir)
  } else {
    Scoring(file.path(InputDirectory,"Output/PosIDed_Fragments.csv"), OutDir = file.path(InputDirectory,"Output"), Mass_col = MZColumn, Retention_col = RTColumn, RowStartForFeatureTableData-1, RT_flagging, KMD_Bin_Window = (PrecursorMassAccuracy*2), upper, lower, RepeatingUnits_dir)
  }
}

if(runNegddMS){
  print(paste("Creating Identifications for Negative Mode", timestamp(), sep = " --> "))
}

if(runNegddMS) { #& !runNegAIF){
  ddMS2directory <- file.path(OutputDirectoryddMSNeg_in,"Confirmed_Compounds")
  Classdirectory <- file.path(OutputDirectoryddMSNegByClass_in,"Confirmed_Compounds")
  AIFdirectory <- "Nothing"
  ForceColumnNames(file.path(InputDirectory,FeatureTable_NEG))
  CreateIDs(file.path(InputDirectory,FeatureTable_NEG), ddMS2directory, Classdirectory, AIFdirectory, ImportLibNEG, OutputDirectory, NegDDLib, NegClassDDLib, NegAIFLib, "Neg")
  #InputDirectory<-"C:/Users/Jeremy Koelmel/Desktop/Desktop/Innovative_Omics/CLIENTS/2024_USGS/JBB_QCs/Annotated/"; CommentColumn<-1; RowStartForFeatureTableData<-2; NegPos = "Neg"; OutDir="/Output/NegIDed_FIN.csv"; OutDirOnlyFrags="/Output/Neg_OnlyIDs_PredFrags.csv"; InputDir_Append="Output/ddMS/NegByClass/Additional_Files"; ID_name="PredictedFrag_IDs"; FragName="Frags"; nFrag="Num_Frags"; fileNames="Files"; ImportTable="/Output/NegIDed_Fragments.csv"
  AppendFrag(CommentColumn, RowStartForFeatureTableData, InputDirectory, NegPos = "Neg", OutDir="Output/NegIDed_Fragments.csv", OutDirOnlyFrags="Output/Neg_OnlyIDs_Fragments.csv", InputDir_Append="Output/ddMS/Neg/Additional_Files", ID_name="Potential_IDs", FragName="Frags", nFrag="Num_Frags", fileNames="Files", ImportTable="Output/NegIDed.csv", firstFragAppender=TRUE)
  ##append by class (EPA MASTER in the case of PFAS as well) - in new columns
  if (Lipid == FALSE) {
    AppendFrag(CommentColumn, RowStartForFeatureTableData, InputDirectory, NegPos = "Neg", OutDir="Output/NegIDed_insilico.csv", OutDirOnlyFrags="Output/Neg_OnlyIDs_PredFrags.csv", InputDir_Append="Output/ddMS/NegByClass/Additional_Files", ID_name="PredictedFrag_IDs", FragName="Frags", nFrag="Num_Frags", fileNames="Files", ImportTable="Output/NegIDed_Fragments.csv", firstFragAppender=FALSE)
    ##Scores all the features, finds homologous series, and sorts data
    Scoring(file.path(InputDirectory,"Output/NegIDed_insilico.csv"), OutDir = file.path(InputDirectory,"Output"), Mass_col = MZColumn, Retention_col = RTColumn, RowStartForFeatureTableData-1, RT_flagging, KMD_Bin_Window = (PrecursorMassAccuracy*2), upper, lower, RepeatingUnits_dir)
  } else {
    Scoring(file.path(InputDirectory,"Output/NegIDed_Fragments.csv"), OutDir = file.path(InputDirectory,"Output"), Mass_col = MZColumn, Retention_col = RTColumn, RowStartForFeatureTableData-1, RT_flagging, KMD_Bin_Window = (PrecursorMassAccuracy*2), upper, lower, RepeatingUnits_dir)
  }
}

##PJS Needs to fix

if(runNegddMS & runPosddMS){
  Neg <- read.csv(file.path(OutputDirectory, "NegIDed_Fragments.csv"), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  Neg <- as.matrix(Neg)
  ncol_Neg <- ncol(Neg)
  Pos <- read.csv(file.path(OutputDirectory, "PosIDed_Fragments.csv"), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  Pos <- as.matrix(Pos)
  ncol_Pos <- ncol(Pos)
  
  # Match Pos and Neg matrix sizes
  if (ncol_Neg > ncol_Pos) {
    extra <- matrix("", nrow(Pos), ncol_Neg - ncol_Pos)
    Pos <- cbind(Pos, extra)
  } else if (ncol_Pos > ncol_Neg) {
    extra <- matrix("", nrow(Neg), ncol_Pos - ncol_Neg)
    Neg <- cbind(Neg, extra)
  }

  # Fix alignment
  if (ncol_Neg > ncol_Pos) {
    Pos[, (ncol_Neg - 8):ncol_Neg] <- Pos[, (ncol_Pos - 8):ncol_Pos]
  } else if (ncol_Pos > ncol_Neg) {
    Neg[, (ncol_Pos - 8):ncol_Pos] <- Neg[, (ncol_Neg - 8):ncol_Neg]
  }

  Data <- rbind(Pos, Neg)
  end_col <- ncol(Data)
  if (!is.na(which(colSums(Data == "") == nrow(Data))[1] - 1)) {
    end_col <- (which(colSums(Data == "") == nrow(Data))[1] - 1)
  }
  Data <- Data[, 1:end_col]
  Data[is.na(Data)] <- 0
  Data <- remove_duplicates(Data, RowStartForFeatureTableData - 1, end_col - 5, end_col - 8, RT_Window, RTColumn)
  if (nrow(Data) != 0) {
    write.table(Data, file.path(OutputDirectory, "CombinedIDed_Fragments.csv"), sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    Scoring(file.path(InputDirectory,"Output/CombinedIDed_Fragments.csv"), OutDir = file.path(InputDirectory,"Output"), Mass_col = MZColumn, Retention_col = RTColumn, RowStartForFeatureTableData-1, RT_flagging, KMD_Bin_Window = (PrecursorMassAccuracy*2), upper, lower, RepeatingUnits_dir)
  }
}

print(paste("Creating final MSMS export", timestamp(), sep = " --> "))

if (runNegddMS) {
  NegPosIDed_dir <- file.path(InputDirectory, "Output/NegIDed.csv")
  FeatureList_in_dir <- file.path(InputDirectory, "Output/NegIDed_FIN.csv")
  final_MSMS_export(NegPosIDed_dir, FeatureList_in_dir, "Neg")
  # write.csv(rawMSMS_df, file.path(InputDirectory, "Output/Neg_rawMSMS.csv"), row.names = FALSE, col.names = TRUE, na = "")
} 

if (runPosddMS) {
  NegPosIDed_dir <- file.path(InputDirectory, "Output/PosIDed.csv")
  FeatureList_in_dir <- file.path(InputDirectory, "Output/PosIDed_FIN.csv")
  final_MSMS_export(NegPosIDed_dir, FeatureList_in_dir, "Pos")
  # write.csv(rawMSMS_df, file.path(InputDirectory, "Output/Pos_rawMSMS.csv"), row.names = FALSE, col.names = TRUE, na = "")
}

options(warn=0)#suppress warning off

if (IMfirst==TRUE) {
  ScoreMSMS_colName<-"Score"
  ScoreIM_colName<-"Score.1"
  Name_or_ClassMSMS_colName<-"Name_or_Class"
  Name_or_ClassIM_colName<-"Name_or_Class.1"
  FormulaMSMS_colName<-"Formula"
  FormulaIM_colName<-"Formula.1"
  SMILESMSMS_colName<-"SMILES"
  SMILESIM_colName<-"SMILES.1"
  #Note the renaming of columns would be much more stable as soon as the feature table is imported to not assume overwrite order
  CombineIMandMSMS(OutputDirectory, ScoreMSMS_colName,ScoreIM_colName,Name_or_ClassMSMS_colName,Name_or_ClassIM_colName,FormulaMSMS_colName, FormulaIM_colName,	SMILESMSMS_colName, SMILESIM_colName)
}


if (IMfirst==TRUE) {
  stop("\r Stopping code here (Code complete), EIC and Mobiligrams, formula prediction, and Kaufmann plots for visualizer from the FluoroMatch IM step")
}

###########################STATS#####################################################


smpls<-c()
if (runNegddMS) {
  if (FLOW || csvInput) {
    if (file.exists(GroupCSVDirectory)&&file.size(GroupCSVDirectory)>2) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Neg")
    }
    else
    {
      statsFilepath<-DetermineFilePath("Neg")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      if (Lipid==FALSE) {
        for (i in 1:length(dataSet$PredictedFrag_IDs)) {
          if (grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\CYP.*", "CYP", dataSet$PredictedFrag_IDs[i])
          }
          if (grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\EC 3.*", "EC3", dataSet$PredictedFrag_IDs[i])
          }
        }}
      output<-dataSet
    }
  }else{
    if (length(GroupCSVDirectory)>0) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Neg")
    }
    else{
      statsFilepath<-DetermineFilePath("Neg")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      if (Lipid==FALSE) {
        for (i in 1:length(dataSet$PredictedFrag_IDs)) {
          if (grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\CYP.*", "CYP", dataSet$PredictedFrag_IDs[i])
          }
          if (grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\EC 3.*", "EC3", dataSet$PredictedFrag_IDs[i])
          }
        }}
      output<-dataSet
    }
  }
  #provide file path for importing and exporting
  outputFile<-DetermineFilePath("Neg")
  write.csv(output,outputFile, row.names = FALSE)
}
if (runPosddMS) {
  if (FLOW || csvInput) {
    if (file.exists(GroupCSVDirectory)&&file.size(GroupCSVDirectory)>2) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Pos")
    }
    else
    {
      statsFilepath<-DetermineFilePath("Pos")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      if (Lipid==FALSE) {
        for (i in 1:length(dataSet$PredictedFrag_IDs)) {
          if (grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\CYP.*", "CYP", dataSet$PredictedFrag_IDs[i])
          }
          if (grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\EC 3.*", "EC3", dataSet$PredictedFrag_IDs[i])
          }
        }}
      output<-dataSet
    }
  }else{
    if (length(GroupCSVDirectory)>0) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Pos")
    }
    else{
      statsFilepath<-DetermineFilePath("Pos")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      if (Lipid==FALSE) {
        for (i in 1:length(dataSet$PredictedFrag_IDs)) {
          if (grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\CYP.*", "CYP", dataSet$PredictedFrag_IDs[i])
          }
          if (grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
            dataSet$PredictedFrag_IDs[i]<-gsub("\\EC 3.*", "EC3", dataSet$PredictedFrag_IDs[i])
          }
        }}
      output<-dataSet
    }
  }
  #provide file path for importing and exporting
  outputFile<-DetermineFilePath("Pos")
  write.csv(output,outputFile, row.names = FALSE)
}
#PJS Needs to Fix
# if (runPosddMS&&runNegddMS) {
#   if (FLOW || csvInput) {
#     if (file.exists(GroupCSVDirectory)&&file.size(GroupCSVDirectory)>2) {
#       output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Combined")
#     }
#     else
#     {
#       statsFilepath<-DetermineFilePath("Combined")
#       data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
#       dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
#       output<-dataSet
#     }
#   }else{
#     if (length(GroupCSVDirectory)>0) {
#       output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Combined")
#     }
#     else{
#       statsFilepath<-DetermineFilePath("Combined")
#       data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
#       dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
#       output<-dataSet
#     }
#   }
#   #provide file path for importing and exporting
#   outputFile<-DetermineFilePath("Combined")
#   write.csv(output,outputFile, row.names = FALSE)
# }

#######################EIC and MS1 R Scripts###################################
#Michael Kummer
#Nandarani Abril
#Eric Siegel
#Jonathan Sparks

print("Beginning generation of EICs and Annotated MS1s")

if(FLOW==TRUE) {
  target_mzXML = paste(dirname(dirname(OutputDirectory)),"/Temp_Work/",sep="")
} else { ##Modular
  target_mzXML = InputDirectory
}

if(TargetEIC_Only==TRUE) {
  mzXML_Files<-list.files(target_mzXML, pattern = ".mzXML|mzML", full.names = FALSE)
  Target_mzXML<-grep(mzXML_Files, pattern = "Target|target|TARGET|Blank|blank|BLANK", invert=FALSE, value=TRUE)
} else {
  Target_mzXML<-list.files(target_mzXML, pattern = ".mzXML|mzML", full.names = FALSE)
}

#mz,rt,id,dt(optional isIM=True),formula
column_names<-c("m.z","Retention.Time","row.ID","Formula")
if(runNegddMS){
  Column_Indices<-get_column_indices(paste0(OutputDirectory,"/NegIDed_FIN.csv"), column_names)
} else {
  Column_Indices<-get_column_indices(paste0(OutputDirectory,"/PosIDed_FIN.csv"), column_names)
}

if (EICcpp==TRUE) {
  arguments = construct_EM_arguments(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Tolerances = c(0.1, 0.5)
    ,DT_Tolerances = c(0.1, 3)
    ,isIM = FALSE 
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = c(Column_Indices[1],Column_Indices[2],Column_Indices[3],Column_Indices[4]) #mz,rt,id,dt(optional isIM=True),formula
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = paste(InputLibrary,"/Scripts/secondary_isotopes.csv",sep="")
    ,path_to_mzXML_Files = target_mzXML
  )
} else {
  arguments = construct_EM_arguments(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Window = RT_Window
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = c(Column_Indices[1],Column_Indices[2],Column_Indices[3],Column_Indices[4]) #mz,rt,id,formula
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = paste(InputLibrary,"/Scripts/secondary_isotopes.csv",sep="")
  )
  arguments$path_to_mzXML_Files = target_mzXML
}

runallEM(Target_mzXML, arguments, runmode = runNegddMS)
runallEM(Target_mzXML, arguments, runmode = runPosddMS, isNeg = FALSE)

#######################Kaufmann###################################
#Jeremy Koelmel

print("Beginning Kaufmann Analysis")

#for testing (source both kaufmann scripts, and Manual_Review)
# InputLibrary<-"C:/SOFTWARE/FluoroMatch_5.4/Flow/LipidMatch_Distribution/LipidMatch_Libraries/"
# outputFile<-"D:/FluoroMatch_Data/2024_12_30_NIST_Revident_Formula_Prediction_Manuscript/NegIDed_FIN_NIST_C_annotated_Dec2024_Kauff.csv"
# RowStartForFeatureTableData <- 2
# library("MASS")
# library("data.table")

#Kaufmann eC, and ratios, values
MD_col_name <- "mass.defect"
MZ_col_name <- "m.z"
C13_col_name <- "13C1"

if (runPosddMS) {
  outputFile<-DetermineFilePath("Pos")
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot information complete, positive mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
  if (Lipid == FALSE && Tween_pos == FALSE) {
    ReferenceKauffData_csv <- paste(InputLibrary,"/Scripts/2020_EPA_MasterList_Kaufmann_Plots.csv",sep="")
    StepSize <- 200
    MZ_col_name <- "MZ_eC"
    MD_col_name <- "MD_eC"
    get_point_density(ReferenceKauffData_csv, outputFile, StepSize, MZ_col_name, MD_col_name)
  }
}
if (runNegddMS) {
  outputFile<-DetermineFilePath("Neg")
  MD_col_name <- "mass.defect"
  MZ_col_name <- "m.z"
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot information complete, negative mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
  if (Lipid == FALSE && Tween_pos == FALSE) {
    ReferenceKauffData_csv <- paste(InputLibrary,"/Scripts/2020_EPA_MasterList_Kaufmann_Plots.csv",sep="")
    StepSize <- 200
    MZ_col_name <- "MZ_eC"
    MD_col_name <- "MD_eC"
    get_point_density(ReferenceKauffData_csv, outputFile, StepSize, MZ_col_name, MD_col_name)
  }
}

if(Lipid == TRUE) {
  if(runNegddMS){
    write.csv(0,paste0(OutputDirectory,'/Neg_Feature_pFormula_MS1s.csv'))
  }
  if(runPosddMS){
    write.csv(0,paste0(OutputDirectory,'/Pos_Feature_pFormula_MS1s.csv'))
  }
  stop("\r Stopping code here (Code complete)")
}


#######################Formula Prediction###################################
#David Schiessel

#for testing (run all parameters first after uploading the csv inputs, then)
#runNegddMS<-TRUE
#OutputDirectory<-"D:/FluoroMatch_Data/2024_06_26_Testing_OUTPUTS/2024_06_30_Shimadzu_NIST_mzML/MZ2/Annotated/Output/"
#q=-1
#Poltxt="-"
#ppmTol=(ppm_Window/2)
#source(paste(InputLibrary,"/Scripts/Formula_Prediction.R",sep=""))
# GLOBALS

print("Beginning formula prediction")

eList1 = c('C','H','N','O','S','F','Br','Cl')
eList2 = c('C','H','N','O','S','F','P')

MF_topN=10  #Choose how many MFs to store

# Formula Prediction can take a significant amount of time, if set to 0 the algorithm will try and predict formula for all potential PFAS (filtered by score)
Override_Predict<-1000 #set to 2 or higher to only predict formula for "n" top ranked features

if (runPosddMS) {
  adducts <- data.frame(Adduct=c('[M+H]+','[M]+'),
                        a=c('H1',NA),
                        d=c(NA,NA))
  fh_Feature_MS1s<-paste0(OutputDirectory,"/Pos_Feature_MS1s.csv")
  fh_Feature_IDList<-paste0(OutputDirectory,"/PosIDed_FIN.csv")
  fh_MS1<-paste0(OutputDirectory,'/Pos_Feature_pFormula_MS1s.csv')
  fh_IDed_FIN<-paste0(OutputDirectory,'/PosIDed_FIN.csv')
  Formula_Prediction(Override_Predict,fh_Feature_MS1s,fh_Feature_IDList,fh_MS1,fh_IDed_FIN,MF_topN,adducts,eList1,eList2,q=1,Poltxt="+",(ppm_Window/2))
}
if (runNegddMS) {
  adducts <- data.frame(Adduct='[M-H]-',
                        a=NA,
                        d='H1')
  fh_Feature_MS1s<-paste0(OutputDirectory,"/Neg_Feature_MS1s.csv")
  fh_Feature_IDList<-paste0(OutputDirectory,"/NegIDed_FIN.csv")
  fh_MS1<-paste0(OutputDirectory,'/Neg_Feature_pFormula_MS1s.csv')
  fh_IDed_FIN<-paste0(OutputDirectory,'/NegIDed_FIN.csv')
  Formula_Prediction(Override_Predict,fh_Feature_MS1s,fh_Feature_IDList,fh_MS1,fh_IDed_FIN,MF_topN,adducts,eList1,eList2,q=-1,Poltxt="-",ppmTol=(ppm_Window/2))
}

#######################Homologous Series Voting###################################
#Michael Kummer

# for testing
# #run get_column_indices function
# library(data.table)
# library(MetaboCoreUtils)
# runNegddMS<-TRUE
# runPosddMS<-FALSE
# OutputDirectory<-"D:/FluoroMatch_Data/2024_12_30_NIST_Revident_Formula_Prediction_Manuscript/"
# source("C:/SOFTWARE/FluoroMatch_5.4/Flow/LipidMatch_Distribution/LipidMatch_Libraries/Scripts/homologousVoting.R")

print("Beginning homologous series voting (which repeating formula are best representative of a series)")

if(runNegddMS){
  fullpath = paste0(OutputDirectory,"/NegIDed_FIN.csv")
  column_names_homovote<-c("SeriesType_Identifier","m.z", "Formula", "zFormula") #Also Formula
  Column_Indices_homovote<-get_column_indices(fullpath, column_names_homovote)
  charge <- -1
  homologous_voting(col_sti = Column_Indices_homovote[1]   #SeriesType_Identifier,
                    , col_mz = Column_Indices_homovote[2]  #m.z. Mass to Charge Ratio
                    , col_mf = Column_Indices_homovote[3] #Formula
                    , col_mfz = Column_Indices_homovote[4] #TopMF (z) Molecular Formula
                    , fn_FT =   fullpath
                    , SeriesFormula_ColName = "SeriesFormula"
                    , SeriesFormula_ColNameZ = "Seriesformula_TopMFz"
                    , SeriesError_ColName = "SeriesError_TopMFz"
                    , charge = charge)
}

if(runPosddMS){
  fullpath = paste0(OutputDirectory,"/PosIDed_FIN.csv")
  column_names_homovote<-c("SeriesType_Identifier","m.z", "Formula", "zFormula") #Also Formula
  Column_Indices_homovote<-get_column_indices(fullpath, column_names_homovote)
  charge <- 1
  homologous_voting(col_sti = Column_Indices_homovote[1]   #SeriesType_Identifier,
                    , col_mz = Column_Indices_homovote[2]  #m.z. Mass to Charge Ratio
                    , col_mf = Column_Indices_homovote[3] #Formula
                    , col_mfz = Column_Indices_homovote[4] #TopMF (z) Molecular Formula
                    , fn_FT =   fullpath
                    , SeriesFormula_ColName = "SeriesFormula"
                    , SeriesFormula_ColNameZ = "Seriesformula_TopMFz"
                    , SeriesError_ColName = "SeriesError_TopMFz"
                    , charge = charge)
}

print("Software Complete!!! Congratulations. Now upload all necessary files into the Visualizer platform for validation and unknown discovery. Please email us with any interesting findings!")


