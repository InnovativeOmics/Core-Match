##Jeremy Koelmel  jeremykoelmel@gmail.com
##Michael Kummer  mk.code.na@gmail.com
##Paul Stelben    paul.stelben@yale.edu
##Nick Kroeger    nkroeger.cs@gmail.com

rm(list = ls()) #remove all R objects from memory, this program can be memory intensive as it is dealing with huge datasets

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
TWeen_pos <- FALSE #PJS: for PolyMatch
FilterAbovePrecursor <- 1 #how far from the precursor should fragment masses be kept (e.g. if precursor is 700, should 702 be considered?)
TargetEIC_Only <- TRUE

#### END Mandatory Parameter to Change ####
#### END "Read Me" Section ####

#Force the directory for packages to be the one that is distributed not the local R directory
if(length(.libPaths())>1){
  R_DistrDir<-.libPaths()[2]
  .libPaths(R_DistrDir)
}

if (FLOW == TRUE) {
  csvInput <- TRUE
  ManuallyInputVariables <- FALSE
}


#Checks for updates, installs packagaes: "installr" "stringr" "sqldf" "gWidgets" "gWidgetstcltk" and "compiler"
# if(!require(installr)) {
#   install.packages("installr"); install.packages("stringr"); require(installr)}
# library(installr)
if (ParallelComputing == TRUE) {
  if("foreach" %in% rownames(installed.packages()) == FALSE) {install.packages("foreach", repos = "http://cran.us.r-project.org")}
  library(foreach)
  if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel", repos = "http://cran.us.r-project.org")}
  library(doParallel)
  DC <- as.numeric(detectCores())
  if (is.na(DC)) {
    DC <- makePSOCKcluster(4)
    registerDoParallel(DC)
  } else if (DC <= 4) {
    DC <- makePSOCKcluster(DC)
    registerDoParallel(DC)
  } else {
    DC <- makePSOCKcluster(DC-2)
    registerDoParallel(DC)
  }
  # Force the directory for parallel computing to be the one that is distributed not the local R directory (doPar error)
  clusterEvalQ(DC, .libPaths()[2])
}

if("tictoc" %in% rownames(installed.packages()) == FALSE) {install.packages("tictoc", repos = "http://cran.us.r-project.org")}
library(tictoc)
tic.clear()
# tic("Main")

if("sqldf" %in% rownames(installed.packages()) == FALSE) {install.packages("sqldf", repos = "http://cran.us.r-project.org")}
if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager", repos = "http://cran.us.r-project.org")


#if("Rdisop" %in% rownames(installed.packages()) == FALSE) {install.packages("Rdisop", repos = "http://cran.us.r-project.org")}

if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages("RSQLite", repos = "http://cran.us.r-project.org")}
# if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets", repos = "http://cran.us.r-project.org")}
# if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk", repos = "http://cran.us.r-project.org")}
if("agricolae" %in% rownames(installed.packages()) == FALSE) {install.packages("agricolae", repos = "http://cran.us.r-project.org")}
if (FLOW == FALSE && csvInput == FALSE && ManuallyInputVariables == FALSE) {
  require(gWidgets)
  require(gWidgetstcltk)
  options(guiToolkit="tcltk")
}

#library(Rdisop)
#BiocManager::install("Rdisop")
library(RSQLite)
library(sqldf)
library(agricolae)
library(Rdisop)
options(warn=-1)#suppress warning on

errorBox <- function(message) {
  window <- gwindow("Confirm")
  group <- ggroup(container = window)
  
  ## A group for the message and buttons
  inner.group <- ggroup(horizontal=FALSE, container = group)
  glabel(message, container=inner.group, expand=TRUE, icon="error")
  
  ## A group to organize the buttons
  button.group <- ggroup(container = inner.group)
  ## Push buttons to right
  addSpring(button.group)
  gbutton("ok", handler=function(h,...) dispose(window), container=button.group)
  return()
}

if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table", repos = "http://cran.us.r-project.org")}
library(data.table)
if("SearchTrees" %in% rownames(installed.packages()) == FALSE) {install.packages("SearchTrees", repos = "http://cran.us.r-project.org")}
library(SearchTrees)
if("comprehenr" %in% rownames(installed.packages()) == FALSE) {install.packages("comprehenr", repos = "http://cran.us.r-project.org")}
library(comprehenr)
if("mzR" %in% rownames(installed.packages()) == FALSE)  {
  BiocManager::install("mzR")}
library(mzR)

#define function to concatenate file paths and files together without issues due to presence of trailing slash in file path.
file.join = function(..., sep = .Platform$file.sep){
    gsub("//", "/", file.path(..., sep = sep))
}

if (ManuallyInputVariables==TRUE){
  
  #Retention Time plus or minus
  RT_Window <- .3 #window of .2 => +/- .1
  
  #parts-per-million window for matching the m/z of fragments obtained in the library to those in experimentally obtained
  ppm_Window <- 10 #window of 10 => +/- 5 ppm error
  
  #Tolerance for mass-to-charge matching at ms1 level (Window)
  PrecursorMassAccuracy<-0.01
  
  #Plus minus range for the mass given after targeting parent ions portrayed in excalibur to match exact mass of Lipid in library
  SelectionAccuracy<-1
  
  #Threshold for determining that the average signal intensity for a given MS/MS ion should be used for confirmation
  intensityCutOff<-1000
  
  #Feature Table information
  CommentColumn <- 1
  MZColumn <- 2
  RTColumn <- 3
  RowStartForFeatureTableData <- 2 #Look at your Feature Table (.csv)...What row do you first see numbers?
  
  
  #ddMS data? set this parameter. If not? Leave it.
  #The minimum number of scans required for the result to be a confirmation
  ScanCutOff<-1
  
  #Have AIF data? set these parameters. If not? Leave them.
  corrMin <-.6
  minNumberOfAIFScans <- 5
  
  # Input Directory for Feature Table and MS2 files ...must have forward slashes but nothing at the end of the directory
  InputDirectory<-"C:/Users/Jeremy Koelmel/Downloads/LipidMatch_Flow_3.5/LipidMatch_Flow_3.5/LipidMatch_Modular/ExampleData"
  
  # Input Directory for Libraries ...must have forward slashes but nothing at the end of the directory
  InputLibrary<-"C:/Users/Jeremy Koelmel/Downloads/LipidMatch_Flow_3.5/LipidMatch_Flow_3.5/LipidMatch_Modular/LipidMatch_Libraries_Acetate"
  
  ## PFAS Specific parameters
  # RT_flagging <- TRUE #JPK: for PFAS analysis
  # Mass defect filtering
  upper <- 0.12 #JPK: Upper limit for PFAS mass defect
  lower <- -0.11 #JPK: Lower limit for PFAS mass defect
  
  #the isotopes that you can add for MS1 isotope labeling and isotope ratio calculations, if you don't differentiate the mass number it will use all from "secondary_isotopes.csv"
  ISOstring <- "13C3;15N;33S;34S;Cl3;18O;Br3;29Si;30Si"
  
  #Check your parameters and see that they are all correct.
  #Then press ctrl+shift+s to execute all code.
  ####################### END MANUALLY INPUTTING VARIABLES SECTION #######################
  
}else if(csvInput == TRUE){
  ####################### Get input from csv VARIABLES SECTION ###############################
  #parametersDirectory
  parametersDir <- "C:/NEW_SOFTWARE/2023_24_UPDATES_FM_LM/RapidTest/"
  parametersFile <- file.join(parametersDir, "PARAMETERS.csv")
  parametersInput_csv <- read.csv(parametersFile, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  parametersInput_csv <- as.matrix(parametersInput_csv)
  ErrorOutput<-0
  #Retention Time plus or minus
  RT_Window <- as.numeric(parametersInput_csv[2,2]) #window of .2 => +/- .1
  
  #parts-per-million window for matching the m/z of fragments obtained in the library to those in experimentally obtained
  ppm_Window <- as.numeric(parametersInput_csv[4,2]) #window of 10 => +/- 5 ppm error
  
  #Tolerance for mass-to-charge matching at ms1 level (Window)
  PrecursorMassAccuracy <- as.numeric(parametersInput_csv[3,2])
  
  #Plus minus range for the mass given after targeting parent ions portrayed in excalibur to match exact mass of Lipid in library
  SelectionAccuracy <- as.numeric(parametersInput_csv[5,2])
  
  #Threshold for determining that the average signal intensity for a given MS/MS ion should be used for confirmation
  intensityCutOff <- as.numeric(parametersInput_csv[7,2])
  
  #Feature Table information
  CommentColumn <- as.numeric(parametersInput_csv[10,2])
  MZColumn <- as.numeric(parametersInput_csv[11,2])
  RTColumn <- as.numeric(parametersInput_csv[12,2])
  RowStartForFeatureTableData <- as.numeric(parametersInput_csv[13,2]) #Look at your Feature Table (.csv)...What row do you first see numbers?
  
  #ddMS data? set this parameter. If not? Leave it.
  #The minimum number of scans required for the result to be a confirmation
  ScanCutOff<-as.numeric(parametersInput_csv[6,2])
  
  #Have AIF data? set these parameters. If not? Leave them.
  corrMin <-as.numeric(parametersInput_csv[9,2])
  minNumberOfAIFScans <- as.numeric(parametersInput_csv[8,2])
  
  # Input Directory for Feature Table and MS2 files ...must have \\ (double backslash) at the end of the directory
  InputDirectory<-as.character(parametersInput_csv[14,2])
  InputDirectory = substring(InputDirectory,1, nchar(InputDirectory)-1)
  
  # Input Directory for Libraries ...must have \\ (double backslash) at the end of the directory
  InputLibrary<-as.character(parametersInput_csv[15,2])
  InputLibrary <- substring(InputLibrary,1, nchar(InputLibrary)-1)
  
  GroupCSVDirectory<-as.character(parametersInput_csv[20,2])
  
  ISOstring <- as.character(parametersInput_csv[23,2])
  
  ## PFAS Specific parameters
  RT_flagging <- TRUE #JPK: for PFAS analysis
  # Mass defect filtering
  upper <- 0.12 #JPK: Upper limit for PFAS mass defect
  lower <- -0.11 #JPK: Lower limit for PFAS mass defect
}else{
  ####################### Pop-up boxes input VARIABLES SECTION ###############################
  
  # check if R version is equal to, or between, version 2.0.3 and 3.3.3, otherwise present pop-up box warning
  # Comment out for now
  # if(!((as.numeric(paste(version$major,version$minor,sep=""))>=20.3) && (as.numeric(paste(version$major,version$minor,sep=""))<=33.3))) {
  # errorBox(message=paste("ERROR: R version must be equal to, or between, 2.0.3 and 3.3.3. Please download 3.3.3. You are using version: ", paste(version$major,version$minor,sep=".")))
  # stop(paste("R version must be equal to, or between, 2.0.3 and 3.3.3. Please download 3.3.3. You are using version: ", paste(version$major,version$minor,sep=".")))
  # }
  ISOstring <- "13C3;15N;33S;34S;Cl3;18O;Br3;29Si;30Si"
  ## Input Directory that holds each folder of .ms2 files & feature table (contains features, peak heights/areas, etc.) (.csv file)
  InputDirectory<-tk_choose.dir(caption="Input Directory of MS2 + Feature Tables")
  if(is.na(InputDirectory)){
    stop()
  }
  foldersToRun <- list.dirs(path=InputDirectory, full.names=FALSE, recursive=FALSE)
  ErrorOutput<-0
  for (i in seq_len(length(foldersToRun))){
    if(foldersToRun[i] == "Output"){
      ErrorOutput<-1
      errorBox(message=paste("ERROR: Remove your 'Output' folder\nfrom the current Input Directory:\n", InputDirectory))
      stop("Warning: Remove your 'Output' folder from the current Input Directory: ", InputDirectory)
    }
  }
  
  ## Input Directory for Libraries
  InputLibrary<-tk_choose.dir(caption="Input Directory of libraries (came with LipidMatch)")
  if(is.na(InputLibrary)){
    stop()
  }
  
  GroupCSVDirectory <- tk_choose.files(caption="If there are no groupings select CANCEL, Otherwise: \nInput file containing groupings \none column file with characters unique to each grouping\neach row as a seperate group",multi=FALSE)
  
  GetInputAndErrorHandle <- function(typeOfVariable, message, title){
    isValidInput <- FALSE
    inputVariable <- ginput(message=message, title=title,icon="question")
    
    while(!isValidInput){
      ##Retention Time plus or minus
      if(inputVariable == "d" || inputVariable == "D"){
        if(typeOfVariable=="RT"){ inputVariable <- 0.3
        }else if(typeOfVariable=="ppm"){ inputVariable <- 10
        }else if(typeOfVariable=="precMassAccuracy"){ inputVariable <- 0.01
        }else if(typeOfVariable=="selectionAccuracy"){ inputVariable <- 1
        }else if(typeOfVariable=="maxInt"){ inputVariable <- 1000
        }else if(typeOfVariable=="scanCutOff"){ inputVariable <- 1
        }else if(typeOfVariable=="FeatRT"){inputVariable<-0.1
        }else if(typeOfVariable=="FeatMZ"){inputVariable<-0.006
        }else if(typeOfVariable=="corrMin"){inputVariable<-0.6
        }else if(typeOfVariable=="minAIFScans"){inputVariable<-5
        }else if(typeOfVariable=="upper"){inputVariable<-.12
        }else if(typeOfVariable=="lower"){inputVariable<-(-0.11)
        }
        isValidInput <- TRUE
      }else if(suppressWarnings(!is.na(as.numeric(inputVariable)))){
        inputVariable <- as.numeric(inputVariable)
        isValidInput <- TRUE
      }else{
        inputVariable <- ginput(message=paste("Error! Invalid Input.\n\n",message), title=title, icon="error")
      }
    }
    return(inputVariable)
  }
  
  RT_Window <- GetInputAndErrorHandle("RT", message="Retention Time Window\n(Window of .3 => +/- .15)\nOr type \"d\" for default: .3", title="RT_Window")
  ppm_Window <- GetInputAndErrorHandle("ppm", message="Parts-per-million window for matching experimental and in-silico fragments m/z \n(Window of 10 => +/- 5)\nOr type \"d\" for default value: 10", title="ppm_Window")
  PrecursorMassAccuracy <- GetInputAndErrorHandle("precMassAccuracy", message="Mass accuracy window for matching experimental and in-silico precursors m.z\nfor full scan (precursor) mass matching \n(Window of .01 Da => +/- .005 Da)\nOr type \"d\" for default: .01", title="PrecMassAccuracyWindow")
  SelectionAccuracy <- GetInputAndErrorHandle("selectionAccuracy", message="MS/MS Isolation Window (Da) \n(For determining MS/MS scans for each feature)\nOr type \"d\" for default: 1", title="SelectionAccuracy")
  intensityCutOff <- GetInputAndErrorHandle("maxInt", message="Threshold for determining what the minimum signal intensity cut off for a\ngiven MS/MS ion should be (used for confirmations)\nOr type \"d\" for default: 1000", title="intensityCutOff")
  CommentColumn <- GetInputAndErrorHandle("CommentCol", message= "Feature Table Info: Comment Column/Row ID \n(look at feature table .csv, first column is 1)\nNote that the column should be the same in negative and positive feature tables \nNo default value.",title="CommentColumn")
  MZColumn <- GetInputAndErrorHandle("MZCol", message= "Feature Table Info: Mass-to-Charge Column \n(look at feature table .csv, first column is 1)\nNote that the column should be the same in negative and positive feature tables \nNo default value.", title="MZColumn")
  RTColumn <- GetInputAndErrorHandle("RTCol", message= "Feature Table Info: Retention Time Column \n(look at feature table .csv, first column is 1)\nNote that the column should be the same in negative and positive feature tables \nNo default value.", title="RTColumn")
  RowStartForFeatureTableData <- GetInputAndErrorHandle("RowStart", message= "Feature Table Info: What row does your numeric data start? \n(First row of .csv starts counting at 1)\nNote that the row should be the same in negative and positive feature tables \nNo default value.", title="RowStartForFeatureTableData")
  upper <- GetInputAndErrorHandle("upper", message= "what is the upper limit for mass defect filtering (scoring, wont remove)? \nOr type \"d\" for default: 0.12 \n(0.12 is the upper 95 percentile of EPA Masterlist mass defects)", title="upper")
  lower <- GetInputAndErrorHandle("lower", message= "what is the lower limit for mass defect filtering (scoring, wont remove)? \nOr type \"d\" for default: -0.11 \n(-0.11 is the lowest 5 percentile of EPA Masterlist mass defects)", title="lower")
  RT_flagging <- TRUE #JPK: will simpy be default, changeable in code for now
  hasAIF <- FALSE
  hasdd <- FALSE
  
  foldersToRun <- list.dirs(path=InputDirectory, full.names=FALSE, recursive=FALSE)
  if(length(foldersToRun)==0){
    lengthFoldersToRun <- 1 #if there are no subfolders, that means you have the faeture table and ms2s in that current directory, therefore, run analysis on those files.
  }else{
    lengthFoldersToRun <- length(foldersToRun)#run analysis on all subfolders
  }
  #checks for output folder, we don't want to overwrite your files
  for (i in seq_len(lengthFoldersToRun)){
    numOfAIFFiles <- vector()
    if(length(foldersToRun)==0){#we're in current (and only) folder that contains feature table and ms2
      fpath <- InputDirectory
      
    }else if(foldersToRun[i] == "Output"){
      fpath <- InputDirectory
      print(paste("Warning: Remove your 'Output' folder from the current Input Directory:", InputDirectory))
    }else{
      fpath <- file.path(InputDirectory, foldersToRun[i])
      
    }
    
    #separate ddMS and AIF
    numOfddFiles <- length(list.files(path=fpath, pattern="[dD][dD]", ignore.case=FALSE))
    numOfAIFFiles <- length(list.files(path=fpath, pattern="[aA][iI][fF]", ignore.case=FALSE))
    if(numOfAIFFiles > 0){ #if # of AIF .ms2 or .ms1 > 0, ask for variable input
      hasAIF <- TRUE
    }else{
      print(paste("Warning: Couldn't find .ms1 or .ms2 files. Make sure you have 'aif' in your .ms1 or .ms2 file name (for data independent analysis)."))
    }
    if(numOfddFiles > 0){
      hasdd <- TRUE
    }else{
      print(paste("Warning: Couldn't find data dependent .ms2 files. Make sure you have 'dd' in your .ms2 file name (for data dependent analysis)."))
    }
  }
  if(hasAIF){
    corrMin <- GetInputAndErrorHandle("corrMin", message= "Minimum adjusted R2 value for AIF confirmation\n(Used in correlating fragment ms1 and ms2 intensities).\nOr type \"d\" for default: 0.6", title="corrMin")
    minNumberOfAIFScans <- GetInputAndErrorHandle("minAIFScans", message= "Minimum number of scans for confirming fragments (AIF)\nOr type \"d\" for default: 5", title="minNumberOfAIFScans")
  }
  if(hasdd){#then has ddMS
    ScanCutOff <- GetInputAndErrorHandle("scanCutOff", message="Minimum number of scans required for the result to be a confirmation (ddMS)\nOr type \"d\" for default: 1", title="ScanCutOff")
  }
  
} ####################### END AUTO INPUTTING VARIABLES SECTION ####################



#Error handling for manual input section... specifically the Input Directory. Checks for Output folder, and stops.
if(ManuallyInputVariables==TRUE || csvInput == TRUE){
  #checks for output folder, we don't want to overwrite your files
  foldersToRun <- list.dirs(path=InputDirectory, full.names=FALSE, recursive=FALSE)
  if(length(foldersToRun)==0){
    lengthFoldersToRun <- 1 #if there are no subfolders, that means you have the faeture table and ms2s in that current directory, therefore, run analysis on those files.
  }else{
    lengthFoldersToRun <- length(foldersToRun)#run analysis on all subfolders
  }
  for (i in seq_len(length(foldersToRun))){
    if(foldersToRun[i] == "Output"){
      errorBox(message=paste("Warning: Remove your 'Output' folder\nfrom the current Input Directory:\n", InputDirectory))
      stop("Warning: Remove your 'Output' folder from the current Input Directory: ", InputDirectory)
    }
    if(length(foldersToRun)==0){#we're in current (and only) folder that contains feature table and ms2
      fpath <- InputDirectory
    }else if(foldersToRun[i] == "Output"){
      fpath <- InputDirectory
      print(paste("Warning: Remove your 'Output' folder from the current Input Directory:", InputDirectory))
    }else{
      fpath <- file.path(InputDirectory, foldersToRun[i])
      
    }
    numOfAIFFiles <- length(list.files(path=fpath, pattern="[AIFaif][AIFaif][AIFaif]", ignore.case=FALSE))
    if(numOfAIFFiles %% 2 == 1){ #if # of AIF .ms2 or .ms1 are odd, STOP!
      stop("Warning: You should have a one-to-one relationship between AIF.ms1 and AIF.ms2 files...We found ", numOfAIFFiles," AIF .ms1 and/or .ms2 files. You should have an even number of AIF files.",sep="")
    }#else, do nothing.
  }
}
#### END VARIABLE DECLARATIONS SECTION ####


#### CODE - DO NOT TOUCH ####

#Library Information
# NAME AND DIRECTORY for exact mass library

ImportLibPOS<-file.path(InputLibrary, "Precursor_Library_POS.csv")
ImportLibNEG<-file.path(InputLibrary, "Precursor_Library_NEG.csv")
# if (Lipid == TRUE || FLOW == TRUE) {
# LibCriteria<- file.path(InputLibrary, "PFAS_ID_CRITERIA.csv")
# } else {
#   LibCriteria<- file.path(InputLibrary, "PFAS_ID_CRITERIA.csv")
# }

##Bernard Brooks 12/15/2022 9:05 AM - Error I found##
# if (FLOW == TRUE || Lipid == TRUE) {
#   LibCriteria<- file.path(InputLibrary, "LIPID_ID_CRITERIA.csv")
# } else {
LibCriteria<- file.path(InputLibrary, "PFAS_ID_CRITERIA.csv")
# }
LibColInfo<-1  #Look at your Library. What column are your IDs? (Columns to retrieve ID information and align with data)
ParentMZcol_in<-2 #Look at your Library. What column are your precursor masses in? (Columns to retrieve ID information and align with data)
ddMS2_Code<-"1"
AIF_Code<-"2"
Class_Code<-"3"
ExactMass_Code<-"4"
NoID_Code<-"5"
PrecursorMassAccuracy<-PrecursorMassAccuracy/2

PPM_CONST <- (10^6 + ppm_Window/2) / 10^6

################### FUNCTIONS ###################
RunAIF <- function(ms1_df, ms2_df, FeatureList, LibraryLipid_self, ParentMZcol, OutputDirectory, ExtraFileNameInfo, ConfirmORcol, ConfirmANDcol, OutputInfo){
  tStart<-Sys.time()
  ColCorrOR <- ConfirmORcol
  ColCorrAND <- ConfirmANDcol
  
  if(!dir.exists(file.path(OutputDirectory,"Confirmed_Compounds"))){
    dir.create(file.path(OutputDirectory,"Confirmed_Compounds"))
  }
  if(!dir.exists(file.path(OutputDirectory,"Additional_Files"))){
    dir.create(file.path(OutputDirectory,"Additional_Files"))
  }
  
  LibraryLipid<-read.csv(LibraryLipid_self, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  libraryLipidMZColumn <- ParentMZcol
  FeatureListMZColumn <- 1
  
  ## Create a Library of Lipids and fragments for those masses in Feature List
  FeatureListHeader <- FeatureList[1,]
  LibraryLipidHeader <- LibraryLipid[1,]
  NewLibraryLipidHeader <- cbind(LibraryLipidHeader, FeatureListHeader)
  
  sqlMerged <- sqldf(paste("select * from LibraryLipid lib inner join FeatureList iList on lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "-", PrecursorMassAccuracy, "<= iList.", colnames(FeatureList)[FeatureListMZColumn], " and iList.", colnames(FeatureList)[FeatureListMZColumn], " <= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "+", PrecursorMassAccuracy, sep = ""))
  # sqlMerged <- sqldf(paste("select * from LibraryLipid lib inner join FeatureList iList on iList.", colnames(FeatureList)[FeatureListMZColumn], " >= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], " -.005 and iList.", colnames(FeatureList)[FeatureListMZColumn], " <= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "+ .005", sep = ""))
  NewLibraryLipid <- as.matrix(sqlMerged)
  # if sqlMerged is not empty the names of the columns need to be adjusted, maybe better to always adjust them
  colnames(NewLibraryLipid) <- colnames(NewLibraryLipidHeader)
  NewLibraryLipid <- rbind(NewLibraryLipidHeader,NewLibraryLipid)
  NewLibraryLipid <- NewLibraryLipid[,-(ncol(LibraryLipid)+1)] #removes mz from Feature list
  
  #Removes 1 rt col and creates RT_min and RT_max cols based on RT_Window
  RT_Col <- as.numeric(NewLibraryLipid[2:nrow(NewLibraryLipid),ncol(NewLibraryLipid)-1])
  NewLibraryLipid <- NewLibraryLipid[,-(ncol(NewLibraryLipid)-1)]
  RT_Window_Vector <- rep(RT_Window/2, each=nrow(NewLibraryLipid)-1)
  
  #Issue found 10/30 Removed Levels function
  RT_min <- RT_Col - as.numeric(RT_Window_Vector)
  RT_max <- RT_Col + as.numeric(RT_Window_Vector)
  
  RT_min <- c("RT_min",RT_min)
  RT_max <- c("RT_max",RT_max)
  Comment <- NewLibraryLipid[,ncol(NewLibraryLipid)]
  NewLibraryLipid <- cbind(NewLibraryLipid[,1:(ncol(NewLibraryLipid)-1)], RT_min, RT_max, Comment)
  colnames(NewLibraryLipid) <- as.character(as.matrix(NewLibraryLipid[1,]))
  NewLibraryLipid <- NewLibraryLipid[-1,]
  NewLibraryLipid <- NewLibraryLipid[!is.na(NewLibraryLipid[,1]),!is.na(NewLibraryLipid[1,])]
  
  startRTCol <- ncol(NewLibraryLipid)-2
  endRTCol <- startRTCol + 1
  LibParentHits_df <- NewLibraryLipid #name change
  nrowLibParentHits <- nrow(LibParentHits_df)
  ncolLibParentHits <- ncol(LibParentHits_df)
  NoMatches_dir<-file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_NoIDs.csv"))
  if(nrowLibParentHits == 0){#No hits between Feature List and Library
    write.table("No Matches Found between Library and Feature List", NoMatches_dir, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }else{#If there was at least one match between Feature List and the Library
    #Output dataframes. Used to build the All Confirmed
    ConfirmedFragments_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    RTMaxIntensity_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    MaxIntensity_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    NumOfScans_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    AverageMZ_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    
    outputHeader <- colnames(LibParentHits_df) #get header from LibParentHits_df
    colnames(ConfirmedFragments_df) <- outputHeader
    colnames(RTMaxIntensity_df) <- outputHeader
    colnames(MaxIntensity_df) <- outputHeader
    colnames(NumOfScans_df) <- outputHeader
    colnames(AverageMZ_df) <- outputHeader
    
    outputLibData <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    ConfirmedFragments_df[,c(1,startRTCol:ncolLibParentHits)] <- outputLibData
    RTMaxIntensity_df[,c(1,startRTCol:ncolLibParentHits)] <- outputLibData
    MaxIntensity_df[,c(1,startRTCol:ncolLibParentHits)] <- outputLibData
    NumOfScans_df[,c(1,startRTCol:ncolLibParentHits)] <- outputLibData
    AverageMZ_df[,c(1,startRTCol:ncolLibParentHits)] <- outputLibData
    
    ########################################################################################
    # loop through lib parent hits M+H column, see                                         #
    #   1) if the precursor(from ms2) falls within this targetAcccuracy mz window          #
    #     2) if so... does this ms2 RT fall within the start and end RT from LibParentHits?#
    #     3) if mz frag from ms2 falls within ppm of LibParentHits                         #
    #   Then fragment is confirmed.                                                        #
    #                                                                                      #
    #   ms2_df[[1, 3]][1, 1]                                                               #
    #         ^matrix ^[row, column] from matrix                                           #
    ########################################################################################
    
    AIFms2SubsettedRTList <- list()
    for(p in seq_len(nrowLibParentHits)){
      RTConditional <- as.numeric(as.character(LibParentHits_df[p,startRTCol])) <= as.numeric(ms2_df[,2]) & as.numeric(ms2_df[,2]) <= as.numeric(as.character(LibParentHits_df[p,endRTCol]))
      AIFms2SubsettedRTList[[p]] <- subset(ms2_df, RTConditional)
    }
    AIFms1SubsettedRTList <- list()
    for(p in seq_len(nrowLibParentHits)){
      RTConditional <- as.numeric(as.character(LibParentHits_df[p,startRTCol])) <= as.numeric(ms1_df[,2]) & as.numeric(ms1_df[,2]) <= as.numeric(as.character(LibParentHits_df[p,endRTCol]))
      AIFms1SubsettedRTList[[p]] <- subset(ms1_df, RTConditional)
    }
    
    ############################################################################
    #Building All_Confirmed data frame for output (ddMS and/or AIF)            #
    # Creates a list of dataframes for each section of the All_Confirmed.csv.  #
    #   Why? To organize our data.                                             #
    # Each list's element corresponds to a row from LibParentHits              #
    # Each list's element holds a fragment found in the scans                  #
    ############################################################################
    #building ALL_Confirmed for AIF
    for(c in ParentMZcol:((startRTCol - ParentMZcol)+1)){ #loop all fragment columns
      subsettedFragList <- list() #has nrow(LibParentHits_df) elements
      for(p in 1:nrowLibParentHits){ #loop through all AIFms2SubsettedRTList elements
        # if(nrow(AIFms2SubsettedRTList[[p]]) != 0){ #don't loop if the list is empty! (for loops run once when the size is 0... :(...)
        subsettedFrag_df <- data.frame(scans=numeric(), rt=numeric(), frags=numeric(), intensity=numeric())
        scans <- vector()
        rt <- vector()
        frags <- vector()
        intensity <- vector()
        nrow_AIFms2SubsettedRTList <- nrow(AIFms2SubsettedRTList[[p]])
        for(d in seq_len(nrow_AIFms2SubsettedRTList)){ #loop through each element in the list
          subsettedFrag <- vector(length=0)
          #Take the theoretical LibParentHits fragment m/z and look within a ppm window through the scan's ms2 fragments
          fragConditional <- as.numeric(as.character(LibParentHits_df[p, c])) - abs(as.numeric(as.character(LibParentHits_df[p, c])) - as.numeric(as.character(LibParentHits_df[p, c]))*PPM_CONST) <= as.numeric(AIFms2SubsettedRTList[[p]][,3][[d]][,1]) & as.numeric(AIFms2SubsettedRTList[[p]][,3][[d]][,1]) <= as.numeric(as.character(LibParentHits_df[p, c])) + abs(as.numeric(as.character(LibParentHits_df[p, c])) - as.numeric(as.character(LibParentHits_df[p, c]))*PPM_CONST)
          subsettedFrag <- subset(AIFms2SubsettedRTList[[p]][,3][[d]], fragConditional)
          nrow_subsettedFrag <- nrow(subsettedFrag)
          if(nrow_subsettedFrag>1){#more than one matching frag
            for(s in 1:nrow_subsettedFrag){
              scans <- append(scans, as.numeric(AIFms2SubsettedRTList[[p]][d,1])) #gets scan number
              rt <- append(rt, as.numeric(AIFms2SubsettedRTList[[p]][d,2])) #get rt
            }
          }else if(nrow_subsettedFrag==1){#only 1 mataching frag
            scans <- append(scans, as.numeric(AIFms2SubsettedRTList[[p]][d,1])) #gets scan number
            rt <- append(rt, as.numeric(AIFms2SubsettedRTList[[p]][d,2])) #get rt
          }#else, length is 0. no fragments found. do nothing.
          frags <- append(frags, as.numeric(subsettedFrag[,1])) #gets fragment from ms2
          intensity <- append(intensity, as.numeric(subsettedFrag[,2])) #gets intensity associated with that fragment
        }
        subsettedFrag_df <- data.frame(scans, rt, frags, intensity)
        subsettedFragList[[p]] <- subsettedFrag_df
        
        #Stores information into respective dataframes for output as All_Confirmed
        if(nrow(subsettedFragList[[p]])>0){#loop over list elements that have values
          #Max Intensity
          maxIntensity <- max(subsettedFragList[[p]][,4])
          MaxIntensity_df[p, c] <- maxIntensity
          
          #Number of Scans
          numScans <- nrow(subsettedFragList[[p]])
          NumOfScans_df[p, c] <- numScans
          
          #Confirming fragments
          if(sum(subsettedFragList[[p]][,4] > intensityCutOff) > 0 & (numScans >= minNumberOfAIFScans)){ #at least one fragment's intensity is above user inputed threshold and number of scans is greater than or equal to user inputed, ScansCutOff
            ConfirmedFragments_df[p, c] <- 1 #1 for yes, 0 for no. (ConfirmedFragments_df was initialized with all 0s)
          }
          
          #RT at max intensity
          maxIntIndex <- which.max(subsettedFragList[[p]][,4])
          RTMaxInt <- subsettedFragList[[p]][maxIntIndex,2]
          RTMaxIntensity_df[p, c] <- RTMaxInt
          
          
          #Average m/z for confirmed fragment
          avgMZ <- mean(subsettedFragList[[p]][,3])
          AverageMZ_df[p, c] <- avgMZ
        }
      }
    }
    
    #Creates df for output "AllConfirmed"
    AllConfirmed_df <- data.frame(matrix("", ncol = ncolLibParentHits*5, nrow = nrowLibParentHits))
    
    blankColumn <- rep("", nrowLibParentHits)
    
    AllConfirmed_df <- cbind(blankColumn, ConfirmedFragments_df, blankColumn, RTMaxIntensity_df, blankColumn, NumOfScans_df, blankColumn, MaxIntensity_df, blankColumn, AverageMZ_df, blankColumn, LibParentHits_df)
    
    colnames(AllConfirmed_df)[1] <- paste("Confirmed fragments if intensity is above", intensityCutOff,"and number of scans is greater than or equal to", minNumberOfAIFScans)
    colnames(AllConfirmed_df)[ncolLibParentHits+2] <- "RT at Max Intensity"
    colnames(AllConfirmed_df)[ncolLibParentHits*2+3] <- "Number of Scans"
    colnames(AllConfirmed_df)[ncolLibParentHits*3+4] <- "Max Intensity"
    colnames(AllConfirmed_df)[ncolLibParentHits*4+5] <- "Average m/z"
    colnames(AllConfirmed_df)[ncolLibParentHits*5+6] <- "Theoretical, Library Parent Hits"
    
    #DON'T WRITE THIS OUT. WE want to write out the final confirmed with Adj R2 and Slope
    # ConfirmedAll_dir<-paste(OutputDirectory,"Additional_Files//", ExtraFileNameInfo,"_All_confirmed.csv",sep="")
    # write.table(AllConfirmed_df, ConfirmedAll_dir, sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    ## Code to reduce the dataset to lipids containing necessary fragments
    ## Start with both true, if no inputs then all rows are retained
    AllConfirmed_matrix<-as.matrix(AllConfirmed_df)
    OR_BinTruth<-1
    AND_BinTruth<-1
    ## For any given confirmed fragment/precursor
    ## add an extra row to avoid error in an apply function encase of a one row matrix
    DataCombined<-rbind(AllConfirmed_matrix,((1:ncol(AllConfirmed_matrix))*0))
    if ((length(ConfirmORcol)>0)&&(is.matrix(DataCombined[1:nrow(DataCombined),ConfirmORcol+1]))) {
      ## Sum all the elements from the binary confirmation table, if any fragment is 1 (above threshold & minimum # of scans) the element is TRUE
      OR_BinTruth<-as.numeric(c((apply(apply(DataCombined[1:nrow(DataCombined),ConfirmORcol+1],2,as.numeric),1,sum)>0)))
    }
    if ((length(ConfirmANDcol)>0)&&(is.matrix(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))) {
      ## Sum all the elements from the binary confirmation table, if ALL fragments are 1 (above threshold & minimum # of scans) the element is TRUE
      AND_BinTruth<-as.numeric(c((apply(apply(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1],2,as.numeric),1,sum)>(length(ConfirmANDcol)-1))))
    }
    if ((length(ConfirmANDcol)>0)&&(is.vector(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))) {
      AND_BinTruth<-as.numeric(c(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))
    }
    ## Reduce the data set to those values which are TRUE for both the OR and AND fragments calculated in the two IF statements above
    DataCombinedConfirmed<-DataCombined[(OR_BinTruth+AND_BinTruth)>1,,drop=FALSE]
    
    NoMatches_dir<-file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_NoIDs.csv"))
    if(nrow(DataCombinedConfirmed)==0){
      write.table("No Confirmed Compounds Found",NoMatches_dir, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }else{
      # ConfirmedReduced_dir <- paste(OutputDirectory,"Confirmed_Lipids//",ExtraFileNameInfo,"_reduced_confirmed.csv",sep="")
      # write.table(DataCombinedConfirmed, ConfirmedReduced_dir, sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      
      rowsToKeep <- as.numeric(row.names(DataCombinedConfirmed))-1 #used to reduce AIFms1SubsettedList and AIFms2SubsettedList
      AIFms1SubsettedRTList <- AIFms1SubsettedRTList[rowsToKeep]
      AIFms2SubsettedRTList <- AIFms2SubsettedRTList[rowsToKeep]
      
      DataCombinedConfirmed_df <- as.data.frame(DataCombinedConfirmed)#reduced confirmed df
      nrowReduced <- nrow(DataCombinedConfirmed_df)
      ncolReduced <- ncol(DataCombinedConfirmed_df)
      
      #gets the original theoretical library that was reduced
      LibInfoFromDataComb_df<- DataCombinedConfirmed_df[,(ncolReduced-ncolLibParentHits+1):ncolReduced]
      
      colsToCorr <- c(ColCorrAND, ColCorrOR)
      numColsForAIF <- 4 + length(colsToCorr) #4 because ID, RTmin, RTmax, Comment
      nameOfMiddleColsForAIF <- vector()
      for(c in colsToCorr){#get "precursor vs fragment" header for output
        nameOfMiddleColsForAIF <- append(nameOfMiddleColsForAIF, paste(colnames(LibParentHits_df[ParentMZcol]), " vs ", colnames(LibParentHits_df[c]), sep=""))
      }
      
      AIFHeader <- c(colnames(LibParentHits_df[1]), nameOfMiddleColsForAIF, colnames(LibParentHits_df[startRTCol:ncolLibParentHits])) #get header from LibParentHits_df
      
      AdjR2_df <- data.frame(matrix(0, ncol = numColsForAIF, nrow = nrowReduced))
      Slope_df <- data.frame(matrix(0, ncol = numColsForAIF, nrow = nrowReduced))
      colnames(AdjR2_df) <- AIFHeader
      colnames(Slope_df) <- AIFHeader
      outputLibData <- LibInfoFromDataComb_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from DataCombinedConfirmed_df
      AdjR2_df[,c(1,(numColsForAIF-2):numColsForAIF)] <- outputLibData
      Slope_df[,c(1,(numColsForAIF-2):numColsForAIF)] <- outputLibData
      
      #correlate precursor from ms1 against predefined columns to correlate against from ms2
      for(c in colsToCorr){ #loop fragment columns to correlate against ms1
        for(p in seq_len(nrowReduced)){#loop all rows of reduced confirmed data frame
          fragCorr_df<-NULL
          ms1ScanNum <- vector()
          ms2ScanNum <- vector()
          ms1Intensity <- vector()
          ms2Intensity <- vector()
          length_AIFms1SubsettedRTList <- length(AIFms1SubsettedRTList[[p]][,1])
          for(d in seq_len(length_AIFms1SubsettedRTList)){#store ms1 scan number and PRECURSOR intensity
            subsettedFrag <- vector(length=0)
            
            #Take the theoretical LibParentHits fragment m/z and look within a ppm window through the scan's ms2 fragments
            AIFms1fragConditional <- (as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol])) - abs(as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol])) - as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol]))*PPM_CONST) <= as.numeric(AIFms1SubsettedRTList[[p]][,3][[d]][,1])) & (as.numeric(AIFms1SubsettedRTList[[p]][,3][[d]][,1]) <= as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol])) + abs(as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol])) - as.numeric(as.character(LibInfoFromDataComb_df[p, ParentMZcol]))*PPM_CONST))
            subsettedFrag <- subset(AIFms1SubsettedRTList[[p]][,3][[d]], AIFms1fragConditional)
            
            #error handling if you find 2 fragments masses within a ppm window. Then what?
            if(nrow(subsettedFrag)>1){
              LibMass<-as.numeric(as.character(LibInfoFromDataComb_df[p, c]))
              fragMass<-as.numeric(subsettedFrag[,1])
              closestMass<-which(abs(LibMass-fragMass)==min(abs(LibMass-fragMass)))
              subsettedFrag<-subsettedFrag[closestMass,]
              subsettedFrag<-matrix(subsettedFrag,1,2)
              if(length(closestMass)>1){#edge case: if the two fragment masses are equidistant from the lib mass... then avg their intensities
                #average the intensities (and masses, but I don't use the mass later, so it doesn't matter)
                subsettedFrag<-matrix(c(mean(as.numeric(subsettedFrag[,1])),mean(as.numeric(subsettedFrag[,2]))), 1, 2)
              }
            }
            
            # print("------------MS1-----------------")
            # print(paste("d: ", d, "subsettedFrag:", subsettedFrag))
            if(nrow(subsettedFrag)==0){
              ms1Intensity <- append(ms1Intensity, 0)
            }else{# if there's 1 fragment
              ms1Intensity <- append(ms1Intensity, subsettedFrag[,2])
            }
            ms1ScanNum <- append(ms1ScanNum, as.numeric(AIFms1SubsettedRTList[[p]][d,1])) #gets scan number
          }
          length_AIFms2SubsettedRTList <- length(AIFms2SubsettedRTList[[p]][,1])
          for(d in seq_len(length_AIFms2SubsettedRTList)){#store ms2 intensity
            subsettedFrag <- vector(length=0)
            
            #Take the theoretical LibParentHits fragment m/z and look within a ppm window through the scan's ms2 fragments
            AIFms2fragConditional <- as.numeric(as.character(LibInfoFromDataComb_df[p, c])) - abs(as.numeric(as.character(LibInfoFromDataComb_df[p, c])) - as.numeric(as.character(LibInfoFromDataComb_df[p, c]))*PPM_CONST) <= as.numeric(AIFms2SubsettedRTList[[p]][,3][[d]][,1]) & as.numeric(AIFms2SubsettedRTList[[p]][,3][[d]][,1]) <= as.numeric(as.character(LibInfoFromDataComb_df[p, c])) + abs(as.numeric(as.character(LibInfoFromDataComb_df[p, c])) - as.numeric(as.character(LibInfoFromDataComb_df[p, c]))*PPM_CONST)
            subsettedFrag <- subset(AIFms2SubsettedRTList[[p]][,3][[d]], AIFms2fragConditional)
            if(nrow(subsettedFrag)>1){
              LibMass<-as.numeric(as.character(LibInfoFromDataComb_df[p, c]))
              fragMass<-as.numeric(subsettedFrag[,1])
              closestMass<-which(abs(LibMass-fragMass)==min(abs(LibMass-fragMass)))
              subsettedFrag<-subsettedFrag[closestMass,]
              subsettedFrag<-matrix(subsettedFrag,1,2)
              if(length(closestMass)>1){#edge case: if the two fragment masses are equidistant from the lib mass... then avg their intensities
                #average the intensities (and masses, but I don't use the mass later, so it doesn't matter)
                subsettedFrag<-matrix(c(mean(as.numeric(subsettedFrag[,1])),mean(as.numeric(subsettedFrag[,2]))), 1, 2)
              }
            }
            
            # print("------------MS2-----------------")
            # print(paste("d: ", d, "subsettedFrag:", subsettedFrag))
            if(nrow(subsettedFrag)==0){
              ms2Intensity <- append(ms2Intensity, 0)
            }else{
              ms2Intensity <- append(ms2Intensity, subsettedFrag[,2])
            }
            ms2ScanNum <- append(ms2ScanNum, as.numeric(AIFms2SubsettedRTList[[p]][d,1])) #gets scan number
          }
          ms1ScanNum<-as.numeric(ms1ScanNum)
          ms2ScanNum<-as.numeric(ms2ScanNum)
          ms1Intensity <- as.numeric(ms1Intensity)
          ms2Intensity <- as.numeric(ms2Intensity)
          ms1IntensityAveraged<-vector()
          ms1ScanNumAveraged <- vector()
          
          #average each ms1's intensity that are next to eachother so we can correlate it against ms2
          length_ms1Intensity <- length(ms1Intensity)
          for(a in seq_len(length_ms1Intensity-1)){#-1 so I don't go out of bounds
            temp <- ms1Intensity[a+1]#next intensity
            averaged <- mean(c(ms1Intensity[a], temp))
            ms1IntensityAveraged <- append(ms1IntensityAveraged, averaged)
          }
          length_ms1ScanNum <- length(ms1ScanNum)
          for(b in seq_len(length_ms1ScanNum-1)){#-1 so I don't go out of bounds
            temp <- ms1ScanNum[b+1]#next intensity
            averaged <- mean(c(ms1ScanNum[b], temp))
            ms1ScanNumAveraged <- append(ms1ScanNumAveraged, averaged)
          }
          
          intensityCol1 <- ms1IntensityAveraged[ms1ScanNumAveraged %in% ms2ScanNum]
          intensityCol2 <- ms2Intensity[ms2ScanNum %in% ms1ScanNumAveraged]
          
          #test for mininum number of couples (found fragments in both ms1 and ms2)
          intensitiesMultiplied <- intensityCol1*intensityCol2
          numOfCouples <- length(subset(intensitiesMultiplied, intensitiesMultiplied != 0))
          
          columnToFill <- match(c, colsToCorr) + 1   # +1 because the ID column in AdjR2_df is column 1.
          #are there min number of user specified couples (Definition: Couples (noun) := fragment found in both ms1 and ms2 at the same scan time)
          if(numOfCouples < minNumberOfAIFScans){#bad, you get no adj R2
            # AdjR2_df[p,columnToFill] <- paste("Not enough ms1 averaged scans & ms2 scans for a linear model. We found: ", numOfCouples, " ms1 and ms2 scan couples. You specified: ", minNumberOfAIFScans, " number of coupled scans.",sep="")
            # Slope_df[p,columnToFill] <- paste("Not enough ms1 averaged scans & ms2 scans for a linear model. We found: ", numOfCouples, " ms1 and ms2 scan couples. You specified: ", minNumberOfAIFScans, " number of coupled scans.",sep="")
            AdjR2_df[p,columnToFill] <- NA
            Slope_df[p,columnToFill] <- NA
          }else{
            fragCorr_df <- data.frame(ms1ScanNumAveraged, intensityCol1, intensityCol2)
            
            intensityFit<-lm(fragCorr_df[,3]~fragCorr_df[,2])
            intensitySummary <- summary(intensityFit)
            intensityAdjR2 <- intensitySummary$adj.r.squared
            intensitySlope <- intensityFit$coefficients[2]
            
            # if(intensityAdjR2 > corrMin){
            AdjR2_df[p,columnToFill] <- round(intensityAdjR2, 3)
            Slope_df[p,columnToFill] <- round(intensitySlope, 5)
            # }else{
            
            #   AdjR2_df[p,columnToFill] <- NA
            #   Slope_df[p,columnToFill] <- NA
            # }
            
            # print(fragCorr_df)
            ##We don't use this. But I'll keep it here. intensityCorr <- cor(fragCorr_df[,2],fragCorr_df[,3],use="complete.obs")##
          }
        }#end row loop
      }#end column loop
      # write the table out.
      
      #remove theoretical library from end of DataCombinedConfirmed_df
      DataCombinedConfirmed_df <- DataCombinedConfirmed_df[,-((ncolReduced-ncolLibParentHits):ncolReduced)] #CAREFUL NICK! THIS MESSES WITH THE HEADER!...but I kinda like the .1 .2 .3. .4
      blankColumn <- rep("", nrowReduced)
      DataCombinedConfirmed_df <- cbind(DataCombinedConfirmed_df, blankColumn, AdjR2_df, blankColumn, Slope_df, blankColumn, LibInfoFromDataComb_df)
      
      colnames(DataCombinedConfirmed_df)[ncolLibParentHits*5+6] <- paste("Adjusted R2, Precursor Intensity vs Fragment Intensity (AIF). AdjR2 > ", corrMin, sep="")
      colnames(DataCombinedConfirmed_df)[ncolLibParentHits*5+6 + numColsForAIF + 1] <- "Slope, Precursor Intensity vs Fragment Intensity (AIF)"
      colnames(DataCombinedConfirmed_df)[ncolLibParentHits*5+6 + numColsForAIF*2 + 2] <- "Theoretical, Library Parent Hits"
      
      ConfirmedReduced_dir <- file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_AdjR2Info.csv"))
      write.table(DataCombinedConfirmed_df, ConfirmedReduced_dir, sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      
      #Reduce matches based on ConfirmAndCol
      if(length(ConfirmANDcol) > 0){
        columnToFill<-vector()
        # AdjR2_df[] <- lapply(AdjR2_df, as.numeric)
        for(i in ConfirmANDcol){
          columnToFill <- append(columnToFill, match(i, ConfirmANDcol) + 1)
        }
        #get R2 rows with text in them, and remove them.
        if(!is.numeric(AdjR2_df[,columnToFill])){#if the adjR2 column is already numeric, (don't change to as.numeric...this gives an error)
          AdjR2_df[,columnToFill] <- lapply(AdjR2_df[,columnToFill], as.numeric)
        }
        if(!is.numeric(Slope_df[,columnToFill])){#if the adjR2 column is already numeric, (don't change to as.numeric...this gives an error)
          Slope_df[,columnToFill] <- lapply(Slope_df[,columnToFill], as.numeric)
        }
        #We want AdjR2 to be > corrMin and slopes to be positive!
        #set non positive slopes to NA
        Slope_df[,columnToFill] <- replace(Slope_df[,columnToFill],Slope_df[,columnToFill]<=0,NA)
        #remove negative slopes
        AdjR2_df <- AdjR2_df[apply(cbind(AdjR2_df[,columnToFill],Slope_df[,columnToFill]), 1, function(x) all(!is.na(x))),]
        AdjR2_df <- AdjR2_df[apply(AdjR2_df[,columnToFill, drop=F], 1, function(x) all(x > corrMin)),]
        rowsToKeep <- as.numeric(row.names(AdjR2_df)) #used to reduce final output table based on the reduced AdjR2_df rows from above ConfirmAND and ConfirmOR columns
        Slope_df <- Slope_df[rowsToKeep,]
        # AdjR2_df <- AdjR2_df[apply(AdjR2_df[columnToFill], 1, function(x) all(x > corrMin & !is.na(x))),]
        
        # AdjR2_df <- AdjR2_df[apply(AdjR2_df[c(2,3)], 1, function(x) all((x > corrMin) & is.numeric(x))),]
      }
      #difference here is in the apply function.... "all" vs "any"
      if(length(ConfirmORcol) > 0){
        columnToFill<-vector()
        # AdjR2_df[] <- lapply(AdjR2_df, as.numeric)
        length_ConfirmANDcol <- length(ConfirmANDcol)
        for(i in ConfirmORcol){
          columnToFill <- append(columnToFill, match(i, ConfirmORcol) + 1 + length_ConfirmANDcol)
        }
        if(nrow(AdjR2_df[,columnToFill]) !=0){#The AND could've removed all rows, so now we have nothing to work with.
          if(!is.numeric(AdjR2_df[,columnToFill])){#if the adjR2 column is already numeric, (don't change to as.numeric...this gives an error)
            AdjR2_df[,columnToFill] <- lapply(AdjR2_df[,columnToFill], as.numeric)
          }
          if(!is.numeric(Slope_df[,columnToFill])){#if the adjR2 column is already numeric, (don't change to as.numeric...this gives an error)
            Slope_df[,columnToFill] <- lapply(Slope_df[,columnToFill], as.numeric)
          }
          
          #We want AdjR2 to be > corrMin and slopes to be positive!
          #set non positive slopes to NA
          Slope_df[,columnToFill] <- replace(Slope_df[,columnToFill],Slope_df[,columnToFill]<=0,NA)
          #remove negative slopes
          AdjR2_df <- AdjR2_df[apply(cbind(AdjR2_df[,columnToFill],Slope_df[,columnToFill]), 1, function(x) all(!is.na(x))),]
          AdjR2_df <- AdjR2_df[apply(AdjR2_df[,columnToFill, drop=F], 1, function(x) any(x > corrMin)),]
        }
      }
      
      rowsToKeep <- as.numeric(row.names(AdjR2_df)) #used to reduce final output table based on the reduced AdjR2_df rows from above ConfirmAND and ConfirmOR columns
      DataCombinedConfirmed_df <- DataCombinedConfirmed_df[rowsToKeep,]
      
      ConfirmedReduced_dir <- file.path(OutputDirectory,"Confirmed_Compounds", paste0(ExtraFileNameInfo,"_IDed.csv"))
      if(nrow(DataCombinedConfirmed_df)>0){#don't output if there's nothing in the reduced file
        write.table(DataCombinedConfirmed_df, ConfirmedReduced_dir, sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
      }
      
      
    }#end else (if there are matches after reduction step)
    
    tEnd<-Sys.time()
    # Capture parameters used for this run and save them
    RunInfo<-matrix("",24,2)
    RunInfo[,1]<-c("Parameters Used","","Run Time","Computation Started:","Computation Ended:","Elapsed Time (largest time unit)","","Files","MS2 file name and Directory:","MS1 file name and Directory","Feature Table name and directory:","Library name and directory:","","MS/MS confirmation criteria","Library columns for confirmation (AND)","Library columns for confirmation (OR)","m/z ppm window:","RT window (min):","minimum scans:","minimum intensity:", "","AIF confirmation criteria","Minimum adjusted R2 value","Minimum scans:")
    RunInfo[,2]<-c("Values","","",as.character(tStart),as.character(tEnd),tEnd-tStart,"","",OutputInfo[1],OutputInfo[4],OutputInfo[3],LibraryLipid_self,"","",paste(ConfirmANDcol,collapse=","),paste(ConfirmORcol,collapse=","),ppm_Window,RT_Window,minNumberOfAIFScans,intensityCutOff,"","",corrMin,minNumberOfAIFScans)
    Info_dir<-file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_Parameters.csv"))
    write.table(RunInfo, Info_dir, sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
  }
}#end RunAIF

# ms2_df <- MS2_df_list[[c]]
# FeatureList <- FeatureList_in
# LibraryLipid_self <- file.path(InputLibrary, LibraryCriteria[i,1])
# ParentMZcol <- ParentMZcol_in
# OutputDirectory <- OutputDirectoryddMSNeg_in
# ExtraFileNameInfo <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="")
# ms2FileName <- ddMS2NEG_in[c]
# NegPos <- "Neg"


#Function for MS/MS ID
RunTargeted <- function(ms2_df, FeatureList, LibraryLipid_self, ParentMZcol, OutputDirectory, ExtraFileNameInfo, ConfirmORcol, ConfirmANDcol, OutputInfo, ms2FileName, NegPos){
  # message("RunTar")
  if(!dir.exists(file.path(OutputDirectory,"Confirmed_Compounds"))){
    dir.create(file.path(OutputDirectory,"Confirmed_Compounds"))
  }
  if(!dir.exists(file.path(OutputDirectory,"Additional_Files"))){
    dir.create(file.path(OutputDirectory,"Additional_Files"))
  }
  #store starting time
  tStart<-Sys.time()
  
  LibraryLipid<-read.csv(LibraryLipid_self, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  libraryLipidMZColumn <- ParentMZcol
  FeatureListMZColumn <- 1
  
  ## Create a Library of Lipids and fragments for those masses in Feature List
  FeatureListHeader <- FeatureList[1,]
  LibraryLipidHeader <- LibraryLipid[1,]
  NewLibraryLipidHeader <- cbind(LibraryLipidHeader, FeatureListHeader)
  
  sqlMerged <- sqldf(paste("select * from LibraryLipid lib inner join FeatureList iList on lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "-", PrecursorMassAccuracy, "<= iList.", colnames(FeatureList)[FeatureListMZColumn], " and iList.", colnames(FeatureList)[FeatureListMZColumn], " <= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "+", PrecursorMassAccuracy, sep = ""))
  # sqlMerged <- sqldf(paste("select * from LibraryLipid lib inner join FeatureList iList on iList.", colnames(FeatureList)[FeatureListMZColumn], " >= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], " -.005 and iList.", colnames(FeatureList)[FeatureListMZColumn], " <= lib.", colnames(LibraryLipid)[libraryLipidMZColumn], "+ .005", sep = ""))
  NewLibraryLipid <- as.matrix(sqlMerged)
  # if sqlMerged is not empty the names of the columns need to be adjusted, maybe better to always adjust them
  colnames(NewLibraryLipid) <- colnames(NewLibraryLipidHeader)
  NewLibraryLipid <- rbind(NewLibraryLipidHeader,NewLibraryLipid)
  NewLibraryLipid <- NewLibraryLipid[,-(ncol(LibraryLipid)+1)] #removes mz from Feature list
  
  #Removes 1 rt col and creates RT_min and RT_max cols based on RT_Window
  RT_Col <- as.numeric(NewLibraryLipid[2:nrow(NewLibraryLipid),ncol(NewLibraryLipid)-1])
  NewLibraryLipid <- NewLibraryLipid[,-(ncol(NewLibraryLipid)-1)]
  RT_Window_Vector <- rep(RT_Window/2, each=nrow(NewLibraryLipid)-1)
  
  #Issue found 10/30 Removed Levels function
  RT_min <- RT_Col - as.numeric(RT_Window_Vector)
  RT_max <- RT_Col + as.numeric(RT_Window_Vector)
  RT_min <- c("RT_min",RT_min)
  RT_max <- c("RT_max",RT_max)
  Comment <- NewLibraryLipid[,ncol(NewLibraryLipid)]
  NewLibraryLipid <- cbind(NewLibraryLipid[,1:(ncol(NewLibraryLipid)-1)], RT_min, RT_max, Comment)
  colnames(NewLibraryLipid) <- as.character(as.matrix(NewLibraryLipid[1,]))
  NewLibraryLipid <- NewLibraryLipid[-1,]
  NewLibraryLipid <- NewLibraryLipid[!is.na(NewLibraryLipid[,1]),!is.na(NewLibraryLipid[1,])]
  
  startRTCol <- ncol(NewLibraryLipid)-2
  endRTCol <- startRTCol + 1
  LibParentHits_df <- NewLibraryLipid #name change
  nrowLibParentHits <- nrow(LibParentHits_df)
  ncolLibParentHits <- ncol(LibParentHits_df)
  NoMatches_dir<-file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_NoIDs.csv"))
  if(nrowLibParentHits == 0){#No hits between Feature List and Library
    write.table("No Matches Found between Library and Feature List", NoMatches_dir, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }else{#If there was at least one match between Feature List and the Library
    #Output dataframes. Used to build the All Confirmed
    ConfirmedFragments_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    RTMaxIntensity_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    MaxIntensity_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    NumOfScans_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    AverageMZ_df <- data.frame(matrix(0, ncol = ncolLibParentHits, nrow = nrowLibParentHits))
    
    colnames(ConfirmedFragments_df) <- colnames(LibParentHits_df) #get header from LibParentHits_df
    colnames(RTMaxIntensity_df) <- colnames(LibParentHits_df) #get header from LibParentHits_df
    colnames(MaxIntensity_df) <- colnames(LibParentHits_df) #get header from LibParentHits_df
    colnames(NumOfScans_df) <- colnames(LibParentHits_df) #get header from LibParentHits_df
    colnames(AverageMZ_df) <- colnames(LibParentHits_df) #get header from LibParentHits_df
    
    ConfirmedFragments_df[,c(1,startRTCol:ncolLibParentHits)] <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    RTMaxIntensity_df[,c(1,startRTCol:ncolLibParentHits)] <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    MaxIntensity_df[,c(1,startRTCol:ncolLibParentHits)] <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    NumOfScans_df[,c(1,startRTCol:ncolLibParentHits)] <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    AverageMZ_df[,c(1,startRTCol:ncolLibParentHits)] <- LibParentHits_df[,c(1,startRTCol:ncolLibParentHits)] #gets data from LibParentHits_df
    
    ########################################################################################
    # loop through lib parent hits M+H column, see                                         #
    #   1) if the precursor(from ms2) falls within this targetAcccuracy mz window          #
    #     2) if so... does this ms2 RT fall within the start and end RT from LibParentHits?#
    #     3) if mz frag from ms2 falls within ppm of LibParentHits                         #
    #   Then fragment is confirmed.                                                        #
    #                                                                                      #
    #   ms2_df[[1, 3]][1, 1]                                                               #
    #         ^matrix ^[row, column] from matrix                                           #
    ########################################################################################
    
    #Create a list (1 element for each LibParentHits row) of subsetted ms2_df based on the
    # mz conditional: LibPrecursor - targetedAccuracy < ms2Precursor & ms2Precursor < LibPrecursor + targetedAccuracy
    # RT conditional: LibRTStart < ms2RT & ms2RT < LibRTEnd
    ms2SubsettedMZList <- list()
    ms2SubsettedRTList <- list()
    simplifiedAndSubsettedList <- list()
    for(p in 1:nrowLibParentHits){
      # MZConditional <- (((as.numeric(LibParentHits_df[p,2]) - abs(as.numeric(LibParentHits_df[p,2]) - as.numeric(LibParentHits_df[p,2])*PPM_CONST)) < as.numeric(ms2_df[,1]))) & (as.numeric(ms2_df[,1]) < (as.numeric(LibParentHits_df[p,2]) + abs(as.numeric(LibParentHits_df[p,2]) - as.numeric(LibParentHits_df[p,2])*PPM_CONST)))
      MZConditional <- (as.numeric(as.character(LibParentHits_df[p,ParentMZcol])) - SelectionAccuracy/2) <= as.numeric(ms2_df[,1]) & as.numeric(ms2_df[,1]) <= (as.numeric(as.character(LibParentHits_df[p,ParentMZcol])) + SelectionAccuracy/2)
      ms2SubsettedMZList[[p]] <- subset(ms2_df, MZConditional)
      RTConditional <- as.numeric(as.character(LibParentHits_df[p,startRTCol])) <= as.numeric(ms2_df[,2]) & as.numeric(ms2_df[,2]) <= as.numeric(as.character(LibParentHits_df[p,endRTCol]))
      ms2SubsettedRTList[[p]] <- subset(ms2_df, RTConditional)
      
      #Combines MZ and RT ms2Subsetted lists
      uniqueMZ <- unique(ms2SubsettedMZList[[p]][,1])
      simplifiedAndSubsettedList[[p]] <- subset(ms2SubsettedRTList[[p]], ms2SubsettedRTList[[p]][,1] %in% uniqueMZ)
    }
    
    ##############################################################
    #Building All_Confirmed data frame for output                #
    # Creates a list of dataframes for each column.              #
    #   Why? To organize our data.                               #
    # Each list's element corresponds to a row from LibParentHits#
    # Each list's element holds a fragment found in the scans    #
    ##############################################################
    for(c in ParentMZcol:((startRTCol - ParentMZcol)+1)){ #loop all fragment columns
      # message(paste("c = ", c, sep = ""))
      subsettedFragList <- list() #has nrow(LibParentHits_df) elements
      for(p in 1:nrowLibParentHits){ #loop through all simplifiedAndSubsettedList elements
        # message(paste("p = ", p, sep = ""))
        if(nrow(simplifiedAndSubsettedList[[p]]) != 0){ #don't loop if the list is empty! (for loops run once when the size is 0... :(...)
          subsettedFrag_df <- data.frame(scans=numeric(), rt=numeric(), frags=numeric(), intensity=numeric())
          scans <- vector()
          rt <- vector()
          frags <- vector()
          intensity <- vector()
          fragConditional <- vector(mode = "list", length = nrow(simplifiedAndSubsettedList[[p]]))
          nrow_simplifiedAndSubsettedList <- nrow(simplifiedAndSubsettedList[[p]])
          for(d in 1:nrow_simplifiedAndSubsettedList){ #loop through each element in the list
            # if (ncol(simplifiedAndSubsettedList[[p]][d,][[3]][[1]]) == 2) {
            #   frags_col <- vector(mode="character", length = nrow(simplifiedAndSubsettedList[[p]][d,][[3]][[1]])) ### NEW ###
            #   simplifiedAndSubsettedList[[p]][d,][[3]][[1]] <- cbind(simplifiedAndSubsettedList[[p]][d,][[3]][[1]], frags_col) ### NEW ###
            # }
            subsettedFrag <- vector(length=0)
            #Take the theoretical LibParentHits fragment m/z and look within a ppm window through the scan's ms2 fragments
            fragConditional[[d]] <- as.numeric(as.character(LibParentHits_df[p, c])) - abs(as.numeric(as.character(LibParentHits_df[p, c])) - as.numeric(as.character(LibParentHits_df[p, c]))*PPM_CONST) <= as.numeric(simplifiedAndSubsettedList[[p]][,3][[d]][,1]) & as.numeric(simplifiedAndSubsettedList[[p]][,3][[d]][,1]) <= as.numeric(as.character(LibParentHits_df[p, c])) + abs(as.numeric(as.character(LibParentHits_df[p, c])) - as.numeric(as.character(LibParentHits_df[p, c]))*PPM_CONST)
            subsettedFrag <- subset(simplifiedAndSubsettedList[[p]][,3][[d]], fragConditional[[d]])
            if(nrow(subsettedFrag)>1){#more than one matching frag
              for(s in 1:nrow(subsettedFrag)){ # Adds the scan and rt twice if 2 rows from the ms2 are matched for example
                scans <- append(scans, as.numeric(row.names(simplifiedAndSubsettedList[[p]][d,]))) #gets scan number
                rt <- append(rt, as.numeric(simplifiedAndSubsettedList[[p]][d,2])) #get rt
              }
            }else if(nrow(subsettedFrag)==1){#only 1 mataching frag
              scans <- append(scans, as.numeric(row.names(simplifiedAndSubsettedList[[p]][d,]))) #gets scan number
              rt <- append(rt, as.numeric(simplifiedAndSubsettedList[[p]][d,2])) #get rt
            }#else, length is 0. no fragments found. do nothing.
            frags <- append(frags, as.numeric(subsettedFrag[,1])) #gets fragment from ms2
            intensity <- append(intensity, as.numeric(subsettedFrag[,2])) #gets intensity associated with that fragment
          }
          subsettedFrag_df <- data.frame(scans, rt, frags, intensity)
          subsettedFragList[[p]] <- subsettedFrag_df
          
          #Stores information into respective dataframes for output as All_Confirmed
          if(nrow(subsettedFragList[[p]])>0){#loop over list elements that have values
            #Max Intensity
            maxIntensity <- max(subsettedFragList[[p]][,4])
            MaxIntensity_df[p, c] <- maxIntensity
            
            #Number of Scans
            numScans <- nrow(subsettedFragList[[p]])
            NumOfScans_df[p, c] <- numScans
            
            #Confirming fragments
            if(sum(subsettedFragList[[p]][,4] > intensityCutOff) > 0 & (numScans >= ScanCutOff)){ #at least one fragment's intensity is above user inputed threshold and number of scans is greater than or equal to user inputed, ScansCutOff
              ConfirmedFragments_df[p, c] <- 1 #1 for yes, 0 for no. (ConfirmedFragments_df was initialized with all 0s)
            }
            
            #RT at max intensity
            maxIntIndex <- which.max(subsettedFragList[[p]][,4])
            RTMaxInt <- subsettedFragList[[p]][maxIntIndex,2]
            RTMaxIntensity_df[p, c] <- RTMaxInt
            
            
            #Average m/z for confirmed fragment
            avgMZ <- mean(subsettedFragList[[p]][,3])
            AverageMZ_df[p, c] <- avgMZ
          }
        }
      }
    }
    
    #Creates df for output "AllConfirmed"
    AllConfirmed_df <- data.frame(matrix("", ncol = ncolLibParentHits*5, nrow = nrowLibParentHits))
    
    blankColumn <- rep("", nrowLibParentHits)
    
    AllConfirmed_df <- cbind(blankColumn, ConfirmedFragments_df, blankColumn, RTMaxIntensity_df, blankColumn, NumOfScans_df, blankColumn, MaxIntensity_df, blankColumn, AverageMZ_df, blankColumn, LibParentHits_df)
    colnames(AllConfirmed_df)[1] <- paste("Confirmed fragments if intensity is above", intensityCutOff,"and number of scans is greater than or equal to", ScanCutOff)
    colnames(AllConfirmed_df)[ncolLibParentHits+2] <- "RT at Max Intensity"
    colnames(AllConfirmed_df)[ncolLibParentHits*2+3] <- "Number of Scans"
    colnames(AllConfirmed_df)[ncolLibParentHits*3+4] <- "Max Intensity"
    colnames(AllConfirmed_df)[ncolLibParentHits*4+5] <- "Average m/z"
    colnames(AllConfirmed_df)[ncolLibParentHits*5+6] <- "Theoretical, Library Parent Hits"
    
    ConfirmedAll_dir<-file.path(OutputDirectory,"Additional_Files", paste0(ExtraFileNameInfo,"_All.csv"))
    write.table(AllConfirmed_df, ConfirmedAll_dir, sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    ## Code to reduce the dataset to lipids containing necessary fragments
    ## Start with both true, if no inputs then all rows are retained
    AllConfirmed_matrix<-as.matrix(AllConfirmed_df)
    OR_BinTruth<-1
    AND_BinTruth<-1
    ## For any given confirmed fragment/precursor
    ## add an extra row to avoid error in an apply function encase of a one row matrix
    DataCombined<-rbind(AllConfirmed_matrix,((1:ncol(AllConfirmed_matrix))*0))
    if ((length(ConfirmORcol)>0)&&(is.matrix(DataCombined[1:nrow(DataCombined),ConfirmORcol+1]))) {
      ## Sum all the elements from the binary confirmation table, if any fragment is 1 (above threshold & minimum # of scans) the element is TRUE
      OR_BinTruth<-as.numeric(c((apply(apply(DataCombined[1:nrow(DataCombined),ConfirmORcol+1],2,as.numeric),1,sum)>0)))
    }
    if ((length(ConfirmANDcol)>0)&&(is.matrix(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))) {
      ## Sum all the elements from the binary confirmation table, if ALL fragments are 1 (above threshold & minimum # of scans) the element is TRUE
      AND_BinTruth<-as.numeric(c((apply(apply(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1],2,as.numeric),1,sum)>(length(ConfirmANDcol)-1))))
    }
    if ((length(ConfirmANDcol)>0)&&(is.vector(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))) {
      AND_BinTruth<-as.numeric(c(DataCombined[1:nrow(DataCombined),ConfirmANDcol+1]))
    }
    ## Reduce the data set to those values which are TRUE for both the OR and AND fragments calculated in the two IF statements above
    DataCombinedConfirmed<-DataCombined[(OR_BinTruth+AND_BinTruth)>1,,drop=FALSE]
    
    NoMatches_dir<-file.path(OutputDirectory,"Additional_Files",paste0(ExtraFileNameInfo,"_NoID.csv"))
    if(nrow(DataCombinedConfirmed)==0){
      write.table("No Confirmed Compounds Found",NoMatches_dir, col.names=FALSE, row.names=FALSE, quote=FALSE)
    }else{
      ConfirmedReduced_dir <- file.path(OutputDirectory,"Confirmed_Compounds",paste0(ExtraFileNameInfo,"_IDed.csv"))
      write.table(DataCombinedConfirmed, ConfirmedReduced_dir, sep=",",col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    }
  }
  tEnd<-Sys.time()
  # Capture parameters used for this run and save them
  RunInfo<-matrix("",19,2)
  RunInfo[,1]<-c("Parameters Used","","Run Time","Computation Started:","Computation Ended:","Elapsed Time (largest time unit)","","Files","MS2 file name and Directory:","Feature Table name and directory:","Library name and directory:","","MS/MS confirmation criteria","Library columns for confirmation (AND)","Library columns for confirmation (OR)","m/z ppm window:","RT window (min):","minimum scans:","minimum intensity:")
  RunInfo[,2]<-c("Values","","",as.character(tStart),as.character(tEnd),tEnd-tStart,"","",OutputInfo[1],OutputInfo[3],LibraryLipid_self,"","",paste(ConfirmANDcol,collapse=","),paste(ConfirmORcol,collapse=","),ppm_Window,RT_Window,ScanCutOff,intensityCutOff)
  Info_dir<-file.path(OutputDirectory,"Additional_Files",paste0(ExtraFileNameInfo,"_Parameters.csv"))
  write.table(RunInfo, Info_dir, sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
}#end of function

ReadFeatureTable <- function(FeatureTable_dir){
  #Converts Feature Table to an "Feature List" which has only MZ, RT, and Comment columns from the Feature table.
  FeatureTable_df<- read.csv(FeatureTable_dir, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE, stringsAsFactors=FALSE)
  
  #Conversion: Feature Table -> Feature List
  IncMZCol <- as.numeric(FeatureTable_df[RowStartForFeatureTableData:nrow(FeatureTable_df),MZColumn])
  IncRTCol <- as.numeric(FeatureTable_df[RowStartForFeatureTableData:nrow(FeatureTable_df),RTColumn])
  IncCommentCol <- as.numeric(FeatureTable_df[RowStartForFeatureTableData:nrow(FeatureTable_df),CommentColumn])
  FeatureList <- matrix(data=c(IncMZCol, IncRTCol, IncCommentCol), nrow=length(IncMZCol), ncol=3)
  FeatureList <- rbind(c("M/Z", "RT", "Comment"), FeatureList)
  FeatureList <- as.data.frame(FeatureList)
  FeatureList<-FeatureList[c(TRUE,(!is.na(as.numeric(FeatureList[2:nrow(FeatureList),1])))),]
  return(FeatureList)
}

createDataFrame <- function (msx_dir){
  #################################
  #Store ms1/2 data into dataframe#
  #################################
  dir_string <- strsplit(msx_dir, split="")[[1]]
  Level <- as.numeric(dir_string[length(dir_string)])
  msx <- scan(file=msx_dir, what='character') #read in .ms1/2 as list of characters
  s_index <- match("S", msx)
  msx <- msx[s_index:length(msx)] #cut off head useless info from .ms1/2
  s_indicies <- which(msx=="S")
  RTime_indicies <- which(msx == "RTime")
  TIC_indicies <- which(msx == "TIC")
  if (length(grep("AIF", msx_dir)) > 0) {
    AIF <- TRUE
  } else {
    AIF <- FALSE
  }
  if (Level == "1" || AIF == TRUE) {
    msx_df <- data.frame(scanNum=numeric(), rt=numeric(), mz_intensity=I(list()))
  }else{
    msx_df <- data.frame(precursor=numeric(), rt=numeric(), mz_intensity=I(list()))
  }
  to_remove <- c()
  
  s_indicies <- append(s_indicies, length(msx) + 1)
  length_s_indicies <- length(s_indicies)
  for (i in 1:(length_s_indicies - 1)){
    if (msx[TIC_indicies[i] + 2] != "Z"){
      if ((s_indicies[i + 1] - TIC_indicies[i]) <= 4){
        to_remove <- append(to_remove, i)
        next
      }
      data <- msx[(TIC_indicies[i] + 2):(s_indicies[i + 1] - 1)]
      mz <- data[c(TRUE, FALSE)]
      intensity <- data[c(FALSE, TRUE)]
    } else {
      if ((s_indicies[i + 1] - TIC_indicies[i]) <= 7){
        to_remove <- append(to_remove, i)
        next
      }
      data <- msx[(TIC_indicies[i] + 5):(s_indicies[i + 1] - 1)]
      mz <- data[c(TRUE, FALSE)]
      intensity <- data[c(FALSE, TRUE)]
    }
    if (Level == "1" || AIF == TRUE) {
      msx_df[i, 1] <- msx[s_indicies[i] + 1]
    }else{
      msx_df[i, 1] <- msx[s_indicies[i] + 3]
    }
    msx_df[i, 2] <- msx[RTime_indicies[i] + 1]
    msx_df[[i, 3]] <- matrix(cbind(mz, intensity), length(mz), 2)
  }
  
  if (length(to_remove) > 0) {
    msx_df <- msx_df[-to_remove,]
  }
  
  return(msx_df)
}

# PeakTableDirectory <- file.path(fpath,FeatureTable_POS)
# ddMS2directory<-ddMS2directory
# ImportLib<-ImportLibPOS
# ddMS2<-PosDDLib
# ddMS2Class<-PosClassDDLib
# AIF<-PosAIFLib
# mode<-"pos"


CreateIDs <- function(PeakTableDirectory, ddMS2directory, Classdirectory, AIFdirectory, ImportLib, OutputDirectory, ddMS2, ddMS2Class, AIF, mode){
  # Generates IDs by matching Comment column from FeatureTable.csv and All_Confirmed.csv files
  # Feature Table Columns:
  # ID
  # Intensity
  # Class
  # Adduct
  # Notation if only 1 class was observed
  
  PeakTable<-read.csv(PeakTableDirectory,sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  LastCol<-ncol(PeakTable)
  IDcolumn<-LastCol+1
  ddMS2_Files<-list.files(ddMS2directory,pattern="*IDed.csv")
  class_Files<-list.files(Classdirectory,pattern="*IDed.csv")
  AIF_Files<-list.files(AIFdirectory,pattern="*IDed.csv")
  
  ##Debug
  # RunTF <- ddMS2
  # Directory <- ddMS2directory
  # Files <- ddMS2_Files
  # Code <- ddMS2_Code
  # isAIF <- FALSE
  
  #i <- 1
  
  getIDCommentAndIntensities <- function (RunTF,Directory,Files,Code,isAIF){
    Compiled<-data.frame(Lipid=character(), Feature=numeric(), SumOfAllFragmentsIntensities=numeric(), class=character())
    if(RunTF==TRUE & length(Files) > 0){
      for(i in seq_len(length(Files))){
        #Read in reduced_confirmed file
        TempFile <- read.csv(file.path(Directory,Files[i]),sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
        if(nrow(TempFile)>1){ #safety net to check if there are more rows than just the header
          CommentCol <- which(TempFile[1,]=="Comment")[1]#receive first comment column
          #get Max Intensity dataframe portion from All_Confirmed (AIF all_confirmed files have specific header)
          if(isAIF){
            MaxIntensityCols <- (which(TempFile[1,]=="Max Intensity")+2):(which(TempFile[1,]=="RT_min.3")-1)
          }else{
            MaxIntensityCols <- (which(TempFile[1,]=="Max Intensity")+2):((which(TempFile[1,]=="RT_min")[4])-1)
            #MaxIntensityCols <- (which(TempFile[1,]=="Max Intensity")+2):(0)
          }
          maxIntensity_df <- TempFile[2:nrow(TempFile),MaxIntensityCols]
          
          summedIntensity<-c()
          nrow_maxIntensity_df <- nrow(maxIntensity_df)
          for(i in seq_len(nrow_maxIntensity_df)){
            summedIntensity <- append(summedIntensity, sum(unique(as.numeric(as.matrix(maxIntensity_df[i,])))))
          }
          
          class<-rep(as.character(TempFile[1,2]), times=(nrow(TempFile)-1))#Gets class and repeats it
          adduct<-rep(as.character(TempFile[1,3]), times=(nrow(TempFile)-1))#Gets adduct and repeats it
          TempFile <- TempFile[2:nrow(TempFile),c(2, CommentCol)] #remove header and keep ID & Comment column
          TempFile <- cbind(TempFile, summedIntensity, class, adduct) #append summed intensity to TempFile
          colnames(TempFile)[1]<-"Lipid"
          colnames(TempFile)[2]<-"Feature"
          colnames(TempFile)[3]<-"SumOfAllFragmentsIntensities"
          colnames(TempFile)[4]<-"Class"
          colnames(TempFile)[5]<-"Adduct"
          
          Compiled <- rbind(Compiled, TempFile)
        }
      }
      Compiled[,1]<-paste(Code,"_",Compiled[,1],sep="")
    }
    return(Compiled)
  }
  compAIF_df <- getIDCommentAndIntensities(AIF,AIFdirectory,AIF_Files,AIF_Code, TRUE)
  
  #Issue found at 1200 - 10/13/2022
  compddMS2_df <- getIDCommentAndIntensities(ddMS2,ddMS2directory,ddMS2_Files,ddMS2_Code, FALSE)
  
  
  compddMS2class_df <- getIDCommentAndIntensities(ddMS2Class,Classdirectory,class_Files,Class_Code, FALSE)
  
  compiled<-rbind(compddMS2_df, compAIF_df, compddMS2class_df) #join all compiled dataframes together (AIF, class, ddms2)
  
  write.table(compiled, file.path(OutputDirectory, paste0(mode,"_OnlyIDs.csv")), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  
  #Create dataframe to append on the end of your feature table
  NOID<-paste(NoID_Code,"_NoID",sep="")
  appendedIDCol<-matrix(NOID, nrow(PeakTable),1)
  if(RowStartForFeatureTableData > 2){
    appendedIDCol[2:(RowStartForFeatureTableData-1)]<-"" #Start "_NoID" where your data/comments/features starts...therefore, don't put "_NoID" if there is no data there.
  }
  appendedIntensityCol<-matrix("", nrow(PeakTable),1)
  appendedClass<-matrix("", nrow(PeakTable),1)
  appendedAdduct<-matrix("", nrow(PeakTable),1)
  appendedOnlyOneClass<-matrix("", nrow(PeakTable),1)
  NewPeakTable<-cbind(PeakTable,appendedIDCol, appendedIntensityCol, appendedClass, appendedAdduct, appendedOnlyOneClass)
  
  #Handy function to convert factor to character so you can manipulate the data in cells
  factorToCharacter <- function(df){
    for(i in which(sapply(df, class) == "factor")) df[[i]] = as.character(df[[i]])
    return(df)
  }
  NewPeakTable[,IDcolumn:(IDcolumn+4)] <- factorToCharacter(as.data.frame(NewPeakTable[,IDcolumn:(IDcolumn+4)]))#factor to character on the last 5 columns
  NewPeakTable[1,IDcolumn:(IDcolumn+4)] <- c("ID_Ranked","Intensity_Ranked","Class_At_Max_Intensity","Adduct_At_Max_Intensity","Only_One_Class")#create header for the last 5 columns
  
  nrow_NewPeakTable <- nrow(NewPeakTable)
  for (i in RowStartForFeatureTableData:nrow_NewPeakTable) {
    TempCompiled <- compiled[(as.numeric(as.character(compiled[,2]))==as.numeric(as.character(NewPeakTable[i,CommentColumn]))),]
    if(nrow(TempCompiled)==1) {
      NewPeakTable[i,IDcolumn]<-TempCompiled[,1]
      NewPeakTable[i,IDcolumn+1] <- paste0(TempCompiled[,3], collapse=" | ") #intensity
      NewPeakTable[i,IDcolumn+2] <- paste0(TempCompiled[1,4], collapse=" | ") #get's class at max intensity
      NewPeakTable[i,IDcolumn+3] <- paste0(TempCompiled[1,5]) #get's adduct at max intensity
      NewPeakTable[i,IDcolumn+4] <- "Yes" #only one class. yes.
    }else if(nrow(TempCompiled)>1) {
      TempCompiled <- TempCompiled[with(TempCompiled, order(-SumOfAllFragmentsIntensities)), ]#sort intensities, max first
      
      rowsToKeep <- match(unique(TempCompiled[,1]),TempCompiled[,1]) #removes duplicates and keeps highest intensity
      TempCompiled <- TempCompiled[rowsToKeep,] #removes duplicates
      
      NewPeakTable[i,IDcolumn]   <- paste0(TempCompiled[,1], collapse=" | ")
      NewPeakTable[i,IDcolumn+1] <- paste0(TempCompiled[,3], collapse=" | ") #intensity
      NewPeakTable[i,IDcolumn+2] <- paste0(TempCompiled[1,4], collapse=" | ") #get's class at max intensity
      NewPeakTable[i,IDcolumn+3] <- paste0(TempCompiled[1,5]) #get's adduct at max intensity
      NewPeakTable[i,IDcolumn+4] <- ifelse(length(unique(TempCompiled[,4]))==1, "Yes", "No") #Only one class?: yes, else? no.
    }
  }
  
  Lib <- read.csv(ImportLib, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
  Lib <- as.matrix(Lib)
  #Exact mass matching for Precursor library to feature table m/zs
  NOID <- paste(NoID_Code,"_NoID",sep="")
  NumLibMZ <- as.numeric(Lib[,ParentMZcol_in])
  NumLibMZ[1] <- 0
  nrow_NewPeakTable <- nrow(NewPeakTable)
  for (i in RowStartForFeatureTableData:nrow_NewPeakTable) {
    NumData <- as.numeric(as.character(NewPeakTable[i,MZColumn]))
    TempID <- Lib[(NumData - PrecursorMassAccuracy < NumLibMZ) & (NumLibMZ < NumData + PrecursorMassAccuracy), LibColInfo]
    # TempID <- Lib[(NumData-abs(NumData-NumData*PPM_CONST) < NumLibMZ) & (NumLibMZ < NumData + (abs(NumData-NumData*PPM_CONST))), LibColInfo]
    if ((length(TempID)==1)&&(NewPeakTable[i,IDcolumn]==NOID)){
      NewPeakTable[i,IDcolumn]<-paste0(ExactMass_Code,"_",TempID,sep="")
    }else if ((length(TempID)>1)&&(NewPeakTable[i,IDcolumn]==NOID)){
      NewPeakTable[i,IDcolumn]<- paste0(paste(ExactMass_Code,"_",TempID, sep=""),collapse=" | ")
    }
  }
  
  #Sorting based on confirmation numbers
  dataRange <- RowStartForFeatureTableData:nrow(NewPeakTable) #defining a constant
  TempSorted <- NewPeakTable[dataRange,]
  ConfirmationNum <- as.numeric(substring(NewPeakTable[dataRange,IDcolumn], 1, 1)) #get's confirmation numbers (first element in string of IDs)
  IDSubstr <- substring(NewPeakTable[dataRange,IDcolumn], 3, 10)
  TempSorted <- cbind(TempSorted,ConfirmationNum)  #add IDs to end of NewPeakTable
  TempSorted <- TempSorted[order(ConfirmationNum, IDSubstr),]  #sort NewPeakTable based on Confirmation numbers
  TempSorted <- TempSorted[,-ncol(TempSorted)]
  NewPeakTable[dataRange,] <- TempSorted
  
  
  write.table(NewPeakTable, file.path(OutputDirectory, paste0(mode,"IDed.csv")), sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
  
  # Export MS/MS scans with assignments into a csv format
  # ID_Ranked_col <- ncol(IDedTable) - 4
  # for (i in RowStartForFeatureTableData:nrow(NewPeakTable)) {
  #   if (substr(NewPeakTable[i, ID_Ranked_col], 1, 2) != "5_") {
  #
  #   }
  # }
}
############function to append fragments to final NegIDed feature table, NegPos is "Neg" or "Pos"################
# Paul Stelben
# Jeremy Koelmel
# 02/13/2020

# NegPos <- "Neg"
# OutDir <- "Output/NegIDed_insilico.csv"
# OutDirOnlyFrags="Output/Neg_OnlyIDs_PredFrags.csv"
# InputDir_Append="Output/ddMS/NegByClass/Additional_Files"
# ID_name="PredictedFrag_IDs"
# FragName="Frags"
# nFrag="Num_Frags"
# fileNames="Files"
# ImportTable="Output/NegIDed_Fragments.csv"
# firstFragAppender <- FALSE

# For Testing:
#InputDirectory<-"C:/Users/jpk72/Desktop/OUT/LipidMatch_Run/"; CommentColumn<-1; RowStartForFeatureTableData<-2; NegPos = "Neg"; OutDir="/Output/NegIDed_FIN.csv"; OutDirOnlyFrags="/Output/Neg_OnlyIDs_PredFrags.csv"; InputDir_Append="Output/ddMS/NegByClass/Additional_Files"; ID_name="PredictedFrag_IDs"; FragName="Frags"; nFrag="Num_Frags"; fileNames="Files"; ImportTable="/Output/NegIDed_Fragments.csv"; firstFragAppender==FALSE
AppendFrag <- function (CommentColumn, RowStartForFeatureTableData, InputDirectory, NegPos, OutDir, OutDirOnlyFrags, InputDir_Append, ID_name, FragName, nFrag, fileNames, ImportTable, firstFragAppender){
  #Directory of tentative identifications
  if (NegPos=="Neg") {
    Dir_Additional_Files <- file.path(InputDirectory,InputDir_Append)
  } else if (NegPos=="Pos") {
    Dir_Additional_Files <- file.path(InputDirectory,InputDir_Append)
  } else {
    stop("error in appending fragments: NegPos input must be either 'Neg' or 'Pos'")
  }
  #List all files with fragment information
  All_Info_Files<-list.files(Dir_Additional_Files, pattern = "_All.csv")
  L <- length(All_Info_Files)
  #matrix to make with appended fragments for comments
  #empty row is included, if a matrix is one dimension it will be converted to a vector
  FragsFilled<-matrix(c("Tentative_IDs", "Fragments", "Comment", "file", "#Fragments"),1,5)
  # timestamp()
  
  
  # start <- 1
  # end <- L
  # All_Info_Files[1:L]
  
  for_loop <- function(start, end, All_Info_Files) {
    # FragsFilled1<-matrix(c("Tentative_IDs",0,"Fragments",0,"Comment",0,"file",0,"#Fragments",0),2,5)
    FragsFilled1<-matrix("", 0, 5)
    #Iteratively search across all files and compile species with fragments into one table
    for (i in 1:length(All_Info_Files)){
      #import file (i = 5 is a good test case)
      # print(i)
      All_Info_Temp<-read.csv(file.path(Dir_Additional_Files,All_Info_Files[i]), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE, check.names = FALSE)
      All_Info_Temp<-as.matrix(All_Info_Temp)
      #get instances of RT_Min in headers to define the fragment columns
      RT_Min_Index<-grep("RT_min", All_Info_Temp[1,], value = F, fixed = T)
      #Index which columns have fragments (with and without the precursor included)
      Frag_Index<-3:(RT_Min_Index[1]-1)
      Frag_Index_No_Adduct<-4:(RT_Min_Index[1]-1)
      ## Count 1's, if no 1's remove rows
      Frags<-which(All_Info_Temp[2:nrow(All_Info_Temp),Frag_Index_No_Adduct]=="1")
      #if there are fragments
      if (length(Frags)>0){
        # store current row
        RowFragsFilled1<-nrow(FragsFilled1)
        # reduce to those with Frags
        nrow_All_Info_Temp <- nrow(All_Info_Temp)
        for (x in 2:nrow_All_Info_Temp) {
          #calculate number of fragments
          Frag_Count<-length(which(All_Info_Temp[x,Frag_Index_No_Adduct]=="1"))
          if(Frag_Count>0) {
            #Combine all fragments with 1's into one line
            Frags_One_Line<-paste(All_Info_Temp[1,which(All_Info_Temp[x,Frag_Index]=="1")+2],collapse=";")
            #Vector of attributes (names, fragments, comment, File), will add file later
            VectorNewTable<-c(All_Info_Temp[x,2],Frags_One_Line,All_Info_Temp[x,RT_Min_Index[1]+2],NA,Frag_Count)
            FragsFilled1<-rbind(FragsFilled1,VectorNewTable)
          }
        }
        if (NegPos=="Neg") {
          FragsFilled1[(RowFragsFilled1+1):nrow(FragsFilled1),4]<-gsub("_Neg.*","",All_Info_Files[i])
        } else {
          FragsFilled1[(RowFragsFilled1+1):nrow(FragsFilled1),4]<-gsub("_Pos.*","",All_Info_Files[i])
        }
      }
      ## 1's will also be used to pull out and append observed fragments
    }
    # FragsFilled <- append(FragsFilled, FragsFilled1)
    return(FragsFilled1)
  }
  
  
  if (ParallelComputing == TRUE) {
    #a <- foreach (i = 1:L, .combine = rbind) %dopar% {
    #for_loop(i, i, All_Info_Files[i:i])
    #}
    a <- for_loop(1, L, All_Info_Files[1:L])
    FragsFilled <- rbind(FragsFilled, a)
  } else {
    a <- for_loop(1, L, All_Info_Files[1:L])
    FragsFilled <- rbind(FragsFilled, a)
  }
  
  #Header
  Header_Frags<-FragsFilled[1,]
  #Remove empty row
  if (nrow(FragsFilled) > 2) {
    #Header
    Header_Frags<-FragsFilled[1,]
    #Remove empty row
    FragsFilled<-FragsFilled[c(-1),]
    #Sort by the number of fragments
    FragsFilledSorted<-FragsFilled[order(as.numeric(FragsFilled[,5]),decreasing=TRUE),]
    #Sort by the comments
    FragsFilledSorted<-FragsFilledSorted[order(FragsFilledSorted[,3]),]
    #Add back headers
    FragsFilledSortedHead<-rbind(Header_Frags,FragsFilledSorted)
  } else if (nrow(FragsFilled) == 2) {
    FragsFilledSorted <- matrix(FragsFilled[2,], 1, ncol(FragsFilled))
    FragsFilledSortedHead <- FragsFilled
  } else {
    FragsFilledSorted <- matrix("", 1, ncol(FragsFilled))
    FragsFilledSortedHead <- matrix(FragsFilled[1,], 1, ncol(FragsFilled))
  }
  #Export total list of fragments and hits by FluoroMatch before appending and reducing
  write.table(FragsFilledSortedHead, file.path(InputDirectory,OutDirOnlyFrags), sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
  #aggregate values
  #FragsSortAggregate<-aggregate(FragsFilledSorted, list(FragsFilledSorted[,3]), function(x) paste0(unique(x)))
  FragsSortAggregate<-aggregate(FragsFilledSorted, by=list(FragsFilledSorted[,3]),function(x) paste(x,sep="|"))
  #Import NegIDed to append new columns
  NegIDed <- read.csv(file.path(InputDirectory,ImportTable), header=FALSE)
  
  #appending to table
  #Add 4 extra columns to be filled with 1) Tentative_IDs	2) Fragments	3) number of fragments 4) Files
  NewCols<-matrix(NA,nrow(NegIDed),4)
  NewCols[1,]<-c(ID_name,FragName,nFrag,fileNames)
  NegIDed <- cbind(NegIDed, NewCols)
  NegIDed <- as.matrix(NegIDed)
  #Break aggregates up and append to NegIDed
  nrow_NegIDed <- nrow(NegIDed)
  for (i in RowStartForFeatureTableData:nrow_NegIDed) {
    if (NegIDed[i,CommentColumn] %in% FragsSortAggregate[,1]) {
      #find the row in the FragSorAggregate table which has a matching identified to NegIDed
      Feature_Position<-match(NegIDed[i,CommentColumn],FragsSortAggregate[,1])
      #Collapse strings of each variable to append and add to NegIDed
      NegIDed[i,(ncol(NegIDed)-3)]<-paste0(FragsSortAggregate[Feature_Position,2][[1]],collapse="|")
      NegIDed[i,(ncol(NegIDed)-2)]<-paste0(FragsSortAggregate[Feature_Position,3][[1]],collapse="|")
      NegIDed[i,(ncol(NegIDed)-1)]<-paste0(FragsSortAggregate[Feature_Position,6][[1]],collapse="|")
      NegIDed[i,ncol(NegIDed)]<-paste0(FragsSortAggregate[Feature_Position,5][[1]],collapse="|")
    }
  }
  #Get the row for which no IDs exist (note it is blank the first round when no hits and NA the second round)
  if (firstFragAppender==TRUE) {
    NoIDs_row<-min(which(nchar(NegIDed[,ncol(NegIDed)-4])<1))
  } else {
    NoIDs_row<-min(which(is.na(NegIDed[,ncol(NegIDed)-4])))
    NegIDed[,ncol(NegIDed)-3]<-substr(NegIDed[,ncol(NegIDed)-3],1,10000)
  }
  # sort from that row down using the first element of number of fragments
  NegIDed[NoIDs_row:nrow(NegIDed),]<-NegIDed[order(NegIDed[NoIDs_row:nrow(NegIDed),ncol(NegIDed)-1],decreasing=TRUE)+NoIDs_row-1,]
  write.table(NegIDed, file.path(InputDirectory,OutDir), sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
}

# For debugging
# data <- Data
# row_start <- RowStartForFeatureTableData - 1
# adduct_col <- end_col - 5
# annotation_col <- end_col - 8
# window <- RT_Window
# RT_col <- RTColumn

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

############function to do KMD analysis, sort, and score for PFAS################
# Paul Stelben
# Jeremy Koelmel
# 09/02/2020

# For debugging
# IDedTable_dir <- file.path(InputDirectory,"Output/PosIDed_Fragments.csv")
# OutDir = file.path(InputDirectory,"Output")
# Mass_col = MZColumn
# Retention_col = RTColumn
# RowStartForFeatureTableData <- RowStartForFeatureTableData-1
# KMD_Bin_Window = (PrecursorMassAccuracy*2)
# NegPos = "Pos"

#Inputs:
# 1) directory after EPA masterlist predicted annotation, 2) output directory
# 3) m/z, 4) RT, 5) the confident annotation column, 6) one or more frag column, 7) same for predicted from EPA master, 8) whether or not RT order should be accounted for
Scoring <- function (IDedTable_dir, OutDir, Mass_col, Retention_col, RowStartForFeatureTableData, RT_flagging, KMD_Bin_Window, upper, lower, RepeatingUnits_dir){
  
  #Import IDedTable
  IDedTable <- read.csv(IDedTable_dir,sep=",")
  IDedTable <- as.matrix(IDedTable)
  
  if (nrow(IDedTable) == 0) {
    message("No MSMS data, stopping scoring function")
    stop()
  }
  
  #Import RepeatingUnits
  RepeatingUnits <- read.csv(RepeatingUnits_dir,sep=",")
  RepeatingUnits <- as.matrix(RepeatingUnits)
  Names <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),1])
  Series <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),2])
  
  # FOR DEBUGGING
  # IDedTable <- IDedTable[1:3000,]
  
  ID_Ranked_col <- which(colnames(IDedTable) == "ID_Ranked")
  Potential_IDs_col <- which(colnames(IDedTable) == "Potential_IDs")
  Frags_col <- ncol(IDedTable)-2
  
  # add columns and calculate mass defect and Kendrick mass defect for Mass_col
  m <- matrix(0, nrow(IDedTable), 2)
  IDedTable <- cbind(IDedTable, m)
  nrow_IDedTable <- nrow(IDedTable)
  ncol_IDedTable <- ncol(IDedTable)
  colnames(IDedTable)[ncol_IDedTable - 1] <- "nominal mass"
  colnames(IDedTable)[ncol_IDedTable] <- "mass defect"
  Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
  nmass <- round(Mass_vec)
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 1] <- as.numeric(nmass)
  massd <- Mass_vec - nmass
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable] <- as.numeric(massd)
  md_col <- ncol_IDedTable
  hseries_count_cols <- c()
  
  # timestamp()
  
  length_Series <- length(Series)
  for(x in 1:length_Series){
    # print(x)
    # x <- 1
    Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
    exactMass <- getMolecule(Series[x], maxisotopes = 1)$exactmass
    Mass <- round(exactMass)
    m <- matrix(0, nrow(IDedTable), 5)
    IDedTable <- cbind(IDedTable, m)
    nrow_IDedTable <- nrow(IDedTable)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable - 4] <- paste(Names[x], "Kendrick_mass", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 3] <- paste(Names[x], "nominal_Kendrick_mass", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 2] <- paste(Names[x], "Kendrick_mass_defect", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 1] <- paste(Names[x], "KMD_Group", sep="_")
    colnames(IDedTable)[ncol_IDedTable] <- paste(Names[x], "homologous_series", sep="_")
    
    KM <- Mass_vec * (Mass/exactMass)
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 4] <- as.numeric(KM)
    nKM <- round(KM)
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 3] <- as.numeric(nKM)
    KMD <- KM - nKM
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 2] <- as.numeric(KMD)
    
    # Sort and group by KMD
    tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
    tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 2])),]
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
    KMD_col <- ncol_IDedTable - 2
    group <- 1
    start <- RowStartForFeatureTableData
    IDedTable[RowStartForFeatureTableData, KMD_col + 1] <- group
    for (i in (RowStartForFeatureTableData + 1):nrow_IDedTable) {
      if (as.numeric(IDedTable[i, KMD_col]) - as.numeric(IDedTable[start, KMD_col]) > KMD_Bin_Window) {
        IDedTable[start:(i-1), KMD_col + 1] <- group
        group <- group + 1
        start <- i
      }
    }
    if (start != nrow_IDedTable) {
      IDedTable[start:nrow_IDedTable, KMD_col + 1] <- group
    }
    
    # Calculate homologous series
    series <- 1
    for (i in RowStartForFeatureTableData:(nrow_IDedTable - 1)) {
      if (is.na(IDedTable[i, KMD_col + 2]) || IDedTable[i, KMD_col + 2] != 0) {
        next
      }
      group <- as.numeric(IDedTable[i, KMD_col + 1])
      grow <- which(IDedTable[, KMD_col + 1] == group)
      mrow <- grow[which(((as.numeric(IDedTable[grow, KMD_col - 1]) - as.numeric(IDedTable[i, KMD_col - 1])) / Mass) %% 1 == 0)]
      members <- as.numeric(IDedTable[mrow, KMD_col - 1])
      if (length(unique(members)) < 2) {
        IDedTable[mrow, KMD_col + 2] <- NA
      } else {
        IDedTable[mrow, KMD_col + 2] <- series
        series <- series + 1
      }
    }
    IDedTable[which(is.na(IDedTable[, KMD_col + 2])), KMD_col + 2] <- 0
    
    # Sort by homologous series
    # IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable])),]
    
    m <- matrix(0, nrow(IDedTable), 2)
    IDedTable <- cbind(IDedTable, m)
    # Phantom columns
    ncol_IDedTable <- ncol(IDedTable) + 8
    colnames(IDedTable)[ncol_IDedTable - 9] <- paste(Names[x], "Homologous_series_count", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 8] <- paste(Names[x], "MD_Filter", sep="_")
    hseries_count_cols <- append(hseries_count_cols, ncol_IDedTable - 9)
    
    # Mass defect filtering
    IDedTable[, ncol_IDedTable - 8] <- (as.numeric(IDedTable[, md_col]) >= lower & as.numeric(IDedTable[, md_col]) <= upper)
    
    # Pre-flagging sorting
    tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
    tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 10]), as.numeric(tmptable[, ncol_IDedTable - 11]), as.numeric(tmptable[, ncol_IDedTable - 13]), as.numeric(tmptable[, Retention_col])),]
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice0.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    # RT Flagging
    # which() function?
    m <- matrix(0, nrow(IDedTable), 2)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol_IDedTable - 1
    colnames(IDedTable)[ncol_IDedTable - 6] <- paste(Names[x], "RT_series", sep="_")
    series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 9]))
    if (series[1] == 0) {
      series <- series[-1]
    }
    for (i in series) {
      curr <- which(IDedTable[, ncol_IDedTable - 9] == i)
      mseries <- as.numeric(unique(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]))
      if (length(mseries) < 2) {
        next
      }
      length_mseries <- length(mseries)
      for (j in 1:length_mseries) {
        mcurr <- curr[which(as.numeric(IDedTable[curr, ncol_IDedTable - 12]) == mseries[j])]
        if (mseries[j] == mseries[1]) {
          nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
          IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
        } else if (mseries[j] == mseries[length(mseries)]) {
          prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
          IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
        } else {
          prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
          IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
          nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
          IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
        }
      }
    }
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 6] == FALSE || IDedTable[i, ncol_IDedTable - 5] == FALSE) {
        IDedTable[i, ncol_IDedTable - 6] <- "RT not ordered"
      } else if (IDedTable[i, ncol_IDedTable - 6] == TRUE || IDedTable[i, ncol_IDedTable - 5] == TRUE) {
        IDedTable[i, ncol_IDedTable - 6] <- "RT ordered"
      } else {
        IDedTable[i, ncol_IDedTable - 6] <- "NA"
      }
    }
    
    IDedTable <- IDedTable[, 1:(ncol_IDedTable - 6)]
    
    ncol_IDedTable <- ncol_IDedTable + 1
    
    ## For debugging
    # write.table(IDedTable, file.path(OutDir, "NegIDed_test_optimized.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    if (RT_flagging == TRUE) {
      # Homologous seris count
      series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
      if (series[1] == 0) {
        series <- series[-1]
      }
      none <- 2
      for (i in series) {
        curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
        stable <- IDedTable[curr, c(ncol_IDedTable - 13, ncol_IDedTable - 7)]
        # Remove not ordered rows
        stable <- stable[!(is.na(stable[, 2]) == FALSE & stable[, 2] != "RT ordered" & stable[, 2] != "RT ordered, RT ordered"),]
        if (is.null(nrow(stable))) {
          none <- 1
        } else if (nrow(stable) == 0 && ncol(stable) == 2) {
          none <- 0
        }
        if (none == 2) {
          svec <- unique(stable[, 1])
          IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
        } else if (none == 1) {
          IDedTable[curr, ncol_IDedTable - 9] <- 1
          none <- 2
        } else if (none == 0) {
          IDedTable[curr, ncol_IDedTable - 9] <- 0
          none <- 2
        }
      }
    } else {
      series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
      if (series[1] == 0) {
        series <- series[-1]
      }
      for (i in series) {
        curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
        svec <- unique(IDedTable[curr, ncol_IDedTable - 13])
        IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
      }
    }
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice1.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  }
  
  # timestamp()
  
  # write.table(IDedTable, file.path(OutDir, "Test.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  
  if (Lipid == FALSE) {
    m <- matrix(0, nrow(IDedTable), 8)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable - 7] <- "Confident_ID"
    colnames(IDedTable)[ncol_IDedTable - 6] <- "2+_homologous_series"
    colnames(IDedTable)[ncol_IDedTable - 5] <- "3+_homologous_series"
    colnames(IDedTable)[ncol_IDedTable - 4] <- "Tentative_ID"
    colnames(IDedTable)[ncol_IDedTable - 3] <- "F_containing"
    colnames(IDedTable)[ncol_IDedTable - 2] <- "Exact_mass_match"
    colnames(IDedTable)[ncol_IDedTable - 1] <- "Score"
    colnames(IDedTable)[ncol_IDedTable] <- "Score_Description"
    
    # Scoring system
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "1_" || substr(IDedTable[i, ID_Ranked_col], 1, 2) == "2_") {
        IDedTable[i, ncol_IDedTable - 7] <- TRUE
      } else {
        IDedTable[i, ncol_IDedTable - 7] <- FALSE
      }
      IDedTable[i, ncol_IDedTable - 6] <- FALSE
      for (hsc_col in hseries_count_cols){
        if (as.numeric(IDedTable[i, hsc_col]) >= 2) {
          IDedTable[i, ncol_IDedTable - 6] <- TRUE
        }
      }
      IDedTable[i, ncol_IDedTable - 5] <- FALSE
      for (hsc_col in hseries_count_cols){
        if (as.numeric(IDedTable[i, hsc_col]) >= 3) {
          IDedTable[i, ncol_IDedTable - 5] <- TRUE
        }
      }
      if (is.na(IDedTable[i, Potential_IDs_col])) {
        IDedTable[i, ncol_IDedTable - 4] <- FALSE
      } else {
        IDedTable[i, ncol_IDedTable - 4] <- TRUE
      }
      Fvec <- strsplit(IDedTable[i, Frags_col], "")
      if (is.na(Fvec) == FALSE) {
        for (j in 1:length(Fvec[[1]])) {
          if (Fvec[[1]][j] == 'F') {
            IDedTable[i, ncol_IDedTable - 3] <- TRUE
            break
          }
        }
      }
      if (IDedTable[i, ncol_IDedTable - 3] != TRUE) {
        IDedTable[i, ncol_IDedTable - 3] <- FALSE
      }
      if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        IDedTable[i, ncol_IDedTable - 2] <- TRUE
      } else {
        IDedTable[i, ncol_IDedTable - 2] <- FALSE
      }
    }
    
    # timestamp()
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 7] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "A"
        IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments observed along with exact mass) and 2+ in homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 7] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "A-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments and exact mass)"
      } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B+"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable[i, ncol_IDedTable - 3] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" || IDedTable[i, ncol_IDedTable - 3] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if ((IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") && IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D+"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: 3+ within homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: Mass defect falling within -0.11 and 0.12 OR exact mass match"
      } else {
        IDedTable[i, ncol_IDedTable - 1] <- "E"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: likely not PFAS"
      }
    }
    
    # timestamp()
    
    # for (i in RowStartForFeatureTableData:nrow_IDedTable) {
    if (ParallelComputing == TRUE) {
      
      IDedTable <- foreach (i = RowStartForFeatureTableData:nrow_IDedTable, .combine = rbind) %dopar% {
        #print(i)
        if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
          for (hsc_col in hseries_count_cols) {
            if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
              tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
              hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
              start <- hsc_rows[1]
              stop <- hsc_rows[length(hsc_rows)]
              if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
                IDedTable[i, ncol_IDedTable - 1] <- "C+"
                IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
              } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
                if (IDedTable[i, hsc_col] >= 3) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                } else if (IDedTable[i, hsc_col] >= 2) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C-"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                }
              }
            }
          }
        }
        return(IDedTable[i,])
      }
    } else {
      for (i in RowStartForFeatureTableData:nrow_IDedTable) {
        #print(i)
        if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
          for (hsc_col in hseries_count_cols) {
            if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
              tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
              hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
              start <- hsc_rows[1]
              stop <- hsc_rows[length(hsc_rows)]
              if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
                IDedTable[i, ncol_IDedTable - 1] <- "C+"
                IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
              } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
                if (IDedTable[i, hsc_col] >= 3) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                } else if (IDedTable[i, hsc_col] >= 2) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C-"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                }
              }
            }
          }
        }
      }
    }
    
    # timestamp()
    
    # Change Ds to D+s when in series with D+s
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      #print(i)
      if (IDedTable[i, ncol_IDedTable - 1] == "D") {
        D_rows <- which(IDedTable[, ncol_IDedTable - 11] == IDedTable[i, ncol_IDedTable - 11])
        start <- D_rows[1]
        stop <- D_rows[length(D_rows)]
        for (j in start:stop) {
          if (IDedTable[j, ncol_IDedTable - 1] == "D+") {
            IDedTable[i, ncol_IDedTable - 1] <- "D+"
            IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
            break
          }
        }
      }
    }
    
    # timestamp()
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    # Next seven lines for debugging
    # rm(list = ls())
    # IDedTable_dir <- "C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/IDedTableIDed_FIN_KMD_scored_final.csv"
    # IDedTable <- read.csv(IDedTable_dir,sep=",")
    # IDedTable <- as.matrix(IDedTable)
    # Retention_col <- 3
    # nrow_IDedTable <- nrow(IDedTable)
    # ncol_IDedTable <- ncol(IDedTable)
    
    # Remove tentative and F containing columns
    IDedTable <- IDedTable[, -((ncol_IDedTable - 4):(ncol_IDedTable - 3))]
    ncol_IDedTable <- ncol(IDedTable)
    
    ######### Sorting
    
    m <- matrix(0, nrow(IDedTable), 4)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol(IDedTable) - 4
    
    # Find max series for each row
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      IDedTable[i, ncol_IDedTable + 2] <- max(which(IDedTable[i, hseries_count_cols] == max(IDedTable[i, hseries_count_cols])))
    }
    
    # Max series series number and count
    m <- matrix(0, nrow_IDedTable, 1)
    IDedTable <- cbind(IDedTable, m)
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      IDedTable[i, ncol_IDedTable + 3] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])] - 1])
      IDedTable[i, ncol_IDedTable + 4] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])]])
    }
    
    # Sort by max series count, max series, and series number
    IDedTable <- IDedTable[order(-as.numeric(IDedTable[, ncol_IDedTable + 4]), -as.numeric(IDedTable[, ncol_IDedTable + 2]), as.numeric(IDedTable[, ncol_IDedTable + 3])),]
    
    # Create rank
    rank <- 1
    IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- rank
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
      if (as.numeric(IDedTable[i, ncol_IDedTable + 4]) == 0) {
        rank <- rank + 1
      } else if (IDedTable[i, ncol_IDedTable + 3] != IDedTable[i - 1, ncol_IDedTable + 3]) {
        rank <- rank + 1
      }
      IDedTable[i, ncol_IDedTable + 5] <- rank
    }
    
    # Score
    for (i in RowStartForFeatureTableData:nrow(IDedTable)) {
      if (IDedTable[i, ncol_IDedTable - 1] == "A" || IDedTable[i, ncol_IDedTable - 1] == "A-") {
        IDedTable[i, ncol_IDedTable + 1] <- 1
      } else if (IDedTable[i, ncol_IDedTable - 1] == "B+" || IDedTable[i, ncol_IDedTable - 1] == "B" || IDedTable[i, ncol_IDedTable - 1] == "B-") {
        IDedTable[i, ncol_IDedTable + 1] <- 2
      } else if (IDedTable[i, ncol_IDedTable - 1] == "E") {
        IDedTable[i, ncol_IDedTable + 1] <- 4
      } else {
        IDedTable[i, ncol_IDedTable + 1] <- 3
      }
    }
    
    # Sort by rank and score
    IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, ncol_IDedTable + 1])),]
    
    # Adjust rank
    old <- IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]
    if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 1) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 10000000
    } else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 2) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 1000000
    } else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 4) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) + 1000000
    }
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
      if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 4) {
        old <- IDedTable[i, ncol_IDedTable + 5]
        IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) + 1000000
      } else if ((IDedTable[i, ncol_IDedTable + 5] != old)) {
        old <- IDedTable[i, ncol_IDedTable + 5]
        if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 1) {
          IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 10000000
        } else if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 2) {
          IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 1000000
        }
      } else {
        old <- IDedTable[i, ncol_IDedTable + 5]
        IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i - 1, ncol_IDedTable + 5])
      }
    }
    
    # Sort by rank and m/z
    IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, Mass_col])),]
    
    ## Is this right? I think so.
    # Final aesthetics
    ncol_IDedTable <- ncol_IDedTable - 1
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 6] == "TRUE" || IDedTable[i, ncol_IDedTable - 6] == " TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "TRUE"
      }
    }
  }
  
  # New formatting
  m <- matrix("", nrow_IDedTable, 5)
  n <- matrix("", nrow_IDedTable, 4)
  IDedTable <- cbind(m, IDedTable[, Mass_col], IDedTable[, Retention_col], n, IDedTable[, -c(Mass_col, Retention_col)])
  colnames(IDedTable)[1] <- "Score"
  colnames(IDedTable)[2] <- "SeriesType_Identifier"
  colnames(IDedTable)[3] <- "Name_or_Class"
  colnames(IDedTable)[4] <- "Formula"
  colnames(IDedTable)[5] <- "SMILES"
  colnames(IDedTable)[6] <- "m/z"
  colnames(IDedTable)[7] <- "Retention Time"
  colnames(IDedTable)[8] <- "Adduct"
  colnames(IDedTable)[9] <- "Unique"
  colnames(IDedTable)[10] <- "Score_Description"
  colnames(IDedTable)[11] <- "Needs_Validation"
  
  ncol_IDedTable <- ncol(IDedTable)
  ID_Ranked_col <- 9 + ID_Ranked_col
  Potential_IDs_col <- 9 + Potential_IDs_col
  Frags_col <- 9 + Frags_col
  
  
  if (Lipid == FALSE && TWeen_pos == FALSE) {
    IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
    IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
    IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- paste(name[3:(which(name == "-")[2] - 1)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
          IDedTable[i, 9] <- IDedTable[i, ID_Ranked_col + 4]
        }
      } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
        primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- paste(name[1:(which(name == "-")[2] - 1)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
        }
        u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      } else if (length(grep("DTXSID", IDedTable[i, Frags_col - 1])) > 0) {
        all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
        primary <- strsplit(all[grep("DTXSID", all)][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          IDedTable[i, 3] <- primary[4]
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[1]
          IDedTable[i, 8] <- "[M-H]-"
        }
        all <- strsplit(all, ";", TRUE)
        u_entries <- c()
        for (j in 1:length(all)) {
          u_entries <- append(u_entries, all[[j]][1])
        }
        u_entries <- unique(u_entries)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(grep(";\\[", IDedTable[i, ID_Ranked_col])) == 0) {
          end <- strsplit(primary[length(primary)], "\\[")[[1]]
          primary[length(primary)] <- end[1]
          primary <- append(primary, paste("[", end[2], sep=""))
        }
        IDedTable[i, 3] <- "NA"
        IDedTable[i, 4] <- primary[length(primary) - 2]
        IDedTable[i, 5] <- "NA"
        IDedTable[i, 8] <- primary[length(primary)]
        all <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], "|", TRUE)[[1]], ";", TRUE)
        u_entries <- c()
        length_all <- length(all)
        for (j in 1:length_all) {
          if (length(grep("\\[", all[[j]][length(all[[j]])])) > 0) {
            end <- strsplit(all[[j]][length(all[[j]])], "\\[")[[1]]
            all[[j]][length(all[[j]])] <- end[1]
            all[[j]] <- append(all[[j]], paste("[", end[2], sep=""))
          }
          u_entries <- append(u_entries, all[[j]][length(all[[j]]) - 1])
        }
        u_entries <- unique(u_entries)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      }
      if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
        IDedTable[i, 11] <- "Yes"
      } else {
        IDedTable[i, 11] <- "No"
      }
    }
  } else if (TWeen_pos == TRUE) {
    IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
    IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
    IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      # Slight name modification
      if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
          IDedTable[i, 9] <- IDedTable[i, ID_Ranked_col + 4]
        }
        # Slight name modification
      } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
        primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- primary[1]
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
        }
        u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
        ### Name modification
      } else if (length(grep("PEG", IDedTable[i, Frags_col - 1])) > 0) {
        all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
        # Take first PEG for name
        # primary <- strsplit(all[grep("PEG", all)][1], ";", TRUE)[[1]]
        primary <- all[grep("PEG", all)]
        if (length(primary) > 1) {
          IDedTable[i, 3] <- primary[1]
          IDedTable[i, 4] <- "NA"
          IDedTable[i, 5] <- "NA"
          # [M+?]+ fill in ? from PEG name
          IDedTable[i, 8] <- paste("[M+", strsplit(strsplit(primary[1], "_", TRUE)[[1]][2], "+", TRUE)[[1]][1], "]+", sep = "")
        }
        # Don't split, just check for unique
        u_entries <- unique(primary)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
        ### Name modification
      } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        primary <- strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]]
        name <- strsplit(primary[1], "")[[1]]
        IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
        IDedTable[i, 4] <- "NA"
        IDedTable[i, 5] <- "NA"
        IDedTable[i, 8] <- "[M+?]+"
        
        u_entries <- unique(primary)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      }
      if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
        IDedTable[i, 11] <- "Yes"
      } else {
        IDedTable[i, 11] <- "No"
      }
    }
  } else if (Lipid == TRUE) {
    ################################## LIPIDMATCH SCORING SYSTEM
    IDedTable[, 1] <- substr(IDedTable[, ID_Ranked_col], 1, 1)
    IDedTable[, 2] <- paste(Names[1], IDedTable[, ncol_IDedTable - 3], sep="_")
    IDedTable[, 3] <- IDedTable[, ID_Ranked_col + 2]
    IDedTable[, 4] <- "CH2"
    IDedTable[, 5] <- "CC"
    
    IDedTable[, 8] <- IDedTable[, ID_Ranked_col + 3]
    IDedTable[which(is.na(IDedTable[, 8])), 8] <- 0
    
    IDedTable[, 9] <- IDedTable[, ID_Ranked_col + 4]
    IDedTable[which(is.na(IDedTable[, 9])), 9] <- "No"
    
    IDedTable[which(IDedTable[, 1] == 1), 10] <- "Annotated using class based rules from standards including by fatty acyl chain constituents with ~5% false positive rate at the level of class and C:DB"
    IDedTable[which(IDedTable[, 1] == 2), 10] <- "Annotated using class based rules from standards and DIA data with ~5% false positive rate at the level of class and C:DB"
    IDedTable[which(IDedTable[, 1] == 3), 10] <- "Annotated using class based rules from standards and DDA data but only annotated by class (fatty acid composition not known)"
    IDedTable[which(IDedTable[, 1] == 4), 10] <- "accurate mass match to database with very high false positive rate (>> 80%) so only use as a starting point for compound identification"
    IDedTable[which(IDedTable[, 1] == 5), 10] <- "no match"
    
    IDedTable[, 11] <- "Yes"
  }
  
  # timestamp()
  
  IDedTable <- IDedTable[, 1:(ncol_IDedTable - 5)]
  
  # Remove whitespace
  IDedTable[IDedTable == ""] <- NA
  
  # Remove "Na" from formulas in column 4
  IDedTable[, 4] <- gsub("Na", "", IDedTable[, 4])
  
  ####################################################################################
  
  if (Lipid == FALSE) {
    # Create a second output file with more info on B scores
    IDedTable_B <- IDedTable
    m <- matrix(0, nrow(IDedTable_B), 5)
    IDedTable_B <- cbind(IDedTable_B, m)
    nrow_IDedTable_B <- nrow(IDedTable_B)
    ncol_IDedTable_B <- ncol(IDedTable_B)
    colnames(IDedTable_B)[ncol_IDedTable_B - 4] <- "Check_Common_Fragment"
    colnames(IDedTable_B)[ncol_IDedTable_B - 3] <- "Number_of_F_Fragments"
    colnames(IDedTable_B)[ncol_IDedTable_B - 2] <- "Formula"
    colnames(IDedTable_B)[ncol_IDedTable_B - 1] <- "MD"
    colnames(IDedTable_B)[ncol_IDedTable_B] <- "Likely_PFAS"
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      count <- 0
      frags_to_check <- c("C2F5", "C3F7", "SO2F", "SO3F", "CO2F", "C2F5O", "CF3O", "C2F3O2", "C3F7O", "C3O2F7", "CF3", "CH3FNSO2", "PO2F", "PO2F2", "SF5")
      for (j in frags_to_check) {
        if (length(grep(j, IDedTable_B[i, Frags_col])) > 0) {
          count <- count + 1
        }
      }
      IDedTable_B[i, ncol_IDedTable_B - 4] <- count
      if (!is.na(IDedTable_B[i, Frags_col + 1])) {
        IDedTable_B[i, ncol_IDedTable_B - 3] <- strsplit(IDedTable_B[i, Frags_col + 1], split="|", fixed=TRUE)[[1]][1]
      }
      if (is.na(IDedTable_B[i, 4])) {
        IDedTable_B[i, ncol_IDedTable_B - 2] <- FALSE
      } else {
        IDedTable_B[i, ncol_IDedTable_B - 2] <- TRUE
      }
      if (as.numeric(IDedTable_B[i, Frags_col + 4]) > -0.25 && as.numeric(IDedTable_B[i, Frags_col + 4]) < 0.1) {
        IDedTable_B[i, ncol_IDedTable_B - 1] <- TRUE
      } else {
        IDedTable_B[i, ncol_IDedTable_B - 1] <- FALSE
      }
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B - 4]) > 0 || as.numeric(IDedTable_B[i, ncol_IDedTable_B - 3]) > 3) {
        IDedTable_B[i, ncol_IDedTable_B] <- TRUE
      } else {
        IDedTable_B[i, ncol_IDedTable_B] <- FALSE
      }
      
      # Assign B--
      if (IDedTable_B[i, ncol_IDedTable_B] == FALSE && (IDedTable_B[i, 1] == "B+" || IDedTable_B[i, 1] == "B" || IDedTable_B[i, 1] == "B-")) {
        IDedTable_B[i, 1] <- "B--"
        IDedTable_B[i, ncol_IDedTable_B - 6] <- "B--"
        IDedTable_B[i, 10] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      }
      
      # Correct B score descriptions that have SMILES
      if (IDedTable_B[i, 1] == "B+" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B-" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B--" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      }
      
    }
    
    ###################### Final Sorting ###########################
    
    m <- matrix(0, nrow(IDedTable_B), 4)
    IDedTable_B <- cbind(IDedTable_B, m)
    ncol_IDedTable_B <- ncol(IDedTable_B) - 4
    hseries_count_cols <- hseries_count_cols + 9
    Mass_col <- Mass_col + 4
    
    # Find max series for each row
    for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
      IDedTable_B[i, ncol_IDedTable_B + 2] <- max(which(IDedTable_B[i, hseries_count_cols] == max(IDedTable_B[i, hseries_count_cols])))
    }
    
    # Max series series number and count
    m <- matrix(0, nrow_IDedTable_B, 1)
    IDedTable_B <- cbind(IDedTable_B, m)
    for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
      IDedTable_B[i, ncol_IDedTable_B + 3] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])] - 1])
      IDedTable_B[i, ncol_IDedTable_B + 4] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])]])
    }
    
    # Sort by max series count, max series, and series number
    IDedTable_B <- IDedTable_B[order(-as.numeric(IDedTable_B[, ncol_IDedTable_B + 4]), -as.numeric(IDedTable_B[, ncol_IDedTable_B + 2]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 3])),]
    
    # Create rank
    rank <- 1
    IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- rank
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 4]) == 0) {
        rank <- rank + 1
      } else if (IDedTable_B[i, ncol_IDedTable_B + 3] != IDedTable_B[i - 1, ncol_IDedTable_B + 3]) {
        rank <- rank + 1
      }
      IDedTable_B[i, ncol_IDedTable_B + 5] <- rank
    }
    
    # Score
    for (i in RowStartForFeatureTableData:nrow(IDedTable_B)) {
      if (IDedTable_B[i, ncol_IDedTable_B - 6] == "A" || IDedTable_B[i, ncol_IDedTable_B - 6] == "A-") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 1
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B+" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B-") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 2
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B--") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 3
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "E") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 5
      } else {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 4
      }
    }
    
    # Sort by rank and score
    IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 1])),]
    
    # Adjust rank
    old <- IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]
    if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 1) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 10000000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 2) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 1000000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 3) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 100000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 5) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) + 1000000
    }
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 5) {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) + 1000000
      } else if ((IDedTable_B[i, ncol_IDedTable_B + 5] != old)) {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 1) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 10000000
        } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 2) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 1000000
        } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 3) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 100000
        }
      } else {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i - 1, ncol_IDedTable_B + 5])
      }
    }
    
    # Sort by rank and m/z
    IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, Mass_col])),]
    
    IDedTable_B <- IDedTable_B[, 1:ncol_IDedTable_B]
  } else if (Lipid == TRUE) {
    m <- matrix(0, nrow(IDedTable), 1)
    IDedTable <- cbind(IDedTable, m)
    nrow_IDedTable <- nrow(IDedTable)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable] <- "Number_of_F_Fragments"
    Num_Frags_col <- which(colnames(IDedTable) == "Num_Frags")
    if (length(Num_Frags_col) > 0) {
      IDedTable[, ncol_IDedTable] <- substr(IDedTable[, Num_Frags_col], 1, 1)
    }
  }
  
  # # Add convenient columns
  # m <- matrix("", nrow(IDedTable_B), 5)
  # IDedTable_B <- cbind(IDedTable_B, m)
  # ncol_IDedTable_B <- ncol(IDedTable_B)
  # colnames(IDedTable_B)[(ncol_IDedTable_B - 4):ncol_IDedTable_B] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
  # 
  # 
  # if (NegPos == "Neg") {
  #   # write.table(IDedTable, file.path(OutDir, "NegIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  #   write.table(IDedTable_B, file.path(OutDir, "NegIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # } else if (NegPos == "Pos") {
  #   # write.table(IDedTable, file.path(OutDir, "PosIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  #   write.table(IDedTable_B, file.path(OutDir, "PosIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # } else if (NegPos == "Combined") {
  #   write.table(IDedTable_B, file.path(OutDir, "CombinedIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # }
  
  # Add convenient columns
  m <- matrix("", nrow(IDedTable), 5)
  IDedTable <- cbind(IDedTable, m)
  ncol_IDedTable <- ncol(IDedTable)
  colnames(IDedTable)[(ncol_IDedTable - 4):ncol_IDedTable] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
  
  
  OutDir_parts <- strsplit(IDedTable_dir, split = "_")[[1]]
  OutDir_parts <- paste(OutDir_parts[-length(OutDir_parts)], collapse = "_")
  OutDir_full <- paste0(OutDir_parts, "_FIN.csv")
  write.table(IDedTable, OutDir_full, sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA", qmethod="double")
  
}

####################################End functions##########################################

####Read in files, create folder structure, and error handle####
if(length(foldersToRun)==0){
  lengthFoldersToRun <- 1 #if there are no subfolders, that means you have the feature table and ms2s in that current directory, therefore, run analysis on those files.
}else{
  lengthFoldersToRun <- length(foldersToRun)#run analysis on all subfolders
}
RepeatingUnits_dir <- file.path(InputLibrary, "REPEATING_UNITS_INPUT.csv")



for(i in seq_len(lengthFoldersToRun)){
  if(length(foldersToRun)==0){#we're in current (and only) folder that contains feature table and ms2
    fpath <- InputDirectory
  }else if(foldersToRun[i] == "Output"){
    fpath <- InputDirectory
    print(paste("Warning: Remove your 'Output' folder from the current Input Directory:", InputDirectory))
  }else{
    fpath <- file.path(InputDirectory, foldersToRun[i])
  }
  fileName <- basename(fpath)
  
  ddMS2NEG_in <- list.files(path=fpath, pattern="[nNgG]\\.ms2", ignore.case=FALSE)
  AIFMS1NEG_in <- list.files(path=fpath, pattern="[nNgG]\\.ms1", ignore.case=FALSE)
  AIFMS2NEG_in <- list.files(path=fpath, pattern="[nNgG]\\.ms2", ignore.case=FALSE)
  ddMS2POS_in <- list.files(path=fpath, pattern="[pPsS]\\.ms2", ignore.case=FALSE)
  AIFMS1POS_in <- list.files(path=fpath, pattern="[pPsS]\\.ms1", ignore.case=FALSE)
  AIFMS2POS_in <- list.files(path=fpath, pattern="[pPsS]\\.ms2", ignore.case=FALSE)
  
  #separate ddMS and AIF
  ddMS2NEG_in <- ddMS2NEG_in[grep("[dD][dD]", ddMS2NEG_in)]
  AIFMS1NEG_in <- AIFMS1NEG_in[grep("[Aa][Ii][Ff]", AIFMS1NEG_in)] #Yang 20180315. Orig: AIFMS1NEG_in <- AIFMS1NEG_in[grep("[AIFaif][AIFaif][AIFaif]", AIFMS1NEG_in)]
  AIFMS2NEG_in <- AIFMS2NEG_in[grep("[Aa][Ii][Ff]", AIFMS2NEG_in)] #Yang 20180315. Orig: AIFMS2NEG_in <- AIFMS2NEG_in[grep("[AIFaif][AIFaif][AIFaif]", AIFMS2NEG_in)]
  ddMS2POS_in <- ddMS2POS_in[grep("[dD][dD]", ddMS2POS_in)]
  AIFMS1POS_in <- AIFMS1POS_in[grep("[Aa][Ii][Ff]", AIFMS1POS_in)] #Yang 20180315. Orig: AIFMS1POS_in <- AIFMS1POS_in[grep("[AIFaif][AIFaif][AIFaif]", AIFMS1POS_in)]
  AIFMS2POS_in <- AIFMS2POS_in[grep("[Aa][Ii][Ff]", AIFMS2POS_in)] #Yang 20180315. Orig: AIFMS2POS_in <- AIFMS2POS_in[grep("[AIFaif][AIFaif][AIFaif]", AIFMS2POS_in)]
  
  #user info outputted for error handling
  if(length(ddMS2POS_in) == 0){
    print(paste("CAUTION: We detected", length(ddMS2POS_in),"positive ddMS .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have 'p', 'P', 'pos', or 'POS' at the end of the file name and you must have a 'dd' within the name. OR Remove the folder: ", fileName))
  }
  if(length(ddMS2NEG_in) == 0){
    print(paste("CAUTION: We detected", length(ddMS2NEG_in),"negative ddMS .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have 'n', 'N', 'neg', or 'NEG' at the end of the file name and you must have a 'dd' within the name. OR Remove the folder: ", fileName))
  }
  if(length(AIFMS1POS_in) == 0){
    print(paste("CAUTION: We detected", length(AIFMS1POS_in),"positive AIF .ms1 files in the folder: ", fileName," ...If this incorrect, check that you have ('AIF') and ('p', 'P', 'pos', or 'POS') at the end of the file name. OR Remove the folder: ", fileName))
  }
  if(length(AIFMS1NEG_in) == 0){
    print(paste("CAUTION: We detected", length(AIFMS1NEG_in),"negative AIF .ms1 files in the folder: ", fileName," ...If this incorrect, check that you have ('AIF') and ('n', 'N', 'neg', or 'NEG') at the end of the file name. OR Remove the folder: ", fileName))
  }
  if(length(AIFMS2POS_in) == 0){
    print(paste("CAUTION: We detected", length(AIFMS2POS_in),"positive AIF .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have ('AIF') and ('p', 'P', 'pos', or 'POS') at the end of the file name. OR Remove the folder: ", fileName))
  }
  if(length(AIFMS2NEG_in) == 0){
    print(paste("CAUTION: We detected", length(AIFMS2NEG_in),"negative AIF .ms2 files in the folder: ", fileName," ...If this incorrect, check that you have ('AIF') and ('n', 'N', 'neg', or 'NEG') at the end of the file name. OR Remove the folder: ", fileName))
  }
  
  FeatureTable_NEG <- list.files(path=fpath, pattern="[nNgG]\\.csv", ignore.case=FALSE)
  if(length(FeatureTable_NEG) > 1){
    stop(paste("ERROR: You should only have 1 Negative mode Feature Table... we detected", length(FeatureTable_NEG)," Feature Tables in the folder:", fileName))
  }else if(length(FeatureTable_NEG) == 0){
    print(paste("CAUTION: Could not find any negative mode Feature Tables... we detected", length(FeatureTable_NEG)," Feature Tables in the folder: ", fileName," ...If this incorrect, check that you have an 'n', 'N', 'neg', or 'NEG' at the end of the file name. OR Remove the folder: ", fileName))
  }
  
  FeatureTable_POS <- list.files(path=fpath, pattern="[PpSs]\\.csv", ignore.case=FALSE)
  if(length(FeatureTable_POS) > 1){
    stop(paste("ERROR: You should only have 1 Positive mode Feature Table... we detected", length(FeatureTable_POS)," Feature Tables in the folder:", fileName))
  }else if(length(FeatureTable_POS) == 0){
    print(paste("CAUTION: Could not find any Positive mode Feature Tables... we detected", length(FeatureTable_POS)," Feature Tables in the folder: ", fileName," ...If this incorrect, check that you have an 'p', 'P', 'pos', or 'POS' at the end of the file name. OR Remove the folder: ", fileName))
  }
  
  #Negative/Positive mode sample names (took .ms2 files and dropped the ".ms2")
  ExtraSampleNameddMSNEG_in <- vector()
  ExtraSampleNameddMSPOS_in <- vector()
  ExtraSampleNameAIFNEG_in <- vector()
  ExtraSampleNameAIFPOS_in <- vector()
  for(j in seq_len(length(ddMS2NEG_in))){   ExtraSampleNameddMSNEG_in[j] <- sub("\\.\\w+", "", ddMS2NEG_in[j])    }
  for(j in seq_len(length(ddMS2POS_in))){   ExtraSampleNameddMSPOS_in[j] <- sub("\\.\\w+", "", ddMS2POS_in[j])    }
  for(j in seq_len(length(AIFMS2NEG_in))){   ExtraSampleNameAIFNEG_in[j] <- sub("\\.\\w+", "", AIFMS2NEG_in[j])    }
  for(j in seq_len(length(AIFMS2POS_in))){   ExtraSampleNameAIFPOS_in[j] <- sub("\\.\\w+", "", AIFMS2POS_in[j])    }
  
  runPosddMS <- FALSE
  runNegddMS <- FALSE
  runPosAIF <- FALSE
  runNegAIF <- FALSE
  #Run Negative mode analysis if there are negative .ms2 and .csv files
  if(length(ddMS2NEG_in) != 0 && length(FeatureTable_NEG) == 1){  runNegddMS <- TRUE  }
  #Run Positive mode analysis if there are positive .ms2 and .csv files
  if(length(ddMS2POS_in) != 0 && length(FeatureTable_POS) == 1) {
    runPosddMS <- TRUE
  }
  if(length(AIFMS2NEG_in) != 0 && length(AIFMS1NEG_in) != 0 && length(FeatureTable_NEG) == 1){  runNegAIF <- TRUE   }
  if(length(AIFMS2POS_in) != 0 && length(AIFMS1POS_in) != 0 && length(FeatureTable_POS) == 1){  runPosAIF <- TRUE   }
  
  #Create output file structure
  #Shrimp
  #--AIF
  #----Neg
  #------Additional_Files
  #------Confirmed_Lipids
  #----Pos
  #------Additional_Files
  #------Confirmed_Lipids
  #--ddMS
  #----Neg
  #------Additional_Files
  #------Confirmed_Lipids
  #----Pos
  #------Additional_Files
  #------Confirmed_Lipids
  #----PosByClass
  #------Additional_Files
  #------Confirmed_Lipids
  
  if(length(foldersToRun)==0){
    #1 root folder
    OutputDirectory<-file.path(InputDirectory, "Output",sep="")
    if(!dir.exists(OutputDirectory)){ dir.create(OutputDirectory) }
    if(runPosAIF || runNegAIF){#AIF
      OutputDirectoryAIF <- file.path(InputDirectory, "Output", "AIF")
      if(!dir.exists(OutputDirectoryAIF)){  dir.create(OutputDirectoryAIF)  }
      if(runPosAIF){#pos AIF
        OutputDirectoryAIFPos_in <- file.path(OutputDirectoryAIF,"Pos")
        if(!dir.exists(OutputDirectoryAIFPos_in)){ dir.create(OutputDirectoryAIFPos_in) }
      }
      if(runNegAIF){#neg AIF
        OutputDirectoryAIFNeg_in <- file.path(OutputDirectoryAIF,"Neg")
        if(!dir.exists(OutputDirectoryAIFNeg_in)){ dir.create(OutputDirectoryAIFNeg_in) }
      }
    }
    if(runPosddMS || runNegddMS){#ddMS
      # OutputDirectoryddMS <- file.path(InputDirectory, "Output", "ddMS") # PJS 11/13/2022
      OutputDirectoryddMS <- file.path(InputDirectory, "Output", "ddMS") # PJS 11/13/2022
      if(!dir.exists(OutputDirectoryddMS)){ dir.create(OutputDirectoryddMS) }
      if(runPosddMS){#pos ddMS
        OutputDirectoryddMSPos_in <- file.path(OutputDirectoryddMS,"Pos")
        if(!dir.exists(OutputDirectoryddMSPos_in)){ dir.create(OutputDirectoryddMSPos_in) }
      }
      if(runPosddMS){#posByClass ddMS
        OutputDirectoryddMSPosByClass_in <- file.path(OutputDirectoryddMS,"PosByClass")
        if(!dir.exists(OutputDirectoryddMSPosByClass_in)){ dir.create(OutputDirectoryddMSPosByClass_in) }
      }
      if(runNegddMS){#negByClass ddMS
        OutputDirectoryddMSNegByClass_in <- file.path(OutputDirectoryddMS,"NegByClass")
        if(!dir.exists(OutputDirectoryddMSNegByClass_in)){ dir.create(OutputDirectoryddMSNegByClass_in) }
      }
      if(runNegddMS){#neg ddMS
        OutputDirectoryddMSNeg_in <- file.path(OutputDirectoryddMS,"Neg")
        if(!dir.exists(OutputDirectoryddMSNeg_in)){ dir.create(OutputDirectoryddMSNeg_in) }
      }
    }
  }else{#more than 1 root folder
    OutputDirectory <- file.path(InputDirectory, foldersToRun[i], "Output")
    if(!dir.exists(OutputDirectory)){ dir.create(OutputDirectory) }
    if(runPosAIF || runNegAIF){#AIF
      OutputDirectoryAIF <- file.path(OutputDirectory, "AIF")
      if(!dir.exists(OutputDirectoryAIF)){  dir.create(OutputDirectoryAIF)  }
      if(runPosAIF){#pos AIF
        OutputDirectoryAIFPos_in <- file.path(OutputDirectoryAIF, "Pos")
        if(!dir.exists(OutputDirectoryAIFPos_in)){ dir.create(OutputDirectoryAIFPos_in) }
      }
      if(runNegAIF){#neg AIF
        OutputDirectoryAIFNeg_in <- file.path(OutputDirectoryAIF, "Neg")
        if(!dir.exists(OutputDirectoryAIFNeg_in)){ dir.create(OutputDirectoryAIFNeg_in) }
      }
    }
    if(runPosddMS || runNegddMS){#ddMS
      OutputDirectoryddMS <- file.path(OutputDirectory, "ddMS")
      if(!dir.exists(OutputDirectoryddMS)){ dir.create(OutputDirectoryddMS) }
      if(runPosddMS){#pos ddMS
        OutputDirectoryddMSPos_in <- file.path(OutputDirectoryddMS, "Pos")
        if(!dir.exists(OutputDirectoryddMSPos_in)){ dir.create(OutputDirectoryddMSPos_in) }
      }
      if(runPosddMS){#posByClass ddMS
        OutputDirectoryddMSPosByClass_in <- file.path(OutputDirectoryddMS, "PosByClass")
        if(!dir.exists(OutputDirectoryddMSPosByClass_in)){ dir.create(OutputDirectoryddMSPosByClass_in) }
      }
      if(runNegddMS){#negByClass ddMS
        OutputDirectoryddMSNegByClass_in <- file.path(OutputDirectoryddMS,"NegByClass")
        if(!dir.exists(OutputDirectoryddMSNegByClass_in)){ dir.create(OutputDirectoryddMSNegByClass_in) }
      }
      if(runNegddMS){#neg ddMS
        OutputDirectoryddMSNeg_in <- file.path(OutputDirectoryddMS, "Neg")
        if(!dir.exists(OutputDirectoryddMSNeg_in)){ dir.create(OutputDirectoryddMSNeg_in) }
      }
    }
  }#end else
  
  
  #### Run the libraries and input data ####
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
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_NEG)
    cat(paste0("Reading in file:\t", FeatureTable_NEG,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
    FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
    # Create list of dataframes with all ms2 data

    length_ddMS2NEG_in <- length(ddMS2NEG_in)
    MS2_df_list <- vector("list", length_ddMS2NEG_in)
    nrow_all_scans <- 0
    for (c in 1:length_ddMS2NEG_in){
      MS2_dir_in <- file.path(fpath, ddMS2NEG_in[c])
      MS2_df_in <- createDataFrame(MS2_dir_in)
      MS2_df_list[[c]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
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
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_NEG)
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
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_POS)
    cat(paste0("Reading in file:\t", FeatureTable_POS,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
    FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
    # Create list of dataframes with all ms2 data

    length_ddMS2POS_in <- length(ddMS2POS_in)
    MS2_df_list <- vector("list", length_ddMS2POS_in)
    nrow_all_scans <- 0
    for (c in 1:length_ddMS2POS_in){
      MS2_dir_in <- file.path(fpath, ddMS2POS_in[c])
      MS2_df_in <- createDataFrame(MS2_dir_in)
      MS2_df_list[[c]] <- MS2_df_in
      nrow_MS2_df_in <- nrow(MS2_df_in)
      for(j in 1:nrow_MS2_df_in){
        nrow_all_scans <- nrow_all_scans + nrow(MS2_df_in[j,][[3]][[1]])
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
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_POS)
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
    print(paste("Finished Posative by class ddMS analysis", timestamp(), sep = " --> "))
    # Export MSMS combined file
    # write.csv(Pos_rawMSMS, file.path(InputDirectory, "Output/Pos_rawMSMS_Example3.csv"), row.names = FALSE, col.names = TRUE, na = "")
  }
  
  ####AIF####
  #Neg AIF
  LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
  LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "NEG",] #subset LibraryCriteria to find negative class libraries
  LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,5]) == "TRUE",] #subset LibraryCriteria to find AIF libraries to run
  if(runNegAIF && nrow(LibraryCriteria)>0){
    NegAIFLib <- TRUE
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_NEG)
    cat(paste0("Reading in file:\t", FeatureTable_NEG,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
    FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
    #sort AIF files
    AIFMS1NEG_in <- AIFMS1NEG_in[order(AIFMS1NEG_in)]
    AIFMS2NEG_in <- AIFMS2NEG_in[order(AIFMS2NEG_in)]
    length_AIFMS1NEG_in <- length(AIFMS1NEG_in)
    for (c in 1:length_AIFMS1NEG_in){
      MS1_dir_in <- file.path(fpath, AIFMS1NEG_in[c])
      cat(paste0("Reading in file:\t", AIFMS1NEG_in[c],"\nFrom Directory:\t\t", MS1_dir_in,"\n"))
      MS1_df_in <- createDataFrame(MS1_dir_in)
      
      MS2_dir_in <- file.path(fpath, AIFMS2NEG_in[c])
      cat(paste0("Reading in file:\t", AIFMS2NEG_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
      MS2_df_in <- createDataFrame(MS2_dir_in)
      
      ExtraSample<-ExtraSampleNameAIFNEG_in[c]
      OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in, MS1_dir_in)
      if (ParallelComputing == TRUE) {
        nrow_LibraryCriteria <- nrow(LibraryCriteria)
        foreach (i = seq_len(nrow_LibraryCriteria), .packages = c("sqldf")) %dopar% {
          LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
          OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
          ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
          ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
          if(length(ConfirmANDCol)==0 || is.na(ConfirmANDCol)){ConfirmANDCol<-NULL}
          if(length(ConfirmORCol)==0 || is.na(ConfirmORCol)){ConfirmORCol<-NULL}
          RunAIF(MS1_df_in, MS2_df_in, FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryAIFNeg_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo)
        }
      } else {
        nrow_LibraryCriteria <- nrow(LibraryCriteria)
        for(i in seq_len(nrow_LibraryCriteria)){
          LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
          OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
          ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
          ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
          if(length(ConfirmANDCol)==0 || is.na(ConfirmANDCol)){ConfirmANDCol<-NULL}
          if(length(ConfirmORCol)==0 || is.na(ConfirmORCol)){ConfirmORCol<-NULL}
          RunAIF(MS1_df_in, MS2_df_in, FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryAIFNeg_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo)
        }
      }
    }
    print(paste("Finished Negative AIF analysis", timestamp(), sep = " --> "))
  }
  
  #Pos AIF
  LibraryCriteria <- read.csv(LibCriteria) #Read-in Library ID criteria (csv) located in the LibrariesReducedAdducts folder
  LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,7]) == "POS",] #subset LibraryCriteria to find positive class libraries
  LibraryCriteria <- LibraryCriteria[toupper(LibraryCriteria[,5]) == "TRUE",] #subset LibraryCriteria to find AIF libraries to run
  if(runPosAIF && nrow(LibraryCriteria)>0){
    PosAIFLib <- TRUE
    FeatureTable_dir_in<-file.path(fpath, FeatureTable_POS)
    cat(paste0("Reading in file:\t", FeatureTable_POS,"\nFrom Directory:\t\t", FeatureTable_dir_in,"\n"))
    FeatureList_in <- ReadFeatureTable(FeatureTable_dir_in)
    # NegByClass_rawMSMS <- FeatureTable(FeatureTable_dir_in)
    AIFMS1POS_in <- AIFMS1POS_in[order(AIFMS1POS_in)]
    AIFMS2POS_in <- AIFMS2POS_in[order(AIFMS2POS_in)]
    length_AIFMS1POS_in <- length(AIFMS1POS_in)
    for (c in 1:length_AIFMS1POS_in){
      MS1_dir_in <- file.path(fpath, AIFMS1POS_in[c])
      cat(paste0("Reading in file:\t", AIFMS1POS_in[c],"\nFrom Directory:\t\t", MS1_dir_in,"\n"))
      MS1_df_in <- createDataFrame(MS1_dir_in)
      
      MS2_dir_in <- file.path(fpath, AIFMS2POS_in[c])
      cat(paste0("Reading in file:\t", AIFMS2POS_in[c],"\nFrom Directory:\t\t", MS2_dir_in,"\n"))
      MS2_df_in <- createDataFrame(MS2_dir_in)
      
      ExtraSample<-ExtraSampleNameAIFPOS_in[c]
      OutputInfo <- c(MS2_dir_in, ExtraSample, FeatureTable_dir_in, MS1_dir_in)
      if (ParallelComputing == TRUE) {
        nrow_LibraryCriteria <- nrow(LibraryCriteria)
        foreach (i = seq_len(nrow_LibraryCriteria), .packages = c("sqldf")) %dopar% {
          LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
          OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
          ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
          ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
          if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
          if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
          RunAIF(MS1_df_in, MS2_df_in, FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryAIFPos_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo)
        }
      } else {
        nrow_LibraryCriteria <- nrow(LibraryCriteria)
        for(i in seq_len(nrow_LibraryCriteria)){
          LibraryFile <- file.path(InputLibrary, LibraryCriteria[i,1]) #create directory/file of each library
          OutputName <- paste(ExtraSample,"_",gsub('.{4}$', '', LibraryCriteria[i,1]), sep="") #get ms2 name and library name
          ConfirmANDCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,2]), ";")))
          ConfirmORCol <- as.numeric(unlist(strsplit(as.character(LibraryCriteria[i,3]), ";")))
          if(length(ConfirmANDCol)==0){ConfirmANDCol<-NULL}
          if(length(ConfirmORCol)==0){ConfirmORCol<-NULL}
          RunAIF(MS1_df_in, MS2_df_in, FeatureList_in, LibraryFile, ParentMZcol_in, OutputDirectoryAIFPos_in, OutputName, ConfirmORCol, ConfirmANDCol, OutputInfo)
        }
      }
    }
    print(paste("Finished Positive AIF analysis", timestamp(), sep = " --> "))
  }
  
  #Compilation/ID code for reduced confirmed files
  if(runPosAIF || runPosddMS){
    print(paste("Creating Identifications for Positive Mode", timestamp(), sep = " --> "))
  }
  if(runPosddMS & !runPosAIF){
    ddMS2directory<-file.path(OutputDirectoryddMSPos_in,"Confirmed_Compounds")
    Classdirectory<-file.path(OutputDirectoryddMSPosByClass_in,"Confirmed_Compounds")
    AIFdirectory<-"Nothing"
    CreateIDs(file.path(fpath,FeatureTable_POS), ddMS2directory, Classdirectory, AIFdirectory, ImportLibPOS, OutputDirectory, PosDDLib, PosClassDDLib, PosAIFLib, "Pos")
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
  
  if(runPosAIF & runPosddMS){
    ddMS2directory<-file.path(OutputDirectoryddMSPos_in,"Confirmed_Compounds")
    Classdirectory<-file.path(OutputDirectoryddMSPosByClass_in,"Confirmed_Compounds")
    AIFdirectory<-file.path(OutputDirectoryAIFPos_in,"Confirmed_Compounds")
    CreateIDs(file.path(fpath,FeatureTable_POS), ddMS2directory, Classdirectory, AIFdirectory, ImportLibPOS, OutputDirectory, PosDDLib, PosClassDDLib, PosAIFLib, "Pos")
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
  
  if(runPosAIF & !runPosddMS){
    ddMS2directory<-"Nothing"
    Classdirectory<-"Nothing"
    AIFdirectory<-file.path(OutputDirectoryAIFPos_in,"Confirmed_Compounds")
    CreateIDs(file.path(fpath,FeatureTable_POS), ddMS2directory, Classdirectory, AIFdirectory, ImportLibPOS, OutputDirectory, PosDDLib, PosClassDDLib, PosAIFLib, "Pos")
  }
  
  if(runNegAIF || runNegddMS){
    print(paste("Creating Identifications for Negative Mode", timestamp(), sep = " --> "))
  }
  
  if(runNegAIF & !runNegddMS){
    ddMS2directory <- "Nothing"
    Classdirectory <- "Nothing"
    AIFdirectory <- file.path(OutputDirectoryAIFNeg_in,"Confirmed_Compounds")
    CreateIDs(file.path(fpath,FeatureTable_NEG), ddMS2directory, Classdirectory, AIFdirectory, ImportLibNEG, OutputDirectory, NegDDLib, NegClassDDLib, NegAIFLib, "Neg")
  }
  
  if(runNegddMS & !runNegAIF){
    ddMS2directory <- file.path(OutputDirectoryddMSNeg_in,"Confirmed_Compounds")
    Classdirectory <- file.path(OutputDirectoryddMSNegByClass_in,"Confirmed_Compounds")
    AIFdirectory <- "Nothing"
    CreateIDs(file.path(fpath,FeatureTable_NEG), ddMS2directory, Classdirectory, AIFdirectory, ImportLibNEG, OutputDirectory, NegDDLib, NegClassDDLib, NegAIFLib, "Neg")
    #InputDirectory<-"C:/Users/jpk72/Desktop/OUT/LipidMatch_Run/"; CommentColumn<-1; RowStartForFeatureTableData<-2; NegPos = "Neg"; OutDir="/Output/NegIDed_FIN.csv"; OutDirOnlyFrags="/Output/Neg_OnlyIDs_PredFrags.csv"; InputDir_Append="Output/ddMS/NegByClass/Additional_Files"; ID_name="PredictedFrag_IDs"; FragName="Frags"; nFrag="Num_Frags"; fileNames="Files"; ImportTable="/Output/NegIDed_Fragments.csv"
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
  
  if(runNegddMS & runNegAIF){
    ddMS2directory <- file.path(OutputDirectoryddMSNeg_in,"Confirmed_Compounds")
    Classdirectory <- file.path(OutputDirectoryddMSNegByClass_in,"Confirmed_Compounds")
    AIFdirectory <- file.path(OutputDirectoryAIFNeg_in,"Confirmed_Compounds")
    CreateIDs(file.path(fpath,FeatureTable_NEG), ddMS2directory, Classdirectory, AIFdirectory, ImportLibNEG, OutputDirectory, NegDDLib, NegClassDDLib, NegAIFLib, "Neg")
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
  
  if(runNegddMS & runPosddMS){
    Neg <- read.csv(file.path(OutputDirectory, "NegIDed_Fragments.csv"), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
    Neg <- as.matrix(Neg)
    ncol_Neg <- ncol(Neg)
    extra <- matrix("", nrow(Neg), 100 - ncol_Neg)
    Neg <- cbind(Neg, extra)
    Pos <- read.csv(file.path(OutputDirectory, "PosIDed_Fragments.csv"), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
    Pos <- as.matrix(Pos)
    ncol_Pos <- ncol(Pos)
    extra <- matrix("", nrow(Pos), 100 - ncol_Pos)
    Pos <- cbind(Pos, extra)
    
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
    MZ_Append <- rev(as.numeric(FeatureList_in[, 6])) # Column should be constant (future error?)
    RT_Append <- rev(as.numeric(FeatureList_in[, 7])) # Column should be constant
    IDs_Append <- rev(as.numeric(FeatureList_in[, 12])) # Column should be constant
    Precursor_rawMSMS <- as.numeric(rawMSMS_df[, PrecursorMZCol_MSMS])
    RT_rawMSMS <- as.numeric(rawMSMS_df[, MSMS_scan_RT_col])
    for (j in 2:nrow_FeatureList_in) { # Row should be constant (future error?)
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
      ppm_error <- ((Precursor_rawMSMS - MZ_Append_current)*10^6)/Precursor_rawMSMS
      MZ_Conditional <- abs(ppm_error)<(ppm_Window/2)
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
  
}#end folder loop

options(warn=0)#suppress warning off

#Rversion<-(paste("ERROR:R version must be equal to, or between, 2.0.3 and 3.3.3. Please download 3.3.3. You are using version: ", paste(version$major,version$minor,sep=".")))
#OutputRemoval<-paste("ERROR: Remove your 'Output' folder from the current Input Directory: ", InputDirectory)

###Code to append fragments and tentative annotations###




#DEBUG CreateIDs
# ddMS2directory<-"Nothing"
# Classdirectory<-"Nothing"
# AIFdirectory<-paste(OutputDirectoryAIFPos_in,"Confirmed_Lipids\\", sep="")
# PeakTableDirectory <- paste(fpath,FeatureTable_POS,sep="")
# ddMS2directory <- ddMS2directory
# Classdirectory <- Classdirectory
# AIFdirectory <- AIFdirectory
# ImportLib <- ImportLibPOS
# OutputDirectory <- OutputDirectory
# ddMS2 <- PosDDLib
# ddMS2Class <- PosClassDDLib
# AIF <- PosAIFLib
# mode <- "Pos"

#CreateIDs(paste(fpath,FeatureTable_POS,sep=""), ddMS2directory, Classdirectory,
#AIFdirectory, ImportLibPOS, OutputDirectory, PosDDLib, PosClassDDLib, PosAIFLib, "Pos")

#Debug AIF
# ms1_df<-MS1_df_in
# ms2_df<-MS2_df_in
# FeatureList<-FeatureList_in
# LibraryLipid_self<-LibraryFile
# ParentMZcol<-ParentMZcol_in
# OutputDirectory<-OutputDirectoryAIFNeg_in
# ExtraFileNameInfo<-OutputName
# ConfirmORcol<-ConfirmORCol
# ConfirmANDcol<-ConfirmANDCol

###########################STATS#####################################################
setwd(paste(InputLibrary,"/Scripts",sep=""))
source(paste(InputLibrary,"/Scripts/Stats.R",sep=""))

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
if (runPosddMS&&runNegddMS) {
  if (FLOW || csvInput) {
    if (file.exists(GroupCSVDirectory)&&file.size(GroupCSVDirectory)>2) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Combined")
    }
    else
    {
      statsFilepath<-DetermineFilePath("Combined")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      output<-dataSet
    }
  }else{
    if (length(GroupCSVDirectory)>0) {
      output <- OutputData(OutputDirectory,GroupCSVDirectory,CommentColumn,"Combined")
    }
    else{
      statsFilepath<-DetermineFilePath("Combined")
      data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
      dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
      output<-dataSet
    }
  }
  #provide file path for importing and exporting
  outputFile<-DetermineFilePath("Combined")
  write.csv(output,outputFile, row.names = FALSE)
}

#######################EIC and MS1 R Scripts###################################
if(FLOW==TRUE) {
  target_mzXML = paste(dirname(dirname(OutputDirectory)),"/Temp_Work/",sep="")
} else { ##Modular
  target_mzXML = InputDirectory
}

if(TargetEIC_Only==TRUE) {
  mzXML_Files<-list.files(target_mzXML, pattern = ".mzXML", full.names = FALSE)
  Target_mzXML<-grep(mzXML_Files, pattern = "Target|target|TARGET|Blank|blank|BLANK", invert=FALSE, value=TRUE)
} else {
  Target_mzXML<-list.files(target_mzXML, pattern = ".mzXML", full.names = FALSE)
}

source(paste(InputLibrary,"/Scripts/EIC_MS1_fns.R",sep=""))
source(paste(InputLibrary,"/Scripts/genEIC.R",sep=""))
source(paste(InputLibrary,"/Scripts/genIsoTable.R",sep=""))
source(paste(InputLibrary,"/Scripts/MS1Spectragen.R",sep=""))
source(paste(InputLibrary,"/Scripts/IsotopePercentages.R",sep=""))
source(paste(InputLibrary,"/Scripts/Kaufmann.R",sep=""))
source(paste(InputLibrary,"/Scripts/Manual_Review.R",sep=""))

arguments = construct_EM_arguments(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Window = RT_Window
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = c(7,8,1,5)
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = paste(InputLibrary,"/Scripts/secondary_isotopes.csv",sep="")
)
arguments$path_to_mzXML_Files = target_mzXML

runallEM(Target_mzXML, arguments, runmode = runNegddMS)
runallEM(Target_mzXML, arguments, runmode = runPosddMS, isNeg = FALSE)

#Kaufmann eC, and ratios, values
MD_col_name <- "mass.defect"
MZ_col_name <- "m.z"
C13_col_name <- "13C1"
if (runPosddMS) {
  outputFile<-DetermineFilePath("Pos")
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot inforation complete, positive mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
}
if (runNegddMS) {
  outputFile<-DetermineFilePath("Neg")
  Kaufmann_eCs(outputFile,MD_col_name,MZ_col_name,C13_col_name)
  print("Kaufmann plot inforation complete, negative mode")
  #Column Headers for Editing and Markup
  Review_Col_Names(outputFile)
}
