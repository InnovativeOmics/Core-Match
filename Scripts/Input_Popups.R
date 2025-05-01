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
  numOfddFiles <- length(list.files(path=fpath, pattern="[dD][dD][Mm][Ss]", ignore.case=FALSE))
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

