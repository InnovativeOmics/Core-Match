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
    #No longer used
    numOfAIFFiles <- length(list.files(path=fpath, pattern="[AIFaif][AIFaif][AIFaif]", ignore.case=FALSE))
    if(numOfAIFFiles %% 2 == 1){ #if # of AIF .ms2 or .ms1 are odd, STOP!
      stop("Warning: You should have a one-to-one relationship between AIF.ms1 and AIF.ms2 files...We found ", numOfAIFFiles," AIF .ms1 and/or .ms2 files. You should have an even number of AIF files.",sep="")
    }#else, do nothing.
  }
}