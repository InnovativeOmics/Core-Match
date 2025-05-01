
# PeakTableDirectory <- file.path(fpath,FeatureTable_NEG)
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