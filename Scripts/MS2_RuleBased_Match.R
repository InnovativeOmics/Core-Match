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
