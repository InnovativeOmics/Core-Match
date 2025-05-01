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
