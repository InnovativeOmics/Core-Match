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
# InputDirectory<-"D:/FluoroMatch_Data/2024_07_01_SarahStow_IM_PNNLtool/3D_DDA_Outputs_Modular/"; CommentColumn<-1; RowStartForFeatureTableData<-2; NegPos = "Neg"; OutDir="/Output/NegIDed_FIN.csv"; OutDirOnlyFrags="/Output/Neg_OnlyIDs_PredFrags.csv"; InputDir_Append="Output/ddMS/NegByClass/Additional_Files"; ID_name="PredictedFrag_IDs"; FragName="Frags"; nFrag="Num_Frags"; fileNames="Files"; ImportTable="/Output/NegIDed_Fragments.csv"; firstFragAppender=FALSE
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
  }
  
  #This will replace all special characters
  # NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "[^[:alnum:]]", "_char_replaced_")
  
  #replace all special characters except certain one commonly needed in SMILES, Adducts, names, etc
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "[^A-Za-z0-9;\\ \\;\\:\\_\\|\\,\\&\\-\\(\\)/\\\\\\[\\]=]", "_char_replaced_")
  #Below are certain "escape characters" which could cause issues
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "\r", "_char_replaced_") #line breaks
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "\n", "_char_replaced_") #line breaks
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "\"", "_char_replaced_") #double quotes (redundant?)
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "\t", "_char_replaced_") #tab delimiter
  NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3]<-str_replace_all(NegIDed[2:nrow(NegIDed),ncol(NegIDed)-3], "^=", "_char_replaced_") #leading = for excel issues
  
  #cutoff the column with predicted fragment results if it is too long for certain other software limits for a cell
  NegIDed[,ncol(NegIDed)-3]<-substr(NegIDed[,ncol(NegIDed)-3],1,10000)
  # sort from that row down using the first element of number of fragments
  NegIDed[NoIDs_row:nrow(NegIDed),]<-NegIDed[order(NegIDed[NoIDs_row:nrow(NegIDed),ncol(NegIDed)-1],decreasing=TRUE)+NoIDs_row-1,]
  write.table(NegIDed, file.path(InputDirectory,OutDir), sep=",",col.names=FALSE, row.names=FALSE, quote=TRUE, na="NA")
}
