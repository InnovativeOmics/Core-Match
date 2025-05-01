#function for maniuplating and grabbing/checking column names

ForceColumnNames<-function(PeakTableDirectory){
  PeakTable<-read.csv(PeakTableDirectory,sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  colnames(PeakTable)[CommentColumn] <- "row.ID"
  write.csv(PeakTable,PeakTableDirectory,sep=",",dec=".",row.names = FALSE)
}


#mz,rt,id,dt(optional isIM=True),formula
# column_names<-c("m.z","Retention.Time","row.ID","Formula")
get_column_indices <- function(FullDirectory, column_names) {
  IDed_Fin_Find_Cols<-read.csv(FullDirectory, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  # Check if all provided column names exist in the data frame
  if (!all(column_names %in% colnames(IDed_Fin_Find_Cols))) {
    stop("Some column names provided do not exist in the data frame.")
  }
  # Get the indices of the column names
  column_indices <- match(column_names, colnames(IDed_Fin_Find_Cols))
  return(column_indices)
}
