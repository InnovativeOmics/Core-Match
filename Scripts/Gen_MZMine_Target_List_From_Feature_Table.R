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