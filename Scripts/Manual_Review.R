#Jeremy Koelmel
#NegIDed_Fin.csv, 
Review_Col_Names <- function(outputFile){
    NegIDed_Fin<-read.csv(outputFile,sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
    NegIDed_Fin<-as.matrix(NegIDed_Fin)
    Review_Cols<-which(NegIDed_Fin[1,]=="Checked_Viz"|NegIDed_Fin[1,]=="TRUE_Viz"|NegIDed_Fin[1,]=="Class_Viz"|NegIDed_Fin[1,]=="Cnumb_Viz"|NegIDed_Fin[1,]=="Comment_Viz")
    if(length(Review_Cols)>0){
      NegIDed_Fin<-NegIDed_Fin[,-Review_Cols]
    }
    NegIDed_Fin<-cbind(NegIDed_Fin,"Checked_Viz","TRUE_Viz","Class_Viz","Cnumb_Viz","Comment_Viz")
    NegIDed_Fin[RowStartForFeatureTableData:nrow(NegIDed_Fin),which(NegIDed_Fin[1,]=="Checked_Viz"|NegIDed_Fin[1,]=="TRUE_Viz"|NegIDed_Fin[1,]=="Class_Viz"|NegIDed_Fin[1,]=="Cnumb_Viz"|NegIDed_Fin[1,]=="Comment_Viz")]<-NA
    write.table(NegIDed_Fin, outputFile, row.names = FALSE, col.names = FALSE, na = "", sep=",")
}