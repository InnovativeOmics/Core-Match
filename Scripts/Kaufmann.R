#Jeremy Koelmel
#David 
#NegIDed_Fin.csv, 
Kaufmann_eCs <- function(outputFile,MD_col_name,MZ_col_name,C13_col_name){
    # Get the symbol from the iso string
    NegIDed_Fin_eCs<-read.csv(outputFile,sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=FALSE)
    NegIDed_Fin_eCs<-as.matrix(NegIDed_Fin_eCs)
    MDcol<-which(NegIDed_Fin_eCs[1,]==MD_col_name)
    MZcol<-which(NegIDed_Fin_eCs[1,]==MZ_col_name)
    C13_col<-which(NegIDed_Fin_eCs[1,]==C13_col_name)
    eC<-as.numeric(NegIDed_Fin_eCs[RowStartForFeatureTableData:nrow(NegIDed_Fin_eCs),C13_col])/1.0816
    MZ_eC<-as.numeric(NegIDed_Fin_eCs[RowStartForFeatureTableData:nrow(NegIDed_Fin_eCs),MZcol])/eC
    MD_eC<-as.numeric(NegIDed_Fin_eCs[RowStartForFeatureTableData:nrow(NegIDed_Fin_eCs),MDcol])/eC
    MZ_eC[!is.finite(MZ_eC)]<-0
    MZ_eC<-append("MZ_eC",MZ_eC)
    MD_eC[!is.finite(MD_eC)]<-0
    MD_eC<-append("MD_eC",MD_eC)
    eC<-append("eC",eC)
    NegIDed_Fin_eCs<-cbind(NegIDed_Fin_eCs,eC,MZ_eC,MD_eC)
    # m <- matrix("", nrow(NegIDed_Fin_eCs), 5)
    # NegIDed_Fin_eCs <- cbind(NegIDed_Fin_eCs, m)
    # ncol_IDedTable <- ncol(IDedTable)
    # NegIDed_Fin_eCs[1,(ncol(NegIDed_Fin_eCs)- 4):ncol(NegIDed_Fin_eCs)] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
    write.table(NegIDed_Fin_eCs, outputFile, row.names = FALSE, col.names = FALSE, na = "", sep=",")
}