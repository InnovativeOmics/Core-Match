
# OutputDirectory<-"D:/FluoroMatch_Data/2024_07_01_SarahStow_IM_PNNLtool/2024_07_14_NIST_C_FIN/IM_Only/Modular/Output/"
# ScoreMSMS_colName<-"Score"
# ScoreIM_colName<-"Score.1"
# Name_or_ClassMSMS_colName<-"Name_or_Class"
# Name_or_ClassIM_colName<-"Name_or_Class.1"
# FormulaMSMS_colName<-"Formula"
# FormulaIM_colName<-"Formula.1"
# SMILESMSMS_colName<-"SMILES"
# SMILESIM_colName<-"SMILES.1"
##Note the renaming of columns would be much more stable as soon as the feature table is imported to not assume overwrite order

CombineIMandMSMS<-function(OutputDirectory, ScoreMSMS_colName,ScoreIM_colName,Name_or_ClassMSMS_colName,Name_or_ClassIM_colName,FormulaMSMS_colName, FormulaIM_colName,	SMILESMSMS_colName, SMILESIM_colName){
  IDed_Fin_IMcomb_Import<-read.csv(paste0(OutputDirectory,"/NegIDed_FIN.csv"), sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  #Insert 5 new columns to add a combined (top) score formula SMILES and name, as well as whether from IM or MSMS
  NewCols_Comb_IM_MSMS<-matrix(NA,nrow(IDed_Fin_IMcomb_Import),5)
  IDed_Fin_IMcomb_Import<-cbind(NewCols_Comb_IM_MSMS,IDed_Fin_IMcomb_Import)
  
  #Get the columns of which are IM and which are MSMS based on name, right now default inputs assume .1 after certain names and an order (for which gets .1 by R by default)
  ScoreMSMS_col<-which(colnames(IDed_Fin_IMcomb_Import) == ScoreMSMS_colName)
  ScoreIM_col<-which(colnames(IDed_Fin_IMcomb_Import) == ScoreIM_colName)
  Name_or_ClassMSMS_col<-which(colnames(IDed_Fin_IMcomb_Import) == Name_or_ClassMSMS_colName)
  Name_or_ClassIM_col<-which(colnames(IDed_Fin_IMcomb_Import) == Name_or_ClassIM_colName)
  FormulaMSMS_col<-which(colnames(IDed_Fin_IMcomb_Import) == FormulaMSMS_colName)
  FormulaIM_col<-which(colnames(IDed_Fin_IMcomb_Import) == FormulaIM_colName)
  SMILESMSMS_col<-which(colnames(IDed_Fin_IMcomb_Import) == SMILESMSMS_colName)
  SMILESIM_col<-which(colnames(IDed_Fin_IMcomb_Import) == SMILESIM_colName)
  
  #rename the columns of which are IM and which are MSMS based on name, right now default inputs assume .1 after certain names and an order
  colnames(IDed_Fin_IMcomb_Import)[ScoreMSMS_col] <- "Score_MSMS"
  colnames(IDed_Fin_IMcomb_Import)[ScoreIM_col] <- "Score_IM"
  colnames(IDed_Fin_IMcomb_Import)[Name_or_ClassMSMS_col] <- "Name_or_Class_MSMS"
  colnames(IDed_Fin_IMcomb_Import)[Name_or_ClassIM_col] <- "Name_or_ClassIM"
  colnames(IDed_Fin_IMcomb_Import)[FormulaMSMS_col] <- "Formula_MSMS"
  colnames(IDed_Fin_IMcomb_Import)[FormulaIM_col] <- "Formula_IM"
  colnames(IDed_Fin_IMcomb_Import)[SMILESMSMS_col] <- "SMILES_IM"
  colnames(IDed_Fin_IMcomb_Import)[SMILESIM_col] <- "SMILES_MSMS"
  colnames(IDed_Fin_IMcomb_Import)[1] <- "Score"
  colnames(IDed_Fin_IMcomb_Import)[2] <- "Name_or_Class"
  colnames(IDed_Fin_IMcomb_Import)[3] <- "Formula"
  colnames(IDed_Fin_IMcomb_Import)[4] <- "SMILES"
  colnames(IDed_Fin_IMcomb_Import)[5] <- "top_is_MSMS"
  
  #rank scores by numbers, to determine whether to use IM or MSMS
  TempRankScoreIM<-as.numeric(str_replace_all(as.character(IDed_Fin_IMcomb_Import[,ScoreIM_col]), c("A\\-" = "2", "A" = "1", "B\\-\\-" = "6","B\\-" = "5", "B\\+" = "3", "B" = "4", "C\\-" = "9", "D\\-" = "12","C\\+" = "7", "C" = "8", "D\\+" = "10","D" = "11","E" = "13")))
  TempRankScoreMSMS<-as.numeric(str_replace_all(as.character(IDed_Fin_IMcomb_Import[,ScoreMSMS_col]), c("A\\-" = "2", "A" = "1", "B\\-\\-" = "6","B\\-" = "5", "B\\+" = "3", "B" = "4", "C\\-" = "9", "D\\-" = "12","C\\+" = "7", "C" = "8", "D\\+" = "10","D" = "11","E" = "13")))
  Logical_MSMS_TRUE<-TempRankScoreMSMS<TempRankScoreIM
  IDed_Fin_IMcomb_Import[Logical_MSMS_TRUE,1:4]<-IDed_Fin_IMcomb_Import[Logical_MSMS_TRUE,c(ScoreMSMS_col,Name_or_ClassMSMS_col,FormulaMSMS_col,SMILESMSMS_col)]
  IDed_Fin_IMcomb_Import[!Logical_MSMS_TRUE,1:4]<-IDed_Fin_IMcomb_Import[!Logical_MSMS_TRUE,c(ScoreIM_col,Name_or_ClassIM_col,FormulaIM_col,SMILESIM_col)]
  IDed_Fin_IMcomb_Import[,5]<-Logical_MSMS_TRUE
  
  write.csv(IDed_Fin_IMcomb_Import, file.path(OutputDirectory, "/NegIDed_FIN.csv"), row.names = FALSE, col.names = TRUE, na = "")
}