manualGrouping <- function(groupCsv){
  data<-read.csv(groupCsv,header = F)
  numOfSamples <- nrow(data)
  count<-nrow(data)
  
  container <- c(rep(" ",count))
  
  for (i in 1:count) {
    container[i]<-data[i,1]
    
  }
  
  return(container)
  
}

getGroupNames<-function(inputData){
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  count <- 0
  
  t<- c(rep(" "),numOfSamples)
  ind = 1;
  
  for (i in 1:length(numOfSamples)) 
  {
    if (grepl("sample",tolower(header[i]),fixed = T)) {
      tmp<- toString(header[i])
      p<-substring(tmp,1,unlist(gregexpr(("_"),tmp))[1])
      
      if ((is.element(p,t)) == FALSE) {
        count <- count+1
        t[ind]<-p
        ind<-ind+1
        
      }
      
      
    }
    
  }
  container <- c(rep(" ",count))
  
  for (i in 1:count) {
    container[i]<-t[i]
    
  }
  
  return(container)

}

averageData <- function(inputData,name, fc){
  t <- c(rep(0.0,nrow(inputData)))
  
  for (i in 1:ncol(inputData)) {
    t <- t + inputData[i]
  }
  
  t<-t/ncol(inputData)
  
  return(t)
  
}

getDataColumns<- function(inputData,name,fc,minimumVal){
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  count <- 0
  
  f <- c(rep(0.0,nrow(inputData)))
  
  for (i in 1:length(numOfSamples)) 
  {
    
    if (grepl(name,header[i],fixed = T)) {
      if (count == 0) {
        f <- data.frame(inputData[i])
      }else{
        f<- cbind(f,inputData[i])
        
      }
      count<- count +1
      
    }
    
  }
  
  set.seed(13)

  for (i in 1:nrow(f)) {
    for (j in 1:ncol(f)) {
      if (f[i,j] == 0) {
        f[i,j] = runif(1,min =(minimumVal/8),max = (minimumVal/2))
      }
    }
  }
  
  return(f)
}

averageNormData <- function(inputData,name, fc){
  numOfSamples <- rownames(inputData)
  header <- rownames(inputData)
  count <- 0
  
  fC <- c(rep(0.0,ncol(inputData)))
  
  for (i in 1:length(numOfSamples)) 
  {
    if (grepl(name,header[i],fixed = T)) {
      
      count <- count+1
      v1 <- c(inputData[i,])
      fc = fc + v1
      
    }
    
  }
  
  lowest <- 100000
  
  for (i in 1:length(nrow(inputData))) {
    if (as.numeric(fc[i])<lowest) {
      lowest <- as.numeric(fc[i])
    }
    
  }
  
  for (i in 1:length(nrow(inputData))) {
    if (fc[i] == 0 || !is.numeric(fc[i])) {
      fc[i] <- (lowest/5)
    }
    
  }
  
  fc <- fc/count
  
  return(fc)
  
}

DetermineFilePath<- function(type){
  statsFilepath <- "tmp"
  
  if (tolower(type) == "neg") {
    statsFilepath<-paste(OutputDirectory,"/NegIDed_FIN.csv",sep="")
  }
  else if (tolower(type) == "pos") {
    statsFilepath<-paste(OutputDirectory,"/PosIDed_FIN.csv",sep="")
  }
  else{
    statsFilepath<-paste(OutputDirectory,"/CombinedIDed_FIN.csv",sep="")
  }
  
  return(statsFilepath)
  
}

groupPCAData<- function(inputData,samples){
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  count <- 0

  for (i in 1:length(numOfSamples)) 
  {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        count <- count+1
        #v1 <- log(inputData[i])
        #print(v1)
      }
      
    }
    
  }
  
  #mean the rows and subtract all values in the row by the mean
  
  oData <- data.frame(matrix(0,nrow(data),count))
  
  ind = 1
  for (i in 1:length(numOfSamples)) 
  {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        v1 <- inputData[i]
        oData[ind] <- v1
        ind <- ind + 1
      }
      
    }
    
  }
  
  #Go through oData and see if its a 0 or a string, and change it to the 1/5 min of that row
  
  
  return(oData)

}

#Volcano's Calculations
{
  twoGroupingGetCalculation <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))

    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    return(tmpData)
  }

  threeGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData

    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    return(tmpData)
    
  }
  
  fourGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData

    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    return(tmpData)
    
  }
  
  fiveGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AE Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #DE calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    return(tmpData)
    
  }
  
  sixGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    g6Data <- data.frame(getDataColumns(tmpData,targets[6],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    g6Avg<- rbind(averageData(g6Data,targets[6],f))

    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AE Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AF Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    #DE calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #EF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    return(tmpData)
    
  }
  
  sevenGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    g6Data <- data.frame(getDataColumns(tmpData,targets[6],f,lowest))
    g7Data <- data.frame(getDataColumns(tmpData,targets[7],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    g6Avg<- rbind(averageData(g6Data,targets[6],f))
    g7Avg<- rbind(averageData(g7Data,targets[7],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AE Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AF Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AG Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }

    #DE calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #EF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #EG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    #FG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    return(tmpData)
    
  }
  
  eightGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    g6Data <- data.frame(getDataColumns(tmpData,targets[6],f,lowest))
    g7Data <- data.frame(getDataColumns(tmpData,targets[7],f,lowest))
    g8Data <- data.frame(getDataColumns(tmpData,targets[8],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    g6Avg<- rbind(averageData(g6Data,targets[6],f))
    g7Avg<- rbind(averageData(g7Data,targets[7],f))
    g8Avg<- rbind(averageData(g8Data,targets[8],f))
    
    for (i in 1:targets) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[1] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AE Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AF Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AG Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AH Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DE calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #EF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #EG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #EH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #FG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #FH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #GH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g7Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g7Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[7],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[7],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[7],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[7],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[7],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[7],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }

    return(tmpData)
    
    
  }
  
  nineGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    g6Data <- data.frame(getDataColumns(tmpData,targets[6],f,lowest))
    g7Data <- data.frame(getDataColumns(tmpData,targets[7],f,lowest))
    g8Data <- data.frame(getDataColumns(tmpData,targets[8],f,lowest))
    g9Data <- data.frame(getDataColumns(tmpData,targets[9],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    g6Avg<- rbind(averageData(g6Data,targets[6],f))
    g7Avg<- rbind(averageData(g7Data,targets[7],f))
    g8Avg<- rbind(averageData(g8Data,targets[8],f))
    g9Avg<- rbind(averageData(g9Data,targets[9],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))
    
    #AB Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g2Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[2],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[2],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[2],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[2],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[2],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[2],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AC calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AD calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AE Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AF Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #AG Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AH Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #AI Calculations
    {
      #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
      #Fold Change
      fc<- as.double(unlist(g1Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g1Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[1],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[1],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[1],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[1],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[1],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[1],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BC Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g3Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[3],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[3],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[3],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[3],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[3],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[3],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BD Calcclations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #BH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #BI calculations
    {
      #Fold Change
      fc<- as.double(unlist(g2Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g2Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[2],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[2],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[2],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[2],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[2],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[2],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CD calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g4Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[4],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[4],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[4],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[4],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[4],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[4],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CE Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #CF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #CI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g3Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g3Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[3],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[3],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[3],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[3],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[3],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[3],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DE calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g5Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[5],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[5],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[5],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[5],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[5],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[5],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DF calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #DH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #DI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g4Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g4Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[4],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[4],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[4],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[4],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[4],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[4],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #EF Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g6Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g6Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[6],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[6],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[6],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[6],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[6],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[6],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
    }
    
    #EG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #EH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #EI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g5Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g5Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[5],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[5],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[5],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[5],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[5],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[5],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #FG Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g7Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g7Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[7],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[7],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[7],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[7],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[7],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[7],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #FH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #FI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g6Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g6Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[6],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[6],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[6],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[6],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[6],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[6],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #GH Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g7Avg))/as.double(unlist(g8Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g7Data[i,],g8Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[7],"/",targets[8],"-FC",sep = "")
      rename2 <- paste(targets[7],"/",targets[8],"-pVal",sep = "")
      rename3 <- paste(targets[7],"/",targets[8],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[7],"/",targets[8],"-Log2FC",sep = "")
      rename5 <- paste(targets[7],"/",targets[8],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[7],"/",targets[8],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #GI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g7Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g7Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[7],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[7],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[7],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[7],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[7],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[7],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    #HI Calculations
    {
      #Fold Change
      fc<- as.double(unlist(g8Avg))/as.double(unlist(g9Avg))
      tmpData$AB_FC <- fc
      
      #Log base 2 of the FC
      log2fc <- log2(fc)
      tmpData$AB_Log2FC <- log2fc
      
      if (length(g1Data) < 3) {
        p_Val <- cbind(rep(NA,nrow(g1Data)))
        
        h <- cbind(rep(NA,nrow(g1Data)))
        
        negpVal <- cbind(rep(NA,nrow(g1Data)))
        
        negfdr <- cbind(rep(NA,nrow(g1Data)))
        
      }
      else
      {
        #P values
        p_Val <- cbind(rep(0.0,nrow(g1Data)))
        for (i in 1:nrow(g1Data)) p_Val[i,1] <- t.test(g8Data[i,],g9Data[i,],paired = F, var.equal = T)$p.value
        
        #Neg Log 10 of pval
        negpVal <- -log10(p_Val)
        
        #fdr Corrected p-value
        h <- p.adjust(c(p_Val),"fdr")
        
        #Neg Log 10 of FDR
        negfdr <- -log10(h)
      }
      
      tmpData$AB_pVal <- p_Val
      tmpData$AB_neg_Log_pVal <- negpVal
      
      tmpData$AB_FDR_pVal <- h
      tmpData$AB_neg_Log_FDR_pVal <- negfdr
      
      rename1 <- paste(targets[8],"/",targets[9],"-FC",sep = "")
      rename2 <- paste(targets[8],"/",targets[9],"-pVal",sep = "")
      rename3 <- paste(targets[8],"/",targets[9],"-FDR_pVal",sep = "")
      rename4 <- paste(targets[8],"/",targets[9],"-Log2FC",sep = "")
      rename5 <- paste(targets[8],"/",targets[9],"-neg_Log_pVal",sep = "")
      rename6 <- paste(targets[8],"/",targets[9],"-neg_Log_FDR_pVal",sep = "")
      
      names(tmpData)[names(tmpData) == "AB_FC"] <- rename1
      names(tmpData)[names(tmpData) == "AB_pVal"] <- rename2
      names(tmpData)[names(tmpData) == "AB_FDR_pVal"] <- rename3
      names(tmpData)[names(tmpData) == "AB_Log2FC"] <- rename4
      names(tmpData)[names(tmpData) == "AB_neg_Log_pVal"] <- rename5
      names(tmpData)[names(tmpData) == "AB_neg_Log_FDR_pVal"] <- rename6
      
    }
    
    return(tmpData)
    
  }
  
  tenGroupingGetCalculations <- function(inputData,targets,normalized){
    #Gets a copy of the current Input Data
    tmpData <- inputData
    
    lowest = getGlobalMin(inputData,targets)
    #Get the names of the grouping
    
    #sets the frame data for the set
    f <- c(rep(0.0,ncol(tmpData)))
    g1Data <- data.frame(getDataColumns(tmpData,targets[1],f,lowest))
    g2Data <- data.frame(getDataColumns(tmpData,targets[2],f,lowest))
    g3Data <- data.frame(getDataColumns(tmpData,targets[3],f,lowest))
    g4Data <- data.frame(getDataColumns(tmpData,targets[4],f,lowest))
    g5Data <- data.frame(getDataColumns(tmpData,targets[5],f,lowest))
    g6Data <- data.frame(getDataColumns(tmpData,targets[6],f,lowest))
    g7Data <- data.frame(getDataColumns(tmpData,targets[7],f,lowest))
    g8Data <- data.frame(getDataColumns(tmpData,targets[8],f,lowest))
    g9Data <- data.frame(getDataColumns(tmpData,targets[9],f,lowest))
    g10Data <- data.frame(getDataColumns(tmpData,targets[10],f,lowest))
    
    g1Avg<- rbind(averageData(g1Data,targets[1],f))
    g2Avg<- rbind(averageData(g2Data,targets[2],f))
    g3Avg<- rbind(averageData(g3Data,targets[3],f))
    g4Avg<- rbind(averageData(g4Data,targets[4],f))
    g5Avg<- rbind(averageData(g5Data,targets[5],f))
    g6Avg<- rbind(averageData(g6Data,targets[6],f))
    g7Avg<- rbind(averageData(g7Data,targets[7],f))
    g8Avg<- rbind(averageData(g8Data,targets[8],f))
    g9Avg<- rbind(averageData(g9Data,targets[9],f))
    g10Avg<- rbind(averageData(g10Data,targets[10],f))
    
    for (i in 1:length(targets)) {
      txt = targets[i]
      rHeader = substring(txt,1, nchar(txt)-1)
      targets[i] <-rHeader
      
    }
    
    p_Val<-cbind(rep(0.0,nrow(g1Data)))
    h<-cbind(rep(0.0,nrow(g1Data)))
    negpVal<-cbind(rep(0.0,nrow(g1Data)))
    negfdr<-cbind(rep(0.0,nrow(g1Data)))

    #Group A/B fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g2Avg))
    tmpData$AB_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AB_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AB_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AB_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AB_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AB_neg_Log_FDR_pVal <- negfdr
    
    #Group A/C fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g3Avg))
    tmpData$AC_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AC_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AC_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AC_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AC_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AC_neg_Log_FDR_pVal <- negfdr
    
    #Group AD fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g4Avg))
    tmpData$AD_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AD_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AD_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AD_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AD_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AD_neg_Log_FDR_pVal <- negfdr
    
    #Group AE fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g5Avg))
    tmpData$AE_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AE_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AE_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AE_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AE_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AE_neg_Log_FDR_pVal <- negfdr
    
    #Group AF fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g6Avg))
    tmpData$AF_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AF_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AF_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AF_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AF_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AF_neg_Log_FDR_pVal <- negfdr
    
    #Group AG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g7Avg))
    tmpData$AG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AG_neg_Log_FDR_pVal <- negfdr
    
    #Group AH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g8Avg))
    tmpData$AH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AH_neg_Log_FDR_pVal <- negfdr
    
    #Group AI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g9Avg))
    tmpData$AI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AI_neg_Log_FDR_pVal <- negfdr
    
    #Group AJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g1Avg))/as.double(unlist(g10Avg))
    tmpData$AJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$AJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$AJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$AJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$AJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$AJ_neg_Log_FDR_pVal <- negfdr
    
    #Group BC fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g3Avg))
    tmpData$BC_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BC_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BC_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BC_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BC_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BC_neg_Log_FDR_pVal <- negfdr
    
    #Group BD fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g4Avg))
    tmpData$BD_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BD_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BD_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BD_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BD_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BD_neg_Log_FDR_pVal <- negfdr
    
    #Group BE fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g5Avg))
    tmpData$BE_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BE_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BE_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BE_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BE_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BE_neg_Log_FDR_pVal <- negfdr
    
    #Group BF fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g6Avg))
    tmpData$BF_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BF_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BF_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BF_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BF_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BF_neg_Log_FDR_pVal <- negfdr
    
    #Group BG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g7Avg))
    tmpData$BG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BG_neg_Log_FDR_pVal <- negfdr
    
    #Group BH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g8Avg))
    tmpData$BH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BH_neg_Log_FDR_pVal <- negfdr
    
    #Group BI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g9Avg))
    tmpData$BI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BI_neg_Log_FDR_pVal <- negfdr
    
    #Group BJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g2Avg))/as.double(unlist(g10Avg))
    tmpData$BJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$BJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$BJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$BJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$BJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$BJ_neg_Log_FDR_pVal <- negfdr
    
    #Group CD fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g4Avg))
    tmpData$CD_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CD_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CD_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CD_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CD_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CD_neg_Log_FDR_pVal <- negfdr
    
    #Group CE fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g5Avg))
    tmpData$CE_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CE_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CE_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CE_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CE_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CE_neg_Log_FDR_pVal <- negfdr
    
    #Group CF fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g6Avg))
    tmpData$CF_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CF_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CF_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CF_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CF_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CF_neg_Log_FDR_pVal <- negfdr
    
    #Group CG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g7Avg))
    tmpData$CG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CG_neg_Log_FDR_pVal <- negfdr
    
    #Group CH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g8Avg))
    tmpData$CH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CH_neg_Log_FDR_pVal <- negfdr
    
    #Group CI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g9Avg))
    tmpData$CI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CI_neg_Log_FDR_pVal <- negfdr
    
    #Group CJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g3Avg))/as.double(unlist(g10Avg))
    tmpData$CJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$CJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$CJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$CJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$CJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$CJ_neg_Log_FDR_pVal <- negfdr
    
    #Group DE fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g5Avg))
    tmpData$DE_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DE_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DE_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DE_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DE_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DE_neg_Log_FDR_pVal <- negfdr
    
    #Group DF fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g6Avg))
    tmpData$DF_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DF_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DF_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DF_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DF_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DF_neg_Log_FDR_pVal <- negfdr
    
    #Group DG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g7Avg))
    tmpData$DG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DG_neg_Log_FDR_pVal <- negfdr
    
    #Group DH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g8Avg))
    tmpData$DH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DH_neg_Log_FDR_pVal <- negfdr
    
    #Group DI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g9Avg))
    tmpData$DH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DI_neg_Log_FDR_pVal <- negfdr
    
    #Group DJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g4Avg))/as.double(unlist(g10Avg))
    tmpData$DJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$DJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$DJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$DJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$DJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$DJ_neg_Log_FDR_pVal <- negfdr
    
    #Group EF fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g5Avg))/as.double(unlist(g6Avg))
    tmpData$EF_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$EF_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$EF_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$EF_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$EF_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$EF_neg_Log_FDR_pVal <- negfdr
    
    #Group EG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g5Avg))/as.double(unlist(g7Avg))
    tmpData$EG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$EG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$EG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$EG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$EG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$EG_neg_Log_FDR_pVal <- negfdr
    
    #Group EH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g5Avg))/as.double(unlist(g8Avg))
    tmpData$EH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$EH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$EH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$EH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$EH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$EH_neg_Log_FDR_pVal <- negfdr
    
    #Group EI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g5Avg))/as.double(unlist(g9Avg))
    tmpData$EI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$EI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$EI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$EI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$EI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$EI_neg_Log_FDR_pVal <- negfdr
    
    #Group EJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g5Avg))/as.double(unlist(g10Avg))
    tmpData$EJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$EJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$EJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$EJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$EJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$EJ_neg_Log_FDR_pVal <- negfdr
    
    #Group FG fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g6Avg))/as.double(unlist(g7Avg))
    tmpData$FG_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$FG_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$FG_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$FG_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$FG_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$FG_neg_Log_FDR_pVal <- negfdr
    
    #Group FH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g6Avg))/as.double(unlist(g8Avg))
    tmpData$FH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$FH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$FH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$FH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$FH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$FH_neg_Log_FDR_pVal <- negfdr
    
    #Group FI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g6Avg))/as.double(unlist(g9Avg))
    tmpData$FI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$FI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$FI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$FI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$FI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$FI_neg_Log_FDR_pVal <- negfdr
    
    #Group FJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g6Avg))/as.double(unlist(g10Avg))
    tmpData$FJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$FJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$FJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$FJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$FJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$FJ_neg_Log_FDR_pVal <- negfdr
    
    #Group GH fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g7Avg))/as.double(unlist(g8Avg))
    tmpData$GH_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$GH_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$GH_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$GH_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$GH_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$GH_neg_Log_FDR_pVal <- negfdr
    
    #Group GI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g7Avg))/as.double(unlist(g9Avg))
    tmpData$GI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$GI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$GI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$GI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$GI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$GI_neg_Log_FDR_pVal <- negfdr
    
    #Group GJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g7Avg))/as.double(unlist(g10Avg))
    tmpData$GJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$GJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$GJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$GJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$GJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$GJ_neg_Log_FDR_pVal <- negfdr
    
    #Group HI fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g8Avg))/as.double(unlist(g9Avg))
    tmpData$HI_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$HI_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$HI_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$HI_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$HI_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$HI_neg_Log_FDR_pVal <- negfdr
    
    #Group HJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g8Avg))/as.double(unlist(g10Avg))
    tmpData$HJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$HJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$HJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$HJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$HJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$HJ_neg_Log_FDR_pVal <- negfdr
    
    #Group IJ fold Change, P-val, fdr Corrected P Values,Log base 2 of the FC, Neg Log 10 of pval,Neg Log 10 of FDR
    #Fold Change
    fc<- as.double(unlist(g9Avg))/as.double(unlist(g10Avg))
    tmpData$IJ_FC <- fc
    
    #P values
    df <- as.integer(length(fc)-1)
    p_Val <- pt(fc,df=df,lower.tail = F)
    tmpData$IJ_pVal <- p_Val
    
    #fdr Corrected p-value
    h <- p.adjust(p_Val,"fdr")
    tmpData$IJ_FDR_pVal <- h
    
    #Log base 2 of the FC
    log2fc <- log2(fc)
    tmpData$IJ_Log2FC <- log2fc
    
    #Neg Log 10 of pval
    negpVal <- -log10(p_Val)
    tmpData$IJ_neg_Log_pVal <- negpVal
    
    #Neg Log 10 of FDR
    negfdr <- -log10(h)
    tmpData$IJ_neg_Log_FDR_pVal <- negfdr
    
    return(tmpData)
    
  }
  
}

calculateTukeyPHA<- function(anovaModel){
  output <- TukeyHSD(anovaModel, conf.level =  .95)
  #FishersLSD add into this
  
  return(output)
  
}

meanCenterData <- function(inputData){
  center = rowMeans(inputData)
  
  inputData - rep(center,rep.int(nrow(inputData),ncol(inputData)))
  
}

getSetsData <- function(inputData,samples){
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  count <- 0
  
  for (i in 1:length(numOfSamples)) 
  {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        count <- count+1
        #v1 <- log(inputData[i])
        #print(v1)
      }
      
    }
    
  }
  
  #mean the rows and subtract all values in the row by the mean
  
  oData <- data.frame(matrix(0,nrow(data),count))
  
  ind = 1
  for (i in 1:length(numOfSamples)) 
  {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        v1 <- inputData[i]
        oData[ind] <- v1
        ind <- ind + 1
      }
      
    }
    
  }
  
  #Go through oData and see if its a 0 or a string, and change it to the 1/5 min of that row
  
  
  return(oData)
  
}

getOutputFile <- function(calData, targets){
  outputData <- rbind(c(rep(NA,ncol(calData))),calData)
  
  numOfSamples <- colnames(calData)
  header <- colnames(calData)
  
  #Loop through container that holds the name here and do this x times
  for(i in targets){
    for (j in 1:length(numOfSamples)) 
    {
      if (grepl(i,header[j],fixed = T)) {
        txt = i
        rHeader = substring(txt,1, nchar(txt)-1)
        
        outputData[1,j] <- rHeader
      }
      
    }
  }
  
  return(outputData)
  
}

getGlobalMin <- function(inputData,samples){
  lowest<-10000000000;
  
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  
  for (i in 1:length(header)) {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        t<-data.frame(inputData[i])
        if ((lowest > min(t[t>0]))) {
          lowest <- min(t[t>0])
        }
      }
    }
  }
  
  return(lowest)
  
}

determineGroupingCalculation<- function(inputData,outdir){
  
  samples <- getGroupNames(inputData)
  
  if (length(samples) == 2) {
    print(samples)
    g <- twoGroupingGetCalculation(inputData,samples,FALSE)
    output<-getOutputFile(g,samples)
    computePCA(inputData,outdir,samples)
    
  }
  else if(length(samples) == 3){
    print(samples)
    output<-getOutputFile(threeGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 4){
    print(samples)
    output<-getOutputFile(fourGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 5){
    print(samples)
    output<-getOutputFile(fiveGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 6){
    print(samples)
    output<-getOutputFile(sixGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 7){
    print(samples)
    output<-getOutputFile(sevenGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 8){
    print(samples)
    output<-getOutputFile(eightGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 9){
    print(samples)
    output<-getOutputFile(nineGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  else if(length(samples) == 10){
    print(samples)
    output<-getOutputFile(tenGroupingGetCalculations(inputData,samples,FALSE),samples)
    computePCA(inputData,outdir,samples)
    #computeTukey(inputData,outdir)
    
  }
  
  if (length(samples) < 2 || length(samples) > 10) {
    print("Warning, no sample groupings were detected. Use sample in the sample names for statistics and have the groupings written before the word sample; for example Control_Sample or ControlSample, would group that file as a control")
    
  }
  else if (length(samples) > 10) {
    print("Warning, sample size is either greater than 10, currently we only support up to ten groupings, go to Stats.R in Libraries/Scripts/")
    
  }
  else{
    return(output)
    
  }
  
  
}

determineManualGroupingCalculation<- function(inputData,outdir,groupCSV){
  samples <- manualGrouping(groupCSV)
  lowest = getGlobalMin(inputData,samples)
  
  numOfSamples <- colnames(inputData)
  header <- colnames(inputData)
  
  tmp <- data.frame(inputData)
  
  for (i in 1:length(header)) {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],header[i],fixed = T)) {
        for (k in 1:nrow(tmp[i])) {
          if (tmp[k,i] == 0 || is.null(tmp[k,i])) {
            tmp[k,i] <- runif(1,min = (lowest/8),max = (lowest/2))
          }
        }
      }
    }
  }
  
  print(samples)
  print("Replacing 0's and Null values with a random number between 1/8th and 1/2th the minimum value of the current dataset")
  
  lowerMin<-paste("1/8th the minimum value: ", (lowest/8))
  UpperMin<-paste("1/2th the minimum value: ", (lowest/2))
  print(lowerMin)
  print(UpperMin)
  
  if (length(samples) == 2) {
    g <- twoGroupingGetCalculation(tmp,samples,FALSE)
    output<-getOutputFile(g,samples)
    
  }
  else if(length(samples) == 3){
    output<-getOutputFile(threeGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 4){
    output<-getOutputFile(fourGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 5){
    output<-getOutputFile(fiveGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 6){
    output<-getOutputFile(sixGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 7){
    output<-getOutputFile(sevenGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 8){
    output<-getOutputFile(eightGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 9){
    output<-getOutputFile(nineGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  else if(length(samples) == 10){
    output<-getOutputFile(tenGroupingGetCalculations(tmp,samples,FALSE),samples)
    
  }
  
  if (length(samples) < 2 || length(samples) > 10) {
    print("Warning, no sample groupings were detected. Use sample in the sample names for statistics and have the groupings written before the word sample; for example Control_Sample or ControlSample, would group that file as a control")
    
  }
  else if (length(samples) > 10) {
    print("Warning, sample size is either greater than 10, currently we only support up to ten groupings, go to Stats.R in Libraries/Scripts/")
    
  }
  else{
    return(output)
    
  }
  
  
}

determineNormalizedCalculations<- function(inputData){
  samples <- getGroupNames(inputData)
  
  s <- getGroupNames(inputData)
  dat<- getSetsData(inputData,s)
  preproc <- qs::qread("preproc.qs")
  dat<- ReplaceMissingByLoD(preproc)
  dat<- logNormTransformation(dat)
  dat<- meanCenterTransformation(dat)
  
  if (length(samples) == 2) {
    print(samples)
    output<-twoGroupingGetCalculation(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 3){
    print(samples)
    output<-threeGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 4){
    print(samples)
    output<-fourGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 5){
    print(samples)
    output<-fiveGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 6){
    print(samples)
    output<-sixGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 7){
    print(samples)
    output<-sevenGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 8){
    print(samples)
    output<-eightGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 9){
    print(samples)
    output<-nineGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  else if(length(samples) == 10){
    print(samples)
    output<-tenGroupingGetCalculations(dat,s,TRUE)
    computePCA(inputData)
    computeTukey(inputData)
    
  }
  
  return(output)
  
}

determineAnovaCalculation<- function(inputData){
  samples <- getGroupNames(inputData)

  #Gets a copy of the current Input Data
  tmpData <- inputData
  
  if (length(samples) == 3) {
    
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    
    calculateThreeGroupAnova(g1Avg,g2Avg,g3Avg, samples)
    
  }
  else if(length(samples) == 4){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    
    calculateFourGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg, samples)
    
  }
  else if(length(samples) == 5){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    
    calculateFiveGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg, samples)
    
    
  }
  else if(length(samples) == 6){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    g6Avg<- rbind(averageData(tmpData,samples[6],f))
    
    calculateSixGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg,g6Avg, samples)
    
  }
  else if(length(samples) == 7){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    g6Avg<- rbind(averageData(tmpData,samples[6],f))
    g7Avg<- rbind(averageData(tmpData,samples[7],f))
    
    calculateSevenGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg,g6Avg,g7Avg, samples)
    
  }
  else if(length(samples) == 8){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    g6Avg<- rbind(averageData(tmpData,samples[6],f))
    g7Avg<- rbind(averageData(tmpData,samples[7],f))
    g8Avg<- rbind(averageData(tmpData,samples[8],f))
    
    calculateEightGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg,g6Avg,g7Avg,g8Avg, samples)
    
    
  }
  else if(length(samples) == 9){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    g6Avg<- rbind(averageData(tmpData,samples[6],f))
    g7Avg<- rbind(averageData(tmpData,samples[7],f))
    g8Avg<- rbind(averageData(tmpData,samples[8],f))
    g9Avg<- rbind(averageData(tmpData,samples[9],f))
    
    calculateNineGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg,g6Avg,g7Avg,g8Avg,g9Avg, samples)
    
    
  }
  else if(length(samples) == 10){
    f <- c(rep(0.0,nrow(tmpData)))
    g1Avg<- rbind(averageData(tmpData,samples[1],f))
    g2Avg<- rbind(averageData(tmpData,samples[2],f))
    g3Avg<- rbind(averageData(tmpData,samples[3],f))
    g4Avg<- rbind(averageData(tmpData,samples[4],f))
    g5Avg<- rbind(averageData(tmpData,samples[5],f))
    g6Avg<- rbind(averageData(tmpData,samples[6],f))
    g7Avg<- rbind(averageData(tmpData,samples[7],f))
    g8Avg<- rbind(averageData(tmpData,samples[8],f))
    g9Avg<- rbind(averageData(tmpData,samples[9],f))
    g10Avg<- rbind(averageData(tmpData,samples[10],f))
    
    calculateTenGroupAnova(g1Avg,g2Avg,g3Avg,g4Avg,g5Avg,g6Avg,g7Avg,g8Avg,g9Avg,g10Avg, samples)
    
    
  }
  else{
    #Throw Error here
  }

}

#PCA Calculations
{
  LogNorm<-function(x, minVal){
    log10((x + sqrt(x^2 + minVal^2))/2)
  }
  
  MeanCenter<-function(x){
    x - mean(x);
  }
  
  logNormTransformation<- function(inputData){
    minVal <- min(abs(inputData[inputData!=0]))/10
    inputData<-apply(inputData,2,LogNorm,minVal)
    
  }
  
  meanCenterTransformation<- function(inputData){
    inputData<-apply(inputData,2,MeanCenter)
    
  }
  
  ReplaceMissingByLoD <- function(inputData){
    inputData <- as.matrix(inputData)
    
    rowNms <- rownames(inputData)
    colNms <- colnames(inputData)
    inputData <- apply(inputData, 2, replaceByLoD)
    rownames(inputData) <- rowNms
    colnames(inputData) <- colNms
    return (inputData)
  }
  
  replaceByLoD <- function(x){
    lod <- min(x[x>0], na.rm=T)/5;
    x[x==0|is.na(x)] <- lod;
    return(x);
  }
  
}

#PCA Handling
computePCA <- function(inputData,outdir,samples,type){
  
  mset<-InitDataObjects("conc","stat",F)
  mset<-Read.TextData(mset,inputData,"colu","disc")
  mset<-SanityCheckData(mset)
  mset<-ReplaceMin(mset)
  mset<-SanityCheckData(mset)
  mset<-FilterVariable(mset,"none","F",25)
  mset<-PreparePrenormData(mset)
  mset<-Normalization(mset,"none","F",25)
  #mset<-PCA.Anal(mset)
  
  pca <- prcomp(mset$dataSet$norm, center=TRUE, scale=F);
  
  # obtain variance explained
  sum.pca <- summary(pca);
  imp.pca <- sum.pca$importance;
  std.pca <- imp.pca[1,]; # standard devietation
  var.pca <- imp.pca[2,]; # variance explained by each PC
  cum.pca <- imp.pca[3,];

  # PCA_Sum <- summary(mset);
  # PCA_Imp <- PCA_Sum$importance;
  # PCA_Std <- PCA_Imp[1,]; # standard devietation
  # PCA_Var <- PCA_Imp[2,]; # variance explained by each PC
  # PCA_Cum <- PCA_Imp[3,]; # cummulated variance explained
  # 
  scores<-paste("/PcaScores_",tolower(type),".csv",sep="")
  loading<-paste("/PcaLoadings_",tolower(type),".csv",sep="")
  
  mset$pca<-append(pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca))
  scoreOut<-paste(outdir,scores,sep="")
  write.csv(signif(mset$pca$x,5), file=scoreOut);
  loadingOut<-paste(outdir,loading,sep="")
  write.csv(signif(mset$pca$rotation,5), file=loadingOut)
  
  #Append from files Heres
  #Scores
  #We need to look at the col[1] row[n] and see if its sample name
  inputData<- read.csv(inputData, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  data <- read.csv(scoreOut)
  
  for (i in 1:length(samples)) {
    txt = samples[i]
    rHeader = substring(txt,1, nchar(txt)-1)
    samples[i] <-rHeader
    
  }
  
  lOut <- cbind(c(rep(NA,nrow(data))),data)
  
  for (i in 1:nrow(lOut)) {
    for (j in 1:length(samples)) {
      if (grepl(samples[j],lOut[i,2],fixed = T)) {
        lOut[i,1]<-samples[j]
      }
      
    }
    
  }
  
  colnames(lOut)[1]<-"Grouping"
  colnames(lOut)[2]<-"Sample_Names"
  write.csv(lOut,scoreOut,row.names = FALSE)
  
  #Loadings
  #append 6 prev columns into the file
  data1 <- read.csv(loadingOut)
  df <- data1[,-1]
  #Grabs Appropiate Columns
  Scores <- inputData$Score[2:length(inputData$Score)]
  SeriesType_Identifer <- inputData$SeriesType_Identifier[2:length(inputData$Score)]
  Name_or_Class <- inputData$Name_or_Class[2:length(inputData$Score)]
  Formula<- inputData$Formula[2:length(inputData$Score)]
  SMILES<-inputData$SMILES[2:length(inputData$Score)]
  row_ID <- inputData$row.ID[2:length(inputData$Score)]
  
  #Appending and Writing
  #sOut<- cbind(row_ID,Scores,SeriesType_Identifer,Name_or_Class,Formula,SMILES,df)
  #write.csv(sOut,loadingOut, row.names = FALSE)
}

#MetaboAnalyst Styled Calculations
{
  parseTukey <- function(tukey, cut.off){
    inx <- tukey$cls[,"p adj"] <= cut.off;
    paste(rownames(tukey$cls)[inx], collapse="; ");
  }
  
  parseFisher <- function(fisher, cut.off){
    inx <- fisher[,"pvalue"] <= cut.off;
    paste(rownames(fisher)[inx], collapse="; ");
  }
  
  FisherLSD <- function(aov.obj, thresh){
    LSD.test(aov.obj,"cls", alpha=thresh)
  }
  
  Ttests.Anal <- function(mSetObj=NA, nonpar=F, threshp=0.05, paired=FALSE, equal.var=TRUE, pvalType="fdr", all_results=FALSE){
    
    res <- GetTtestRes(mSetObj, paired, equal.var, nonpar);
    
    t.stat <- res[,1];
    p.value <- res[,2];
    names(t.stat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
    
    p.log <- -log10(p.value);
    fdr.p <- p.adjust(p.value, "fdr");
    
    if(all_results==TRUE){
      
      all.mat <- data.frame(signif(t.stat,5), signif(p.value,5), signif(p.log,5), signif(fdr.p,5));
      
      if(nonpar){
        tt.nm = "Wilcoxon Rank Test";  
        file.nm <- "wilcox_rank_all.csv"
        colnames(all.mat) <- c("V", "p.value", "-log10(p)", "FDR");
      }else{
        tt.nm = "T-Tests";
        file.nm <- "t_test_all.csv";
        colnames(all.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
      }
      write.csv(all.mat, file=file.nm);
    }
    
    if(pvalType=="fdr"){
      inx.imp <- fdr.p <= threshp;
    }else{
      inx.imp <- p.value <= threshp;
    }
    
    sig.num <- sum(inx.imp, na.rm = TRUE);
    #print(paste("A total of", sig.num, "significant features were found."));
    
    if(sig.num > 0){
      sig.t <- t.stat[inx.imp];
      sig.p <- p.value[inx.imp];
      lod<- -log10(sig.p);
      sig.q <-fdr.p[inx.imp];
      
      sig.mat <- cbind(sig.t, sig.p, lod, sig.q);
      colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
      ord.inx <- order(sig.p);
      sig.mat <- sig.mat[ord.inx,,drop=F];
      sig.mat <- signif(sig.mat, 5);
      
      if(nonpar){
        tt.nm = "Wilcoxon Rank Test";  
        file.nm <- "wilcox_rank.csv"
        colnames(sig.mat) <- c("V", "p.value", "-log10(p)", "FDR");
      }else{
        tt.nm = "T-Tests";
        file.nm <- "t_test.csv";
        colnames(sig.mat) <- c("t.stat", "p.value", "-log10(p)", "FDR");
      }
      write.csv(sig.mat, file=file.nm);
      
      tt <- list (
        tt.nm = tt.nm,
        sig.nm = file.nm,
        sig.num = sig.num,
        paired = paired,
        pval.type = pvalType,
        raw.thresh = threshp,
        t.score = t.stat,
        p.value = p.value,
        p.log = p.log,
        thresh = -log10(threshp), # only used for plot threshold line
        inx.imp = inx.imp,
        sig.mat = sig.mat
      );
    }else{
      tt <- list (
        sig.num = sig.num,
        paired = paired,
        pval.type = pvalType,
        raw.thresh = threshp,
        t.score = t.stat,
        p.value = p.value,
        p.log = p.log,
        thresh = -log10(threshp), # only used for plot threshold line
        inx.imp = inx.imp
      );
    }
    
    mSetObj$analSet$tt <- tt;
    
    
    return(mSetObj);
  }
  
  fast.write.csv <- function(dat, file, row.names=TRUE){
    tryCatch(
      {
        if(is.data.frame(dat)){
          # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
          data.table::fwrite(dat, file, row.names=row.names);
        }else{
          write.csv(dat, file, row.names=row.names);  
        }
      }, error=function(e){
        #print(e);
        write.csv(dat, file, row.names=row.names);   
      }, warning=function(w){
        #print(w);
        write.csv(dat, file, row.names=row.names); 
      });
  }
  
  Volcano.Anal <- function(mSetObj=NA, paired=FALSE, fcthresh, cmpType, nonpar=F, threshp, equal.var=TRUE, pval.type="raw"){
    
    # Note, volcano is based on t-tests and fold change analysis
    #### t-tests and p values. If performed and identical parameters
    mSetObj <- Ttests.Anal(mSetObj, nonpar, threshp, paired, equal.var, pval.type);
    p.value <- mSetObj$analSet$tt$p.value;
    
    if(pval.type == "fdr"){
      p.value <- p.adjust(p.value, "fdr");
    }   
    
    inx.p <- p.value <= threshp;
    p.log <- -log10(p.value);
    
    #### fc analysis
    mSetObj <- FC.Anal(mSetObj, fcthresh, cmpType, paired);
    
    fcthresh = ifelse(fcthresh>1, fcthresh, 1/fcthresh);
    max.xthresh <- log2(fcthresh);
    min.xthresh <- log2(1/fcthresh);
    
    fc.log <- mSetObj$analSet$fc$fc.log;
    fc.all <- mSetObj$analSet$fc$fc.all;
    
    inx.up <- mSetObj$analSet$fc$inx.up;
    inx.down <- mSetObj$analSet$fc$inx.down;
    
    # subset inx.p to inx.up/down
    keep.inx <- names(inx.p) %in% names(inx.up)
    inx.p <- inx.p[keep.inx]
    p.log <- p.log[keep.inx]
    
    # create named sig table for display
    inx.imp <- (inx.up | inx.down) & inx.p;
    sig.var <- cbind(fc.all[inx.imp,drop=F], fc.log[inx.imp,drop=F], p.value[inx.imp,drop=F], p.log[inx.imp,drop=F]);
    
    if(pval.type == "fdr"){
      colnames(sig.var) <- c("FC", "log2(FC)", "p.ajusted", "-log10(p)");
    }else{
      colnames(sig.var) <- c("FC", "log2(FC)", "raw.pval", "-log10(p)");
    }
    
    # first order by log(p), then by log(FC)
    ord.inx <- order(sig.var[,4], abs(sig.var[,2]), decreasing=T);
    sig.var <- sig.var[ord.inx,,drop=F];
    sig.var <- signif(sig.var,5);
    
    fileName <- "volcano.csv";
    write.csv(signif(sig.var,5), file=fileName);
    
    volcano <- list (
      pval.type = pval.type,
      raw.threshx = fcthresh,
      raw.threshy = threshp,
      paired = paired,
      max.xthresh = max.xthresh,
      min.xthresh = min.xthresh,
      thresh.y = -log10(threshp),
      fc.all = fc.all,
      fc.log = fc.log,
      inx.up = inx.up,
      inx.down = inx.down,
      p.log = p.log,
      inx.p = inx.p,
      sig.mat = sig.var
    );
    
    mSetObj$analSet$volcano <- volcano;
    return(mSetObj);
  }
  
  ANOVA.Anal<-function(mSetObj=NA, nonpar=FALSE, thresh=0.05, post.hoc="fisher", all_results=FALSE,outDir,type) {
    
    sig.num <- 0;
    if(nonpar){
      aov.nm <- "Kruskal Wallis Test";
      anova.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, kwtest, cls=mSetObj$dataSet$cls);
      
      #extract all p values
      res <- unlist(lapply(anova.res, function(x) {c(x$statistic, x$p.value)}));
      res <- data.frame(matrix(res, nrow=length(anova.res), byrow=T), stringsAsFactors=FALSE);
      
      fstat <- res[,1];
      p.value <- res[,2];
      
      names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
      fdr.p <- p.adjust(p.value, "fdr");
      
      #inx.imp <- p.value <= thresh;
      inx.imp <- fdr.p <= thresh;
      sig.num <- sum(inx.imp);
      
      if(sig.num > 0){ 
        sig.f <- fstat[inx.imp];
        sig.p <- p.value[inx.imp];
        fdr.p <- fdr.p[inx.imp];
        
        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), 'NA');
        rownames(sig.mat) <- names(sig.p);
        colnames(sig.mat) <- c("chi.squared", "p.value", "-log10(p)", "FDR", "Post-Hoc");
        
        # order the result simultaneously
        ord.inx <- order(sig.p, decreasing = FALSE);
        sig.mat <- sig.mat[ord.inx,,drop=F];
        
        fileName <- "kw_posthoc.csv";
        my.mat <- sig.mat[,1:4];
        colnames(my.mat) <- c("chi_squared", "pval_KW", "-log10(p)", "FDR");
      }
    }
    else{
      aov.nm <- "One-way ANOVA";
      
      aov.res <- apply(as.matrix(mSetObj$dataSet$norm), 2, aof, cls=mSetObj$dataSet$cls);
      anova.res <- lapply(aov.res, anova);
      
      #extract all p values
      res <- unlist(lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}));
      res <- data.frame(matrix(res, nrow=length(aov.res), byrow=T), stringsAsFactors=FALSE);
      
      
      fstat <- res[,1];
      p.value <- res[,2];
      names(fstat) <- names(p.value) <- colnames(mSetObj$dataSet$norm);
      
      fdr.p <- p.adjust(p.value, "fdr");
      
      if(all_results==TRUE){
        all.mat <- data.frame(signif(p.value,5), signif(-log10(p.value),5), signif(fdr.p,5));
        rownames(all.mat) <- names(p.value);
        colnames(all.mat) <- c("p_value", "neg_log10(p)", "FDR");
        write.csv(all.mat, paste(OutputDirectory,"/anova_all_results.csv",sep=""))
      }
      
      # do post-hoc only for signficant entries
      # inx.imp <- p.value <= thresh;
      inx.imp <- fdr.p <= thresh;
      sig.num <- sum(inx.imp);
      if(sig.num > 0){
        # note aov obj is not avaible using fast version
        # need to recompute using slower version for the sig ones
        aov.imp <- aov.res[inx.imp];
        sig.f <- fstat[inx.imp];
        sig.p <- p.value[inx.imp];
        fdr.p <- fdr.p[inx.imp];
        cmp.res <- NULL;
        post.nm <- NULL;
        if(post.hoc=="tukey"){
          tukey.res<-lapply(aov.imp, TukeyHSD, conf.level=1-thresh);
          cmp.res <- unlist(lapply(tukey.res, parseTukey, cut.off=thresh));
          post.nm = "Tukey's HSD";
        }else{
          fisher.res<-lapply(aov.imp, FisherLSD, thresh);
          cmp.res <- unlist(lapply(fisher.res, parseFisher, cut.off=thresh));
          post.nm = "Fisher's LSD";
        }
        
        # create the result dataframe,
        # note, the last column is string, not double
        
        sig.mat <- data.frame(signif(sig.f,5), signif(sig.p,5), signif(-log10(sig.p),5), signif(fdr.p,5), cmp.res);
        rownames(sig.mat) <- names(sig.p);
        colnames(sig.mat) <- c("f_value", "p_value", "neg_log10(p)", "FDR", post.nm);
        
        # order the result simultaneously
        ord.inx <- order(sig.p, decreasing = FALSE);
        sig.mat <- sig.mat[ord.inx,,drop=F];
        fileName <- "anova_posthoc.csv";
      }
    }
    
    #print(paste(c("A total of", sum(inx.imp), "significant features were found."), collapse=" "));
    posthoc <- paste(OutputDirectory,"/anova_posthoc_",tolower(type),".csv",sep="")
    
    if(sig.num> 0){
      res <- 1;
      write.csv(sig.mat,file = posthoc)
      #fast.write.csv(sig.mat,file=fileName);
      aov<-list (
        aov.nm = aov.nm,
        sig.num = sig.num,
        sig.nm = fileName,
        raw.thresh = thresh,
        thresh = -log10(thresh), # only used for plot threshold line
        p.value = p.value,
        p.log = -log10(p.value),
        inx.imp = inx.imp,
        post.hoc = post.hoc,
        sig.mat = sig.mat
      );
    }
    else{
      res <- 0;
      aov<-list (
        aov.nm = aov.nm,
        sig.num = sig.num,
        raw.thresh = thresh,
        thresh = -log10(thresh), # only used for plot threshold line
        p.value = p.value,
        p.log = -log10(p.value),
        inx.imp = inx.imp
      );
    }
    mSetObj$analSet$aov <- aov;
    return(mSetObj);
  }
  
  PCA.Anal <- function(mSetObj=NA){
    pca <- prcomp(mSetObj$dataSet$norm, center=TRUE, scale=F);
    
    # obtain variance explained
    sum.pca <- summary(pca);
    imp.pca <- sum.pca$importance;
    std.pca <- imp.pca[1,]; # standard devietation
    var.pca <- imp.pca[2,]; # variance explained by each PC
    cum.pca <- imp.pca[3,]; # cummulated variance explained
    
    # store the item to the pca object
    mSetObj$analSet$pca<-append(pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
    fast.write.csv(signif(mSetObj$analSet$pca$x,5), file="pca_score.csv");
    fast.write.csv(signif(mSetObj$analSet$pca$rotation,5), file="pca_loadings.csv");
    mSetObj$analSet$pca$loading.type <- "all";
    mSetObj$custom.cmpds <- c();
    return(mSetObj);
  }
  
  Normalization <- function(mSetObj=NA, rowNorm, transNorm, scaleNorm, ref=NULL, ratio=FALSE, ratioNum=20){
    
    # PreparePrenormData() called already
    data <- qs::qread("prenorm.qs");
    
    cls <- mSetObj$dataSet$prenorm.cls;
    
    # note, setup time factor
    if(substring(mSetObj$dataSet$format,4,5)=="ts"){
      if(is.null(mSetObj$dataSet$prenorm.facA)){
        nfacA <- mSetObj$dataSet$facA;
        nfacB <- mSetObj$dataSet$facB;
      }else{
        nfacA <- mSetObj$dataSet$prenorm.facA;
        nfacB <- mSetObj$dataSet$prenorm.facB;
      }
      
      mSetObj$dataSet$facA <- nfacA;
      mSetObj$dataSet$facB <- nfacB;
      if(mSetObj$dataSet$design.type =="time" | mSetObj$dataSet$design.type =="time0"){
        # determine time factor and should order first by subject then by each time points
        if(tolower(mSetObj$dataSet$facA.lbl) == "time"){ 
          time.fac <- nfacA;
          exp.fac <- nfacB;
        }
        else{
          time.fac <- nfacB;
          exp.fac <- nfacA;
        }
        # now make sure time fac is ordered
        lvls <- levels(time.fac);
        time.points <- as.numeric(as.character(lvls));
        ord.lvls <- lvls[order(time.points)];
        time.fac <- ordered(time.fac, levels = ord.lvls);
        mSetObj$dataSet$time.fac <- time.fac;
        mSetObj$dataSet$exp.fac <- exp.fac;
      }
    }
    
    colNames <- colnames(data);
    rowNames <- rownames(data);
    
    # row-wise normalization
    if(rowNorm=="QuantileNorm"){
      data<-QuantileNormalize(data);
      # this can introduce constant variables if a variable is 
      # at the same rank across all samples (replaced by its average across all)
      
      varCol <- apply(data, 2, var, na.rm=T);
      constCol <- (varCol == 0 | is.na(varCol));
      constNum <- sum(constCol, na.rm=T);
      if(constNum > 0){
        #print(paste("After quantile normalization", constNum, "features with a constant value were found and deleted."));
        data <- data[,!constCol, drop=FALSE];
        colNames <- colnames(data);
        rowNames <- rownames(data);
      }
      rownm<-"Quantile Normalization";
    }
    else if(rowNorm=="GroupPQN"){
      grp.inx <- cls == ref;
      ref.smpl <- apply(data[grp.inx, , drop=FALSE], 2, mean);
      data<-t(apply(data, 1, ProbNorm, ref.smpl));
      rownm<-"Probabilistic Quotient Normalization by a reference group";
    }
    else if(rowNorm=="SamplePQN"){
      ref.smpl <- data[ref, , drop=FALSE];
      data<-t(apply(data, 1, ProbNorm, ref.smpl));
      rownm<-"Probabilistic Quotient Normalization by a reference sample";
    }
    else if(rowNorm=="CompNorm"){
      data<-t(apply(data, 1, CompNorm, ref));
      rownm<-"Normalization by a reference feature";
    }
    else if(rowNorm=="SumNorm"){
      data<-t(apply(data, 1, SumNorm));
      rownm<-"Normalization to constant sum";
    }
    else if(rowNorm=="MedianNorm"){
      data<-t(apply(data, 1, MedianNorm));
      rownm<-"Normalization to sample median";
    }
    else if(rowNorm=="SpecNorm"){
      if(!exists("norm.vec")){
        norm.vec <- rep(1,nrow(data)); # default all same weight vec to prevent error
        #print("No sample specific information were given, all set to 1.0");
      }
      rownm<-"Normalization by sample-specific factor";
      data<-data/norm.vec;
    }
    else{
      # nothing to do
      rownm<-"N/A";
    }
    
    # use apply will lose dimension info (i.e. row names and colnames)
    rownames(data)<-rowNames;
    colnames(data)<-colNames;
    
    # if the reference by feature, the feature column should be removed, since it is all 1
    if(rowNorm=="CompNorm" && !is.null(ref)){
      inx<-match(ref, colnames(data));
      data<-data[,-inx, drop=FALSE];
      colNames <- colNames[-inx];
    }
    
    # record row-normed data for fold change analysis (b/c not applicable for mean-centered data)
    row.norm <- as.data.frame(CleanData(data, T, T)); #moved below ratio 
    qs::qsave(row.norm, file="row_norm.qs");
    # this is for biomarker analysis only (for compound concentration data)
    if(ratio){
      min.val <- min(abs(data[data!=0]))/2;
      norm.data <- log2((data + sqrt(data^2 + min.val))/2);
      transnm<-"Log2 Normalization";
      ratio.mat <- CalculatePairwiseDiff(norm.data);
      
      fstats <- Get.Fstat(ratio.mat, cls);
      hit.inx <- rank(-fstats) < ratioNum;  # get top n
      
      ratio.mat <- ratio.mat[, hit.inx, drop=FALSE];
      
      data <- cbind(norm.data, ratio.mat);
      
      colNames <- colnames(data);
      rowNames <- rownames(data);
      mSetObj$dataSet$use.ratio <- TRUE;
      mSetObj$dataSet$proc.ratio <- data;
      
    }else{
      mSetObj$dataSet$use.ratio <- FALSE;
      # transformation
      # may not be able to deal with 0 or negative values
      if(transNorm=='LogNorm'){
        min.val <- min(abs(data[data!=0]))/10;
        data<-apply(data, 2, LogNorm, min.val);
        transnm<-"Log10 Normalization";
      }else if(transNorm=='SrNorm'){
        min.val <- min(abs(data[data!=0]))/10;
        data<-apply(data, 2, SquareRootNorm, min.val);
        transnm<-"Square Root Transformation";
      }else if(transNorm=='CrNorm'){
        norm.data <- abs(data)^(1/3);
        norm.data[data<0] <- - norm.data[data<0];
        data <- norm.data;
        transnm<-"Cubic Root Transformation";
      }else{
        transnm<-"N/A";
      }
    }
    
    # scaling
    if(scaleNorm=='MeanCenter'){
      data<-apply(data, 2, MeanCenter);
      scalenm<-"Mean Centering";
    }else if(scaleNorm=='AutoNorm'){
      data<-apply(data, 2, AutoNorm);
      scalenm<-"Autoscaling";
    }else if(scaleNorm=='ParetoNorm'){
      data<-apply(data, 2, ParetoNorm);
      scalenm<-"Pareto Scaling";
    }else if(scaleNorm=='RangeNorm'){
      data<-apply(data, 2, RangeNorm);
      scalenm<-"Range Scaling";
    }else{
      scalenm<-"N/A";
    }
    
    # note after using "apply" function, all the attribute lost, need to add back
    rownames(data)<-rowNames;
    colnames(data)<-colNames;
    
    # need to do some sanity check, for log there may be Inf values introduced
    data <- CleanData(data, T, F);
    
    if(ratio){
      mSetObj$dataSet$ratio <- CleanData(ratio.mat, T, F)
    }
    
    mSetObj$dataSet$norm <- as.data.frame(data);
    if(substring(mSetObj$dataSet$format,4,5)=="ts"){
      if(rownames(mSetObj$dataSet$norm) != rownames(mSetObj$dataSet$meta.info)){
        #print("Metadata and data norm are not synchronized.")
      }
      mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[rownames(data),]  
    }
    
    qs::qsave(mSetObj$dataSet$norm, file="complete_norm.qs");
    mSetObj$dataSet$cls <- cls;
    
    mSetObj$dataSet$rownorm.method <- rownm;
    mSetObj$dataSet$trans.method <- transnm;
    mSetObj$dataSet$scale.method <- scalenm;
    mSetObj$dataSet$combined.method <- FALSE;
    mSetObj$dataSet$norm.all <- NULL; # this is only for biomarker ROC analysis
    
    return(mSetObj);
  }
  
  FilterVariable <- function(mSetObj=NA, filter, qcFilter, rsd){
    
    #Reset to default
    mSetObj$dataSet$filt <- NULL;
    
    if(is.null(mSetObj$dataSet$proc)){
      int.mat <- as.matrix(qs::qread("data_proc.qs"));
    }else{
      int.mat <- as.matrix(mSetObj$dataSet$proc);
    }
    
    cls <- mSetObj$dataSet$proc.cls;
    
    # save a copy
    mSetObj$dataSet$filt.cls <- cls;
    if(substring(mSetObj$dataSet$format,4,5)=="ts"){
      mSetObj$dataSet$filt.facA <- mSetObj$dataSet$proc.facA; 
      mSetObj$dataSet$filt.facB <- mSetObj$dataSet$proc.facB; 
    }
    
    msg <- "";
    if(qcFilter == "T"){
      rsd <- rsd/100;
      # need to check if QC exists
      qc.hits <- tolower(as.character(cls)) %in% "qc";
      if(sum(qc.hits) > 1){ # require at least 2 QC for RSD
        qc.mat <- int.mat[qc.hits,];
        sds <- apply(qc.mat, 2, sd, na.rm=T);
        mns <- apply(qc.mat, 2, mean, na.rm=T);
        rsd.vals <- abs(sds/mns);  
        gd.inx <- rsd.vals < rsd;
        int.mat <- int.mat[,gd.inx];
        if(mSetObj$analSet$type == "mummichog"){
          msg <- paste("Removed ", sum(!gd.inx), " features based on QC RSD values. QC samples are excluded from downstream functional analysis.");
        }else{
          msg <- paste("Removed ", sum(!gd.inx), " features based on QC RSD values. QC samples are still kept. You can remove them later.");
        }
      }
      else if(sum(qc.hits) > 0){
        #print("RSD requires at least 2 QC samples, and only non-QC based filtering can be applied.");
        return(0);
      }
      else{
        #print("No QC Samples (with class label: QC) found.  Please use non-QC based filtering.");
        return(0);
      }
    }
    filt.res <- PerformFeatureFilter(int.mat, filter, NULL, mSetObj$analSet$type);
    mSetObj$dataSet$filt <- filt.res$data;
    msg <- paste(msg, filt.res$msg);
    print(msg);
    mSetObj$msgSet$filter.msg <- msg;
    return(mSetObj);
  }
  
  SanityCheckData <- function(mSetObj=NA){
    anal.type<-"stat"
    if(file.exists("data_orig.qs")){  
      orig.data <- qs::qread("data_orig.qs");
    } else {
      return(0);
    }  
    msg <- NULL;
    cls <- mSetObj$dataSet$orig.cls;
    mSetObj$dataSet$small.smpl.size <- 0;
    
    # check class info
    if(mSetObj$dataSet$cls.type == "disc"){
      if(substring(mSetObj$dataSet$format,4,5)=="ts"){
        metadata <- mSetObj$dataSet$meta.info
        if(mSetObj$dataSet$design.type =="time"){
          msg<-c(msg, "The data is time-series data.");
        }else if(mSetObj$dataSet$design.type =="time0"){
          msg<-c(msg, "The data is time-series only data.");
        }else{
          msg<-c(msg, "The data is not time-series data.");
        }
        clsA.num <- length(levels(metadata[,1]));
        clsB.num <- length(levels(metadata[,2]));
        if(ncol(metadata) == 2){
          msg<-c(msg, paste(clsA.num, "groups were detected in samples for factor", colnames(metadata)[1]));
          msg<-c(msg, paste(clsB.num, "groups were detected in samples for factor", colnames(metadata)[2]));
        }else{
          msg <- c(msg, paste0(clsA.num, " groups were detected from primary meta-data factor: ", colnames(mSetObj$dataSet$meta.info)[1], "."));
        }
      }else{
        if(mSetObj$dataSet$paired){
          msg<-c(msg,"Samples are paired.");
          # need to first set up pair information if not csv file
          if(!(mSetObj$dataSet$type=="conc" | mSetObj$dataSet$type=="specbin" | mSetObj$dataSet$type=="pktable" )){
            pairs <- ReadPairFile();
            # check if they are of the right length
            if(length(pairs)!=length(mSetObj$dataSet$url.smp.nms)){
              #print("Error: the total paired names are not equal to sample names.");
              return(0);
            }else{
              # matching the names of the files
              inx<-match(rownames(orig.data), names(pairs));
              #check if all matched exactly
              if(sum(is.na(inx))>0){
                #print("Error: some paired names not match the sample names.");
                return(0);
              }else{
                mSetObj$dataSet$pairs <- pairs[inx];
              }
            }
          }
          
          pairs <- mSetObj$dataSet$pairs;
          
          # check if QC samples are present
          qc.hits <- tolower(as.character(cls)) %in% "qc";
          if(sum(qc.hits) > 0){
            #print("<font color='red'>Error: QC samples not supported in paired analysis mode.</font>");
            #print("You can perform QC filtering using regular two-group labels.");
            #print("Then re-upload your data (without QC samples) for paired analysis.");
            return(0);
          }else{
            pairs <- as.numeric(pairs);
          }
          
          label <- as.numeric(pairs);
          cls <- as.factor(ifelse(label>0,1,0));
          mSetObj$dataSet$pairs <- label;
          
          lev <- unique(pairs);
          uni.cl <- length(lev);
          uni.cl.abs <- uni.cl/2;             
          sorted.pairs <- sort(pairs,index=TRUE);
          
          if(!all(sorted.pairs$x==c(-uni.cl.abs:-1,1:uni.cl.abs))){
            #print("There are some problems in paired sample labels! ");
            if(uni.cl.abs != round(uni.cl.abs)){
              duplicates <- pairs[duplicated(pairs)]
              dup.msg <- paste0("Duplicated labels:", duplicates)
              #print(paste("The total samples must be of even number!", dup.msg));
            }else{
              #print(paste("And class labels between ",-uni.cl.abs,
              #" and 1, and between 1 and ",uni.cl.abs,".",sep=""));
            }
            return(0);
          } else {  
            msg <- c(msg,"The labels of paired samples passed sanity check.");
            msg <- c(msg, paste("A total of", uni.cl.abs, "pairs were detected."));
            # make sure paired samples are sorted 1:n/2 and -1:-n/2
            
            x<-sorted.pairs$ix[(uni.cl.abs+1):uni.cl]
            y<-sorted.pairs$ix[uni.cl.abs:1]
            index<-as.vector(cbind(x,y));
            cls<-cls[index];
            pairs <- pairs[index];
            mSetObj$dataSet$pairs <- pairs;
            mSetObj$dataSet$orig.cls <- cls;
            orig.data<- orig.data[index,];
            qs::qsave(orig.data, file="data_orig.qs");
          }
        } else {
          
          # check for class labels at least two replicates per class but QC and BLANK
          
          cls.lbl <- mSetObj$dataSet$orig.cls;
          qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
          if(sum(qb.inx) > 0){
            cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
          } else {
            cls.Clean <- cls.lbl;
          }
          
          # allow it pass to sanity check and correct there
          if(anal.type != "network"){ # add exception for DSPC correlation network 
            if(min(table(cls.Clean)) < 3 | length(levels(cls.Clean)) < 2){
              #print(paste ("A total of", length(levels(cls.Clean)), "groups found with", length(cls.Clean), "samples."));
              #print("<font color='red'>At least <b>two</b> groups and <b>three replicates</b> per group are required for analysis</font>!");
              #print("You can click the <b>Edit Groups</b> button below to see the group labels for each sample and make corrections.");
              return(-1);
            }
          }
          
          if("NMDR_id" %in% names(mSetObj$dataSet)){
            msg <- c(msg, paste("Study", mSetObj$dataSet$NMDR_id, "was successfully downloaded from the Metabolomics Workbench!"))
          }
          msg <- c(msg,"Samples are not paired.");
        }
        
        # checking if too many groups but a few samples in each group
        cls.lbl <- mSetObj$dataSet$orig.cls;
        # need to exclude QC or blank
        qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
        if(sum(qb.inx) > 0){
          cls.lbl <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
        }
        min.grp.size <- min(table(cls.lbl));
        cls.num <- length(levels(cls.lbl));
        if(cls.num/min.grp.size > 3){
          mSetObj$dataSet$small.smpl.size <- 1;
          msg <- c(msg, "<font color='red'>Too many groups with very small number of replicates!</font>");
          msg <- c(msg, "<font color='red'>Only a subset of methods will be available for analysis!</font>");
        }
        
        # if(is.null(mSetObj$dataSet$meta.info)){
        if(mSetObj$analSet$type == "ts"){
          msg <- c(msg, paste0(cls.num, " groups were detected from primary meta-data factor: ", colnames(mSetObj$dataSet$meta.info)[1], "."));
          #msg <- c(msg, paste0(length(colnames(mSetObj$dataSet$meta.info)), " meta-data factors were detected: ", paste0(colnames(mSetObj$dataSet$meta.info), collapse=", "), "."));
          cls.vec <- vector()
          meta.info  <- mSetObj$dataSet$meta.info
          meta.types <- mSetObj$dataSet$meta.types
          for(i in 1:ncol(meta.info)){
            if(meta.types[i] == "disc"){
              cls.lbl <- meta.info[,i];
              qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
              if(sum(qb.inx) > 0){
                cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
              } else {
                cls.Clean <- cls.lbl;
              }
              meta.name <- colnames(meta.info)[i]
              if(min(table(cls.Clean)) < 3 | length(levels(cls.Clean)) < 2){
                cls.vec <- c(cls.vec, meta.name)
              }
            }
          }
        }else{
          msg <- c(msg, paste(cls.num, "groups were detected in samples."));
        }
        
        if("NMDR_id" %in% names(mSetObj$dataSet)){
          msg <- c(msg, paste("Study", mSetObj$dataSet$NMDR_id, "group labels:", paste0(unique(cls.lbl), collapse = ", ")))
        }
        
        mSetObj$dataSet$cls.num <- cls.num;
        mSetObj$dataSet$min.grp.size <- min.grp.size;
      }
      
      #samples may not be sorted properly, need to do some sorting at the beginning 
      if(substring(mSetObj$dataSet$format,4,5)=="ts"){
        metadata <- mSetObj$dataSet$meta.info
        nfacA <- metadata[,1];
        nfacB <- metadata[,2];
        if(mSetObj$dataSet$design.type =="time" | mSetObj$dataSet$design.type =="time0"){
          # determine time factor and should order first by subject then by each time points
          if(tolower(colnames(metadata)[1]) == "time"){ 
            time.fac <- nfacA;
            exp.fac <- nfacB;
          }else{
            time.fac <- nfacB;
            exp.fac <- nfacA;
          }
          # update with new index
          ord.inx <- order(exp.fac);
        }else{
          ord.inx <- order(nfacA);
        }
        mSetObj$dataSet$orig.cls <- mSetObj$dataSet$orig.cls[ord.inx];
        mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[ord.inx];
        mSetObj$dataSet$facA <- mSetObj$dataSet$orig.facA <- metadata[,1][ord.inx];
        mSetObj$dataSet$facB <- mSetObj$dataSet$orig.facB <- metadata[,2][ord.inx];
        orig.data <- orig.data[ord.inx,];
        mSetObj$dataSet$meta.info <- mSetObj$dataSet$meta.info[rownames(orig.data), ];
        qs::qsave(orig.data, file="data_orig.qs");
      }else{
        ord.inx <- order(mSetObj$dataSet$orig.cls);
        mSetObj$dataSet$orig.cls <- cls[ord.inx];
        mSetObj$dataSet$url.smp.nms <- mSetObj$dataSet$url.smp.nms[ord.inx];
        orig.data <- orig.data[ord.inx, , drop=FALSE];
        qs::qsave(orig.data, file="data_orig.qs");
        if(mSetObj$dataSet$paired){
          mSetObj$dataSet$pairs <- mSetObj$dataSet$pairs[ord.inx];
        }
      }
    }
    msg<-c(msg,"Only English letters, numbers, underscore, hyphen and forward slash (/) are allowed.");
    msg<-c(msg,"<font color=\"orange\">Other special characters or punctuations (if any) will be stripped off.</font>");
    
    int.mat <- orig.data;
    
    if(ncol(int.mat)==1){
      if(anal.type=="roc"){
        mSetObj$dataSet$roc_cols <- 1;
      } else {
        #print("<font color='red'>One-column data is only supported for biomarker analysis.</font>");
        return(0);
      }
    } else {
      mSetObj$dataSet$roc_cols <- 2;
    }
    
    # check numerical matrix
    rowNms <- rownames(int.mat);
    colNms <- colnames(int.mat);
    naNms <- sum(is.na(int.mat));
    
    for(c in 1:ncol(int.mat)) {
      if(class(int.mat[,c]) == "integer64"){
        int.mat[,c] <- as.double(int.mat[,c]);
      }
    }
    
    num.mat <- apply(int.mat, 2, as.numeric)
    
    if(sum(is.na(num.mat)) > naNms){
      # try to remove "," in thousand seperator if it is the cause
      num.mat <- apply(int.mat,2,function(x) as.numeric(gsub(",", "", x)));
      if(sum(is.na(num.mat)) > naNms){
        msg<-c(msg,"<font color=\"red\">Non-numeric values were found and replaced by NA.</font>");
      }else{
        msg<-c(msg,"All data values are numeric.");
      }
    }else{
      msg<-c(msg,"All data values are numeric.");
    }
    
    int.mat <- num.mat;
    rownames(int.mat) <- rowNms;
    colnames(int.mat)<- colNms;
    
    # check for columns with all constant (var =0)
    varCol <- apply(int.mat, 2, var, na.rm=T);
    
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      msg<-c(msg, paste("<font color=\"red\">", constNum, "features with a constant or single value across samples were found and deleted.</font>"));
      int.mat <- int.mat[,!constCol, drop=FALSE];
    }
    
    # check zero, NA values
    totalCount <- nrow(int.mat)*ncol(int.mat);
    naCount <- sum(is.na(int.mat));
    naPercent <- round(100*naCount/totalCount,1)
    #  print(naCount)
    mSetObj$dataSet$missingCount <- naCount;
    
    msg<-c(msg, paste("A total of ", naCount, " (", naPercent, "%) missing values were detected.", sep=""));
    msg<-c(msg, "<u>By default, missing values will be replaced by 1/5 of min positive values of their corresponding variables</u>",
           "Click the <b>Proceed</b> button if you accept the default practice;",
           "Or click the <b>Missing Values</b> button to use other methods.");
    
    qs::qsave(as.data.frame(int.mat), "preproc.qs");
    mSetObj$dataSet$proc.cls <- mSetObj$dataSet$cls <- mSetObj$dataSet$orig.cls;
    if(is.null(mSetObj$dataSet$meta.info)){
      mSetObj$dataSet$meta.info <- data.frame(mSetObj$dataSet$cls);
      colnames(mSetObj$dataSet$meta.info) = "Class";
    }
    
    if(substring(mSetObj$dataSet$format,4,5)=="ts"){
      mSetObj$dataSet$proc.facA <- mSetObj$dataSet$orig.facA;
      mSetObj$dataSet$proc.facB <- mSetObj$dataSet$orig.facB;
    }
    
    mSetObj$msgSet$check.msg <- c(mSetObj$msgSet$read.msg, msg);
    #print(c("Successfully passed sanity check!", msg))
    
    return(mSetObj);
  }
  
  computeAnova<-function(inputData, outDir){
    s <- getGroupNames(inputData)
    dat<- getSetsData(inputData,s)
    preproc <- qs::qread("preproc.qs")
    dat<- ReplaceMissingByLoD(preproc)
    
    return(dat)
  }
  
  InitDataObjects <- function(data.type, anal.type, paired=FALSE){
    if(exists("mSet")){
      mSetObj <- .get.mSet(mSet);
      mSetObj$dataSet$type <- data.type;
      mSetObj$dataSet$paired <- paired;
      mSetObj$analSet$type <- anal.type;
      mSetObj<-CleanDataObjects(mSetObj, anal.type);
      return(mSetObj);
    }
    
    dataSet <- list();
    dataSet$type <- data.type;
    dataSet$design.type <- "regular"; # one factor to two factor
    dataSet$cls.type <- "disc"; # default until specified otherwise
    dataSet$format <- "rowu";
    dataSet$paired <- paired;
    analSet <- list();
    analSet$type <- anal.type;
    Sys.setenv("OMP_NUM_THREADS" = 2); # to control parallel computing for some packages
    Sys.setenv("OPENBLAS_NUM_THREADS" = 2);
    mSetObj <- list();
    mSetObj$dataSet <- dataSet;
    mSetObj$analSet <- analSet;
    mSetObj$imgSet <- list();
    mSetObj$msgSet <- list(); # store various message during data processing
    mSetObj$msgSet$msg.vec <- vector(mode="character");     # store error messages
    mSetObj$cmdSet <- vector(mode="character"); # store R command
    
    if (anal.type == "mummichog") {
      # Define this parameter set to avoid global variable
      # Author: Zhiqiang
      mSetObj$paramSet$mumRT <- NA;
      mSetObj$paramSet$mumRT.type <- NA;
      mSetObj$paramSet$version <- NA;
      mSetObj$paramSet$mumDataContainsPval <- 1;
      mSetObj$paramSet$mode <- NA;
      mSetObj$paramSet$adducts <- NA;
      mSetObj$paramSet$peakFormat <- "mpt";
    } 
    else if (anal.type == "metapaths") {
      # Define this parameter set to avoid global variable
      # Author: Zhiqiang
      paramSet <- list();
      paramSet$mumRT <- NA;
      paramSet$mumRT.type <- NA;
      paramSet$version <- NA;
      paramSet$mumDataContainsPval <- 1;
      paramSet$mode <- NA;
      paramSet$adducts <- NA;
      paramSet$peakFormat <- "mpt";
      paramSet$metaNum <- 0;
      mSetObj$paramSet <- paramSet;
      # This is an empty paramSet, and will be copied for multiple datasets
      dataNMs <- names(mSetObj)[grepl("MetaData",names(mSetObj))];
      if(length(dataNMs)>0){
        for(n in dataNMs){
          mSetObj[[n]] <- NULL;
        }
      }
    }
    
    #print("MetaboAnalyst R objects initialized ...");
    return(mSetObj);
  }
  
  Read.TextData <- function(mSetObj=NA, filePath, format="rowu", lbl.type="disc", nmdr = FALSE){
    mSetObj$dataSet$cls.type <- lbl.type;
    mSetObj$dataSet$format <- format;
    anal.type <- "stat"
    
    if(nmdr){
      dat <- qs::qread("nmdr_study.qs")
    }else{
      dat <- .readDataTable(filePath);
    }
    
    if(class(dat) == "try-error" || ncol(dat) == 1){
      #print("Data format error. Failed to read in the data!");
      #print("Make sure the data table is saved as comma separated values (.csv) format!");
      #print("Please also check the followings: ");
      #print("Either sample or feature names must in UTF-8 encoding; Latin, Greek letters are not allowed.");
      #print("We recommend to use a combination of English letters, underscore, and numbers for naming purpose.");
      #print("Make sure sample names and feature (peak, compound) names are unique.");
      #print("Missing values should be blank or NA without quote.");
      #print("Make sure the file delimeters are commas.");
      return(0);
    }
    
    msg <- NULL;
    
    if(substring(format,4,5)=="ts"){
      # two factor time series data
      if(substring(format,1,3)=="row"){ # sample in row
        msg<-c(msg, "Samples are in rows and features in columns");
        smpl.nms <-dat[,1];
        all.nms <- colnames(dat);
        facA.lbl <- all.nms[2];
        cls.lbl <- facA <- dat[,2]; # default assign facA to cls.lbl in order for one-factor analysis
        facB.lbl <- all.nms[3];
        facB <- dat[,3];
        conc <- dat[,-c(1:3)];
        var.nms <- colnames(conc);
      }else{ # sample in col
        msg<-c(msg, "Samples are in columns and features in rows.");
        all.nms <- dat[,1];
        facA.lbl <- all.nms[1];
        cls.lbl <- facA <- dat[1,-1];
        facB.lbl <- all.nms[2];
        facB <- dat[2,-1];
        var.nms <- dat[-c(1:2),1];
        conc<-t(dat[-c(1:2),-1]);
        smpl.nms <- rownames(conc);
      }
      
      metadata <- data.frame(metaA=as.factor(as.character(facA)), metaB=as.factor(as.character(facB)), stringsAsFactors=T);
      colnames(metadata) <- c(facA.lbl, facB.lbl)
      mSetObj$dataSet$meta.info <- metadata
      mSetObj$dataSet$types.cls.lbl <- sapply(metadata, function(x) class(x) ) 
      mSetObj$dataSet$meta.types <- c("disc", "disc");
      names(mSetObj$dataSet$types.cls.lbl) <- c(facA.lbl, facB.lbl)
      names(mSetObj$dataSet$meta.types) <- c(facA.lbl, facB.lbl)
      
    }else{
      
      if(substring(format,1,3)=="row"){ # sample in row
        msg <- c(msg, "Samples are in rows and features in columns");
        smpl.nms <-dat[,1];
        dat[,1] <- NULL; #remove sample names
        if(lbl.type == "qc"){
          rownames(dat) <- smpl.nms;
          #mSetObj$dataSet$orig <- dat;
          qs::qsave(dat, file="data_orig.qs");
          mSetObj$dataSet$cmpd <- colnames(dat);
          return(1);
        }
        cls.lbl <- dat[,1];
        conc <- dat[,-1, drop=FALSE];
        var.nms <- colnames(conc);
        if(lbl.type == "no"){ #no class label
          cls.lbl <- rep(1, nrow(dat));
          conc <- dat[,, drop=FALSE]; 
          var.nms <- colnames(conc);
        }
      }else{ # sample in col
        msg<-c(msg, "Samples are in columns and features in rows.");
        if(lbl.type == "no"){
          cls.lbl <- rep(1, ncol(dat));
          conc <- t(dat[,-1]);
          var.nms <- dat[,1];
          smpl.nms <- colnames(dat[,-1]);
        }else{
          var.nms <- dat[-1,1];
          dat[,1] <- NULL;
          smpl.nms <- colnames(dat);
          cls.lbl <- dat[1,];
          conc <- t(dat[-1,]);
        }
      }
    }
    
    mSetObj$dataSet$type.cls.lbl <- class(cls.lbl);
    
    msg <- c(msg, "The uploaded file is in comma separated values (.csv) format.");
    
    # try to remove empty line if present
    # identified if no sample names provided
    
    empty.inx <- is.na(smpl.nms) | smpl.nms == ""
    if(sum(empty.inx) > 0){
      msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty rows</font> were detected and excluded from your data."));
      smpl.nms <- smpl.nms[!empty.inx];
      cls.lbl <-  cls.lbl[!empty.inx];
      conc <- conc[!empty.inx, ];
    }
    
    # try to check & remove empty lines if class label is empty
    # Added by B. Han
    empty.inx <- is.na(cls.lbl) | cls.lbl == ""
    if(sum(empty.inx) > 0){
      if(mSetObj$analSet$type != "roc"){
        msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty labels</font> were detected and excluded from your data."));
        smpl.nms <- smpl.nms[!empty.inx];
        cls.lbl <-  cls.lbl[!empty.inx];
        conc <- conc[!empty.inx, ];
      }else{
        # force all NA to empty string, otherwise NA will become "NA" class label
        cls.lbl[is.na(cls.lbl)] <- "";
        msg <- c(msg, paste("<font color=\"orange\">", sum(empty.inx), "new samples</font> were detected from your data."));
      }
    }
    
    if(anal.type == "roc"){
      if(length(unique(cls.lbl[!empty.inx])) > 2){
        print("ROC analysis is only defined for two-group comparisions!");
        return(0);
      }
    }
    
    # check for uniqueness of dimension name
    if(length(unique(smpl.nms))!=length(smpl.nms)){
      dup.nm <- paste(smpl.nms[duplicated(smpl.nms)], collapse=" ");
      print("Duplicate sample names are not allowed!");
      print(dup.nm);
      return(0);
    }
    
    # try to remove check & remove empty line if feature name is empty
    empty.inx <- is.na(var.nms) | var.nms == "";
    if(sum(empty.inx) > 0){
      msg <- c(msg, paste("<font color=\"red\">", sum(empty.inx), "empty features</font> were detected and excluded from your data."));
      var.nms <- var.nms[!empty.inx];
      conc <- conc[,!empty.inx];
    }
    
    if(length(unique(var.nms))!=length(var.nms)){
      dup.nm <- paste(var.nms[duplicated(var.nms)], collapse=" ");
      print("Duplicate feature names are not allowed!");
      print(dup.nm);
      return(0);
    }
    
    if(anal.type == "mummichog"){
      is.rt <- mSetObj$paramSet$mumRT;
      if(!is.rt){
        mzs <- as.numeric(var.nms);
        if(sum(is.na(mzs) > 0)){
          print("Make sure that feature names are numeric values (mass or m/z)!");
          return(0);
        }
      }
    }
    
    # now check for special characters in the data labels
    if(sum(is.na(iconv(smpl.nms)))>0){
      na.inx <- is.na(iconv(smpl.nms));
      nms <- paste(smpl.nms[na.inx], collapse="; ");
      print(paste("No special letters (i.e. Latin, Greek) are allowed in sample names!", nms, collapse=" "));
      return(0);
    }
    
    if(sum(is.na(iconv(var.nms)))>0){
      na.inx <- is.na(iconv(var.nms));
      nms <- paste(var.nms[na.inx], collapse="; ");
      print(paste("No special letters (i.e. Latin, Greek) are allowed in feature names!", nms, collapse=" "));
      return(0);
    }
    
    # black slash could kill Rserve immediately
    smpl.nms <- gsub("\\\\", "-", smpl.nms);
    url.smp.nms <- CleanNames(smpl.nms);
    names(url.smp.nms) <- smpl.nms;
    
    var.nms <- gsub("\\\\", "-", var.nms);
    url.var.nms <- CleanNames(var.nms); # allow space, comma and period
    names(url.var.nms) <- var.nms;
    
    cls.lbl <- ClearStrings(as.vector(cls.lbl));
    
    # now assgin the dimension names
    rownames(conc) <- smpl.nms;
    colnames(conc) <- var.nms;
    
    # check if paired or not
    if(mSetObj$dataSet$paired){
      # save as it is and process in sanity check step
      mSetObj$dataSet$orig.cls <- mSetObj$dataSet$pairs <- cls.lbl;
    } else {
      if(lbl.type == "disc"){
        mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- as.factor(as.character(cls.lbl));
        
        if(substring(format,4,5)=="ts"){
          
          mSetObj$dataSet$facA.type <- is.numeric(facA);
          mSetObj$dataSet$orig.facA <- mSetObj$dataSet$facA <- as.factor(as.character(facA));
          mSetObj$dataSet$facA.lbl <- facA.lbl;
          
          mSetObj$dataSet$facB.type <- is.numeric(facB);
          mSetObj$dataSet$orig.facB <- mSetObj$dataSet$facB <- as.factor(as.character(facB));
          mSetObj$dataSet$facB.lbl <- facB.lbl;
        }
        
      } else { # continuous
        
        mSetObj$dataSet$orig.cls <- mSetObj$dataSet$cls <- tryCatch({
          as.numeric(cls.lbl);
        },warning=function(na) {
          print("Class labels must be numeric and continuous!");
          return(0);
        })
        
        if(mSetObj$dataSet$cls == 0){
          print("Class labels must be numeric and continuous!");
          return(0)
        }
      }
    }
    
    # for the current being to support MSEA and MetPA
    if(mSetObj$dataSet$type == "conc"){
      mSetObj$dataSet$cmpd <- var.nms;
    }
    
    mSetObj$dataSet$mumType <- "table";
    mSetObj$dataSet$url.var.nms <- url.var.nms;
    mSetObj$dataSet$url.smp.nms <- url.smp.nms;
    #mSetObj$dataSet$orig <- conc; # copy to be processed in the downstream
    qs::qsave(conc, file="data_orig.qs");
    mSetObj$msgSet$read.msg <- c(msg, paste("The uploaded data file contains ", nrow(conc),
                                            " (samples) by ", ncol(conc), " (", tolower(GetVariableLabel(mSetObj$dataSet$type)), ") data matrix.", sep=""));
    
    return(mSetObj);
  }
  
  CleanNames <- function(query){
    query <- gsub("[^[:alnum:].@_-]", "", query);
    return(make.unique(query));
  }
  
  ClearStrings<-function(query){
    # kill multiple white space
    query <- gsub(" +"," ",query);
    
    # black slash escape sign could kill Rserve immediately
    query <- gsub("\\\\", "-", query);
    
    # remove leading and trailing space
    query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
    return(query);
  }
  
  GetVariableLabel<-function(data.type){
    if(data.type=="conc"){
      return("Compounds");
    }else if(data.type=="specbin"){
      return("Spectra Bins");
    }else if(data.type=="nmrpeak"){
      return("Peaks (ppm)");
    }else if(data.type=="mspeak"){
      return("Peaks (mass)");
    }else{
      return("Peaks(mz/rt)");
    }
  }
  
  .readDataTable <- function(fileName){
    
    dat <- tryCatch(
      {
        data.table::fread(fileName, header=TRUE, check.names=FALSE, 
                          blank.lines.skip=TRUE, data.table=FALSE, integer64 = "numeric");
      }, error=function(e){
        #print(e);
        return(.my.slowreaders(fileName));    
      }, warning=function(w){
        #print(w);
        return(.my.slowreaders(fileName));
      });
    
    if(any(dim(dat) == 0)){
      dat <- .my.slowreaders(fileName);
    }
    return(dat);
  }
  
  .my.slowreaders <- function(fileName){
    #print("Using slower file reader ...");
    formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
    if(formatStr == "txt"){
      dat <- try(read.table(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
    }else{ # note, read.csv is more than read.table with sep=","
      dat <- try(read.csv(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
    }  
    return(dat);
  }
  
  ReplaceMin <- function(mSetObj=NA){
    
    #Reset to default
    mSetObj$dataSet$filt <- mSetObj$dataSet$edit <- NULL;
    
    # replace zero and missing values using Detection Limit for each variable 
    preproc <- qs::qread("preproc.qs");
    int.mat <- ReplaceMissingByLoD(preproc);  
    
    # note, this is last step of processing, also save to proc
    #mSetObj$dataSet$proc <- as.data.frame(int.mat);
    mSetObj$dataSet$proc.feat.num <- ncol(int.mat);
    qs::qsave(as.data.frame(int.mat), file="data_proc.qs");
    
    mSetObj$msgSet$replace.msg <- paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable.");
    invisible(gc()); # suppress gc output
    
    return(mSetObj);
  }
  
  PerformFeatureFilter <- function(int.mat, filter, remain.num = NULL, anal.type = NULL){
    feat.num <- ncol(int.mat);
    feat.nms <- colnames(int.mat);
    nm <- NULL;
    msg <- "";
    if(filter == "none" && feat.num < 5000){ # only allow for less than 4000
      remain <- rep(TRUE, feat.num);
      #msg <- paste(msg, "No filtering was applied");
    }else{
      if (filter == "rsd"){
        sds <- apply(int.mat, 2, sd, na.rm=T);
        mns <- apply(int.mat, 2, mean, na.rm=T);
        filter.val <- abs(sds/mns);
        nm <- "Relative standard deviation";
      }else if (filter == "nrsd" ){
        mads <- apply(int.mat, 2, mad, na.rm=T);
        meds <- apply(int.mat, 2, median, na.rm=T);
        filter.val <- abs(mads/meds);
        nm <- "Non-paramatric relative standard deviation";
      }else if (filter == "mean"){
        filter.val <- apply(int.mat, 2, mean, na.rm=T);
        nm <- "mean";
      }else if (filter == "sd"){
        filter.val <- apply(int.mat, 2, sd, na.rm=T);
        nm <- "standard deviation";
      }else if (filter == "mad"){
        filter.val <- apply(int.mat, 2, mad, na.rm=T);
        nm <- "Median absolute deviation";
      }else if (filter == "median"){
        filter.val <- apply(int.mat, 2, median, na.rm=T);
        nm <- "median";
      }else{ # iqr
        filter.val <- apply(int.mat, 2, IQR, na.rm=T);
        nm <- "Interquantile Range";
      }
      
      # get the rank of the filtered variables
      rk <- rank(-filter.val, ties.method='random');
      
      if(is.null(remain.num)){ # apply empirical filtering based on data size
        if(feat.num < 250){ # reduce 5%
          remain <- rk < feat.num*0.95;
          #msg <- paste(msg, "Further feature filtering based on", nm);
        }else if(feat.num < 500){ # reduce 10%
          remain <- rk < feat.num*0.9;
          #msg <- paste(msg, "Further feature filtering based on", nm);
        }else if(feat.num < 1000){ # reduce 25%
          remain <- rk < feat.num*0.75;
          #msg <- paste(msg, "Further feature filtering based on", nm);
        }else{ # reduce 40%, if still over 5000, then only use top 5000
          remain <- rk < feat.num*0.6;
          #msg <- paste(msg, "Further feature filtering based on", nm);
          
          if(anal.type == "mummichog"){
            max.allow <- 7500;
          }else if(anal.type == "power" || anal.type == "ts"){
            max.allow <- 5000;
          }else{
            max.allow <- 2500;
          }
          
          if(sum(remain) > max.allow){
            remain <- rk < max.allow;
            #msg <- paste(msg, paste("Reduced to", max.allow, "features based on", nm));
          }
        }
      }else{
        remain <- rk < remain.num;
      }
    }
    #print(msg);
    return(list(data=int.mat[, remain], msg=msg));
  }
  
  PreparePrenormData <- function(mSetObj=NA){
    
    if(!is.null(mSetObj$dataSet$edit)){
      mydata <- mSetObj$dataSet$edit;
      if(!is.null(mSetObj$dataSet$filt)){
        # some features could be removed
        hit.inx <- colnames(mydata) %in% colnames(mSetObj$dataSet$filt);
        mydata <- mydata[,hit.inx, drop=FALSE];
      }
      prenorm <- mydata;
      mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$edit.cls;
      if(substring(mSetObj$dataSet$format,4,5) == "ts"){
        mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$edit.facA;
        mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$edit.facB;
      }
    }else if(!is.null(mSetObj$dataSet$filt)){
      prenorm <- mSetObj$dataSet$filt;
      mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$filt.cls;
      if(substring(mSetObj$dataSet$format,4,5)=="ts"){
        mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$filt.facA;
        mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$filt.facB;
      }
    }else{
      prenorm <- qs::qread("data_proc.qs");
      mSetObj$dataSet$prenorm.cls <- mSetObj$dataSet$proc.cls;
      if(substring(mSetObj$dataSet$format,4,5) == "ts"){
        mSetObj$dataSet$prenorm.facA <- mSetObj$dataSet$proc.facA;
        mSetObj$dataSet$prenorm.facB <- mSetObj$dataSet$proc.facB;
      }
    }
    qs::qsave(prenorm, "prenorm.qs");
    mSetObj$dataSet$prenorm.smpl.nms <- rownames(prenorm);
    mSetObj$dataSet$prenorm.feat.nms <- colnames(prenorm);
    return(mSetObj)
  }
  
  CleanData <-function(bdata, removeNA=T, removeNeg=T, removeConst=T){
    
    if(sum(bdata==Inf, na.rm=TRUE)>0){
      inx <- bdata == Inf;
      bdata[inx] <- NA;
      bdata[inx] <- max(bdata, na.rm=T)*2
    }
    if(sum(bdata==-Inf, na.rm=TRUE)>0){
      inx <- bdata == -Inf;
      bdata[inx] <- NA;
      bdata[inx] <- min(bdata, na.rm=T)/2
    }
    
    if(removeNA){
      if(sum(is.na(bdata))>0){
        bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
      }
    }
    if(removeNeg){
      if(sum(as.numeric(bdata<=0)) > 0){
        inx <- bdata <= 0;
        bdata[inx] <- NA;
        bdata[inx] <- min(bdata, na.rm=T)/2
      }
    }
    if(removeConst){
      varCol <- apply(data.frame(bdata), 2, var, na.rm=T); # getting an error of dim(X) must have a positive length, fixed by data.frame
      constCol <- (varCol == 0 | is.na(varCol));
      constNum <- sum(constCol, na.rm=T);
      
    }
    bdata;
  }
  
  aof <- function(x, cls) {
    aov(x ~ cls);
  }
  
  LSD.test <- function(y, trt, alpha = 0.05){
    clase<-c("aov","lm")
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    if("aov"%in%class(y) | "lm"%in%class(y)){
      A<-y$model
      DFerror<-df.residual(y)
      MSerror<-deviance(y)/DFerror
      y<-A[,1]
      ipch<-pmatch(trt,names(A))
      name.t <-names(A)[ipch]
      trt<-A[,ipch]
      name.y <- names(A)[1]
    }
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    means <- tapply.stat(junto[, 1], junto[, 2], stat="mean") #change
    sds <- tapply.stat(junto[, 1], junto[, 2], stat="sd")     #change
    nn <- tapply.stat(junto[, 1], junto[, 2], stat="length")  #change
    std.err <- sds[, 2]/sqrt(nn[, 2])
    Tprob <- qt(1 - alpha/2, DFerror)
    LCL <- means[,2]-Tprob*std.err
    UCL <- means[,2]+Tprob*std.err
    means <- data.frame(means, std.err, replication = nn[, 2], LCL, UCL)
    names(means)[1:2] <- c(name.t, name.y)
    #row.names(means) <- means[, 1]
    ntr <- nrow(means)
    nk <- choose(ntr, 2)
    nr <- unique(nn[, 2])
    
    comb <- combn(ntr, 2)
    nn <- ncol(comb)
    dif <- rep(0, nn)
    LCL1<-dif
    UCL1<-dif
    sig<-NULL
    pvalue <- rep(0, nn)
    for (k in 1:nn) {
      i <- comb[1, k]
      j <- comb[2, k]
      if (means[i, 2] < means[j, 2]){
        comb[1, k]<-j
        comb[2, k]<-i
      }
      dif[k] <- abs(means[i, 2] - means[j, 2])
      sdtdif <- sqrt(MSerror * (1/means[i, 4] + 1/means[j,4]))
      pvalue[k] <- 2 * (1 - pt(dif[k]/sdtdif, DFerror));
      pvalue[k] <- round(pvalue[k],6);
      LCL1[k] <- dif[k] - Tprob*sdtdif
      UCL1[k] <- dif[k] + Tprob*sdtdif
      sig[k]<-" "
      if (pvalue[k] <= 0.001) sig[k]<-"***"
      else  if (pvalue[k] <= 0.01) sig[k]<-"**"
      else  if (pvalue[k] <= 0.05) sig[k]<-"*"
      else  if (pvalue[k] <= 0.1) sig[k]<-"."
    }
    tr.i <- means[comb[1, ],1]
    tr.j <- means[comb[2, ],1]
    output<-data.frame("Difference" = dif, pvalue = pvalue,sig,LCL=LCL1,UCL=UCL1)
    rownames(output)<-paste(tr.i,tr.j,sep=" - ");
    output;
  }
  
}

strippedData <- function(inputData, endPoint){
  fC <- c(rep(0.0,ncol(inputData)))
  count <- 0
  
  for (i in 1:endPoint-1) 
  {
    if (count == 0) {
      fc <- data.frame(inputData[i,])
    }
    else{
      fc<- rbind(fc,inputData[i,])
      
    }
    count<- count +1
    
  }
  
  return(fc)
  
}

getCalculationsData <- function(inputData, endPoint){
  target<-c(1:endPoint-1)
  data<- inputData[-target,]
  
  return(data)
  
}

OutputData<-function(OutputDirectory,GroupCSVDirectory,CommentColumn,type){
  statsFilepath <- "tmp"
  
  if (tolower(type) == "neg") {
    statsFilepath<-paste(OutputDirectory,"/NegIDed_FIN.csv",sep="")
  }else if (tolower(type) == "pos") {
    statsFilepath<-paste(OutputDirectory,"/PosIDed_FIN.csv",sep="")
  }else{
    statsFilepath<-paste(OutputDirectory,"/CombinedIDed_FIN.csv",sep="")
  }
  
  data<- read.csv(statsFilepath, sep=",", na.strings="NA", dec=".", strip.white=TRUE,header=TRUE)
  dataSet <-data[,c(CommentColumn+11,1:(ncol(data)-1))]
  
  for (i in 1:length(dataSet$PredictedFrag_IDs)) {
    if (isTRUE(grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) && grepl("CYP", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
      dataSet$PredictedFrag_IDs[i]<-gsub("\\CYP.*", "CYP", dataSet$PredictedFrag_IDs[i])
    }
    
    if (isTRUE(grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE))&&grepl("EC 3", dataSet$PredictedFrag_IDs[i], fixed = TRUE)) {
      dataSet$PredictedFrag_IDs[i]<-gsub("\\EC 3.*", "EC3", dataSet$PredictedFrag_IDs[i])
    }
  }
  
  smpls<-c()
  
  if (file.exists(GroupCSVDirectory)&&file.size(GroupCSVDirectory)>2) {
    smpls<-manualGrouping(GroupCSVDirectory)
  }else{
    print("Warning, no sample groupings were detected. Use sample in the sample names for statistics and have the groupings written before the word sample; for example Control_Sample or ControlSample, would group that file as a control")
    return()
  }
  
  output<-determineManualGroupingCalculation(dataSet,OutputDirectory,GroupCSVDirectory)
  
  write.csv(output,statsFilepath)
  
  computePCA(statsFilepath,OutputDirectory,smpls,type)

  if (length(smpls)>2) {
    mset<-InitDataObjects("conc","stat",F)

    mset<-Read.TextData(mset,statsFilepath,"colu","disc")

    mset<-SanityCheckData(mset)

    mset<-ReplaceMin(mset)

    mset<-SanityCheckData(mset)

    mset<-FilterVariable(mset,"none","F",25)

    mset<-PreparePrenormData(mset)

    mset<-Normalization(mset,"none","F",25)

    mset<-ANOVA.Anal(mset,"F",0.05,"fisher",FALSE,OutputDirectory,type)

  }

  return(output)
  
}
