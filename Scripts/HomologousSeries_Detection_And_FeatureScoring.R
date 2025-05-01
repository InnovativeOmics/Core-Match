

############function to do KMD analysis, sort, and score for PFAS################
# Paul Stelben
# Jeremy Koelmel
# 09/02/2020

# For debugging
# IDedTable_dir <- file.path(InputDirectory,"Output/PosIDed_insilico.csv")
# OutDir = file.path(InputDirectory,"Output")
# Mass_col = MZColumn
# Retention_col = RTColumn
# RowStartForFeatureTableData <- RowStartForFeatureTableData-1
# KMD_Bin_Window = (PrecursorMassAccuracy*2)
# NegPos = "Pos"

#Inputs:
# 1) directory after EPA masterlist predicted annotation, 2) output directory
# 3) m/z, 4) RT, 5) the confident annotation column, 6) one or more frag column, 7) same for predicted from EPA master, 8) whether or not RT order should be accounted for

#Testing for Neg and Pos
# IDedTable_dir=file.path(InputDirectory,"Output/CombinedIDed_Fragments.csv")
# OutDir = file.path(InputDirectory,"Output")
# Mass_col = MZColumn
# Retention_col = RTColumn
# RowStartForFeatureTableData=RowStartForFeatureTableData-1
# KMD_Bin_Window = (PrecursorMassAccuracy*2)
# 
# #Testing for Neg
# IDedTable_dir=file.path(InputDirectory,"Output/NegIDed_Fragments.csv")
# OutDir = file.path(InputDirectory,"Output")
# Mass_col = MZColumn
# Retention_col = RTColumn
# RowStartForFeatureTableData=RowStartForFeatureTableData-1
# KMD_Bin_Window = (PrecursorMassAccuracy*2)


Scoring <- function (IDedTable_dir, OutDir, Mass_col, Retention_col, RowStartForFeatureTableData, RT_flagging, KMD_Bin_Window, upper, lower, RepeatingUnits_dir){
  
  #Import IDedTable
  IDedTable <- read.csv(IDedTable_dir,sep=",")
  IDedTable <- as.matrix(IDedTable)
  
  if (nrow(IDedTable) == 0) {
    message("No MSMS data, stopping scoring function")
    stop()
  }
  
  #Import RepeatingUnits
  RepeatingUnits <- read.csv(RepeatingUnits_dir,sep=",")
  RepeatingUnits <- as.matrix(RepeatingUnits)
  Names <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),1])
  Series <- rev(RepeatingUnits[which(RepeatingUnits[,3] == TRUE | RepeatingUnits[,3] == "TRUE" | RepeatingUnits[,3] == " TRUE"),2])
  
  # FOR DEBUGGING
  # IDedTable <- IDedTable[1:3000,]
  if (IMfirst==TRUE) {
    ID_Ranked_col <- which(colnames(IDedTable) == "ID_Ranked.1")
  } else {
    ID_Ranked_col <- which(colnames(IDedTable) == "ID_Ranked")
  }
  Potential_IDs_col <- which(colnames(IDedTable) == "Potential_IDs") ## PJS 4/7/2025
  Frags_col <- which(colnames(IDedTable) == "Frags") ## PJS 4/7/2025
  
  # add columns and calculate mass defect and Kendrick mass defect for Mass_col
  m <- matrix(0, nrow(IDedTable), 2)
  IDedTable <- cbind(IDedTable, m)
  nrow_IDedTable <- nrow(IDedTable)
  ncol_IDedTable <- ncol(IDedTable)
  colnames(IDedTable)[ncol_IDedTable - 1] <- "nominal mass"
  colnames(IDedTable)[ncol_IDedTable] <- "mass defect"
  Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
  nmass <- round(Mass_vec)
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 1] <- as.numeric(nmass)
  massd <- Mass_vec - nmass
  IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable] <- as.numeric(massd)
  md_col <- ncol_IDedTable
  hseries_count_cols <- c()
  
  # timestamp()
  
  length_Series <- length(Series)
  for(x in 1:length_Series){
    # print(x)
    # x <- 1
    Mass_vec <- as.numeric(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, Mass_col])
    exactMass <- Rdisop::getMolecule(Series[x], maxisotopes = 1)$exactmass ##added Rdisop:: as Rcdk has the same function
    Mass <- round(exactMass)
    m <- matrix(0, nrow(IDedTable), 5)
    IDedTable <- cbind(IDedTable, m)
    nrow_IDedTable <- nrow(IDedTable)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable - 4] <- paste(Names[x], "Kendrick_mass", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 3] <- paste(Names[x], "nominal_Kendrick_mass", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 2] <- paste(Names[x], "Kendrick_mass_defect", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 1] <- paste(Names[x], "KMD_Group", sep="_")
    colnames(IDedTable)[ncol_IDedTable] <- paste(Names[x], "homologous_series", sep="_")
    
    KM <- Mass_vec * (Mass/exactMass)
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 4] <- as.numeric(KM)
    nKM <- round(KM)
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 3] <- as.numeric(nKM)
    KMD <- KM - nKM
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 2] <- as.numeric(KMD)
    
    # Sort and group by KMD
    tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
    tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 2])),]
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
    KMD_col <- ncol_IDedTable - 2
    group <- 1
    start <- RowStartForFeatureTableData
    IDedTable[RowStartForFeatureTableData, KMD_col + 1] <- group
    for (i in (RowStartForFeatureTableData + 1):nrow_IDedTable) {
      if (as.numeric(IDedTable[i, KMD_col]) - as.numeric(IDedTable[start, KMD_col]) > KMD_Bin_Window) {
        IDedTable[start:(i-1), KMD_col + 1] <- group
        group <- group + 1
        start <- i
      }
    }
    if (start != nrow_IDedTable) {
      IDedTable[start:nrow_IDedTable, KMD_col + 1] <- group
    }
    
    # Calculate homologous series
    series <- 1
    for (i in RowStartForFeatureTableData:(nrow_IDedTable - 1)) {
      if (is.na(IDedTable[i, KMD_col + 2]) || IDedTable[i, KMD_col + 2] != 0) {
        next
      }
      group <- as.numeric(IDedTable[i, KMD_col + 1])
      grow <- which(IDedTable[, KMD_col + 1] == group)
      mrow <- grow[which(((as.numeric(IDedTable[grow, KMD_col - 1]) - as.numeric(IDedTable[i, KMD_col - 1])) / Mass) %% 1 == 0)]
      members <- as.numeric(IDedTable[mrow, KMD_col - 1])
      if (length(unique(members)) < 2) {
        IDedTable[mrow, KMD_col + 2] <- NA
      } else {
        IDedTable[mrow, KMD_col + 2] <- series
        series <- series + 1
      }
    }
    IDedTable[which(is.na(IDedTable[, KMD_col + 2])), KMD_col + 2] <- 0
    
    # Sort by homologous series
    # IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable])),]
    
    m <- matrix(0, nrow(IDedTable), 2)
    IDedTable <- cbind(IDedTable, m)
    # Phantom columns
    ncol_IDedTable <- ncol(IDedTable) + 8
    colnames(IDedTable)[ncol_IDedTable - 9] <- paste(Names[x], "Homologous_series_count", sep="_")
    colnames(IDedTable)[ncol_IDedTable - 8] <- paste(Names[x], "MD_Filter", sep="_")
    hseries_count_cols <- append(hseries_count_cols, ncol_IDedTable - 9)
    
    # Mass defect filtering
    IDedTable[, ncol_IDedTable - 8] <- (as.numeric(IDedTable[, md_col]) >= lower & as.numeric(IDedTable[, md_col]) <= upper)
    
    # Pre-flagging sorting
    tmptable <- IDedTable[RowStartForFeatureTableData:nrow_IDedTable,]
    tmptable <- tmptable[order(as.numeric(tmptable[, ncol_IDedTable - 10]), as.numeric(tmptable[, ncol_IDedTable - 11]), as.numeric(tmptable[, ncol_IDedTable - 13]), as.numeric(tmptable[, Retention_col])),]
    IDedTable[RowStartForFeatureTableData:nrow_IDedTable,] <- tmptable
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice0.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    # RT Flagging
    # which() function?
    m <- matrix(0, nrow(IDedTable), 2)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol_IDedTable - 1
    colnames(IDedTable)[ncol_IDedTable - 6] <- paste(Names[x], "RT_series", sep="_")
    series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 9]))
    if (series[1] == 0) {
      series <- series[-1]
    }
    for (i in series) {
      curr <- which(IDedTable[, ncol_IDedTable - 9] == i)
      mseries <- as.numeric(unique(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]))
      if (length(mseries) < 2) {
        next
      }
      length_mseries <- length(mseries)
      for (j in 1:length_mseries) {
        mcurr <- curr[which(as.numeric(IDedTable[curr, ncol_IDedTable - 12]) == mseries[j])]
        if (mseries[j] == mseries[1]) {
          nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
          IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
        } else if (mseries[j] == mseries[length(mseries)]) {
          prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
          IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
        } else {
          prev <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j - 1])]
          IDedTable[mcurr, ncol_IDedTable - 6] <- IDedTable[mcurr, Retention_col] > IDedTable[prev[1], Retention_col]
          nex <- curr[which(as.numeric(IDedTable[curr[1]:curr[length(curr)], ncol_IDedTable - 12]) == mseries[j + 1])]
          IDedTable[mcurr, ncol_IDedTable - 5] <- IDedTable[mcurr, Retention_col] < IDedTable[nex[length(nex)], Retention_col]
        }
      }
    }
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 6] == FALSE || IDedTable[i, ncol_IDedTable - 5] == FALSE) {
        IDedTable[i, ncol_IDedTable - 6] <- "RT not ordered"
      } else if (IDedTable[i, ncol_IDedTable - 6] == TRUE || IDedTable[i, ncol_IDedTable - 5] == TRUE) {
        IDedTable[i, ncol_IDedTable - 6] <- "RT ordered"
      } else {
        IDedTable[i, ncol_IDedTable - 6] <- NA ## PJS NA 4/8/2025
      }
    }
    
    IDedTable <- IDedTable[, 1:(ncol_IDedTable - 6)]
    
    ncol_IDedTable <- ncol_IDedTable + 1
    
    ## For debugging
    # write.table(IDedTable, file.path(OutDir, "NegIDed_test_optimized.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    if (RT_flagging == TRUE) {
      # Homologous seris count
      series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
      if (series[1] == 0) {
        series <- series[-1]
      }
      none <- 2
      for (i in series) {
        curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
        stable <- IDedTable[curr, c(ncol_IDedTable - 13, ncol_IDedTable - 7)]
        # Remove not ordered rows
        stable <- stable[!(is.na(stable[, 2]) == FALSE & stable[, 2] != "RT ordered" & stable[, 2] != "RT ordered, RT ordered"),]
        if (is.null(nrow(stable))) {
          none <- 1
        } else if (nrow(stable) == 0 && ncol(stable) == 2) {
          none <- 0
        }
        if (none == 2) {
          svec <- unique(stable[, 1])
          IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
        } else if (none == 1) {
          IDedTable[curr, ncol_IDedTable - 9] <- 1
          none <- 2
        } else if (none == 0) {
          IDedTable[curr, ncol_IDedTable - 9] <- 0
          none <- 2
        }
      }
    } else {
      series <- as.numeric(unique(IDedTable[RowStartForFeatureTableData:nrow_IDedTable, ncol_IDedTable - 10]))
      if (series[1] == 0) {
        series <- series[-1]
      }
      for (i in series) {
        curr <- which(IDedTable[, ncol_IDedTable - 10] == i)
        svec <- unique(IDedTable[curr, ncol_IDedTable - 13])
        IDedTable[curr, ncol_IDedTable - 9] <- length(svec)
      }
    }
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final_output_practice1.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  }
  
  # timestamp()
  
  # write.table(IDedTable, file.path(OutDir, "Test.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  
  if (Lipid == FALSE) {
    m <- matrix(0, nrow(IDedTable), 8)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable - 7] <- "Confident_ID"
    colnames(IDedTable)[ncol_IDedTable - 6] <- "2+_homologous_series"
    colnames(IDedTable)[ncol_IDedTable - 5] <- "3+_homologous_series"
    colnames(IDedTable)[ncol_IDedTable - 4] <- "Tentative_ID"
    colnames(IDedTable)[ncol_IDedTable - 3] <- "F_containing"
    colnames(IDedTable)[ncol_IDedTable - 2] <- "Exact_mass_match"
    colnames(IDedTable)[ncol_IDedTable - 1] <- "Score"
    colnames(IDedTable)[ncol_IDedTable] <- "Score_Description"
    
    # Scoring system
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "1_" || substr(IDedTable[i, ID_Ranked_col], 1, 2) == "2_") {
        IDedTable[i, ncol_IDedTable - 7] <- TRUE
      } else {
        IDedTable[i, ncol_IDedTable - 7] <- FALSE
      }
      IDedTable[i, ncol_IDedTable - 6] <- FALSE
      for (hsc_col in hseries_count_cols){
        if (as.numeric(IDedTable[i, hsc_col]) >= 2) {
          IDedTable[i, ncol_IDedTable - 6] <- TRUE
        }
      }
      IDedTable[i, ncol_IDedTable - 5] <- FALSE
      for (hsc_col in hseries_count_cols){
        if (as.numeric(IDedTable[i, hsc_col]) >= 3) {
          IDedTable[i, ncol_IDedTable - 5] <- TRUE
        }
      }
      if (is.na(IDedTable[i, Potential_IDs_col])) {
        IDedTable[i, ncol_IDedTable - 4] <- FALSE
      } else {
        IDedTable[i, ncol_IDedTable - 4] <- TRUE
      }
      Fvec <- strsplit(IDedTable[i, Frags_col], "")
      if (is.na(Fvec) == FALSE) {
        for (j in 1:length(Fvec[[1]])) {
          if (Fvec[[1]][j] == 'F') {
            IDedTable[i, ncol_IDedTable - 3] <- TRUE
            break
          }
        }
      }
      if (IDedTable[i, ncol_IDedTable - 3] != TRUE) {
        IDedTable[i, ncol_IDedTable - 3] <- FALSE
      }
      if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        IDedTable[i, ncol_IDedTable - 2] <- TRUE
      } else {
        IDedTable[i, ncol_IDedTable - 2] <- FALSE
      }
    }
    
    # timestamp()
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 7] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "A"
        IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments observed along with exact mass) and 2+ in homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 7] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "A-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 2 Schymanski: Confident ID (class specific dominant fragments and exact mass)"
      } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B+"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable[i, ncol_IDedTable - 3] == "TRUE" && IDedTable[i, ncol_IDedTable - 6] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable[i, ncol_IDedTable - 4] == "TRUE" || IDedTable[i, ncol_IDedTable - 3] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "B-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if ((IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") && IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D+"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 5] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: 3+ within homologous series"
      } else if (IDedTable[i, ncol_IDedTable - 9] == "TRUE" || IDedTable[i, ncol_IDedTable - 2] == "TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "D-"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS: Mass defect falling within -0.11 and 0.12 OR exact mass match"
      } else {
        IDedTable[i, ncol_IDedTable - 1] <- "E"
        IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: likely not PFAS"
      }
    }
    
    # timestamp()
    
    # for (i in RowStartForFeatureTableData:nrow_IDedTable) {
    if (ParallelComputing == TRUE) {
      
      IDedTable <- foreach (i = RowStartForFeatureTableData:nrow_IDedTable, .combine = rbind) %dopar% {
        #print(i)
        if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
          for (hsc_col in hseries_count_cols) {
            if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
              tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
              hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
              start <- hsc_rows[1]
              stop <- hsc_rows[length(hsc_rows)]
              if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
                IDedTable[i, ncol_IDedTable - 1] <- "C+"
                IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
              } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
                if (IDedTable[i, hsc_col] >= 3) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                } else if (IDedTable[i, hsc_col] >= 2) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C-"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                }
              }
            }
          }
        }
        return(IDedTable[i,])
      }
    } else {
      for (i in RowStartForFeatureTableData:nrow_IDedTable) {
        #print(i)
        if (IDedTable[i, ncol_IDedTable - 1] == "D+" || IDedTable[i, ncol_IDedTable - 1] == "D" || IDedTable[i, ncol_IDedTable - 1] == "D-") {
          for (hsc_col in hseries_count_cols) {
            if (as.numeric(IDedTable[i, hsc_col]) > 1) { ## TEST
              tmptable <- IDedTable[order(as.numeric(IDedTable[, hsc_col - 1])),]
              hsc_rows <- which(tmptable[, hsc_col - 1] == IDedTable[i, hsc_col - 1])
              start <- hsc_rows[1]
              stop <- hsc_rows[length(hsc_rows)]
              if (length(grep("A", tmptable[start:stop, ncol_IDedTable - 1])) > 0) {
                IDedTable[i, ncol_IDedTable - 1] <- "C+"
                IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 2+ in homologous series, and at least one confident PFAS identification within homologous series (A- or higher grade)"
              } else if (length(grep("B", tmptable[start:stop, ncol_IDedTable - 1])) > 0 && IDedTable[i, ncol_IDedTable - 1] != "C+" && IDedTable[i, ncol_IDedTable - 1] != "C") {
                if (IDedTable[i, hsc_col] >= 3) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: 3+ in homologous series, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                } else if (IDedTable[i, hsc_col] >= 2) {
                  IDedTable[i, ncol_IDedTable - 1] <- "C-"
                  IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: Possible ID, possible PFAS: homologous series 2+ within, and at least one highly likely PFAS identified within homologous series (B- or higher grade)"
                }
              }
            }
          }
        }
      }
    }
    
    # timestamp()
    
    # Change Ds to D+s when in series with D+s
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      #print(i)
      if (IDedTable[i, ncol_IDedTable - 1] == "D") {
        D_rows <- which(IDedTable[, ncol_IDedTable - 11] == IDedTable[i, ncol_IDedTable - 11])
        start <- D_rows[1]
        stop <- D_rows[length(D_rows)]
        for (j in start:stop) {
          if (IDedTable[j, ncol_IDedTable - 1] == "D+") {
            IDedTable[i, ncol_IDedTable - 1] <- "D+"
            IDedTable[i, ncol_IDedTable] <- "Lvl 5 Schymanski: No ID, possible PFAS:  Mass defect falling within -0.11 and 0.12 OR exact mass match, and 3+ within homologous series"
            break
          }
        }
      }
    }
    
    # timestamp()
    
    # write.table(IDedTable, paste("C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/", "IDedTableIDed_FIN_KMD_scored_final.csv",sep=""), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
    
    # Next seven lines for debugging
    # rm(list = ls())
    # IDedTable_dir <- "C:/Users/pstel/Documents/PFAS Lab/Agilent/SMILES2MSMS/IDedTableIDed_FIN_KMD_scored_final.csv"
    # IDedTable <- read.csv(IDedTable_dir,sep=",")
    # IDedTable <- as.matrix(IDedTable)
    # Retention_col <- 3
    # nrow_IDedTable <- nrow(IDedTable)
    # ncol_IDedTable <- ncol(IDedTable)
    
    # Remove tentative and F containing columns
    IDedTable <- IDedTable[, -((ncol_IDedTable - 4):(ncol_IDedTable - 3))]
    ncol_IDedTable <- ncol(IDedTable)
    
    ######### Sorting
    
    m <- matrix(0, nrow(IDedTable), 4)
    IDedTable <- cbind(IDedTable, m)
    ncol_IDedTable <- ncol(IDedTable) - 4
    
    # Find max series for each row
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      IDedTable[i, ncol_IDedTable + 2] <- max(which(IDedTable[i, hseries_count_cols] == max(IDedTable[i, hseries_count_cols])))
    }
    
    # Max series series number and count
    m <- matrix(0, nrow_IDedTable, 1)
    IDedTable <- cbind(IDedTable, m)
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      IDedTable[i, ncol_IDedTable + 3] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])] - 1])
      IDedTable[i, ncol_IDedTable + 4] <- as.numeric(IDedTable[i, hseries_count_cols[as.numeric(IDedTable[i, ncol_IDedTable + 2])]])
    }
    
    # Sort by max series count, max series, and series number
    IDedTable <- IDedTable[order(-as.numeric(IDedTable[, ncol_IDedTable + 4]), -as.numeric(IDedTable[, ncol_IDedTable + 2]), as.numeric(IDedTable[, ncol_IDedTable + 3])),]
    
    # Create rank
    rank <- 1
    IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- rank
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
      if (as.numeric(IDedTable[i, ncol_IDedTable + 4]) == 0) {
        rank <- rank + 1
      } else if (IDedTable[i, ncol_IDedTable + 3] != IDedTable[i - 1, ncol_IDedTable + 3]) {
        rank <- rank + 1
      }
      IDedTable[i, ncol_IDedTable + 5] <- rank
    }
    
    # Score
    for (i in RowStartForFeatureTableData:nrow(IDedTable)) {
      if (IDedTable[i, ncol_IDedTable - 1] == "A" || IDedTable[i, ncol_IDedTable - 1] == "A-") {
        IDedTable[i, ncol_IDedTable + 1] <- 1
      } else if (IDedTable[i, ncol_IDedTable - 1] == "B+" || IDedTable[i, ncol_IDedTable - 1] == "B" || IDedTable[i, ncol_IDedTable - 1] == "B-") {
        IDedTable[i, ncol_IDedTable + 1] <- 2
      } else if (IDedTable[i, ncol_IDedTable - 1] == "E") {
        IDedTable[i, ncol_IDedTable + 1] <- 4
      } else {
        IDedTable[i, ncol_IDedTable + 1] <- 3
      }
    }
    
    # Sort by rank and score
    IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, ncol_IDedTable + 1])),]
    
    # Adjust rank
    old <- IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]
    if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 1) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 10000000
    } else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 2) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) - 1000000
    } else if (IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 1] == 4) {
      IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5] <- as.numeric(IDedTable[RowStartForFeatureTableData, ncol_IDedTable + 5]) + 1000000
    }
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable)) {
      if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 4) {
        old <- IDedTable[i, ncol_IDedTable + 5]
        IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) + 1000000
      } else if ((IDedTable[i, ncol_IDedTable + 5] != old)) {
        old <- IDedTable[i, ncol_IDedTable + 5]
        if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 1) {
          IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 10000000
        } else if (as.numeric(IDedTable[i, ncol_IDedTable + 1]) == 2) {
          IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i, ncol_IDedTable + 5]) - 1000000
        }
      } else {
        old <- IDedTable[i, ncol_IDedTable + 5]
        IDedTable[i, ncol_IDedTable + 5] <- as.numeric(IDedTable[i - 1, ncol_IDedTable + 5])
      }
    }
    
    # Sort by rank and m/z
    IDedTable <- IDedTable[order(as.numeric(IDedTable[, ncol_IDedTable + 5]), as.numeric(IDedTable[, Mass_col])),]
    
    ## Is this right? I think so.
    # Final aesthetics
    ncol_IDedTable <- ncol_IDedTable - 1
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, ncol_IDedTable - 6] == "TRUE" || IDedTable[i, ncol_IDedTable - 6] == " TRUE") {
        IDedTable[i, ncol_IDedTable - 1] <- "TRUE"
      }
    }
  }
  
  # New formatting
  m <- matrix("", nrow_IDedTable, 5)
  n <- matrix("", nrow_IDedTable, 4)
  IDedTable <- cbind(m, IDedTable[, Mass_col], IDedTable[, Retention_col], n, IDedTable[, -c(Mass_col, Retention_col)])
  colnames(IDedTable)[1] <- "Score"
  colnames(IDedTable)[2] <- "SeriesType_Identifier"
  colnames(IDedTable)[3] <- "Name_or_Class"
  colnames(IDedTable)[4] <- "Formula"
  colnames(IDedTable)[5] <- "SMILES"
  colnames(IDedTable)[6] <- "m/z"
  colnames(IDedTable)[7] <- "Retention Time"
  colnames(IDedTable)[8] <- "Adduct"
  colnames(IDedTable)[9] <- "Unique"
  colnames(IDedTable)[10] <- "Score_Description"
  colnames(IDedTable)[11] <- "Needs_Validation"
  
  ncol_IDedTable <- ncol(IDedTable)
  ID_Ranked_col <- 9 + ID_Ranked_col
  Potential_IDs_col <- 9 + Potential_IDs_col
  Frags_col <- 9 + Frags_col
  
  
  if (Lipid == FALSE && Tween_pos == FALSE) {
    IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
    IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
    IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          #Not the best parsing function: This one looks for the second dash to pull out the PFAS name, as some PFAS have one - in the name and some don't it sometimes pulls in literature and sometimes not
          IDedTable[i, 3] <- paste(name[3:(which(name == "-")[2] - 1)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
          IDedTable[i, 9] <- IDedTable[i, ID_Ranked_col + 4]
        }
      } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
        primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- paste(name[1:(which(name == "-")[2] - 1)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
        }
        u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      } else if (length(grep("DTXSID", IDedTable[i, Frags_col - 1])) > 0) {
        all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
        primary <- strsplit(all[grep("DTXSID", all)][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          IDedTable[i, 3] <- primary[4]
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[1]
          IDedTable[i, 8] <- "[M-H]-"
        }
        all <- strsplit(all, ";", TRUE)
        u_entries <- c()
        for (j in 1:length(all)) {
          u_entries <- append(u_entries, all[[j]][1])
        }
        u_entries <- unique(u_entries)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(grep(";\\[", IDedTable[i, ID_Ranked_col])) == 0) {
          end <- strsplit(primary[length(primary)], "\\[")[[1]]
          primary[length(primary)] <- end[1]
          primary <- append(primary, paste("[", end[2], sep=""))
        }
        IDedTable[i, 3] <- NA ## PJS NA 4/8/2025
        IDedTable[i, 4] <- primary[length(primary) - 2]
        IDedTable[i, 5] <- NA ## PJS NA 4/8/2025
        IDedTable[i, 8] <- primary[length(primary)]
        all <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], "|", TRUE)[[1]], ";", TRUE)
        u_entries <- c()
        length_all <- length(all)
        for (j in 1:length_all) {
          if (length(grep("\\[", all[[j]][length(all[[j]])])) > 0) {
            end <- strsplit(all[[j]][length(all[[j]])], "\\[")[[1]]
            all[[j]][length(all[[j]])] <- end[1]
            all[[j]] <- append(all[[j]], paste("[", end[2], sep=""))
          }
          u_entries <- append(u_entries, all[[j]][length(all[[j]]) - 1])
        }
        u_entries <- unique(u_entries)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      }
      if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
        IDedTable[i, 11] <- "Yes"
      } else {
        IDedTable[i, 11] <- "No"
      }
    }
  } else if (Tween_pos == TRUE) {
    IDedTable[, 1] <- IDedTable[, ncol_IDedTable - 6]
    IDedTable[, 2] <- paste(Names[as.numeric(IDedTable[, ncol_IDedTable - 3])], IDedTable[, ncol_IDedTable - 2], sep="_") ## Error??
    IDedTable[, 10] <- IDedTable[, ncol_IDedTable - 5]
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      # Slight name modification
      if (IDedTable[i, 1] == "A" || IDedTable[i, 1] == "A-") {
        primary <- strsplit(strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
          IDedTable[i, 9] <- IDedTable[i, ID_Ranked_col + 4]
        }
        # Slight name modification
      } else if (!is.na(IDedTable[i, Potential_IDs_col])) {
        primary <- strsplit(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]][1], ";", TRUE)[[1]]
        if (length(primary) > 1) {
          name <- strsplit(primary[1], "")[[1]]
          IDedTable[i, 3] <- primary[1]
          IDedTable[i, 4] <- primary[2]
          IDedTable[i, 5] <- primary[3]
          IDedTable[i, 8] <- primary[4]
        }
        u_entries <- unique(strsplit(IDedTable[i, Potential_IDs_col], "|", TRUE)[[1]])
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
        ### Name modification
      } else if (length(grep("PEG", IDedTable[i, Frags_col - 1])) > 0) {
        all <- strsplit(IDedTable[i, Frags_col - 1], "|", TRUE)[[1]]
        # Take first PEG for name
        # primary <- strsplit(all[grep("PEG", all)][1], ";", TRUE)[[1]]
        primary <- all[grep("PEG", all)]
        if (length(primary) > 1) {
          IDedTable[i, 3] <- primary[1]
          IDedTable[i, 4] <- NA ## PJS NA 4/8/2025
          IDedTable[i, 5] <- NA ## PJS NA 4/8/2025
          # [M+?]+ fill in ? from PEG name
          IDedTable[i, 8] <- paste("[M+", strsplit(strsplit(primary[1], "_", TRUE)[[1]][2], "+", TRUE)[[1]][1], "]+", sep = "")
        }
        # Don't split, just check for unique
        u_entries <- unique(primary)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
        ### Name modification
      } else if (substr(IDedTable[i, ID_Ranked_col], 1, 2) == "4_") {
        primary <- strsplit(IDedTable[i, ID_Ranked_col], " | ", TRUE)[[1]]
        name <- strsplit(primary[1], "")[[1]]
        IDedTable[i, 3] <- paste(name[3:length(name)], collapse = "")
        IDedTable[i, 4] <- NA ## PJS NA 4/8/2025
        IDedTable[i, 5] <- NA ## PJS NA 4/8/2025
        IDedTable[i, 8] <- "[M+?]+"
        
        u_entries <- unique(primary)
        if (length(u_entries) == 1) {
          IDedTable[i, 9] <- "Yes"
        } else {
          IDedTable[i, 9] <- "No"
        }
      }
      if (IDedTable[i, 1] != "A" & IDedTable[i, 1] != "A-") {
        IDedTable[i, 11] <- "Yes"
      } else {
        IDedTable[i, 11] <- "No"
      }
    }
  } else if (Lipid == TRUE) {
    ################################## LIPIDMATCH SCORING SYSTEM
    IDedTable[, 1] <- substr(IDedTable[, ID_Ranked_col], 1, 1)
    IDedTable[, 2] <- paste(Names[1], IDedTable[, ncol_IDedTable - 3], sep="_")
    IDedTable[, 3] <- IDedTable[, ID_Ranked_col + 2]
    IDedTable[, 4] <- "CH2"
    IDedTable[, 5] <- "CC"
    
    IDedTable[, 8] <- IDedTable[, ID_Ranked_col + 3]
    IDedTable[which(is.na(IDedTable[, 8])), 8] <- 0
    
    IDedTable[, 9] <- IDedTable[, ID_Ranked_col + 4]
    IDedTable[which(is.na(IDedTable[, 9])), 9] <- "No"
    
    IDedTable[which(IDedTable[, 1] == 1), 10] <- "Annotated using class based rules from standards including by fatty acyl chain constituents with ~5% false positive rate at the level of class and C:DB"
    IDedTable[which(IDedTable[, 1] == 2), 10] <- "Annotated using class based rules from standards and DIA data with ~5% false positive rate at the level of class and C:DB"
    IDedTable[which(IDedTable[, 1] == 3), 10] <- "Annotated using class based rules from standards and DDA data but only annotated by class (fatty acid composition not known)"
    IDedTable[which(IDedTable[, 1] == 4), 10] <- "accurate mass match to database with very high false positive rate (>> 80%) so only use as a starting point for compound identification"
    IDedTable[which(IDedTable[, 1] == 5), 10] <- "no match"
    
    IDedTable[, 11] <- "Yes"
  }
  
  # timestamp()
  
  IDedTable <- IDedTable[, 1:(ncol_IDedTable - 5)]
  
  # Remove whitespace
  IDedTable[IDedTable == ""] <- NA
  IDedTable[IDedTable == "NA"] <- NA ## PJS 4/8/2025
  
  # Remove "Na" from formulas in column 4
  # IDedTable[, 4] <- gsub("NA", "", IDedTable[, 4]) ## PJS 4/8/2025 - this doesn't work, don't think it's necessary
  
  ####################################################################################
  
  if (Lipid == FALSE) {
    # Create a second output file with more info on B scores
    IDedTable_B <- IDedTable
    m <- matrix(0, nrow(IDedTable_B), 5)
    IDedTable_B <- cbind(IDedTable_B, m)
    nrow_IDedTable_B <- nrow(IDedTable_B)
    ncol_IDedTable_B <- ncol(IDedTable_B)
    colnames(IDedTable_B)[ncol_IDedTable_B - 4] <- "Check_Common_Fragment"
    colnames(IDedTable_B)[ncol_IDedTable_B - 3] <- "Number_of_F_Fragments"
    colnames(IDedTable_B)[ncol_IDedTable_B - 2] <- "Formula"
    colnames(IDedTable_B)[ncol_IDedTable_B - 1] <- "MD"
    colnames(IDedTable_B)[ncol_IDedTable_B] <- "Likely_PFAS"
    
    for (i in RowStartForFeatureTableData:nrow_IDedTable) {
      count <- 0
      frags_to_check <- c("C2F5", "C3F7", "SO2F", "SO3F", "CO2F", "C2F5O", "CF3O", "C2F3O2", "C3F7O", "C3O2F7", "CF3", "CH3FNSO2", "PO2F", "PO2F2", "SF5")
      for (j in frags_to_check) {
        if (length(grep(j, IDedTable_B[i, Frags_col])) > 0) {
          count <- count + 1
        }
      }
      IDedTable_B[i, ncol_IDedTable_B - 4] <- count
      if (!is.na(IDedTable_B[i, Frags_col + 1])) {
        IDedTable_B[i, ncol_IDedTable_B - 3] <- strsplit(IDedTable_B[i, Frags_col + 1], split="|", fixed=TRUE)[[1]][1]
      }
      if (is.na(IDedTable_B[i, 4])) { 
        IDedTable_B[i, ncol_IDedTable_B - 2] <- FALSE
      } else {
        IDedTable_B[i, ncol_IDedTable_B - 2] <- TRUE
      }
      
      md_col <- which(colnames(IDedTable_B) == "mass defect") ## PJS 4/7/2025 - replaced "Frags_col + 4" with md_col
      
      #*** JPK Edited 03/14/2025 with a tryCatch (was erroring in PosNeg run). Note this value should technically take the user input for mass defect filtering...
      tryCatch({
        if (as.numeric(IDedTable_B[i, md_col]) > -0.25 && as.numeric(IDedTable_B[i, md_col]) < 0.1) {
          IDedTable_B[i, ncol_IDedTable_B - 1] <- TRUE
        } else {
          IDedTable_B[i, ncol_IDedTable_B - 1] <- FALSE
        }
      },error=function(e){cat("Could not decipher mass defect value =",IDedTable_B[i, md_col],"Error Message:", conditionMessage(e), "\n")})
      
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B - 4]) > 0 || as.numeric(IDedTable_B[i, ncol_IDedTable_B - 3]) > 3) {
        IDedTable_B[i, ncol_IDedTable_B] <- TRUE
      } else {
        IDedTable_B[i, ncol_IDedTable_B] <- FALSE
      }
      
      # Assign B--
      if (IDedTable_B[i, ncol_IDedTable_B] == FALSE && (IDedTable_B[i, 1] == "B+" || IDedTable_B[i, 1] == "B" || IDedTable_B[i, 1] == "B-")) {
        IDedTable_B[i, 1] <- "B--"
        IDedTable_B[i, ncol_IDedTable_B - 6] <- "B--"
        IDedTable_B[i, 10] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 5 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 3/2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      }
      
      # Correct B score descriptions that have SMILES
      if (IDedTable_B[i, 1] == "B+" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS common fragment (F containing) and exact mass) and 2+ in homologous series. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B-" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, highly likely PFAS (1+ PFAS fragments from standards and exact mass OR 1+ common PFAS fragment (F containing), and exact mass). Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      } else if (IDedTable_B[i, 1] == "B--" && !is.na(IDedTable_B[i, 5])) {
        IDedTable_B[i, 10] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
        IDedTable_B[i, ncol_IDedTable_B - 5] <- "Lvl 3 Schymanski: Tentative ID, possible PFAS (1+ PFAS fragment (F containing) and exact mass). Note the fragment(s) observed may have other fragment formula assignments not containing F. Can be assigned a Lvl 2 after manual review of fragment evidence provided, homologous series evidence provided, EIC, and spectra"
      }
      
    }
    
    ###################### Final Sorting ###########################
    
    m <- matrix(0, nrow(IDedTable_B), 4)
    IDedTable_B <- cbind(IDedTable_B, m)
    ncol_IDedTable_B <- ncol(IDedTable_B) - 4
    hseries_count_cols <- hseries_count_cols + 9
    Mass_col <- Mass_col + 4
    
    # Find max series for each row
    for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
      IDedTable_B[i, ncol_IDedTable_B + 2] <- max(which(IDedTable_B[i, hseries_count_cols] == max(IDedTable_B[i, hseries_count_cols])))
    }
    
    # Max series series number and count
    m <- matrix(0, nrow_IDedTable_B, 1)
    IDedTable_B <- cbind(IDedTable_B, m)
    for (i in RowStartForFeatureTableData:nrow_IDedTable_B) {
      IDedTable_B[i, ncol_IDedTable_B + 3] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])] - 1])
      IDedTable_B[i, ncol_IDedTable_B + 4] <- as.numeric(IDedTable_B[i, hseries_count_cols[as.numeric(IDedTable_B[i, ncol_IDedTable_B + 2])]])
    }
    
    # Sort by max series count, max series, and series number
    IDedTable_B <- IDedTable_B[order(-as.numeric(IDedTable_B[, ncol_IDedTable_B + 4]), -as.numeric(IDedTable_B[, ncol_IDedTable_B + 2]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 3])),]
    
    # Create rank
    rank <- 1
    IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- rank
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 4]) == 0) {
        rank <- rank + 1
      } else if (IDedTable_B[i, ncol_IDedTable_B + 3] != IDedTable_B[i - 1, ncol_IDedTable_B + 3]) {
        rank <- rank + 1
      }
      IDedTable_B[i, ncol_IDedTable_B + 5] <- rank
    }
    
    # Score
    for (i in RowStartForFeatureTableData:nrow(IDedTable_B)) {
      if (IDedTable_B[i, ncol_IDedTable_B - 6] == "A" || IDedTable_B[i, ncol_IDedTable_B - 6] == "A-") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 1
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B+" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B" || IDedTable_B[i, ncol_IDedTable_B - 6] == "B-") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 2
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "B--") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 3
      } else if (IDedTable_B[i, ncol_IDedTable_B - 6] == "E") {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 5
      } else {
        IDedTable_B[i, ncol_IDedTable_B + 1] <- 4
      }
    }
    
    # Sort by rank and score
    IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, ncol_IDedTable_B + 1])),]
    
    # Adjust rank
    old <- IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]
    if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 1) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 10000000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 2) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 1000000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 3) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) - 100000
    } else if (IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 1] == 5) {
      IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[RowStartForFeatureTableData, ncol_IDedTable_B + 5]) + 1000000
    }
    for (i in (RowStartForFeatureTableData + 1):nrow(IDedTable_B)) {
      if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 5) {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) + 1000000
      } else if ((IDedTable_B[i, ncol_IDedTable_B + 5] != old)) {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 1) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 10000000
        } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 2) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 1000000
        } else if (as.numeric(IDedTable_B[i, ncol_IDedTable_B + 1]) == 3) {
          IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i, ncol_IDedTable_B + 5]) - 100000
        }
      } else {
        old <- IDedTable_B[i, ncol_IDedTable_B + 5]
        IDedTable_B[i, ncol_IDedTable_B + 5] <- as.numeric(IDedTable_B[i - 1, ncol_IDedTable_B + 5])
      }
    }
    
    # Sort by rank and m/z
    IDedTable_B <- IDedTable_B[order(as.numeric(IDedTable_B[, ncol_IDedTable_B + 5]), as.numeric(IDedTable_B[, Mass_col])),]
    
    IDedTable_B <- IDedTable_B[, 1:ncol_IDedTable_B]
  } else if (Lipid == TRUE) {
    m <- matrix(0, nrow(IDedTable), 1)
    IDedTable <- cbind(IDedTable, m)
    nrow_IDedTable <- nrow(IDedTable)
    ncol_IDedTable <- ncol(IDedTable)
    colnames(IDedTable)[ncol_IDedTable] <- "Number_of_F_Fragments"
    Num_Frags_col <- which(colnames(IDedTable) == "Num_Frags")
    if (length(Num_Frags_col) > 0) {
      IDedTable[, ncol_IDedTable] <- substr(IDedTable[, Num_Frags_col], 1, 1)
    }
  }
  
  # # Add convenient columns
  # m <- matrix("", nrow(IDedTable_B), 5)
  # IDedTable_B <- cbind(IDedTable_B, m)
  # ncol_IDedTable_B <- ncol(IDedTable_B)
  # colnames(IDedTable_B)[(ncol_IDedTable_B - 4):ncol_IDedTable_B] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
  # 
  # 
  # if (NegPos == "Neg") {
  #   # write.table(IDedTable, file.path(OutDir, "NegIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  #   write.table(IDedTable_B, file.path(OutDir, "NegIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # } else if (NegPos == "Pos") {
  #   # write.table(IDedTable, file.path(OutDir, "PosIDed_FIN_KMD_scored.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  #   write.table(IDedTable_B, file.path(OutDir, "PosIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # } else if (NegPos == "Combined") {
  #   write.table(IDedTable_B, file.path(OutDir, "CombinedIDed_FIN.csv"), sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA")
  # }
  
  # Add convenient columns
  m <- matrix("", nrow(IDedTable), 5)
  IDedTable <- cbind(IDedTable, m)
  ncol_IDedTable <- ncol(IDedTable)
  colnames(IDedTable)[(ncol_IDedTable - 4):ncol_IDedTable] <- c("Checked_Viz", "TRUE_Viz", "Class_Viz", "Cnumb_Viz", "Comment_Viz")
  
  
  OutDir_parts <- strsplit(IDedTable_dir, split = "_")[[1]]
  OutDir_parts <- paste(OutDir_parts[-length(OutDir_parts)], collapse = "_")
  OutDir_full <- paste0(OutDir_parts, "_FIN.csv")
  write.table(IDedTable, OutDir_full, sep=",", col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA", qmethod="double")
  
}
