#!/usr/bin/env Rscript

# MZcol <- 3 # where is your m/z column in the feature table?
# RTcol <- 4 # where is your retention time column in the feature table?
# rowID <- 1 # where is your column containing numeric identifiers in the feature table?
# RT_Tolerance <- .11 #Tolerance in minutes for how far from the peak center to include scans for correlation. So .2 minutes would be +/- 0.2 minutes
# PrecursorMassAccuracy <- 0.01 #Window in Daltons for full-scan data (and All-Ions MS2 level data). So 0.01 Da would be +/- 0.005 Da
# corThreshold <- 0.7 #Correlation threshold for determining if a fragment belongs to a precursor ion based on co-eluting EIC profiles
#                     #(EIC profiles should have the same shape between fragments and precursors)
# Min_EIC_Len <- 6 #Minimum number of scans for full-scan and fragment EICs in order to proceed with correlation and retaining fragments
# min_intensity <- 40 #ignore intensities under this threshold
#
# library(data.table)
# library(SearchTrees)
# library(comprehenr)
# library(mzR)

# InputDirectory <- "../../../4_Jeremy_AllIons/Explorer list/Output/"

# source(paste("","./AIF.R", sep=""))


# FeatureTables = c("Unique to A cpds group.csv", "Unique to B cpds group.csv", "Unique to C cpds group.csv")
# MZMLs = c("AllIons NIST A.mzML", "AllIons NIST B.mzML", "AllIons NIST C.mzML")
# MZML_folders = c("../../../4_Jeremy_AllIons/mzML/", "../../../AllIons mzML Converted files from Qual/")

IonDecon_AIF <- function(
    InputDirectory, FeatureTable, MZMLs,
    MZcol, RTcol, rowID, 
    RT_Tolerance, PrecursorMassAccuracy, 
    corThreshold, Min_EIC_Len, min_intensity ){
  
  args = construct_EMAQ_arguments(
      PrecursorMassAccuracy = PrecursorMassAccuracy
      ,RT_Tolerance = RT_Tolerance
      ,OutputDirectory = InputDirectory
      ,FT_Cols = c(MZcol,RTcol,rowID)
      ,GroupCSVDirectory = c(1)
  )
  args$min_intensity = min_intensity
  args$corThreshold = corThreshold
  args$UseAgg = FALSE
  args$Min_EIC_Len = Min_EIC_Len #min Min, is joke
  
  args$fn_FT = FeatureTable
  if(length(MZMLs)>0){
    for(i in 1:length(MZMLs)){
      args$fn_mzml = MZMLs[i]
      mzML_Base <- gsub('.mzML','', args$fn_mzml, ignore.case = TRUE)
      mzML_Base <- gsub('All_Ions','', mzML_Base, ignore.case = TRUE)
      args$fn_ms2_output = paste("ddMS2_",mzML_Base,".ms2",sep="")
      args$fn_corr_output = paste("Cor_",mzML_Base,".csv",sep="")
      extract_AIF(args)
    }
  }
}

# IonDecon_AIF()