#!/usr/bin/env Rscript
source(paste("","../Scripts/EIC_MS1_fns.R",   sep=""))
source(paste("","../Scripts/genIsoTable.R",   sep=""))
source(paste("","../Scripts/MS1Spectragen.R", sep=""))
Rcpp::sourceCpp(paste("","../Scripts/MS1im.cpp", sep=""))

test_MS1 <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Window = 1/6
        ,OutputDirectory = "../../../datafiles/"
        ,FeatureID_Cols = c(6,7,12,13,4)
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "../Scripts/secondary_isotopes.csv"
    )
    args$dttol = c(0.3, 4.0)
    args$rttol = c(0.1, 1.0)
    args$cols = c("rt", "mz", "intensity", "dt")
    args$fn_FeatureID = "CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_FIN.csv"
    args$fn_mzxml = "CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
    args$fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv"
    args$path_to_mzXML_Files = "../../../datafiles/"
    extract_MS1(args)
}

test_MS1()