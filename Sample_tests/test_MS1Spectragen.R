#!/usr/bin/env Rscript
source(paste("","../Scripts/EIC_MS1_fns.R",   sep=""))
source(paste("","../Scripts/genIsoTable.R",   sep=""))
source(paste("","../Scripts/MS1Spectragen.R", sep=""))
Rcpp::sourceCpp(paste("","../Scripts/MS1im.cpp", sep=""))

test_MS1 <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Window = 1/6
        ,OutputDirectory = "../RapidTestModular_FM/Output"
        ,FeatureID_Cols = c(6,7,12,4)+1 #mz, rt, rowID, formula
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "../Scripts/secondary_isotopes.csv"
    )
    args$dttol = 3
    args$cols = c("rt", "mz", "intensity")
    args$fn_FeatureID = "NegIDed_FIN.csv"
    args$fn_mzxml = "AFFF3_Target_Neg.mzXML"
    args$fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv"
    args$path_to_mzXML_Files = "../RapidTestModular_FM/Input"
    extract_MS1(args)
}

test_MS1()