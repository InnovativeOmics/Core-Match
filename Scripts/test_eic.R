#!/usr/bin/env Rscript
source(paste("","./EIC_MS1_fns.R",  sep=""))
source(paste("","./genEIC.R",       sep=""))

test_eic <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Window = 1/6
        ,OutputDirectory = "../OUTPUT"
        ,FeatureID_Cols = c(6,7,12)+1
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "./secondary_isotopes.csv"
    )
    args$fn_FeatureID = "NegIDed_FIN.csv"
    args$fn_mzxml = "20230807_Target_Standard2_Neg.mzXML"
    args$fn_eic_output = "EXAMPLE_EIC_OUTPUT.csv"
    args$path_to_mzXML_Files = "../mzxml"
    extract_EICs(args)
}

test_eic()
