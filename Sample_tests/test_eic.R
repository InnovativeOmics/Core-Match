#!/usr/bin/env Rscript
source(paste("","EIC_MS1_fns.R",  sep=""))
source(paste("","genEIC.R",       sep=""))

test_eic <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Tolerances = c(0.1, 0.5)
        ,DT_Tolerances = c(0.1, 3)
        ,OutputDirectory = "../../Output"
        ,FeatureID_Cols = c(6,7,12)+1
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "secondary_isotopes.csv"
        ,isIM = FALSE
    )
    args$fn_FeatureID = "NegIDed_FIN.csv"
    args$fn_mzxml = "AFFF3_Target_Neg.mzXML"
    args$fn_eic_output = "EXAMPLE_EIC_OUTPUT_NEW.csv"
    args$path_to_mzXML_Files = "../../INPUT"
    extract_EICs(args)
}

test_eic()
