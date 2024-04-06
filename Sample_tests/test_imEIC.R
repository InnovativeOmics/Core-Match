#!/usr/bin/env Rscript
source(paste("","../Scripts/EIC_MS1_fns.R",  sep=""))
source(paste("","../Scripts/genEIC.R",       sep=""))
Rcpp::sourceCpp(paste("","../Scripts/EICim.cpp",sep=""))

test_eic_im <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Tolerances = c(0.1, 0.5)
        ,DT_Tolerances = c(0.1, 3)
        ,OutputDirectory = "../../../datafiles/"
        ,FeatureID_Cols = c(6,7,12,13)
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "secondary_isotopes.csv"
        ,isIM = TRUE
    )
    args$fn_FeatureID = "CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_FIN.csv"
    args$fn_mzxml = "CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
    args$fn_eic_output = "EXAMPLE_EIC_IM_OUTPUT_NEW-modified.csv"
    args$path_to_mzXML_Files = "../../../datafiles/"
    extract_EICs(args)
}
test_eic_im()
