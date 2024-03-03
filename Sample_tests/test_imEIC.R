#!/usr/bin/env Rscript
source(paste("","../Scripts/EIC_MS1_fns.R",  sep=""))
source(paste("","../Scripts/genEIC.R",       sep=""))
Rcpp::sourceCpp(paste("","../Scripts/EICim.cpp",sep=""))

test_eic_im <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Window = 1/6
        ,OutputDirectory = "../../../datafiles/"
        ,FeatureID_Cols = c(6,7,12,13)
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "../Scripts/secondary_isotopes.csv"
    )
    args$dttol = c(0.3, 4.0)
    args$rttol = c(0.1, 1.0)
    args$cols = c("rt", "mz", "intensity", "dt")
    # args$colheader = c('Feature', 'RT', 'Intensity', 'mz', 'File', "DriftTime")
    args$colheader = c("Feature","MZ","RT","DT","Intensity,File")
    args$fn_FeatureID = "CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_FIN.csv"
    args$fn_mzxml = "CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
    args$fn_eic_output = "EXAMPLE_EIC_IM_OUTPUT.csv"
    args$path_to_mzXML_Files = "../../../datafiles/"
    extract_EICs(args)
}
test_eic_im()
