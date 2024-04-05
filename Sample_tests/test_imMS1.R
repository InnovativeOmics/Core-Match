#!/usr/bin/env Rscript
source(paste("","../Scripts/EIC_MS1_fns.R",   sep=""))
source(paste("","../Scripts/genIsoTable.R",   sep=""))
source(paste("","../Scripts/MS1Spectragen.R", sep=""))
Rcpp::sourceCpp(paste("","../Scripts/MS1im.cpp", sep=""))

test_MS1 <- function(){
    args = construct_EM_arguments(
        PrecursorMassAccuracy = 0.01
        ,RT_Tolerances = c(0.1, 0.5)
        ,DT_Tolerances = c(0.1, 3)
        ,OutputDirectory = "../../../datafiles/"
        ,FeatureID_Cols = c(6,7,12,13,4)
        ,GroupCSVDirectory = c()
        ,isostring = "13C3;N;S;Cl2;18O;Br2"
        ,isotable = "../Scripts/secondary_isotopes.csv"
        ,isIM = TRUE
    )
    args$fn_FeatureID = "CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_FIN.csv"
    # args$fn_FeatureID = "CM_PFAS_Level6_IM4_bit_DI3_d_DeMP_FIN_small.csv"

    fn_mzxmls = list("CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
        , "NoParameters_Uncompressed_CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
        , "PeakPicking_CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
        , "PeakPicking_CWT_CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
        , "PeakPickingUncompressed_CM PFAS Level 6 IM 4 bit_DI3.d.DeMP.mzML"
    )

    # args$fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv"
    for(i in 1:1){
    # for(i in 1:length(fn_mzxmls)){
        args$path_to_mzXML_Files = if(i == 1) "../../../datafiles/" else "../../../datafiles/newones"
        args$fn_mzxml = fn_mzxmls[[i]]
        args$fn_MS1_output = paste("MS1_Ion_Mobility_", args$fn_mzxml, "NEW-modified.csv", sep="")
        print(args$fn_mzxml)
        print(args$fn_MS1_output)
        extract_MS1(args)
    }
}

test_MS1()