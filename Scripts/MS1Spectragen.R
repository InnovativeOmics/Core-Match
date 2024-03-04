#!/usr/bin/env Rscript

dfIsotopes <- function(arguments){
    # Get isotopes from NIST data based on the isotope string
    df = genIsoTable(arguments)
    iso = list(df$RelativeMass, df$Symbol, df$IsotopicComposition)
    isoIndices = order( iso[[1]] )
    iso[[1]] = iso[[1]][isoIndices]
    iso[[2]] = iso[[2]][isoIndices]
    iso[[3]] = iso[[3]][isoIndices]
    return( iso )
}

printMS1header <- function(fn_MS1_output){
    # Write the header for the MS1 file
    cat(paste("Feature","m/z","Intensity","Isotope","Formula","PredAbundance"
        ,"File","Zoom", sep=","), file = fn_MS1_output, append = FALSE, sep = "\n")
}

extract_MS1 <- function( arguments ){
    # Convert an mzxml file with multiple collision energies to ms2 with the fragments
    arguments$fn_mzxml_wpath = file.path(arguments$path_to_mzXML_Files, arguments$fn_mzxml)
    AllIons <- openMSfile(arguments$fn_mzxml_wpath)
    CEids = getCEids(AllIons)
    minid = which.min(CEids$uCE)

    arguments$fn_FeatureID_wpath = file.path(arguments$path_to_output_folder, arguments$fn_FeatureID)
    df_FeatureID <- readFeatureTable(arguments$fn_FeatureID_wpath, arguments$FeatureID_Cols, arguments$FeatureID_SkipRows)
    df_FeatureID_mzrtid = df_FeatureID[,1:(ncol(df_FeatureID)-1)]
    df_FeatureID_formula = df_FeatureID[,ncol(df_FeatureID)]
    class(df_FeatureID_mzrtid) <- "numeric"

    mzrtdt_ids = order(df_FeatureID_mzrtid[,2])
    df_FeatureID_mzrtid = df_FeatureID_mzrtid[mzrtdt_ids,]
    df_FeatureID_formula = df_FeatureID_formula[mzrtdt_ids]
    df_FeatureID = list(df_FeatureID_mzrtid, df_FeatureID_formula)

    arguments$fn_MS1_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_MS1_output)
    if (!file.exists(arguments$fn_MS1_output_wpath)) {
        printMS1header(arguments$fn_MS1_output_wpath)
    }

    arguments$Iso = dfIsotopes(arguments)
    arguments$mzZoomWindowLow  = -2
    arguments$mzZoomWindowHigh = max(arguments$Iso[[1]]) + 1

    MS1s = getAllSpectras(arguments, AllIons, CEids$IDs[[minid]], isMS1=TRUE)
    RTs = sort( unique( header(AllIons,CEids$IDs[[minid]])$retentionTime/60 ) )
    RTw = (RTs[2] - RTs[1])/2

    hasdt = FALSE
    DTw = 0
    if("dt" %in% arguments$cols){
        hasdt = TRUE
        DTs = sort( unique( header(AllIons,CEids$IDs[[minid]])$ionMobilityDriftTime ) )
        DTw = (DTs[2] - DTs[1])/2
    }

    save_MS1s_to_file(MS1s, df_FeatureID, arguments$Iso, RTw, DTw
    , arguments$fn_MS1_output_wpath, arguments$fn_MS1_output
    , arguments$mzZoomWindowLow, arguments$mzZoomWindowHigh
    , arguments$mztol, hasdt
    )

}
