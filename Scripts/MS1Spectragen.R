#!/usr/bin/env Rscript

generateMS1forFeatureID <- function(arguments, MS1s, RTs, df_FeatureID, FeatureID_row){
    # Generate the MS1 for a row of FeatureID
    mzZoomLow  = arguments$mzZoomWindowLow
    mzZoomHigh = arguments$mzZoomWindowHigh
    tmz = df_FeatureID[[1]][FeatureID_row,1]
    trt = df_FeatureID[[1]][FeatureID_row,2]
    MS1index = which.min(abs(RTs - trt))
    df = MS1s[[MS1index]][,-c(4,5)] #remove the extra pMZ and snum columns from the shared getAllSpectras function
    MS1_Zoom = df[ tmz + mzZoomLow < df[,2]
                  & df[,2] < tmz + mzZoomHigh,, drop = FALSE]

    # print( sprintf("%d, %d", FeatureID_row, nrow(MS1_Zoom)) )
    if(nrow(MS1_Zoom) > 0){
        mz_index = which.min(abs(MS1_Zoom[,2] - tmz))
        # this is just the closest, it's not necessarily close
        # maybe set a flag with mztol?
        if(mz_index < nrow(MS1_Zoom) ){
            dmz_scan = MS1_Zoom[(mz_index+1):nrow(MS1_Zoom),2] - MS1_Zoom[mz_index,2]
        }
    }

    empty_cols <- matrix("", nrow = nrow(MS1_Zoom), ncol = 3)
    MS1_Zoom <- cbind(MS1_Zoom, empty_cols)

    rowid_col <- rep(df_FeatureID[[1]][FeatureID_row,3], nrow(MS1_Zoom))
    MS1_Zoom = cbind(MS1_Zoom, rowid_col)

    fn_col <- rep(arguments$fn_mzxml, nrow(MS1_Zoom))
    MS1_Zoom = cbind(MS1_Zoom, fn_col)

    zoom_col <- rep(1, nrow(MS1_Zoom))
    MS1_Zoom = cbind(MS1_Zoom, zoom_col)

    if(nrow(MS1_Zoom) > 0){
        MS1_Zoom[mz_index,4] = "M"
        f = df_FeatureID[[2]][[FeatureID_row]]
        MS1_Zoom[mz_index,5] = f

        if(mz_index < nrow(MS1_Zoom)){
            MS1_Zoom[(mz_index+1):(mz_index+length(dmz_scan)), 4] = get_iso_strings(arguments, dmz_scan)
        }
    }
    return( MS1_Zoom )
}

get_iso_strings <- function(arguments, dmz_scan){
    # Generate the iso string for each of the rows.
    IsoStrings = c()
    for (dmz_i in 1:length(dmz_scan)){
        isoProx = abs(arguments$Iso[[1]] - dmz_scan[dmz_i])
        isoMatches = which(isoProx < arguments$mztol)
        s = sort( isoProx[isoMatches], index.return = TRUE )
        isoMatchesSorted = isoMatches[ s$ix ]
        IsoStrings = append(IsoStrings, iso_string_genlist(arguments, isoProx, isoMatchesSorted) )
    }
    return( IsoStrings )
}

iso_string_genlist <- function(arguments, isoProx, isoMatches){
    # Generate the iso string for a single row.
    # Example output
    # 81Br(97.51%;0.00174Da);37Cl(32.40%;0.00264Da);34S(4.43%;0.00389Da);18O(0.20%;0.00455Da)
    IsoStrings = c()
    for (match_i in isoMatches){
        IsoStrings = append( IsoStrings, iso_string_gen(arguments$Iso, isoProx, match_i) )
    }
    IsoString = paste(IsoStrings, collapse="_")
    return(IsoString)
}

iso_string_gen <- function(isotopes, isoProx, match_i){
    # For a particular index, compose the isotopic string.
    return( sprintf("%s(%.2f%%;%.5fDa)", isotopes[[2]][[match_i]], isotopes[[3]][[match_i]], isoProx[[match_i]]) )
}

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

printMS1s <- function(arguments, MS1_FeatureID, df_FeatureID, FeatureID_row){
    # Print the data to the MS1 output file.
    colnames(MS1_FeatureID) = c("RT", "mz", "Intensity", "Isotope","Formula","PredAbundance", "Feature", "File", "Zoom")
    write.table(MS1_FeatureID[,c("Feature", "mz", "Intensity", "Isotope","Formula","PredAbundance", "File", "Zoom"),drop=FALSE]
        , file = arguments$fn_MS1_output_wpath, row.names = FALSE, quote = FALSE
        , col.names = FALSE, append = TRUE, sep = ",")
}

calculateMS1s <- function(arguments, MS1s, RTs, df_FeatureID){
    # Generate the FeatureID MS1s
    # Write MS1 data to file
    for(FeatureID_row in 1:nrow(df_FeatureID[[1]])){
        MS1_FeatureID = generateMS1forFeatureID(arguments, MS1s, RTs, df_FeatureID, FeatureID_row)
        printMS1s(arguments, MS1_FeatureID, df_FeatureID, FeatureID_row)
    }

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
    print(df_FeatureID_mzrtid[1:10,])

    arguments$fn_MS1_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_MS1_output)
    if (!file.exists(arguments$fn_MS1_output_wpath)) {
        printMS1header(arguments$fn_MS1_output_wpath)
    }

    arguments$Iso = dfIsotopes(arguments)
    arguments$mzZoomWindowLow  = -2
    arguments$mzZoomWindowHigh = max(arguments$Iso[[1]]) + 1
    # print(arguments$Iso)

    # MS1s = getAllSpectras(arguments, AllIons, CEids$IDs[[minid]], isMS1=TRUE)
    MS1s = getAllSpectras(arguments, AllIons, CEids$IDs[[minid]][1:1000], isMS1=TRUE)
    RTs = sort( unique( header(AllIons,CEids$IDs[[minid]])$retentionTime/60 ) )
    RTw = (RTs[2] - RTs[1])/2
    DTs = sort( unique( header(AllIons,CEids$IDs[[minid]])$ionMobilityDriftTime ) )
    DTw = (DTs[2] - DTs[1])/2

    # print(MS1s[[647]])

    hasdt = if("dt" %in% arguments$cols) TRUE else FALSE
    save_MS1s_to_file(MS1s, df_FeatureID, arguments$Iso, RTw, DTw
    , arguments$fn_MS1_output_wpath, arguments$fn_MS1_output
    , arguments$mzZoomWindowLow, arguments$mzZoomWindowHigh
    , arguments$mztol, hasdt
    )

    # calculateMS1s(arguments, MS1s, RTs, df_FeatureID)
}
