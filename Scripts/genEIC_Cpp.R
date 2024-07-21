#!/usr/bin/env Rscript

printEICs <- function(arguments, EIC_FeatureID, df_FeatureID, FeatureID_row){
    # rowsum(EIC_FeatureID, EIC_FeatureID[,1])
    # Total Ion Current vs Base Peak Intensity
    if(arguments$UseAgg){
        EIC_FeatureID = aggregate(EIC_FeatureID, by = list(EIC_FeatureID[,1])
                        , FUN = arguments$AggFUN)[,arguments$cols]
    }
    EIC_FeatureID[,4] = df_FeatureID[FeatureID_row,3]
    EIC_FeatureID[,5] = arguments$fn_mzxml
    EIC_FeatureID = cbind(EIC_FeatureID, 1)

    colnames(EIC_FeatureID) = c("RT", "MZ", "Intensity", "Feature", "File")
    write.table(EIC_FeatureID[,arguments$colheader]
        , file = arguments$fn_eic_output_wpath, row.names = FALSE, quote = FALSE
        , col.names = FALSE, append = TRUE, sep = ",")
}

calculateEICs <- function(arguments, CEs, RTs, trees, df_FeatureID){
    # generate the FeatureID EICs
    # write EIC data to file
    for(FeatureID_row in 1:nrow(df_FeatureID)){
        # print(FeatureID_row)
        result = generateEICforFeatureID(CEs, RTs, trees, df_FeatureID, FeatureID_row, arguments$mztol, arguments$rttol[2])
        EIC_FeatureID = result[[2]]

        # when Min_EIC_Len is 1, skips printing when there is only one row
        if (length(EIC_FeatureID) > length(arguments$cols) * arguments$Min_EIC_Len){
            printEICs(arguments, EIC_FeatureID, df_FeatureID, FeatureID_row)
        }
    }
}

calculateEICsIM <- function(arguments, CEs, df_FeatureID){
    # generate the FeatureID EICs
    # write EIC data to file
    M = CEs[order(CEs[,2]),]
    ft = df_FeatureID[order(df_FeatureID[,1]),]
    # ids = 1:50
    # print(paste(arguments$mztol, arguments$rttol, arguments$dttol
    #     , arguments$fn_eic_output_wpath, arguments$fn_mzxml, sep=","))

    # print(arguments$fndt)
    # print(arguments$fnrt)

    write_eics(
        M[,2], M[,1], M[,4], M[,3]
        , ft[,1], ft[,2], ft[,3], ft[,4]
        # , ft[ids,1], ft[ids,2], ft[ids,3], ft[ids,4]
        , arguments$mztol, arguments$rttol[1], arguments$dttol[2]
        , arguments$fndt, arguments$fn_mzxml, FALSE )

    write_eics(
        M[,2], M[,1], M[,4], M[,3]
        , ft[,1], ft[,2], ft[,3], ft[,4]
        # , ft[ids,1], ft[ids,2], ft[ids,3], ft[ids,4]
        , arguments$mztol, arguments$rttol[2], arguments$dttol[1]
        , arguments$fnrt, arguments$fn_mzxml, TRUE )

}

printCSVheader <- function(fn_csv_output){
    # Write the header for the csv file
    cat(paste("rt", "mz", "intensity", sep=","), file = fn_csv_output, append = FALSE, sep = "\n")
}

printEICheader <- function(fn_eic_output, colheader){
    # Write the header for the eic file
    cat(paste(colheader, collapse=","), file = fn_eic_output, append = FALSE, sep = "\n")
}

extract_EICs <- function( arguments ){
    # Convert an mzxml file with multiple collision energies to ms2 with the fragments
    arguments$fn_mzxml_wpath = file.path(arguments$path_to_mzXML_Files, arguments$fn_mzxml)
    arguments$fn_eic_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_eic_output)
    arguments$fn_FeatureID_wpath = file.path(arguments$path_to_output_folder, arguments$fn_FeatureID)

    AllIons <- openMSfile(arguments$fn_mzxml_wpath)  
    CEids = getCEids(AllIons)
    minid = which.min(CEids$uCE)

    df_FeatureID <- readFeatureTable(arguments$fn_FeatureID_wpath, arguments$FeatureID_Cols, arguments$FeatureID_SkipRows)
    class(df_FeatureID) <- "numeric"

    CEs = getAllSpectras(arguments, AllIons, CEids$IDs[[minid]])

    if("dt" %in% arguments$cols){
        arguments$fnrt = gsub('_IM_', '_IMrt_', arguments$fn_eic_output_wpath)
        if (!file.exists(arguments$fnrt)) {
            printEICheader(arguments$fnrt, arguments$colheader)
        }
        arguments$fndt = gsub('_IM_', '_IMdt_', arguments$fn_eic_output_wpath)
        if (!file.exists(arguments$fndt)) {
            printEICheader(arguments$fndt, arguments$colheader)
        }
        calculateEICsIM(arguments, CEs, df_FeatureID)
    } else{
        if (!file.exists(arguments$fn_eic_output_wpath)) {
            printEICheader(arguments$fn_eic_output_wpath, arguments$colheader)
        }
        CEs[,1] = round(CEs[,1],arguments$precision_rt) #retention time
        CEs[,2] = round(CEs[,2],arguments$precision_mz) #mass to charge ratio
        CEs[,3] = round(CEs[,3],arguments$precision_i)  #intensity
        RTs = unique(CEs[,1])
        trees = createTree(CEs, treeType = "quad", dataType = "point", columns = 1:2)
        calculateEICs(arguments, CEs, RTs, trees, df_FeatureID)
    }
}