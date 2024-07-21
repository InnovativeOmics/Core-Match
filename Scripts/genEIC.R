#!/usr/bin/env Rscript

printEICs <- function(arguments, EIC_FeatureID, df_FeatureID, FeatureID_row){
    # when Min_EIC_Len is 1, skips printing when there is only one row
    if (length(EIC_FeatureID) > length(arguments$cols) * arguments$Min_EIC_Len){
        # rowsum(EIC_FeatureID, EIC_FeatureID[,1])
        # Total Ion Current vs Base Peak Intensity
        if(arguments$UseAgg){
            EIC_FeatureID = aggregate(EIC_FeatureID, by = list(EIC_FeatureID[,1])
                            , FUN = arguments$AggFUN)[,arguments$cols]
        }
        EIC_FeatureID[,4] = df_FeatureID[FeatureID_row,3]
        EIC_FeatureID[,5] = arguments$fn_mzxml
        EIC_FeatureID = cbind(EIC_FeatureID, 1)

        colnames(EIC_FeatureID) = c("RT", "mz", "Intensity", "Feature", "File", "Zoom")
        write.table(EIC_FeatureID[,c("Feature", "RT", "Intensity", "mz", "File", "Zoom")]
            , file = arguments$fn_eic_output_wpath, row.names = FALSE, quote = FALSE
            , col.names = FALSE, append = TRUE, sep = ",")
    }
}

calculateEICs <- function(arguments, CEs, RTs, trees, df_FeatureID){
    # generate the FeatureID EICs
    # write EIC data to file
    for(FeatureID_row in 1:nrow(df_FeatureID)){
        # print(FeatureID_row)
        result = generateEICforFeatureID(CEs, RTs, trees, df_FeatureID, FeatureID_row, arguments$mztol, arguments$rttol)
        EIC_FeatureID = result[[2]]

        printEICs(arguments, EIC_FeatureID, df_FeatureID, FeatureID_row)
    }
}

printCSVheader <- function(fn_csv_output){
    # Write the header for the csv file
    cat(paste("rt", "mz", "intensity", sep=","), file = fn_csv_output, append = FALSE, sep = "\n")
}

printEICheader <- function(fn_eic_output){
    # Write the header for the eic file
    cat(paste('Feature', 'RT', 'Intensity', 'mz', 'File', 'Zoom', sep=","), file = fn_eic_output, append = FALSE, sep = "\n")
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

    if (!file.exists(arguments$fn_eic_output_wpath)) {
        printEICheader(arguments$fn_eic_output_wpath)
    }

    # CEs = to_list(for(i in 1:length(CEids$uCE)) getAllSpectras(arguments, AllIons, CEids$IDs[[i]]) )
    # RTs = to_list(for(i in 1:length(CEids$uCE)) unique(CEs[[i]][,1]))
    # trees = to_list(for(i in 1:length(CEids$uCE)) createTree(CEs[[i]], treeType = "quad", dataType = "point", columns = 1:2))

    CEs = getAllSpectras(arguments, AllIons, CEids$IDs[[minid]])
    RTs = unique(CEs[,1])
    trees = createTree(CEs, treeType = "quad", dataType = "point", columns = 1:2)


    CEs[,1] = round(CEs[,1],arguments$precision_rt) #retention time
    CEs[,2] = round(CEs[,2],arguments$precision_mz) #mass to charge ratio
    CEs[,3] = round(CEs[,3],arguments$precision_i)  #intensity

    calculateEICs(arguments, CEs, RTs, trees, df_FeatureID)
}