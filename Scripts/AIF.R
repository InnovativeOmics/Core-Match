#!/usr/bin/env Rscript
source(paste(DirectoryScript,"EIC_MS1_AIF_QRAI_fns.R",sep=""))

printToCorFile <- function(arguments, df_FT, FT_row, mzFrag, corrFrag, CEsx, ce){
    ft_mz = df_FT[FT_row,1]
    ft_rt = df_FT[FT_row,2]
    ft_id = df_FT[FT_row,3]
    output = sprintf("%f,%f,%d,%f,%f,%d", ft_mz, ft_rt, ft_id, mzFrag, corrFrag, CEsx[ce])
    cat(output, file = arguments$fn_corr_output, append = TRUE, sep = "\n")
}

fragMatrix <- function(arguments, CEX_EICs, EIC_FT, MS2_Scan, df_FT, FT_row, CEsx){
    # this processes the EICs to get fragments above a correlation threshold
    # and adds those fragments to the matrix to print out in MS2


    M = matrix(nrow=0, ncol=4)
    colnames(M) = c("mz", "intensity", "CE", "correlation")

    if(arguments$UseAgg){
        EIC_FT = aggregate(EIC_FT, by = list(EIC_FT[,1])
                        , FUN = arguments$AggFUN)[,arguments$cols]
    }

    for(CEi in 1:length(CEX_EICs)){
        # print(CEi)
        for(EICi in 1:length(CEX_EICs[[CEi]])){
            c = getCorrelation(arguments, CEX_EICs[[CEi]], EICi, EIC_FT)
            if( !is.na(c) & (arguments$corThreshold < c) ){  ##### user input for threshhold
                printToCorFile(arguments, df_FT, FT_row, CEX_EICs[[CEi]][[EICi]][1,2], c, CEsx, CEi)
                # M = rbind(M, c(mz,intensity,4))
                mz = MS2_Scan[[CEi]][EICi,2]
                intensity = MS2_Scan[[CEi]][EICi,3]
                M = rbind( M, c( mz, intensity, CEi, c)  )
            }
        }
    }
    return(M)
}

ProcessMS2 <- function(arguments, CEs, RTs, trees, df_FT, EIC_FT, FT_row, scan_num, trt_CE0, num_eics, CEsx){
    if( length(EIC_FT) > arguments$Min_EIC_Len * arguments$ncols){

        MS2_Scan = to_list(for(CEi in 2:length(CEs)) genMS2_Scan(CEs, CEi, RTs, trees, trt_CE0))
        CEX_EICs = to_list(for(CEi in 2:length(CEs)) genCEXEICs(CEs, CEi, trees, MS2_Scan[[CEi-1]], arguments$mztol, arguments$rttol))

        for(CEi in 1:length(CEX_EICs)){
            num_eics = num_eics + length(CEX_EICs[[CEi]])
        }

        M = fragMatrix(arguments, CEX_EICs, EIC_FT, MS2_Scan, df_FT, FT_row, CEsx)
        printMS2data(arguments, M, df_FT, FT_row, scan_num)
    }
    return( num_eics )
}

calculateFragmentCorrelations <- function(arguments, CEs, RTs, trees, df_FT, CEsx){
    # generate the FT EICs
    # generate the MS2 Scan data at the target retention time
    # generate the Fragment Matrix
    # print MS2 data to file fn_output
    num_eics = 0
    for(FT_row in 1:nrow(df_FT)){
        # print("----------------")
        result = generateEICforFT(CEs, RTs, trees, df_FT, FT_row, arguments$mztol, arguments$rttol)
        trt_CE0 = result[[1]]
        EIC_FT = result[[2]]

        closest = which.min(abs(RTs[[1]] - trt_CE0))
        scan_num = CEs[[1]][closest,ncol(CEs[[1]])]
        # rowID_FT = df_FT[FT_row,3]
        if (FT_row %% 100 == 0){ 
          print(  paste(FT_row, length(EIC_FT), df_FT[FT_row, 1], sep=",")  )
        }
        num_eics = num_eics + 1
        num_eics = ProcessMS2(arguments, CEs, RTs, trees, df_FT, EIC_FT, FT_row, scan_num, trt_CE0, num_eics, CEsx)
    }
    print(num_eics)
}

printCorheader <- function(fn_output){
    output = print("ft_mz, ft_rt, ft_id, mzFrag, corrFrag, ce")
    cat(output, file = fn_output, sep = "\n")

}

printCSVs <- function(CEs, CEids){
    # Write the header for the csv file

    for(i in 1:length(CEids$uCE)){
        fn_csv_output = sprintf("EXAMPLE_CSV_OUTPUT_CE%d.csv",CEids$uCE[i])
        cat(paste("rt", "mz", "intensity", sep=","), file = fn_csv_output, append = FALSE, sep = "\n")
        write.table(CEs[[i]][,c("rt", "mz", "intensity")], file = fn_csv_output
        , row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    }
}

extract_AIF <- function( arguments ){
    # Convert an mzxml file with multiple collision energies to ms2 with the fragments
    arguments$fn_mzxml_wpath = file.path(arguments$path_to_input_folder, arguments$fn_mzxml)
    arguments$fn_corr_output = file.path(arguments$path_to_output_folder, arguments$fn_corr_output)
    arguments$fn_ms2_output =  file.path(arguments$path_to_output_folder, arguments$fn_ms2_output)
    arguments$fn_FT_wpath = file.path(arguments$path_to_input_folder, arguments$fn_FT)
    AIF <- openMSfile(arguments$fn_mzxml_wpath)
    print("loaded AllIons...")
    CEids = getCEids(AIF)
    print(CEids$uCE)
    CEsx = CEids$uCE[-1]

    df_FT <- readFeatureTable(arguments$fn_FT_wpath, arguments$FT_Cols, arguments$FT_SkipRows)
    printMS2header(arguments$fn_ms2_output, arguments$fn_mzxml)
    printCorheader(arguments$fn_corr_output)

    # print(df_FT)

    CEs = to_list(for(i in 1:length(CEids$uCE)) getAllSpectras(arguments, AIF, CEids$IDs[[i]]))
    RTs = to_list(for(i in 1:length(CEids$uCE)) unique(CEs[[i]][,1]))
    trees = to_list(for(i in 1:length(CEids$uCE)) createTree(CEs[[i]], treeType = "quad", dataType = "point", columns = 1:2))

    # printCSVs(CEs, CEids)

    calculateFragmentCorrelations(arguments, CEs, RTs, trees, df_FT, CEsx)
}