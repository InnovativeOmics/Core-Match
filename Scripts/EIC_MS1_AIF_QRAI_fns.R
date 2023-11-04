# options(warn=2) # to check for errors/warnings

construct_EMAQ_arguments <- function(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Tolerance = RT_Tolerance
    ,OutputDirectory = OutputDirectory
    ,FT_Cols = FT_Cols
    ,GroupCSVDirectory = GroupCSVDirectory
){
    # Commented items need to be defined for those algorithms.
    arguments = list(
        #USER
            path_to_output_folder = OutputDirectory
            # ,fn_mzxml = ""
            # ,fn_FT = "NegIDed_FIN.csv"

        # ALL
            ,mztol = PrecursorMassAccuracy/2
            ,precision_mz = 5
            ,precision_rt = 2
            ,precision_i = 3
            ,min_intensity = 40
            ,use_min_i_filter = TRUE
            ,cols = c("rt", "mz", "intensity", "pMz", "snum")
            ,FT_SkipRows = c()

        # Some
            ,ncols = 5
            ,toggle_frag_row = TRUE
            ,FT_Cols = FT_Cols # c(6,7,12) #mz, rt, rowID
            ,FT_SkipRows = c() #c() means don't skip, c(1) means skip the first row
            ,Min_EIC_Len = 6
            ,corThreshold = 0.8
            ,AggFUN = max #min, mean, or max: sum will sum the mz and rt would need more code
            ,UseAgg = FALSE 
            # ,fn_ms2_output = "EXAMPLE_MS2_OUTPUT.ms2" # this is where the MS2 is stored using feature table

        # MS1
            ,mzZoomWindow = 5
            # ,FT_Cols = c(6,7,12,4) #mz, rt, rowID, formula
            # ,fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv" # this is where the MS1 is stored using feature table

        # Not MS1
            ,rttol = RT_Tolerance

        # EIC
            # ,fn_eic_output = "EXAMPLE_EIC_OUTPUT.csv" # this is where the EIC is stored using feature table

        # EIC AIF
            # ,fn_csv_output = "EXAMPLE_CSV_OUTPUT.csv" # this is converted from mzxml data

        # AIF
            # ,fn_corr_output = "EXAMPLE_COR_OUTPUT.csv" # this is where the MS2 is stored using feature table

        # QRAI
        #     ,fn_csv = "QRAI_Example/DIA_list.csv"
    )
    if (length(GroupCSVDirectory)>0) { 
        arguments$FT_SkipRows = c(1)
    }
    return(arguments)
}


readFeatureTable <- function(fn, columns, skip){
    # Read in Feature Table (Typically NegID) file with columns
    # Allow for skipping extra info column headers
    data <- fread(fn)
    cols <- data[, ..columns]
    if(length(skip) > 0){ 
        cols = cols[-skip,]
    }
    cols <- as.matrix(cols)
    return(cols)
}

getSpectra <- function(arguments, AllIons, PEAKS, data, IDs, rt, pMz, snum, i){
    # makes a mini matrix of the spectra
    IDi = IDs[i]
    p = PEAKS[[i]]
    if( arguments$use_min_i_filter ){
        p = p[ arguments$min_intensity < p[,2], ]
    }

    m = matrix(nrow = nrow(p), ncol = length(arguments$cols))
    m[,1] = rt[i]
    m[,2:3] = p
    m[,4] = pMz[i]
    m[,5] = snum[i]
    # If needed, set extra column data
    colnames(m) = arguments$cols
    return(m)
}

getAllSpectras <- function(arguments, AllIons, IDs, isMS1=FALSE){
    # creates matrix for all scans for a particular CE
    end = length(AllIons)
    PEAKS = lapply(seq_along(IDs), function(i) peaks(AllIons, IDs[i]))
    data = header(AllIons,IDs)
    rt = data$retentionTime/60 
    pMz = data$precursorMZ #get precursorMZ
    snum = data$acquisitionNum #get acquisition number (scan num)

    print(paste("length of CE processing: ",length(IDs)))
    MS = lapply(seq_along(IDs), function(i) getSpectra(arguments, AllIons, PEAKS, data, IDs, rt, pMz, snum, i))
    M = do.call(rbind, MS)
    if(isMS1){return(MS)}
    return(M)
}

getNearestRt <- function(rts, rt){
    # return the nearest retention time 
    closest = which.min(abs(rts - rt))
    return(rts[closest])
}

EIC_Tree <- function(tree, tmz, trt, mztol, rttol){
    # Return the indicies within a bounding box (mz, rt)
    v = rectLookup(tree, xlim = c(trt-rttol, trt+rttol)
                       , ylim = c(tmz-mztol, tmz+mztol))
    return(v)
}

generateEICforFT <- function(CEs, RTs, trees, df_FT, FT_row, mztol, rttol){
    # Generate the EIC for a row of FT
    tmz = df_FT[FT_row,1]
    trt = df_FT[FT_row,2]
    trt_CE0 = getNearestRt(RTs[[1]], trt)
    EIC_Indices = EIC_Tree(trees[[1]], tmz, trt, mztol, rttol)
    EIC_FT = CEs[[1]][EIC_Indices,]
    return( list(trt_CE0, EIC_FT) )
}

genCEXEICs <- function(CEs, CEi, trees, MS2_Scan, mztol, rttol){
    # Generate the EIC for each CEi mz fragment near target retention time
    trt_CEX = MS2_Scan[1,1]
    mz_fragments = MS2_Scan[,2]
    CEX_EICs = to_list(for(i in 1:length(mz_fragments)) CEs[[CEi]]
        [ EIC_Tree( trees[[CEi]], mz_fragments[i], trt_CEX, mztol, rttol), ])
    return(CEX_EICs)
}

getCEids <- function(DATA){
    # Get the unique collision energies and their indices
    # Set NA to 0 which corresponds to Full Scan
    data = header(DATA)
    ce = data$collisionEnergy
    ce[is.na(ce)] = 0
    uv = unique(ce)
    L = to_list(for(i in 1:length(uv)) which(ce == ce[i]))
    return(list( "uCE" = uv, "IDs" = L))
}

genMS2_Scan <- function(CEs, CEi, RTs, trees, trt_CE0){
    # get nearest retention time
    # return ms2 scan at that retention time (all rt, mz, intensity, pMZ)
    trt_CEX = getNearestRt(RTs[[CEi]], trt_CE0)
    highestmz = 1000000000 #an arbitrarily high value so that all mz are included
    # get max explicitly
    v = EIC_Tree(trees[[CEi]], highestmz, trt_CEX, highestmz, 0)
    return( CEs[[CEi]][v,] )
}

getCorrelation <- function(arguments, CEX_EICs, EICi, EIC_FT){
    # intensities at the FT retention times
    s = arguments$Min_EIC_Len # the minimum number of points in the FT EIC
    b1 = s < length( CEX_EICs[[EICi]] ) / ncol(EIC_FT)
    b2 = s < nrow( EIC_FT ) 
    b3 = b1 & b2
    if(b3){
        if(arguments$UseAgg){
            CEX_EICs[[EICi]] = aggregate(CEX_EICs[[EICi]], by = list(CEX_EICs[[EICi]][,1])
                        , FUN = arguments$AggFUN)[,arguments$cols]
        }

        x = CEX_EICs[[EICi]][,1] # rt
        y = CEX_EICs[[EICi]][,3] # intensity
        if(!arguments$UseAgg){
            dups = duplicated(x)
            x = x[!dups]
            y = y[!dups]
        }
        EIC_frag_intensity = approx(x, y, xout=EIC_FT[,1])$y
        return( cor(EIC_FT[,3], EIC_frag_intensity) )
    }
    #-2 indicates an error since correlation is between -1 and 1
    # This allows it to continue without crashing, but might need addressing
    return(-2) 
}

printMS2data <- function(arguments, M, df_FT, FT_row, scan_num){
    # Print MS2 data
    # BPI: Base peak intensity
    # BPM: Base peak mass
    # TIC: Total intensity
    # Z: Charge
    fn_output = arguments$fn_ms2_output
    nr = nrow(M)
    if( 0 < nr ){
        S11mz = paste("S", scan_num, scan_num, df_FT[FT_row,1], sep="\t")
        cat(S11mz, file = arguments$fn_ms2_output, append = TRUE, sep = "\n")
     
        # This orders by the mz
        M[1:nr,] = M[order(M[,1], decreasing=FALSE),]
        idx = which.max(M[,2]) # Max mz index
        BPI = M[idx,2]
        BPM = M[idx,1]
        TIC = sum(M[,2])
        cat(paste("I\tRTime  ", df_FT[FT_row,2]), file = fn_output, append = TRUE, sep = "\n")
        cat(paste("I\tBPI", BPI, sep="\t"), file = fn_output, append = TRUE, sep = "\n")
        cat(paste("I\tBPM", BPM, sep="\t"), file = fn_output, append = TRUE, sep = "\n")
        cat(paste("I\tTIC", TIC, sep="\t"), file = fn_output, append = TRUE, sep = "\n")
        cat(paste("Z\t1", df_FT[FT_row,1], sep='\t'), file = fn_output, append = TRUE, sep = "\n")
        if(nr == 1 & arguments$toggle_frag_row){
            cat(paste(0.1,0.1), file = fn_output, append = TRUE, sep = "\n")
        }
        for(i in 1:nr){
            mz = format(round(M[i,1], arguments$precision_mz), nsmall = arguments$precision_mz)
            intensity = format(round(M[i,2], arguments$precision_i), nsmall = arguments$precision_i)
            cat(paste(mz, intensity), file = fn_output, append = TRUE, sep = "\n")
        }
    }
}
printMS2header <- function(fn_output, fn_mzxml){
    # Write the header for the ms2 file
    cat(paste("H	CreationDate ", format(Sys.time(), "%a %b %d %H:%M:%S %Y")), file = fn_output, sep = "\n")
    cat("H	Extractor	ProteoWizard", file = fn_output, append = TRUE, sep = "\n")
    cat("H	Extractor version	Xcalibur", file = fn_output, append = TRUE, sep = "\n")
    cat(paste("H	Source file	", fn_mzxml), file = fn_output, append = TRUE, sep = "\n")
}

runallEMAQ <- function(Target_mzXML, args, runmode = TRUE, isNeg = TRUE){
    if(isNeg) { prefix = "Neg"; pattern = "_Neg|_neg|_NEG" } 
    else      { prefix = "Pos"; pattern = "_Pos|_pos|_POS" }
    suppressWarnings(if (runmode) {
        args$fn_FT         <- paste(prefix, "IDed_FIN.csv",      sep="")
        args$fn_eic_output <- paste(prefix, "_Feature_EICs.csv", sep="")
        args$fn_MS1_output <- paste(prefix, "_Feature_MS1s.csv", sep="")
        NPTarget_mzXML <- grep(Target_mzXML, pattern = pattern, invert=FALSE, value=TRUE)
        for (i in 1:length(NPTarget_mzXML)) {
            args$fn_mzxml = NPTarget_mzXML[i] #CHECK
            extract_EICs(args)
            extract_MS1(args)
        }
    })
}

# usage

# arguments = construct_EMAQ_arguments <- function(
#     PrecursorMassAccuracy = PrecursorMassAccuracy
#     ,RT_Tolerance = RT_Tolerance
#     ,OutputDirectory = OutputDirectory
#     ,FT_Cols = c(7,8,1,5)
#     ,GroupCSVDirectory = GroupCSVDirectory
# )

# runallEMAQ(Target_mzXML, arguments, runmode = runNegddMS)
# runallEMAQ(Target_mzXML, arguments, runmode = runPosddMS, isNeg = FALSE)
