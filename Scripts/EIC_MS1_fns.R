library(data.table)
library(SearchTrees)
library(comprehenr)
library(mzR)

# options(warn=2) # to check for errors/warnings

construct_EM_arguments <- function(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Window = RT_Window
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = FeatureID_Cols
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = isotable 
){
    # Commented items need to be defined for those algorithms.
    arguments = list(
        #USER
            path_to_output_folder = OutputDirectory
            # ,fn_mzxml = ""
            # ,fn_FeatureID = "NegIDed_FIN.csv"
            # ,isostring = "13C5;N;S;Cl4;18O;Br2"
            ,isostring = isostring
            ,isotable = isotable

        # ALL
            ,mztol = PrecursorMassAccuracy/2
            ,precision_mz = 5
            ,precision_rt = 2
            ,precision_i = 3
            ,min_intensity = 1
            ,use_min_i_filter = FALSE
            ,cols = c("rt", "mz", "intensity", "pMz", "snum")
            ,FeatureID_SkipRows = c()

        # Some
            ,ncols = 5
            ,toggle_frag_row = TRUE
            ,FeatureID_Cols = FeatureID_Cols # c(6,7,12) #mz, rt, rowID
            ,FeatureID_SkipRows = c() #c() means don't skip, c(1) means skip the first row
            ,Min_EIC_Len = 5
            ,AggFUN = max #min, mean, or max: sum will sum the mz and rt would need more code
            ,UseAgg = TRUE 
            # ,fn_ms2_output = "EXAMPLE_MS2_OUTPUT.ms2" # this is where the MS2 is stored using feature table

        # MS1
            # ,mzZoomWindow = 5 #now based on max RelativeMass
            # ,FeatureID_Cols = c(6,7,12,4) #mz, rt, rowID, formula
            # ,fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv" # this is where the MS1 is stored using feature table

        # EIC
            ,rttol = (RT_Window/2) * 6
            # ,fn_eic_output = "EXAMPLE_EIC_OUTPUT.csv" # this is where the EIC is stored using feature table
            # ,fn_csv_output = "EXAMPLE_CSV_OUTPUT.csv" # this is converted from mzxml data
    )
    if (length(GroupCSVDirectory)>0) { 
        arguments$FeatureID_SkipRows = c(1)
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

getSpectra <- function(arguments, PEAKS, IDs, rt, pMz, snum, imdt, i){
    # makes a mini matrix of the spectra
    IDi = IDs[i]
    p = PEAKS[[i]]
    if( arguments$use_min_i_filter ){
        p = p[ arguments$min_intensity < p[,2], ]
    }

    m = matrix(nrow = nrow(p), ncol = length(arguments$cols))
    m[,1] = rt[i]
    m[,2:3] = p
    # m[,4] = pMz[i]
    # m[,5] = snum[i]
    if("dt" %in% arguments$cols){ m[,4] = imdt[i] } # Ion Mobility Drift Time
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
    if("dt" %in% arguments$cols){ 
        # Ion Mobility Drift Time
        imdt = data$ionMobilityDriftTime
    } else{
        imdt = NULL
    }

    print(paste("length of CE processing: ",length(IDs)))
    MS = lapply(seq_along(IDs), function(i) getSpectra(arguments, PEAKS, IDs, rt, pMz, snum, imdt, i))
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

generateEICforFeatureID <- function(CEs, RTs, trees, df_FeatureID, FeatureID_row, mztol, rttol){
    # Generate the EIC for a row of FeatureID
    tmz = df_FeatureID[FeatureID_row,1]
    trt = df_FeatureID[FeatureID_row,2]
    trt_CE0 = getNearestRt(RTs, trt)
    EIC_Indices = EIC_Tree(trees, tmz, trt, mztol, rttol)
    EIC_FeatureID = CEs[EIC_Indices,]
    return( list(trt_CE0, EIC_FeatureID) )
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

runallEM <- function(Target_mzXML, arguments, runmode = TRUE, isNeg = TRUE){
    if(isNeg) { prefix = "Neg"; pattern = "_Neg|_neg|_NEG" } 
    else      { prefix = "Pos"; pattern = "_Pos|_pos|_POS" }
    suppressWarnings(if (runmode) {
      arguments$fn_FeatureID <- paste(prefix, "IDed_FIN.csv", sep="")
      arguments$fn_eic_output <- paste(prefix, "_Feature_EICs.csv", sep="")
      arguments$fn_MS1_output <- paste(prefix, "_Feature_MS1s.csv", sep="")
        NPTarget_mzXML <- grep(Target_mzXML, pattern = pattern, invert=FALSE, value=TRUE)
        for (i in 1:length(NPTarget_mzXML)) {
            arguments$fn_mzxml = NPTarget_mzXML[i]
            extract_EICs(arguments)
            extract_MS1(arguments)
        }
        append_percentages(arguments)
    })
}
