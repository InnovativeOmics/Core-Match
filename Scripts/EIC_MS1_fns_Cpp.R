library(data.table)
library(SearchTrees)
library(comprehenr)
library(mzR)

# options(warn=2) # to check for errors/warnings

construct_EM_arguments <- function(
    PrecursorMassAccuracy = PrecursorMassAccuracy
    ,RT_Tolerances = RT_Tolerances
    ,DT_Tolerances = DT_Tolerances
    ,OutputDirectory = OutputDirectory
    ,FeatureID_Cols = FeatureID_Cols
    ,GroupCSVDirectory = GroupCSVDirectory
    ,isostring = ISOstring
    ,isotable = isotable 
    ,isIM = isIM
    ,path_to_mzXML_Files
){
    # Commented items need to be defined for those algorithms.
    arguments = list(
        #USER
      path_to_mzXML_Files = path_to_mzXML_Files
            ,path_to_output_folder = OutputDirectory
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
            ,FeatureID_SkipRows = c()
            ,rttol = RT_Tolerances
            ,dttol = DT_Tolerances

        # Some=
            ,toggle_frag_row = TRUE
            ,FeatureID_Cols = FeatureID_Cols # c(6,7,12) #mz, rt, rowID
            ,FeatureID_SkipRows = c() #c() means don't skip, c(1) means skip the first row
            ,Min_EIC_Len = 5
            ,AggFUN = max #min, mean, or max: sum will sum the mz and rt would need more code
            ,UseAgg = TRUE 
            ,isIM = isIM
    )
    if (length(GroupCSVDirectory)>0) { 
        arguments$FeatureID_SkipRows = c(1)
    }
    if(isIM){
        arguments$cols = c("rt", "mz", "intensity", "dt")
        arguments$colheader = c("Feature","mz","RT","DT","Intensity","File")
    } else {
        arguments$cols = c("rt", "mz", "intensity")
        arguments$colheader = c("Feature", "mz", "RT", "Intensity", "File")
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
    print("Getting feature table should be (mz, rt, id, (dt), formula) and got:")
    print(colnames(cols))
    return(cols)
}

getSpectra <- function(arguments, PEAKS, IDs, rt, imdt, i){
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
    PEAKS = lapply(seq_along(IDs), function(i) peaks(AllIons, IDs[i]))
    data = header(AllIons,IDs)
    rt = data$retentionTime/60 # retention times in minutes
    # pMz = data$precursorMZ #get precursorMZ
    # snum = data$acquisitionNum #get acquisition number (scan num)
    if("dt" %in% arguments$cols){ 
        # Ion Mobility Drift Time
        imdt = data$ionMobilityDriftTime
    } else{
        imdt = NULL
    }

    print(paste("length of CE processing: ",length(IDs)))
    MS = lapply(seq_along(IDs), function(i) getSpectra(arguments, PEAKS, IDs, rt, imdt, i))
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
