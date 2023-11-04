#!/usr/bin/env Rscript

get_ref_row <- function(group){
    # Get the reference row, sum the intensities over multiple samples
    rows = group[ group$Isotope == "M", ]
    rows[1, "Intensity"] = sum(rows[,"Intensity"])
    rows[1, "Feature"] = group[1,"Feature"] #in case M isn't found
    return(rows[1,])
}

get_iso_rows <- function(group){
    # get rows which have isotope entries from MS1 data
    return( group[ endsWith(group$Isotope,")"), ] )
}

split_iso_string <- function(iso_str){
    # Example Input string
    # 18O(0.20%;0.00262Da)_81Br(97.51%;0.00367Da)_37Cl(32.40%;0.00457Da)

    # Split the string into fields
    fields <- strsplit(iso_str, "_", fixed = TRUE)[[1]]

    # Extract symbol, composition, and distance
    symbols <- gsub("\\(.*", "", fields)
    abundances <- as.numeric(gsub(".*\\(|%.*", "", fields))
    distances <- as.numeric(gsub(".*%;|Da)", "", fields))

    df <- data.frame(Symbol = symbols, Distance = distances, theoreticalAbundance = abundances)
    return(df)
}

mini_iso_df <- function(rows, i, ref){
    # split the isostring for an MS1 row and turn it into a data frame
    s = split_iso_string(rows[i,4])
    s$m1Distance = rows[i,2] - ref$m.z
    s$relativeAbundance = 100 * rows[i,3] / ref$Intensity
    s$predictedCount = round( s$relativeAbundance / s$theoreticalAbundance )
    s$mz = rows[i,2]
    s$Feature = ref$Feature
    return(s)
}

process_group <- function(group){
    # Make a data frame of all the isotopic information for one feature
    ref = get_ref_row(group)
    rows = get_iso_rows(group)
    ISOs = lapply(seq_len(nrow(rows)), function(i) mini_iso_df(rows, i, ref))
    df = do.call(rbind, ISOs)
    return(df)
}

append_to_FeatureID <- function(arguments, featureID, groups){
    # Append empty columns in the appropriate place
    # For each of the rows of the feature table
    # Check if there's data from the MS1 groups
    # For each MS1 group row
    # Check if it's one of the selected elements
    # Accumulate that abundance to the appropriate symbol

    # Get the index of the column before which you want to insert the new columns
    insert_index <- which(names(featureID) == "Checked_Viz")
    # Create empty vectors or matrices with zeros
    empty_columns <- matrix(0, nrow = nrow(featureID), ncol = length(arguments$selected_isotopes))
    # Create the updated dataframe by combining the existing columns and new columns
    updated_featureID <- cbind(featureID[, 1:(insert_index - 1)], empty_columns, featureID[, insert_index:length(featureID)])
    # Update the column names if necessary
    colnames(updated_featureID) <- c(names(featureID)[1:(insert_index - 1)], arguments$selected_isotopes, names(featureID)[insert_index:length(featureID)])
    featureID = updated_featureID

    featureID[,arguments$selected_isotopes] = 0
    for (i in 1:nrow(featureID)) {
        rowID = paste(featureID[i,"row.ID"])
        if( !is.null(groups[[rowID]]) ){
            for ( j in 1:nrow(groups[[rowID]]) ) {
                group_row = groups[[rowID]][j,]
                if(is.element(group_row$Symbol, arguments$selected_isotopes)){
                    featureID[i, c(group_row$Symbol)] = featureID[i, c(group_row$Symbol)] + group_row$relativeAbundance
                }
            }
        }
    }
    featureID[,arguments$selected_isotopes] = round(featureID[,arguments$selected_isotopes],6)

    return(featureID)
}

printWarning <- function(arguments, df){
    items_not_found <- arguments$selected_isotopes[!(arguments$selected_isotopes %in% unique(df$Symbol))]
    if( length(items_not_found) > 0){
        items_not_found_message = paste("Isotope(s)", paste(items_not_found, collapse = ", "), "not found in any of the provided samples. Could be: due to typo, not existing in our NIST table, M peak not found, or precision issues.")
        print(items_not_found_message)
    }
}

append_percentages <- function( arguments ){
    options(width = 160)

    arguments$selected_isotopes = genIsoTable(arguments)$Symbol

    arguments$fn_MS1_input = file.path(arguments$path_to_output_folder, arguments$fn_MS1_output)
    #arguments$fn_FeatureID_output = file.path(arguments$path_to_output_folder, paste("isotopes", arguments$fn_FeatureID, sep="_"))
    arguments$fn_FeatureID_wpath = file.path(arguments$path_to_output_folder, arguments$fn_FeatureID)
    ms1_data = read.csv(arguments$fn_MS1_input)
    # print(ms1_data)

    # load the ms1 data and remove filenames which contain blank, then split into groups by feature
    ms1_data = ms1_data[!grepl("blank", ms1_data$File, ignore.case=TRUE), ]
    groups = split(ms1_data, ms1_data$Feature)
    # print(groups[[1]])

    # Create a dataframe which has all of the isotopes and relative abundances
    MS = lapply(seq_along(groups), function(i) process_group(groups[[i]]))
    M = do.call(rbind, MS)
    df = M[,c("Feature", "Symbol", "relativeAbundance")]
    df$relativeAbundance[is.infinite(df$relativeAbundance)] <- NA  #set infinities to NA

    printWarning(arguments, df)
    groups = split(df, df$Feature)

    # append the data to columns in the Feature Table
    featureID = read.csv(arguments$fn_FeatureID_wpath)
    featureID = append_to_FeatureID(arguments, featureID, groups)

    write.table(featureID, file=arguments$fn_FeatureID_wpath, row.names = FALSE, sep = ",")
}