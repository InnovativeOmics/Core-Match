#!/usr/bin/env Rscript
# library(data.table)
# library(MetaboCoreUtils)

readFeatureTable <- function(fn, columns, skip){
    # Read in Feature Table (Typically NegID) file with columns
    # Allow for skipping extra info column headers
    data <- fread(fn)
    cols <- data[, ..columns]
    if(length(skip) > 0){ 
        cols = cols[-skip,]
    }
    cols <- as.matrix(cols)
    # print("Getting feature table should be (sti, mz, topmf) and got:")
    # print(colnames(cols))
    return(cols)
}


get_series_ranges <- function(data){
    # gets the adjacent groupings of the homologous series from the FT
    # This works because they're ordered by series already
    # output is a flatten list of [[start, end]]
    # to access, you'll have index properly
    # insight found on SO here:
    # https://stackoverflow.com/questions/32489277/find-index-of-change-in-a-column
    reps = data[,1]
    changes = which(reps[-1] != reps[-length(reps)])

    changes = sort(c(1, changes, changes+1))
    return(changes)
}
# example usage where data's first column is the series:
# changes = get_series_ranges(data)


get_rep <- function(s){
    # Get repeating unit from string
    # "[CF2]n_100" -> "CF2"
    # https://stackoverflow.com/questions/39086400/extracting-a-string-between-other-two-strings-in-r
    pattern <- "\\[\\s*(.*?)\\s*\\]"
    result <- regmatches(s, regexec(pattern, s))
    return(result[[1]][2])
    # return(countElements(result[[1]][2]))
}
# example usage: 
# print(get_rep("[CF2]n_100"))


get_remainder <- function(formula, rep){
    # repeatedly multiply the rep unit until it can't be subtracted
    # report last valid remainder
    if( !containsElements(formula, rep) || formula == "NA"){
        # Catch the case where rep is not in formula
        return( list("", 0) )
    } 
    n = 1
    while (TRUE) {
        rep_n = multiplyElements( rep, n ) # n must be a positive integer
        rep_n1 = multiplyElements( rep, n + 1 )
        # get the last instance where rep is in the formula
        # get the remainder: formula - rep_n
        if(     containsElements(formula, rep_n ) 
            && !containsElements(formula, rep_n1) ){
            return( list(subtractElements(formula, rep_n), n) )
            # return( subtractElements(formula, rep_n) )
        }
        n = n + 1
    }
}
# example usage: 
# print( get_remainder("C3HO3SF7", "CF2"))


get_remainders_reps <- function(formulas, rep){
    # get the remainders and number of repeats of the repeating unit
    # input is: 
    #    a vector of formulas, ex: ("HF", "H2FN", ...)
    #    a repeating unit, ex: "H"
    # output is:
    #    the list of remaining atoms after subtracting out rep: ("F", "FN", ...)
    #    the list of the number of times the rep unit can be subtracted: (1, 2, ...)

    remainders = c()
    reps = c()
    # print(rep)
    for( i in 1:length(formulas)){
        # print(formulas[i])
        if(!is.na(formulas[i])){
            re = get_remainder(formulas[i], rep)
            # print(paste("[", rep, "]", re[[2]], re[[1]], sep=""))
        } else {
            re = list("", 0)
        }
        remainders = c(remainders, re[[1]])
        reps = c(reps, re[[2]])
        # print(paste("[", rep, "]", re[[2]], re[[1]], sep=""))
        # print(rep_mass*re[[2]] + calculateMass(re[[1]]) )
    }
    return(list(remainders, reps))
}
# usage:
# get_remainders_reps(formulas, rep)

process_series <- function(series, charge, formula_column, calculate_error){
    # print( series )
    # print( series[1,1] )
    # print( nrow(series) )
    series_number = as.numeric( strsplit(series[1,1], split="_")[[1]][2] ) # "series_123" -> 123
    formulas = series[,formula_column]
    unique_formulas = unique(formulas[!is.na(formulas)])
    unique_formulas = unique_formulas[unique_formulas != ""]
    # unique_adducts9 = unique(adducts[!is.na(adducts)])
    if(series_number < 1 || length(unique_formulas) < 2){ return(series) }
    # print(unique_formulas)

    rep = get_rep( series[1,1] )
    rep_mass = calculateMass(rep)[[1]]

    result = get_remainders_reps(unique_formulas, rep)
    remainders = result[[1]]
    reps = result[[2]]

    # print(remainders)
    # This is a table of subclasses and the number of votes in each for the unique formulas
    t = table(remainders[nzchar(remainders)])
    # print(t)
    max_votes = max(t)
    # print(paste("max_votes = ", max_votes, " |  length(t) = ", length(t), sep=""))
    if(max_votes <= 1){ return(series) }
    max_vote_ids = grep(max_votes, t) # get columns of max votes from table, could tie (2,2)

    for(vid in 1:length(max_vote_ids)){
        mode_remainder = names(t[max_vote_ids[vid]])
        mode_remainder_mass = as.numeric(calculateMass(mode_remainder)) #- adductmass
        # print(paste(vid, "mode_remainder", mode_remainder, mode_remainder_mass, sep=", "))

        # print("mz, v_reps, v_rep_mass, mode_remainder_mass, theoretical_mass, mz_err_ppm")
        for( i in 1:nrow(series)){
            mz = as.numeric( series[i,2] )
            v_reps = round(  (mz - mode_remainder_mass) / rep_mass  )
            v_rep_mass = rep_mass * v_reps

            if( calculate_error ){
                # doing ZFormula
                electron = (0.00055 * charge * -1)
                theoretical_mass = mode_remainder_mass + v_rep_mass + electron
                # charge = round(theoretical_mass / mz)
                # mz = mz * charge
                # Assuming singly charged adducts

                mz_err_ppm = (  (mz - theoretical_mass) / theoretical_mass  ) * 10^6
                mz_err_ppm = round(mz_err_ppm, digits=4)
                # print(paste(mz, v_reps, v_rep_mass, mode_remainder_mass, theoretical_mass, mz_err_ppm, sep=", "))
            }

            hs = paste("[", rep, "]", v_reps, mode_remainder, sep="")
            if( calculate_error ) {
                if(series[i,6] == ""){
                    series[i,6] = hs
                    series[i,7] = mz_err_ppm
                } else{
                    series[i,6] = paste(series[i,6], hs, sep=";")
                    series[i,7] = paste(series[i,7], mz_err_ppm, sep=";")
                }
            } else {
                if(series[i,5] == ""){
                    series[i,5] = hs
                } else{
                    series[i,5] = paste(series[i,5], hs, sep=";")
                }
            }
        }
    }
    # print(series[,c(2, 3, 4)])
    return(series)
}
# example usage: 
# process_series( data[1:29,], 0.01 )
# process_series( data[22:25,], 0.01 )  #example tie


homologous_voting <- function(
    col_sti = 3   #SeriesType_Identifier
    , col_mz = 7  #m.z. Mass to Charge Ratio
    , col_mf = 5 #Formula Molecular Formula
    # , col_mf = 13 #TopMF Molecular Formula
    , col_mfz = 15 #zFormula Molecular Formula
    , fn_FT = 'NegIDed_FIN.csv'
    , SeriesFormula_ColName = "HSformula"
    , SeriesFormula_ColNameZ = "HSZformula"
    , SeriesError_ColNameZ = "HSZppm_error"
    , charge = -1 #polarity, negmode = -1, posmode = 1, assumes 1 charge
    ){
    skip = c() #for use in skipping rows for reading feature table

    # repunits = readFeatureTable("REPEATING_UNITS_INPUT.csv", c(1,2,3), skip)
    columns = c(col_sti, col_mz, col_mf, col_mfz)
    data = readFeatureTable(fn_FT, columns, skip)
    data[,2] = as.numeric(data[,2]) 
    data = cbind(data, "")
    data = cbind(data, "")
    data = cbind(data, "")
    # colnames(data) = c("SeriesType_Identifier", "m.z", "TopMF", SeriesFormula_ColName, SeriesError_ColName)
    colnames(data) = c("SeriesType_Identifier", "m.z", "Formula", "zFormula", SeriesFormula_ColName, SeriesFormula_ColNameZ, SeriesError_ColNameZ)
    changes = get_series_ranges(data)

    i = 0
    # while(i*2 + 2 < 5){
    while(i*2 + 2 < length(changes)){
        low = changes[2*i+1]
        high = changes[2*i+2]
        # print(paste(low, high, sep=" - "))
        if( 0 < high - low){
            data[ low : high , ] = process_series( 
                data[ low : high , ], charge, 4, TRUE ) #matches col_mfz in columns
            data[ low : high , ] = process_series( 
                data[ low : high , ], charge, 3, FALSE ) #matches col_mf in columns
            # print(data[ low : high , ])
        }
        i = i+1
    }

    FT = fread(fn_FT)

    updated_FT <- cbind( FT[, 1:(col_mf)]
                        , data[,c(5)]
                        , FT[, (col_mf+1):col_mfz]
                        , data[,c(6,7)]
                        , FT[, (col_mfz+1):length(FT)]  )

    colnames(updated_FT) <- c(  names(FT)[1:(col_mf)]
                                , c(SeriesFormula_ColName)
                                , names(FT)[(col_mf+1):col_mfz]
                                , c(SeriesFormula_ColNameZ, SeriesError_ColNameZ)
                                , names(FT)[(col_mfz+1):length(FT)]
                                )

    # write.table(updated_FT, file="output_test.csv", row.names = FALSE, sep = ",")
    write.table(updated_FT, file=fn_FT, row.names = FALSE, sep = ",")
}
