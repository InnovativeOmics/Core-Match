#!/usr/bin/env Rscript

getSymbol <- function(isostring){
    # Get the symbol from the iso string
    symbol = strsplit(isostring, "\\d+")[[1]]
    return(symbol[length(symbol)])
}

getIso <- function(isostring){
    sym = getSymbol(isostring)
    # split about our Atomic Symbol
    nums = as.numeric(strsplit(isostring, sym)[[1]])

    # get the mass number from the string, or set to 0
    massnum = nums[1]
    if(is.na(massnum)){
        massnum = 0
    }

    # get the count, if any, or set to 1
    if(length(nums) == 1){
        nums = c(nums[1], 1)
    }
    count = nums[2]

    # construct a data frame row for each of these variables
    df <- data.frame(Symbol=sym, MassNumber=massnum, Count=count)
    return(df)
}

parseIsoString <- function(isostring){
    # replace spaces with empty string
    # split about commas
    isostring = gsub(" ", "", isostring)
    isos = strsplit(isostring,";")[[1]]
    
    # get mini data frames (rows) 13C5 -> 13 C 5
    # then combine them
    df = lapply(seq_len(length(isos)), function(i) getIso(isos[i]))
    df = do.call(rbind, df)
    return(df)
}

getHarmonics <- function(df, i){
    # Duplicate the table up to count, 
    # 13C5 -> 5 copies
    # Mulitply relative mass (13C + 13C) - (12C + 12C)
    # Exponentiate abundance 0.0108^i
    # Convert to percentage 
    h = data.frame(df)
    h$RelativeMass = h$RelativeMass * i
    h$IsotopicComposition = 100 * h$IsotopicComposition ^ i
    h$Symbol = paste(h$MassNumber, h$AtomicSymbol, i, sep="")
    return(h[,c("Symbol", "RelativeMass", "IsotopicComposition")])
}

getIsoData <- function(iso, data){
    # Select the rows from the NIST table for the symbol and mass number
    if(is.element(iso$Symbol, data$AtomicSymbol)){
        df = data[
                data$AtomicSymbol == iso$Symbol 
                & (iso$MassNumber == 0 
                | data$MassNumber == iso$MassNumber )
            , ]

        # Construct a df for multiple copies of the isotope and combine them together
        M = lapply(seq_len(iso$Count), function(i) getHarmonics(df, i))
        M = do.call(rbind, M)
    }else{
        M = data.frame()
    }
    return(M)
}

genIsoTable <- function(arguments){
    # Construct an isotope data frame from the input string (Symbol, MassNumber, Count)
    # Get the data from the NIST table for each isotope
    isos = parseIsoString(arguments$isostring)
    nistdf = read.csv(arguments$isotable)
    df = lapply(seq_len(nrow(isos)), function(i) getIsoData(isos[i,], nistdf))
    df = do.call(rbind, df)
    # print(df)
    # write.table(df, row.names = FALSE, sep = ",")
    return(df)
}

example <- function(){
    arguments = list()
    arguments$isostring = "13C3;N;S;Cl2;18O;Br2;  NotAnElement  ; Ti"
    genIsoTable(arguments)
}
# example()