#!/usr/bin/env Rscript

# Credit to cjp : https://stackoverflow.com/questions/16225530/contours-of-percentiles-on-level-plot

# packages <- c("MASS", "data.table");

readFeatureTable <- function(fn, columns, skip){
    # Read in Feature Table (Typically NegID) file with columns
    # Allow for skipping extra info column headers
    data <- fread(fn)
    cols <- data[, ..columns]
    if(length(skip) > 0){ 
        cols = cols[-skip,]
    }
    cols <- as.matrix(cols)
    # cols <- as.numeric(cols)
    return(cols)
}

load_kauf <- function(kauf_Table = "2020_EPA_MasterList_Kaufmann_Plots.csv"){
    columns = c(10, 12)
    df = readFeatureTable(kauf_Table, columns, skip=c())

    df = cbind(df, 1) # EPA data
    df = as.matrix(df)
    class(df) <- "numeric"
    # colnames(df) <- c("MZ_eC", "MD_eC", "group")
    colnames(df) <- c("x", "y", "group")
    return(df)
}

get_density_prob2 <- function(test_points, dens, levels) {
    # Compute x and y indices for all test points
    x_idx <- findInterval(test_points[, 1], dens$x)
    y_idx <- findInterval(test_points[, 2], dens$y)

    # Filter points within bounds
    in_bounds <- x_idx > 0 & x_idx <= length(dens$x) & y_idx > 0 & y_idx <= length(dens$y)
    density <- rep(0, nrow(test_points))  # Initialize densities with 0
    density[in_bounds] <- dens$z[cbind(x_idx[in_bounds], y_idx[in_bounds])]

    # Map densities to probabilities
    level_idx <- findInterval(density, levels)
    probability <- pmax(0, 1 - level_idx / 100)  # Ensure probabilities are non-negative

    cbind(density, probability)
}

get_point_density <- function(ReferenceKauffData_csv, outputFile, StepSize, MZ_col_name, MD_col_name) {
    # Load data
    data <- load_kauf(ReferenceKauffData_csv)
    
    # Extract x and y
    x <- data[, 1]
    y <- data[, 2]
    
    # Specify grid limits
    #       If the EPA dataset changes
    #           plot the points and make a new bounding box
    xlim <- c(10, 120)
    ylim <- c(-0.04, 0.015)
    
    # Perform 2D KDE
    dens <- kde2d(x, y, n = StepSize, lims = c(xlim, ylim));
    dx <- diff(dens$x[1:2])
    dy <- diff(dens$y[1:2])

    # Compute cumulative density and levels
    sz <- sort(dens$z)
    c1 <- cumsum(sz) * dx * dy
    c1_unique <- unique(c1)
    sz_unique <- sz[!duplicated(c1)]
    
    # Probability for contour levels
    prob <- seq(1, 0, -0.01)
    levels <- sapply(prob, function(p) { approx(c1_unique, sz_unique, xout = 1 - p)$y})
    levels[1] <- 0
    levels[length(levels)] <- 1.1 * levels[length(levels)-1]
 
    # load feature table, grab slice from it, calculate (density, probability)
    FT <- fread(outputFile);
    #grab md and mz columns
    FT_matrix<-as.matrix(FT)
    col_MD_eC<-which(FT_matrix[1,]==MD_col_name)
    if (length(col_MD_eC)==0) {
      col_MD_eC<-which(colnames(FT_matrix)==MD_col_name)
    }
    col_MZ_eC<-which(FT_matrix[1,]==MZ_col_name)
    if (length(col_MZ_eC)==0) {
      col_MZ_eC<-which(colnames(FT_matrix)==MZ_col_name)
    }
    MZMD <- as.matrix(FT[, .SD, .SDcols = c(col_MZ_eC, col_MD_eC)])
    density_prob <- get_density_prob2(MZMD, dens, levels)

    # splice the data in, correct column names, write to an output file
    updated_FT <- cbind( FT[, 1:col_MD_eC]
                        , density_prob[, 1]
                        , density_prob[, 2]
                        , FT[, (col_MD_eC+1):length(FT)]  )
    colnames(updated_FT) <- c(  names(FT)[1:(col_MD_eC)]
                                , c("Kauff_Density"), c("Kauff_Probability")
                                , names(FT)[(col_MD_eC+1):length(FT)]
                                )
    write.table(updated_FT, outputFile, row.names = FALSE, col.names = TRUE, na = "", sep=",")

}
