
#Force the directory for packages to be the one that is distributed not the local R directory
if(length(.libPaths())>1){
  R_DistrDir<-.libPaths()[2]
  .libPaths(R_DistrDir)
}

#Or just force .libPaths()
# .libPaths("C:/NEW_SOFTWARE/2023_24_UPDATES_FM_LM/FluoroMatch-4.4_GitHub/Flow/R-4.2.1/library/")

if (FLOW == TRUE) {
  csvInput <- TRUE
  ManuallyInputVariables <- FALSE
}


#Checks for updates, installs packagaes: "installr" "stringr" "sqldf" "gWidgets" "gWidgetstcltk" and "compiler"
# if(!require(stringr)) {
# install.packages("stringr"); require(installr)}
# library(installr)
if (!require("stringr", quietly = TRUE))install.packages("stringr", repos = "http://cran.us.r-project.org")
library("stringr")

if (EICcpp==TRUE) {
  #These are distributed by us (thanks Jonathan)
  #cannot installed online
  if("Rcpp" %in% (.packages())){
    detach("package:mzR", unload=TRUE) 
    detach("package:Rcpp", unload=TRUE) 
  }
  library("EICim")
  library("MS1im")
}

if("MASS" %in% rownames(installed.packages()) == FALSE) {install.packages("MASS", repos = "http://cran.us.r-project.org")}
library("MASS")

if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (ParallelComputing == TRUE) {
  if("foreach" %in% rownames(installed.packages()) == FALSE) {install.packages("foreach", repos = "http://cran.us.r-project.org")}
  library(foreach)
  if("doParallel" %in% rownames(installed.packages()) == FALSE) {install.packages("doParallel", repos = "http://cran.us.r-project.org")}
  library(doParallel)
  if("parallelly" %in% rownames(installed.packages()) == FALSE) {install.packages("parallelly", repos = "http://cran.us.r-project.org")}
  library("parallelly")
  # tic("Main")
  ###NOTE not sure if any of this is necessary but was causing error
  DC <- as.numeric(availableCores(constraints = "connections"))
  print(paste0(freeConnections()," connections are free and will be used of ",availableConnections()," R hard limited connections; the following numbers of cores exist on your system: ",as.numeric(detectCores())))
  if (is.na(DC)) {
    DC <- makePSOCKcluster(4)
    registerDoParallel(DC)
  } else if (DC <= 4) {
    DC <- makePSOCKcluster(DC)
    registerDoParallel(DC)
  } else {
    DC <- makePSOCKcluster(DC-2)
    registerDoParallel(DC)
  }
  # Force the directory for parallel computing to be the one that is distributed not the local R directory (doPar error)
  invisible(clusterEvalQ(DC, .libPaths()[2]))
}

if("tictoc" %in% rownames(installed.packages()) == FALSE) {install.packages("tictoc", repos = "http://cran.us.r-project.org")}
library(tictoc)
tic.clear()
# tic("Main")

if("remotes" %in% rownames(installed.packages()) == FALSE) {install.packages("remotes", repos = "http://cran.us.r-project.org")}
library("remotes")

if("MassTools" %in% rownames(installed.packages()) == FALSE) {remotes::install_github("mjhelf/MassTools")}
library("MassTools")

#library(MetaboCoreUtils)
if("enviPat" %in% rownames(installed.packages()) == FALSE) {install.packages("enviPat", repos = "http://cran.us.r-project.org")}
library(enviPat, quietly=T)

if("RMassBank" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("RMassBank")}
library("RMassBank")

if("MetaboCoreUtils" %in% rownames(installed.packages()) == FALSE) {BiocManager::install("MetaboCoreUtils")}
library("MetaboCoreUtils")

if("sqldf" %in% rownames(installed.packages()) == FALSE) {install.packages("sqldf", repos = "http://cran.us.r-project.org")}
library("sqldf")

if("Rdisop" %in% rownames(installed.packages()) == FALSE) {install.packages("Rdisop", repos = "http://cran.us.r-project.org")}
library("Rdisop")

if("RSQLite" %in% rownames(installed.packages()) == FALSE) {install.packages("RSQLite", repos = "http://cran.us.r-project.org")}
library("RSQLite")

if("agricolae" %in% rownames(installed.packages()) == FALSE) {install.packages("agricolae", repos = "http://cran.us.r-project.org")}
library("agricolae")

if("gWidgets" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgets", repos = "http://cran.us.r-project.org")}
if("gWidgetstcltk" %in% rownames(installed.packages()) == FALSE) {install.packages("gWidgetstcltk", repos = "http://cran.us.r-project.org")}
# if (FLOW == FALSE && csvInput == FALSE && ManuallyInputVariables == FALSE) {
require(gWidgets)
require(gWidgetstcltk)
options(guiToolkit="tcltk")
# }

#library(Rdisop)
#BiocManager::install("Rdisop")
library(RSQLite)
library(sqldf)
library(agricolae)
library(Rdisop)
options(warn=-1)#suppress warning on

errorBox <- function(message) {
  window <- gwindow("Confirm")
  group <- ggroup(container = window)
  
  ## A group for the message and buttons
  inner.group <- ggroup(horizontal=FALSE, container = group)
  glabel(message, container=inner.group, expand=TRUE, icon="error")
  
  ## A group to organize the buttons
  button.group <- ggroup(container = inner.group)
  ## Push buttons to right
  addSpring(button.group)
  gbutton("ok", handler=function(h,...) dispose(window), container=button.group)
  return()
}

if("data.table" %in% rownames(installed.packages()) == FALSE) {install.packages("data.table", repos = "http://cran.us.r-project.org")}
library(data.table)
if("SearchTrees" %in% rownames(installed.packages()) == FALSE) {install.packages("SearchTrees", repos = "http://cran.us.r-project.org")}
library(SearchTrees)
if("comprehenr" %in% rownames(installed.packages()) == FALSE) {install.packages("comprehenr", repos = "http://cran.us.r-project.org")}
library(comprehenr)
if("mzR" %in% rownames(installed.packages()) == FALSE)  {
  BiocManager::install("mzR")}
library(mzR)
