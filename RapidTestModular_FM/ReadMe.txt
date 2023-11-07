Use Modular from GitHub Core: 
https://github.com/InnovativeOmics/Core-Match/Modular.R

download project
unzip
    test files

To run full project
Copy the contents of the root Scripts folder into:
    RapidTest/LipidMatch_Libraries/Scripts/


In the script set the following:
FLOW <- FALSE
csvInput <- TRUE
ManuallyInputVariables <- FALSE
RT_flagging <- TRUE #JPK: for PFAS analysis
ParallelComputing <- TRUE
Lipid <- FALSE
TWeen_pos <- FALSE #PJS: for PolyMatch
FilterAbovePrecursor <- 1 #how far from the precursor should fragment masses be kept (e.g. if precursor is 700, should 702 be considered?)
TargetEIC_Only <- TRUE

AND change the following directory (currently line 171 but that changes) to the one with the PARAMETERS.csv file in the test directory:
parametersDir <- "C:/NEW_SOFTWARE/2023_24_UPDATES_FM_LM/FluoroMatch-4.3/Flow/LipidMatch_Distribution/"

AND change the directories in the PARAMTERS.csv file (two that came with the RapidTest folder), make sure to have the directories with forward slashes and ending in a forward slash:
Directory of Inputs (INPUT folder)
Directory of Libraries (LipidMatch_Libraries folder)