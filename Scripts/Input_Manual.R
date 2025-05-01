
#Retention Time plus or minus
RT_Window <- .3 #window of .2 => +/- .1

#parts-per-million window for matching the m/z of fragments obtained in the library to those in experimentally obtained
ppm_Window <- 10 #window of 10 => +/- 5 ppm error

#Tolerance for mass-to-charge matching at ms1 level
PrecursorMassAccuracy<-0.005

#Plus minus range for the mass given after targeting parent ions portrayed in excalibur to match exact mass of Lipid in library
SelectionAccuracy<-1

#Threshold for determining that the average signal intensity for a given MS/MS ion should be used for confirmation
intensityCutOff<-1000

#Feature Table information
CommentColumn <- 1
MZColumn <- 2
RTColumn <- 3
RowStartForFeatureTableData <- 2 #Look at your Feature Table (.csv)...What row do you first see numbers?


#ddMS data? set this parameter. If not? Leave it.
#The minimum number of scans required for the result to be a confirmation
ScanCutOff<-1

#Have AIF data? set these parameters. If not? Leave them.
corrMin <-.6
minNumberOfAIFScans <- 5

# Input Directory for Feature Table and MS2 files ...must have forward slashes but nothing at the end of the directory
InputDirectory<-"C:/Users/Jeremy Koelmel/Downloads/LipidMatch_Flow_3.5/LipidMatch_Flow_3.5/LipidMatch_Modular/ExampleData"


## PFAS Specific parameters
# RT_flagging <- TRUE #JPK: for PFAS analysis
# Mass defect filtering
upper <- 0.12 #JPK: Upper limit for PFAS mass defect
lower <- -0.11 #JPK: Lower limit for PFAS mass defect

#the isotopes that you can add for MS1 isotope labeling and isotope ratio calculations, if you don't differentiate the mass number it will use all from "secondary_isotopes.csv"
ISOstring <- "13C3;15N;33S;34S;Cl3;18O;Br3;29Si;30Si"

#Check your parameters and see that they are all correct.
#Then press ctrl+shift+s to execute all code.