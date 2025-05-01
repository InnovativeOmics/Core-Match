
ErrorOutput<-0
#Retention Time plus or minus
RT_Window <- as.numeric(parametersInput_csv[2,2]) #window of .2 => +/- .1

#parts-per-million window for matching the m/z of fragments obtained in the library to those in experimentally obtained
ppm_Window <- as.numeric(parametersInput_csv[4,2]) #window of 10 => +/- 5 ppm error

#Tolerance for mass-to-charge matching at ms1 level (Window)
PrecursorMassAccuracy <- as.numeric(parametersInput_csv[3,2])

#Plus minus range for the mass given after targeting parent ions portrayed in excalibur to match exact mass of Lipid in library
SelectionAccuracy <- as.numeric(parametersInput_csv[5,2])

#Threshold for determining that the average signal intensity for a given MS/MS ion should be used for confirmation
intensityCutOff <- as.numeric(parametersInput_csv[7,2])

#Feature Table information
CommentColumn <- as.numeric(parametersInput_csv[10,2])
MZColumn <- as.numeric(parametersInput_csv[11,2])
RTColumn <- as.numeric(parametersInput_csv[12,2])
RowStartForFeatureTableData <- as.numeric(parametersInput_csv[13,2]) #Look at your Feature Table (.csv)...What row do you first see numbers?

#ddMS data? set this parameter. If not? Leave it.
#The minimum number of scans required for the result to be a confirmation
ScanCutOff<-as.numeric(parametersInput_csv[6,2])

#Have AIF data? set these parameters. If not? Leave them.
corrMin <-as.numeric(parametersInput_csv[9,2])
minNumberOfAIFScans <- as.numeric(parametersInput_csv[8,2])

# Input Directory for Feature Table and MS2 files ...must have \\ (double backslash) at the end of the directory
InputDirectory<-as.character(parametersInput_csv[14,2])
InputDirectory = substring(InputDirectory,1, nchar(InputDirectory)-1)

GroupCSVDirectory<-as.character(parametersInput_csv[20,2])

ISOstring <- as.character(parametersInput_csv[23,2])

## PFAS Specific parameters
RT_flagging <- TRUE #JPK: for PFAS analysis
# Mass defect filtering
upper <- 0.12 #JPK: Upper limit for PFAS mass defect
lower <- -0.11 #JPK: Lower limit for PFAS mass defect