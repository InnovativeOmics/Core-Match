# Core-Match
Core scripts for LipidMatch, PolyMatch, and FluoroMatch
For LipidMatch, FluoroMatch, and PolyMatch **users should directly download from innovativeomics.com/software** for the latest stable release.

**Developers can edit the main algorithms here on github** as a team. Then the code needs to be integrated by placing the edited code in the correct directory from the downloaded distribution for innovativeomics.com, please contact us as well (jeremykoelmel@gmail.com) with major changes or project ideas, or discuss here.

Instructions for integrating the code into the software framework:
Install FluoroMatch or LipidMatch from InnovativeOmics.com

If edited, put Modular.r into two directories (replace existing code):

`LipidMatch-4.2\Flow\LipidMatch_Distribution`

`LipidMatch-4.2\FluoroMatch_Modular`

Place any edited script files (genEIC.r, MS1Spectragen.r, Stats.R) into:

`LipidMatch-4.2\Flow\LipidMatch_Distribution\LipidMatch_Libraries\Scripts`

For the Modular version set

`FLOW <- FALSE`

`csvInput <- FALSE #Alternatively you can keep this true and use csv inputs`

`ManuallyInputVariables <- FALSE`

For the Flow version set

`FLOW <- TRUE`

Make sure to toggle the following parameters depending on your application, if both are FALSE it defaults to PFAS analysis:

`Lipid <- TRUE`

`TWeen_pos <- FALSE`

Make sure to follow all instructions for installing dependencies and package installation:

LipidMatch-4.2\FluoroMatch_Flow_Manual.docx

LipidMatch-4.2\FluoroMatch_Modular\FluoroMatch_Manual.docx

LipidMatch-4.2\FluoroMatch_Modular\Trouble_Shooting_Common_Issues.docx

If you have issues please contact Jeremy Koelmel at jeremykoelmel@gmail.com or open them up here on github

## To test the modular project
 * Download or clone project
 * Copy the contents of the root Scripts folder into:
   * RapidTest/LipidMatch_Libraries/Scripts/
