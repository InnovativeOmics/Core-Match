
createDataFrame <- function (msx_dir){
  #################################
  #Store ms1/2 data into dataframe#
  #################################
  dir_string <- strsplit(msx_dir, split="")[[1]]
  Level <- as.numeric(dir_string[length(dir_string)])
  msx <- scan(file=msx_dir, what='character') #read in .ms1/2 as list of characters
  s_index <- match("S", msx)
  msx <- msx[s_index:length(msx)] #cut off head useless info from .ms1/2
  s_indicies <- which(msx=="S")
  RTime_indicies <- which(msx == "RTime")
  TIC_indicies <- which(msx == "TIC")
  if (length(grep("AIF", msx_dir)) > 0) {
    AIF <- TRUE
  } else {
    AIF <- FALSE
  }
  if (Level == "1" || AIF == TRUE) {
    msx_df <- data.frame(scanNum=numeric(), rt=numeric(), mz_intensity=I(list()))
  }else{
    msx_df <- data.frame(precursor=numeric(), rt=numeric(), mz_intensity=I(list()))
  }
  to_remove <- c()
  
  s_indicies <- append(s_indicies, length(msx) + 1)
  length_s_indicies <- length(s_indicies)
  for (i in 1:(length_s_indicies - 1)){
    if (msx[TIC_indicies[i] + 2] != "Z"){
      if ((s_indicies[i + 1] - TIC_indicies[i]) <= 4){
        to_remove <- append(to_remove, i)
        next
      }
      data <- msx[(TIC_indicies[i] + 2):(s_indicies[i + 1] - 1)]
      mz <- data[c(TRUE, FALSE)]
      intensity <- data[c(FALSE, TRUE)]
    } else {
      if ((s_indicies[i + 1] - TIC_indicies[i]) <= 7){
        to_remove <- append(to_remove, i)
        next
      }
      data <- msx[(TIC_indicies[i] + 5):(s_indicies[i + 1] - 1)]
      mz <- data[c(TRUE, FALSE)]
      intensity <- data[c(FALSE, TRUE)]
    }
    if (Level == "1" || AIF == TRUE) {
      msx_df[i, 1] <- msx[s_indicies[i] + 1]
    }else{
      msx_df[i, 1] <- msx[s_indicies[i] + 3]
    }
    msx_df[i, 2] <- msx[RTime_indicies[i] + 1]
    msx_df[[i, 3]] <- matrix(cbind(mz, intensity), length(mz), 2)
  }
  
  if (length(to_remove) > 0) {
    msx_df <- msx_df[-to_remove,]
  }
  
  return(msx_df)
}
