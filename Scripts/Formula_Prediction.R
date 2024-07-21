Formula_Prediction <- function(Override_Predict,fh_Feature_MS1s,fh_Feature_IDList,fh_MS1,fh_IDed_FIN,MF_topN,adducts,eList1,eList2,q,Poltxt,ppmTol) {
  data("isotopes")  #requires enviPat
  mf.stor=T  #For debugging purposes (stores times to run)
  z = abs(q)  # z is charge, may need to feed this variable on per feature basis
  pks_MS1_Spectra <- read.csv(fh_Feature_MS1s)
  pks_FeatureIDs <- read.csv(fh_Feature_IDList)
  pks_FeatureIDs$Adduct[which(pks_FeatureIDs$Adduct=='[M-H]- ')] <- '[M-H]-'  #Workaround to fix extra space in Adducts
  # Load Peaks - Spectra defined as M
  PeakList <- pks_MS1_Spectra[which(pks_MS1_Spectra$Isotope=='M'),]
  MS1_pks <- merge(aggregate(Intensity~Feature,PeakList,FUN=max), PeakList) #Find Spectra with M and max intensity
  MS1_pks <- MS1_pks[order(-MS1_pks$Intensity),]
  
  
  
  ## Calculate Scores for MFs in FeatureID List
  Feature_MF_List <- pks_FeatureIDs[pks_FeatureIDs$Score %in% c('A','A-','B+','B','B-','C+','C'),]
  ## Note I removed Hg from predictions... because it was causing an error
  Feature_MF_List <- Feature_MF_List[which(Feature_MF_List$Formula != '' & !is.na(Feature_MF_List$Formula) & !grepl("Hg", Feature_MF_List$Formula)),] #MF cannot be blank to score it
  
  ## Start MF Reference Scoring for High Conf hits || BM = 1.6 sec per known feature
  FormulaMatches <- data.frame()
  n_FeaturesMF <- nrow(Feature_MF_List)
  #if ((Override_Predict>1)&&(Override_Predict<n_FeaturesMF)){
  #  n_FeaturesMF<-Override_Predict
  #}
  # Override n_FeaturesMF=1000.   DO NOT DO THIS HERE DS 3/19/24 - you want to score all MFs from upstream.  Only limit Denovo below
  for(fmf in 1:n_FeaturesMF){
    
    ft <- Feature_MF_List[fmf,]
    ft.MF = ft$Formula
    ft.id <- ft$row.ID
    pks <- PeakList[which(PeakList$Feature==ft.id),]  #Potentially grab all Files for a feature
    if(nrow(pks) > 0){
      
      idx_pk_spec <- which(pks$Intensity==max(pks$Intensity))
      pk <- pks[idx_pk_spec,]  #Choose file with highest int
      
      s.unk <- pks_MS1_Spectra[which(pks_MS1_Spectra$Feature==ft.id & 
                                       pks_MS1_Spectra$File==pk$File & 
                                       (pks_MS1_Spectra$m.z - pk$m.z > -2.5) &
                                       (pks_MS1_Spectra$m.z - pk$m.z <  4.5) ), 
                               c('m.z','Intensity',"Isotope")] ## Lookup spectra with that FileID and Feature
      s.mdf <- mzSpec_MDFilter2(spectrum=s.unk[,1:2], Q1 = pk$m.z, mdTol=ppmTol*4, z=1, splot=F, plotBG=F, minIntN = 1)  # Deconvolute spectrum
      s.cln <- s.mdf$Spectrum[which(s.mdf$Spectrum$Class=='_Target'),c('mz','Int','Int_N','mzIP')]
      n_mz <- nrow(s.cln)
      
      # Setup Reference Mass Spectrum
      ion_op <- merge(ft, adducts)[,c('a','d')]
      MF.ion <- ionizeMF(m=ft.MF, a=ion_op$a, d=ion_op$d)
      
      # Generate Ref Spec
      s.ref <- MFtoSpectrum(formula=MF.ion$Mion, q=q, Res=F, ppmTol, patOnly=T)
      colnames(s.ref) <- c('mz','Int')  #rename columns
      s.ref <- cbind(s.ref, mzIP=seq(1:nrow(s.ref))) #Gets mzIP
      
      s.unk.ol <- s.cln[which(s.cln$mz > (min(s.ref$mz)-0.1) & s.cln$mz < (max(s.ref$mz)+0.1)),] #Later need to refine spectrum for charge
      s.unk.ol$Int_N <- 100*s.unk.ol$Int_N/max(s.unk.ol$Int_N[which(s.unk.ol$mzIP==0)])  #Must renormalize
      
      # Old check for ppm diff
      #idx_rows <- 1:min(nrow(s.unk.ol),nrow(s.ref))
      #check.ppm <- length(which(abs(s.unk.ol[idx_rows,'mz'] - s.ref[idx_rows,'mz'])*1e6/ft$m.z < ppmTol*2)) == min(nrow(s.unk.ol),nrow(s.ref))
      
      #Check overlap extent.  Is spec "close enough" to compare
      checkOL <- abs(s.ref$mz[which(s.ref$Int==100)] - 
                       s.unk.ol$mz[which(s.unk.ol$Int_N==100 & s.unk.ol$mzIP<=2)])
      
      if(!is.na((checkOL[1] < 1) && (checkOL[1] > 0))){ 
        #SPS <- SpectrumSimilarity_custom(s.ref, s.unk.ol, dppm = z*ppmTol, b = 0.1, int.prec=2,
        #                                 bottom.label = paste0(round(ft$m.z,4),' @',round(Feature_MF_List$Retention.Time[fmf], 2),' min'), 
        #                                 top.label = paste('Simulated',MF.ion$Mion,Poltxt), 
        #                                 xlim = c(min(s.ref$mz[which(s.ref$Int>1)])-0.2, max(s.ref$mz[which(s.ref$Int>1)])+0.2), 
        #                                 x.threshold = 0, print.alignment = F, print.graphic = T, output.list = T)
        #MFScore1 <- round(SPS$similarity.score,5)
        
        ft.MFScores <- SpectrumSimilarity_custom2(obs=s.unk.ol, the=s.ref, dppm=ppmTol, int_prec = 0.3, rnd_prec = 5)
        
        
        rnf <- max(s.unk.ol$Int) #renormalization factor
        rpf <- nrow(ft.MFScores$Alignment)# rows per feature
        FormulaMatches <- rbind(FormulaMatches, cbind(fmf, Feature =pk$Feature,  
                                                      MZ = ft.MFScores$Alignment$m_obs, Abundance=ft.MFScores$Alignment$i_obs*rnf,
                                                      Isotope = seq(1:rpf), Formula = MF.ion$MFinput, zFormula=MF.ion$Mion,
                                                      PredMZ = ft.MFScores$Alignment$m_ref, PredAbundance=ft.MFScores$Alignment$i_ref*rnf,
                                                      File = pk$File,
                                                      Score1 = ft.MFScores$Score1, Score2=ft.MFScores$Score2, Comment='', Rank='DB_Match' ) )
        
        #MFScore2 <- InterpretMSSpectrum::mScore(obs = s.unk.ol, the = s.ref)
        
      }else{
        FormulaMatches <- rbind(FormulaMatches, cbind(fmf, Feature = pk$Feature, 
                                                      MZ = ft$m.z, Abundance=NA, 
                                                      Isotope='M', Formula = MF.ion$MFinput, zFormula=MF.ion$Mion,
                                                      PredMZ=NA, PredAbundance = NA, File=pk$File,
                                                      Score1=NA, Score2=NA, Comment='Spectral overlap too low', Rank=NA) )
      } #endif check overlap
    }else{ #endif spec found  
      FormulaMatches <- rbind(FormulaMatches, cbind(fmf, Feature = ft$row.ID, 
                                                    MZ = ft$m.z, Abundance=NA, 
                                                    Isotope=NA, Formula = NA, zFormula=NA,
                                                    PredMZ=NA, PredAbundance = NA, File=ft$Files.1[1],
                                                    Score1=NA, Score2=NA, Comment='MS1 Spectra not found for feature', Rank=NA) )
    }
  }
  
  # write.csv(FormulaMatches, file=paste0(OutputDirectory,'Neg_Feature_MS1s_MFScores.csv'), row.names = F)
  rm(n_FeaturesMF, fmf, ft, ft.MF, ft.id,pks,idx_pk_spec,pk,s.unk,s.mdf,s.cln,n_mz,ion_op,MF.ion, s.ref,s.unk.ol,checkOL,
     ft.MFScores, rnf, rpf)
  
  
  ## Prioritize prediction on D+ or higher scores
  Feature_MF_List <- pks_FeatureIDs[pks_FeatureIDs$Score %in% c('A','A-','B+','B','B-','C+','C','C-','D+'),]
  
  
  ## START DeNovo MFs for features
  ss_MS1_pks <- MS1_pks[MS1_pks$Feature %in% Feature_MF_List$row.ID,]
  row.names(ss_MS1_pks) <- NULL
  n_DN_MFs <- nrow(ss_MS1_pks)
  ss_MS1_spectra <- pks_MS1_Spectra[pks_MS1_Spectra$Feature %in% ss_MS1_pks$Feature[1:n_DN_MFs] & pks_MS1_Spectra$File %in% ss_MS1_pks$File[1:n_DN_MFs], ]
  
  ## This is where we should implement OverridePredict limits (DeNovo) DS 3/19/24
  if ((Override_Predict > 1)&&(Override_Predict < n_DN_MFs)){
    n_DN_MFs <- Override_Predict
  }
  
  
  ## Careful ***
  pks_MS1_MF <- data.frame()
  pks_Top_MFs <- data.frame()
  dn.FormulaMatches <- data.frame() # DF to be added to FormulaMatches
  MF_Denovo <- data.frame()
  ## Careful ***
  
  for(f in 1:n_DN_MFs){
    
    print(paste0('Processing feature ',f))
    ft.comment <- ''
    #Initialize
    decIso <- T # Decompose Isotopes/spectrum (vs decomposeMass only)
    MFSolved <- F # MF List made successfully
    MFEscape <- F # Used at end of while loop, Must be reset to F
    
    
    pk <- ss_MS1_pks[f,] 
    ft <- pks_FeatureIDs[which(pks_FeatureIDs$row.ID==pk$Feature),]
    
    if(nrow(ft)>0 & pk$m.z <= 1000){
      
      s.unk <- ss_MS1_spectra[which(ss_MS1_spectra$Feature==ss_MS1_pks$Feature[f] & 
                                      ss_MS1_spectra$File==ss_MS1_pks$File[f] &
                                      (ss_MS1_spectra$m.z - pk$m.z > -2.5) &
                                      (ss_MS1_spectra$m.z - pk$m.z <  4.5)  ),c('m.z','Intensity',"Isotope")]
      MI <- pk$m.z
      
      s.mdf <- mzSpec_MDFilter2(spectrum=s.unk[,1:2], Q1 = MI, mdTol=ppmTol*2.8, z=1, 
                                splot=F, minIntN = 0, forceQ1=T, labelBG=F, plotBG=F)  # Deconvolute spectrum
      # s.mdf$Plot
      
      s.cln <- s.mdf$Spectrum[which(s.mdf$Spectrum$Class=='_Target' & s.mdf$Spectrum$mzIP>=0),c('mz','Int','Int_N','mzIP')]
      s.cln <- merge(s.cln, aggregate(Int_N ~ round(mzIP), s.cln, FUN=max) )
      s.cln <- s.cln[order(s.cln$mzIP),c('mz','Int','Int_N','mzIP')]
      s.cln <- s.cln[!duplicated(s.cln[,c('Int_N','mzIP')]),]
      n_mz <- nrow(s.cln)
      
      # Preprocessing - Check for spectral fidelity issues, mzIP = 1 position 
      if(n_mz>1){
        # check for spectral discontinuity
        mzIPdc <- which(diff(s.cln$mzIP) != 1)  # removed +1
        txt_mzIPdc <- paste0(mzIPdc, collapse = ',')
        
        if(any(s.cln$Int_N[s.cln$mzIP==2] < 80 & s.cln$mzIP>2 & mzIPdc>2 & s.cln$Int_N>80)){
          s.cln <- s.cln[-which(s.cln$mzIP>2 & s.cln$Int_N>80),]
        }
        
        #Method 2 - % RSD in spectra is high > 1
        s.RSD <- sd(s.cln$Int_N)/mean(s.cln$Int_N)
        if(s.RSD < 0.9 & MI < 800){decIso <- F}
        
        #Method 3 - mzIP position 1 intensity must make sense for mz(MI)
        M1Cutoff <- MI/(12*0.9259)+30
        idx_mzIP1 <- s.cln$mzIP %in% 1
        if(any(idx_mzIP1)){
          if(any(s.cln$Int_N[idx_mzIP1] > M1Cutoff)){
            paste0(ft.comment,'| Noisy spectrum ')
            decIso <- F
          }
        }
      }else # switch off deciso if M+1 position does not make sense
      {
        ft.comment <- paste0(ft.comment,'| Single mz spectrum ')
        decIso <- F
        MFSolved <- F
        txt_mzIPdc <- NA
      } #switch off decIso
      
      
      ## Ideal case where n of mzs at least two
      if(n_mz > 1 && decIso==T){
        #####-----------------------------------------------------------------------------------
        ### DeNovo MF List generation from mass spectrum decomposition
        
        # Check MI Mass
        if(MI != s.mdf$MI){
          MI <- s.mdf$MI
          ft.comment <- paste0(ft.comment, '| MI Updated ')
        }
        
        #Print warnings
        if(length(mzIPdc >0)){print(paste0('Warning: mzIP discontinuity detected at mzIP = ', txt_mzIPdc))}
        if(s.mdf$SPC<0.5){print(paste0('Warning: Spectral clarity = ',round(s.mdf$SPC,2)))}
        
        nMFs.st <- Sys.time()
        
        ## Generate MF using spectrum to support element isotope pattern
        MFList1 <- SpectrumToMFList2(s.unk = s.cln, Q1ref=MI, q = -1, eList1 = eList1, ppm=ppmTol - MI/100, maxCounts = F)
        if(MI <= 800){
          MFList2 <- SpectrumToMFList2(s.unk = s.cln, Q1ref=MI, q = -1, eList1 = eList2, ppm=ppmTol, maxCounts = F) 
          MFList <- rbind(MFList1,MFList2)
          rm(MFList1,MFList2) #Trash
        }else{ MFList <- MFList1}
        
        if(length(MFList)==0) {MFList <- NULL} else if (nrow(MFList)==0) {MFList <- NULL}
        
        # Deduplication, Senior Filter, NRule Filter
        while(!is.null(MFList) && MFEscape==F){
          MFList <- MFList[!duplicated(MFList$MF),]  #remove duplicates
          numMFs <- nrow(MFList)
          
          # Logic for choosing Senior3 cutoff
          filt.Senior <-  round(MI/20) #default
          idx_Senior <- which(MFList$SENIOR3 <= filt.Senior & MFList$SENIOR3>=0)
          if(length(idx_Senior)>0){
            MFList <- MFList[idx_Senior,]  #pre-filter so you don't NRule validate Junk based on Senior's Rule
          }
          
          MFList <- MFnRuleValidate(MFList, Type=1)  # Only Validate MFs for nitrogen rule on charged species
          MFList <- MFList[which(MFList$nrule=='Valid'),]  #Apply NRule
          MFEscape <- T
        }
        MFEscape <- F #reset
        
        # Continue to filter using elemental heuristics, and scoring
        if(!is.null(MFList)){
          
          ## Heuristic based on literature
          ## https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-105/tables/2
          MFList <- MFList[which(  (MFList$HFClBrItoC <= 6 | is.na(MFList$HFClBrItoC)) & 
                                     (MFList$FtoC <= 6 | is.na(MFList$FtoC)) &
                                     (MFList$CltoC <= 2 | is.na(MFList$CltoC)) &
                                     (MFList$NtoC <= 0.5 | is.na(MFList$NtoC)) &
                                     (MFList$OtoC <= 1 | is.na(MFList$OtoC) | is.na(MFList$OCS)) &
                                     (MFList$OCS <= 1 | is.na(MFList$OCS)) &
                                     (MFList$PtoC <= 0.3 | is.na(MFList$PtoC)) &  #Exception
                                     (MFList$StoC <=3 | is.na(MFList$StoC)) &
                                     (MFList$OStoP <= 3 | is.na(MFList$OStoP))
          ),]   # Keep only non-Invalids
          
          ## Use 99.7% heuristics if MFList is large
          if(nrow(MFList)>2000){
            MFList <- MFList[which(  (MFList$HFClBrItoC <= 3.1 | is.na(MFList$HFClBrItoC)) & 
                                       (MFList$HFClBrItoC >= 0.2 | is.na(MFList$HFClBrItoC)) & 
                                       (MFList$FtoC <= 3.5 | is.na(MFList$FtoC)) &   ## Modified 99.7%
                                       (MFList$CltoC <= 0.8 | is.na(MFList$CltoC)) &
                                       (MFList$NtoC <= 0.5 | is.na(MFList$NtoC)) &
                                       (MFList$OtoC <= 1 | is.na(MFList$OtoC) | is.na(MFList$OCS)) &
                                       (MFList$OCS <= 1 | is.na(MFList$OCS)) &
                                       (MFList$PtoC <= 0.15 | is.na(MFList$PtoC)) & 
                                       (MFList$StoC <=0.8 | is.na(MFList$StoC)) &
                                       (MFList$OStoP <= 0.7 | is.na(MFList$OStoP))
            ),] } # Higher filter if n_MFs > 100
          
          
          ##----------------------------------------------- 
          #### Cycle through MFList and score results
          # Revisit this (sorted??)
          n.MFs.to.score <- min(nrow(MFList), 50) #Choose how many MFs to calculate scores
          numMFs <- n.MFs.to.score
          if(n.MFs.to.score>0){
            
            ft.FormulaMatches <- data.frame()
            for(mf in 1:n.MFs.to.score){
              ft.comment <- ''
              
              s <- s.cln
              # Setup Reference Mass Spectrum
              if(is.na(ft$Adduct) || ft$Adduct==''){RefIon='[M-H]-'}else{RefIon=ft$Adduct}
              ion_op <- merge(RefIon, adducts, by.x='x', by='Adduct')[,c('a','d')]
              MF.ion <- ionizeMF(m=MFList$MF[mf], a=ion_op$a, d=ion_op$d)
              
              #if(mf<=MF_topN){topNplot=T}else{topNplot=F}
              
              if(!is.na(MF.ion$Mion)){
                # Generate Ref Spec
                # MF.ion$Mion <- 'C18H33O2'
                s.ref <- MFtoSpectrum(formula=MF.ion$Mion, q=q, Res=F, ppmTol, patOnly=T)
                colnames(s.ref) <- c('mz','Int')  #rename columns
                s.ref <- cbind(s.ref, mzIP=round(s.ref$mz - s.ref$mz[1])) #Sets mzIP 
                s.ref <- merge(aggregate(Int ~ mzIP, s.ref, FUN=max), s.ref)
                s.ref <- s.ref[,c('mz','Int','mzIP')]
                
                
                s.unk.ol <- s.cln[which(s.cln$mz > (min(s.ref$mz)-0.1) & s.cln$mz < (max(s.ref$mz)+0.1)),] #Later need to refine spectrum for charge
                # s.unk.ol$Int_N <- 100*s.unk.ol$Int_N/max(s.unk.ol$Int_N)
                
                #Check overlap extent.  Is spec "close enough" to compare
                idx_rows <- 1:min(nrow(s.unk.ol),nrow(s.ref))
                check.ppm <- length(which(abs(s.unk.ol[idx_rows,'mz'] - s.ref[idx_rows,'mz'])*1e6/MI < ppmTol*2)) == min(nrow(s.unk.ol),nrow(s.ref))
                
                checkOL <- abs(s.ref$mz[s.ref$mzIP==0] - s.unk.ol$mz[s.unk.ol$mzIP==0])
                
                if(check.ppm == T || (checkOL[1] < 1) && (checkOL[1] > 0)){ 
                  #SPS <- SpectrumSimilarity_custom(s.ref, s.unk.ol, dppm=ppmTol*z*2, b = 0.1, int.prec=2,
                  #                                 bottom.label = paste0(round(pk$m.z,4),' @',round(ft$Retention.Time, 2),' min'), 
                  #                                 top.label = paste('Simulated',MF.ion$Mion,Poltxt), 
                  #                                 xlim = c(min(s.ref$mz[which(s.ref$Int>1)])-0.2, max(s.ref$mz[which(s.ref$Int>1)])+0.2), 
                  #                                 x.threshold = 0, print.alignment = F, print.graphic = T, output.list = T)
                  
                  
                  ft.MFScores <- SpectrumSimilarity_custom2(obs=s.unk.ol, the=s.ref, dppm=ppmTol, int_prec = 0.3, rnd_prec = 5)
                  
                  rnf <- max(s.unk.ol$Int) #renormalization factor
                  rpf <- nrow(ft.MFScores$Alignment)# rows per feature
                  ft.FormulaMatches <- rbind(ft.FormulaMatches, cbind(mf, f, Feature = pk$Feature,  
                                                                      MZ = ft.MFScores$Alignment$m_obs, Abundance=ft.MFScores$Alignment$i_obs*rnf,
                                                                      Isotope = seq(1:rpf), Formula = MF.ion$MFinput, zFormula=MF.ion$Mion,
                                                                      PredMZ = ft.MFScores$Alignment$m_ref, PredAbundance=ft.MFScores$Alignment$i_ref*rnf,
                                                                      File = pk$File,
                                                                      Score1 = ft.MFScores$Score1, Score2=ft.MFScores$Score2, Comment=ft.comment) )
                  
                }else{
                  ft.comment <- paste0(ft.comment,'| Spectrum did not overlap Ref ')
                  ft.FormulaMatches <- rbind(ft.FormulaMatches, cbind(mf, f, Feature = pk$Feature,  
                                                                      MZ = ft.MFScores$Alignment$m_obs, Abundance=ft.MFScores$Alignment$i_obs*rnf,
                                                                      Isotope = seq(1:rpf), Formula = MF.ion$MFinput, zFormula=MF.ion$Mion,
                                                                      PredMZ = NA, PredAbundance=NA,
                                                                      File = pk$File,
                                                                      Score1 = '0', Score2='0', Comment=ft.comment) )
                } ## If overlap not present
              }else{
                ft.comment <- paste0(ft.comment,'| Cannot ionize MF ')
                ft.FormulaMatches <- rbind(ft.FormulaMatches, cbind(mf, f, Feature = pk$Feature,  
                                                                    MZ = ft.MFScores$Alignment$m_obs, Abundance=ft.MFScores$Alignment$i_obs*rnf,
                                                                    Isotope = seq(1:rpf), Formula = MF.ion$MFinput, zFormula=NA,
                                                                    PredMZ = NA, PredAbundance=NA,
                                                                    File = pk$File,
                                                                    Score1 = '0', Score2='0', Comment=ft.comment) )
              }  ## If Ionization failed
            } # end for loop
            
            #Add 10% penalty for Phosphorus
            idx_MFs_P <- which(grepl('P',ft.FormulaMatches$Formula))
            if(any(idx_MFs_P)){ft.FormulaMatches$Score2[idx_MFs_P] <- as.numeric(ft.FormulaMatches$Score2[idx_MFs_P]) * 0.9 }
            
            ft.FormulaMatches <- ft.FormulaMatches[order(-as.numeric(ft.FormulaMatches$Score2)),]
            ## Re-rank based on Score2
            ft.FormulaMatches <- merge(ft.FormulaMatches, data.frame(mf=unique(ft.FormulaMatches$mf), 
                                                                     Rank=seq(1:mf)), by='mf')
            #Keep TopN only
            ft.FormulaMatches <- ft.FormulaMatches[ft.FormulaMatches$Rank<=MF_topN, ]
            ft.FormulaMatches <- ft.FormulaMatches[order(ft.FormulaMatches$Rank),]
            
            nMFs.et <- Sys.time()
            nMFs.tt <- round(as.numeric (nMFs.et - nMFs.st, units = "secs"), 1)
            MFSolved <- T
            
          }else{
            ft.comment <- paste0(ft.comment,'| No valid MFs found for spectrum ')
            decIso <- F # Try decomposeMass
            MFSolved <- F
          }
          
        }else{
          ft.comment <- paste0(ft.comment,'| Unable to predict MF on spectrum ')
          decIso <- F # Try decomposeMass
          MFSolved <- F
        } #end if nSP >=1
      } #end decompose isotopes
      
      ## If there is only 1 mz or decIso got switched to false, we can only do decompose Mass
      if((n_mz==1 || decIso==F) && MFSolved==F){
        
        nMFs.st <- Sys.time()
        # Case where only 1 mz peak - Can't have Br or Cl
        if(n_mz==1){
          eList_t <- eList1[-which(eList1=='Br' | eList1=='Cl')]
        }else{
          eList_t <- eList1
        }
        
        #Predict MF on neutral, single mz
        MFList <- calcMF2(mz = MI + q*5.48579909070e-4 - q*1.007276, z = 0, ppm = ppmTol, top = NULL, elements = c(Rdisop::initializeElements(eList_t)), maxCounts = F,
                          SeniorRule = T, HCratio = T, moreRatios = T, elementHeuristic = T, 
                          Filters = list(maxCounts =F, SENIOR3 = F, HCratio = F, 
                                         moreRatios = F, elementHeuristic = F), summarize = F, BPPARAM = NULL)
        
        filt.Senior <-  round(MI/20) #default
        idx_Senior <- which(MFList$SENIOR3 <= filt.Senior & MFList$SENIOR3>=0)
        #pre-filter so you don't NRule validate Junk based on Senior's Rule
        if(length(idx_Senior)>0){MFList <- MFList[idx_Senior,]}
        
        #Apply Filter BEFORE NRule Validate
        MFList <- MFList[which( MFList$SENIOR3 < filt.Senior & MFList$SENIOR3 >= 0 & MFList$unsat > -10 &
                                  (MFList$HFClBrItoC <= 6 | is.na(MFList$HFClBrItoC)) & 
                                  (MFList$FtoC <= 6 | is.na(MFList$FtoC)) &
                                  (MFList$CltoC <= 2 | is.na(MFList$CltoC)) &
                                  (MFList$NtoC <= 0.5 | is.na(MFList$NtoC)) &
                                  (MFList$OtoC <= 1 | is.na(MFList$OtoC) | is.na(MFList$OCS)) &
                                  (MFList$OCS <= 1 | is.na(MFList$OCS)) &
                                  (MFList$PtoC <= 2 | is.na(MFList$PtoC)) & 
                                  (MFList$StoC <=3 | is.na(MFList$StoC)) &
                                  (MFList$OStoP <= 3 | is.na(MFList$OStoP))
        ),] 
        
        if(!is.null(MFList) && nrow(MFList)>0){
          MFList <- MFnRuleValidate(MFList, Type=1) 
          MFList <- MFList[which(MFList$nrule=='Valid'),]
        }
        
        while((!is.null(MFList) && nrow(MFList)>0) && MFEscape==F){
          
          #Simple heuristic to choose MF with RDBE closest to 0, take first one(lowest ppm) if two exist
          MFList <- MFList[order(abs(MFList$unsat), abs(MFList$ppm)),]
          numMFs <- nrow(MFList) # Get this number before storing TopN
          
          #Final MFList
          MFList <- MFList[1:(min(MF_topN,numMFs)),]
          
          #Get MF ion formulas
          ion_op <- merge(RefIon, adducts, by.x='x', by.y='Adduct')[,c('a','d')]
          MF.ion <- lapply(X=MFList$MF, FUN=ionizeMF,  a=ion_op$a, d=ion_op$d)
          MF.ion <- as.data.frame(t(do.call(cbind, MF.ion)))
          
          nMFs.et <- Sys.time()
          nMFs.tt <- round(as.numeric (nMFs.et - nMFs.st, units = "secs"), 1)
          MFEscape <- T
        }
        MFEscape <- F #reset
        
        
        # Basic final score assigned only on mass error and if present
        if(!is.null(MFList) && nrow(MFList)>0){
          mzScore <- 1 - abs(MFList$ppm)/ppmTol
          MFSolved <- T
          ft.FormulaMatches <- data.frame(mf = seq(1:nrow(MFList)), f, Feature = ft$row.ID, 
                                          MZ = s.cln[1,'mz'], Abundance=s.cln[1,'Int'],
                                          Isotope = 1, Formula = MFList$MF, zFormula=unlist(MF.ion$Mion),
                                          PredMZ = MFList$mz, PredAbundance=s.cln[1,'Int'],
                                          File = pk$File,Score1 = 0, Score2=mzScore, 
                                          Comment=ft.comment, Rank = seq(1:nrow(MFList)))
        }else{ #end final score and store of mz=1
          ft.comment <- 'No MF from spectrum or single mz'
          MFSolved <- F
        }
        
      }# End if 1 mz value in spec
      
    }else{ #End Feature in Master table present
      MFSolved <- F
      
      if(nrow(ft)==0){ft.comment <- 'Spectrum not linked to a feature'}
      if(pk$m.z > 1000){ft.comment <- 'High mz MF not predicted'}
    } #else no spectrum found
    
    
    # Write MFList based on variables
    if(MFSolved==T){
      
      dn.FormulaMatches <- rbind(dn.FormulaMatches, ft.FormulaMatches)
      
      if(mf.stor==T){
        idx_ft_toprow <- which(ft.FormulaMatches$Isotope==1)
        pks_MS1_MF <- rbind(pks_MS1_MF, cbind(Feature=pk$Feature, Src='DN', MF=ft.FormulaMatches$Formula[idx_ft_toprow], MFscore=ft.FormulaMatches$Score2[idx_ft_toprow], MFRank=ft.FormulaMatches$Rank[idx_ft_toprow], nMFs=numMFs))
        pks_Top_MFs <- rbind(pks_Top_MFs, data.frame(Feature=pk$Feature, TopMF=ft.FormulaMatches$Formula[1], TopMFScore=ft.FormulaMatches$Score2[1],  SPC=s.mdf$SPC, numMFs, nMFs.tt, mzIPdc=txt_mzIPdc) )
      }
    }
    
    if(MFSolved==F){
      dn.FormulaMatches <- rbind(dn.FormulaMatches, cbind(mf=NA, f, Feature = ft$row.ID, 
                                                          MZ = ft$m.z, Abundance=NA, 
                                                          Isotope=NA, Formula = NA, zFormula=NA,
                                                          PredMZ=NA, PredAbundance = NA, File=ft$Files.1[1],
                                                          Score1=NA, Score2=NA, Comment=ft.comment, Rank=NA) ) 
      
      if(mf.stor==T){
        pks_MS1_MF <- rbind(pks_MS1_MF, data.frame(Feature=pk$Feature, Src='DN', MF=NA, MFscore=NA, MFRank=NA, nMFs=0 ))
        pks_Top_MFs <- rbind(pks_Top_MFs, data.frame(Feature=pk$Feature, TopMF=NA, TopMFScore=NA,  SPC=s.mdf$SPC, numMFs=0, nMFs.tt=NA,mzIPdc=txt_mzIPdc) )
      } 
    }
  } #end Feature List cycle
  
  
  ###------------  Start reshaping the outputs to insert in NegIDed_FIN.csv ----------------###
  
  
  ## Make Neg_Feature_pFormula_MS1s.csv
  
  temp_MF_Denovo <- merge(pks_Top_MFs, Feature_MF_List[,1:15], by.x='Feature', by.y='row.ID')
  ## Evaluate Match Rates
  MF.match <- NA
  for(m in 1:nrow(temp_MF_Denovo)){
    MF.match[m] <- compareMFs(c(temp_MF_Denovo$TopMF[m], temp_MF_Denovo$Formula[m]))
  }
  
  temp_MF_Denovo <- cbind(isMatch=MF.match, temp_MF_Denovo)
  
  MF_Denovo <- rbind(MF_Denovo, temp_MF_Denovo)
  MF_Denovo$isMatch[which(MF_Denovo$Formula=='')] <- NA
  MF_Denovo <- MF_Denovo[which(!duplicated(MF_Denovo)),]
  # match.rate <- length(which(MF_Denovo$isMatch==T))/length(which(MF_Denovo$isMatch==T | MF_Denovo$isMatch==F))
  write.csv(dn.FormulaMatches, file=fh_MS1, row.names = T)
  
  
  ## Create DF that will be inserted into NegIDed_FIN.csv
  
  n.f <- nrow(pks_FeatureIDs)
  idx_insert_colA <- which(colnames(pks_FeatureIDs)=='Needs_Validation')
  idx_insert_colB <- which(colnames(pks_FeatureIDs)=='row.ID.1')  #added DIS 3/19/24
  
  ft.ncol <- ncol(pks_FeatureIDs)
  MS1_Features_MF <- data.frame()
  for(ft.i in 1:n.f){
    ft <- pks_FeatureIDs[ft.i,]
    ft.dnMFs.all <- dn.FormulaMatches[dn.FormulaMatches$Feature==ft$row.ID ,]
    ft.dnMFs.top <- ft.dnMFs.all[ft.dnMFs.all$Rank==1,]
    
    if(nrow(ft.dnMFs.all)>0){
      dl.mz <- paste0(round(as.numeric(ft.dnMFs.top$PredMZ),4), sep=';', collapse='')
      dl.int <- paste0(round(as.numeric(ft.dnMFs.top$PredAbundance)), sep=';', collapse='')
      MFList <- ft.dnMFs.all[,c('Formula','Score1','Score2','Rank')]
      MFList <- MFList[!duplicated(MFList),]
      
      # Check to see if pMF is correct
      dn.MFMatch <- NA
      if(!is.na(ft$Formula)){  
        for(m in 1:nrow(MFList)){
          if(compareMFs(c(MFList$Formula[m], ft$Formula))){
            dn.MFMatch = as.numeric(MFList$Rank[m])
            break 
          }
          else{dn.MFMatch <- NA}
        }
      }
      
      dl.MFlist <- paste0(MFList$Formula, sep=';', collapse='')
      dl.MFScore2 <- paste0(MFList$Score2, sep=';', collapse='')
      #dl.Rank <- paste0(MFList$Rank, sep=';', collapse='')
      dl.Rank <- nrow(MFList)
      MS1_Features_MF <- rbind(MS1_Features_MF, 
                               cbind(TopMF=ft.dnMFs.top$Formula[1],
                                     dn_MFMatch = dn.MFMatch,
                                     zFormula = ft.dnMFs.top$zFormula[1],
                                     #Pred_mzs = dl.mz,
                                     #Pred_ints = dl.int,
                                     MFscore = ft.dnMFs.top$Score2[1],
                                     MFComment = ft.dnMFs.top$Comment[1],
                                     dn_MFRanks = dl.Rank,
                                     dn_MFs = dl.MFlist,
                                     dn_MFScores = dl.MFScore2
                               ))
    }else{
      MS1_Features_MF <- rbind(MS1_Features_MF,  
                               cbind( TopMF=NA,
                                      dn_MFMatch = NA,
                                      zFormula = NA,
                                      #Pred_mzs = NA,
                                      #Pred_ints = NA,
                                      MFscore = NA,
                                      MFComment = 'No MFs predicted',
                                      dn_MFRanks = NA,
                                      dn_MFs = NA,
                                      dn_MFScores = NA
                               ))
    }
    
    if(round(ft.i,-2)==ft.i){   ## Show processed peaks in 100
      print(paste('Processing Feature',ft.i,'of',n.f))
    }
    
  }
  
  # MS1_Features_MF gets sandwiched in between "Needs_Validation" and "row.ID.1" columns
  final_MS1_df <- cbind(pks_FeatureIDs[,1:idx_insert_colA],  MS1_Features_MF, pks_FeatureIDs[,(idx_insert_colB):ft.ncol] )
  write.csv(final_MS1_df, file=fh_IDed_FIN, row.names = F)
  
} ## End of main MFPredict function

# Apply Senior Rule to MF
getSenior3.MFobject2 <- function (MF){
  sum(MF * getOption("MassTools.maxValences")) - 2 * (sum(MF) -  1)
}


# MD calculator
mzdefect <- function(x, c=0, p=3){
  d = x - floor(x+c)
  d = round(d, digits=p)
  return(d)
}

# Estimate # carbons based on mzIP=1
M1Int2NC <- function(I){
  nC <- round(I*0.92379 + 0.002309, digits=2)
  return(nC)
}


mzSpec_MDFilter2 <- function(spectrum, Q1, mdTol=NULL, minIntN=1, IntPrec=0, IntNPrec=2, z=1, splot=T, forceQ1=T, labelBG=F, plotBG=F){
  
  mzIntList <- spectrum
  colnames(mzIntList) <- c('mz','Int')  #rename columns
  mzIntList <- mzIntList[order(mzIntList[,'mz']),]  #order by Q1
  
  #Prefilter
  mzIntList <- mzIntList[which(abs(mzIntList[,'mz'] - Q1)*z < 5),]  #Upper Q1 limit of M+5
  
  if(is.null(mdTol)){mdTol = 0.25*1e6*z/(Q1*sqrt(Q1))}    # empirically derived
  # Note MD ppm tol should be > Mass ppm error to reduce over filtering
  
  mzIP <- round((mzIntList[,'mz'] - Q1),1) #mz integer position
  MD <- mzdefect(mzIntList[,'mz']) #mz defect
  MDz <- mzdefect(mzIntList[,'mz'] * z) #mass defect
  dMDz <- MDz - mzdefect(Q1*z) #delta mass defect in daltons
  idx1 <- which(dMDz < -0.5) #left side
  idx2 <- which(dMDz >  0.5) #right side
  if(length(idx1)>0){dMDz[idx1] <- dMDz[idx1] + 1} #Shift mzIP ?
  if(length(idx2)>0){dMDz[idx2] <- dMDz[idx2] - 1}
  dMDz <- round(dMDz*1e6/Q1)
  adMDz <- abs(dMDz)
  
  mzIntList <- cbind(mzIntList, mzIP, MD, MDz, dMDz, adMDz)
  
  #idx_MI <- min(which(mzIntList$mzIP==0 & mzIntList$adMDz==min(adMDz))) #Original Idx of MI
  idx_MI <- which(mzIntList$Int == max(mzIntList$Int[which(mzIntList$adMDz<=mdTol & abs(mzIntList$mzIP)<=1)]))
  
  #Find potential MI related to FFT distortion.
  if(forceQ1==F){
    idx_MI_dist <- which(mzIntList$Int > 10*mzIntList$Int[idx_MI] & mzIntList$dMDz < 100)  
    if(length(idx_MI_dist)>0){
      print(paste0('FFT mz distortion detected for  ',Q1))
      idx_MI <- min(idx_MI_dist)
    }
  }
  
  #Re-center mzIP based on MI
  mzIntList$mzIP <- mzIntList$mzIP - mzIntList$mzIP[idx_MI]
  mzIntList$dMDz <- mzIntList$dMDz - mzIntList$dMDz[idx_MI]
  mzIntList$adMDz <- abs(mzIntList$dMDz)
  
  
  idx <- which(abs(mzIntList$adMDz) < 40 & mzIntList$mzIP*z>=-2 & mzIntList$mzIP*z<=4)
  Target <- mzIntList[idx,]  #Apply MD Filter
  Background <- mzIntList[-idx,]
  
  #Calculate spectral clarity
  SPC <-  sum(Target$Int)/(sum(Target$Int) + sum(Background$Int))
  
  Q1new <- mzIntList[,'mz'][idx_MI]
  
  Target <- cbind(Target, Class=rep('_Target',length(Target$mz)))
  Background <- cbind(Background, Class=rep('Background',length(Background$mz)))
  
  specOL <- data.frame()
  specOL <- rbind(Background,Target)
  
  ##  Normalized Intensity Filter
  Int_N = specOL$Int*100/max(Target$Int[which(Target$mzIP<2)])  #Calc
  Int_N = round(Int_N, IntNPrec)  # Round
  specOL <- cbind(specOL, Int_N)   # Join
  specOL <- specOL[which(Int_N>minIntN),]    #Apply
  
  label <- as.character(round(specOL$mz,4))
  specOL <- cbind(specOL,label)
  
  if(splot==T){
    xlims <- c(round(min(specOL$mz)-3/z),round(max(specOL$mz)+3/z))
    ylims <- c(0,1.03*max(Target$Int))
    
    title=paste('MS1 Scan of ',round(Q1,4))
    
    lbl_specOL <- specOL
    if(labelBG==F){lbl_specOL <- subset(specOL, Class %in% '_Target')}
    if(plotBG==F){specOL <- lbl_specOL}
    
    gg <- ggplot()
    gg <- gg + labs(title=title, x='m/z', y='Response') + lims(x=xlims,y=ylims)
    gg <- gg + geom_segment(specOL, mapping = aes(x=mz, y=Int, xend = mz, yend = 0, colour = Class), linewidth=0.75)
    if(plotBG==T){
      gg <- gg + geom_text(lbl_specOL, mapping = aes(x=mz, y=Int + 0.02*max(Int), colour = Class, label=label), size=4)
    }else{
      gg <- gg + geom_text(lbl_specOL, mapping = aes(x=mz, y=Int + 0.02*max(Int), colour = Class, label=label), size=4)
    }
    return(list(Q1in=Q1, MI=Q1new, Spectrum=specOL, SPC=SPC, Plot=gg))
  }else{
    return(list(Q1in=Q1, MI=Q1new, Spectrum=specOL, SPC=SPC))
  }
  
}


SpectrumToMFList2 <- function(s.unk, Q1ref, q = -1, ppm = 10,top = NULL, nC=NULL, noC=F,
                              eList1 = c('C','H','N','O','P','S'), eList2=NULL,
                              maxCounts = T,
                              SeniorRule = T,
                              HCratio = T,
                              moreRatios = T,
                              elementHeuristic = T,
                              Filters = list(DBErange = c(-5,40),
                                             minElements = NULL,
                                             maxElements = NULL,
                                             parity = NULL,
                                             maxCounts = T,
                                             SENIOR3 = 0,
                                             HCratio = F,
                                             moreRatios = F,
                                             elementHeuristic = F),
                              summarize = F,
                              BPPARAM = NULL
){
  
  z=abs(q)
  
  s.unk.t <- s.unk
  #Charge transform all m/z, mzIP, Q1ref
  s.unk.t[,'mz'] <- s.unk.t[,'mz']*z
  s.unk.t$mzIP <- s.unk.t$mzIP*z  
  Q1ref.t = Q1ref*z  
  
  s.unk.t <- s.unk.t[which(s.unk.t[,'mz']>(Q1ref.t-5) & s.unk.t[,'mz']<(Q1ref.t+5)),]  #Refilter
  s.unk.t$Int_N  <-  s.unk.t$Int_N*100/max(s.unk.t$Int_N)
  
  idx_mzPos1 <- which(s.unk.t$mzIP==1)
  idx_mzPos1 <- idx_mzPos1[which(s.unk.t$Int_N[idx_mzPos1] == max(s.unk.t$Int_N[idx_mzPos1]))]
  
  ## The first section of this builds an element list (eList) based on isotopic abundances in the spectrum. 
  ## It first looks for carbon in the mzIP = 1 position, then uses a special function to detect Cl, Br, and S in mzIP=2 position
  ## Added a check for silicon (siloxanes).  Occurs when M+1 is <<< expected
  
  #s.unk <- cbind(s.unk,Q1_dist)  #Attach to DF
  if(noC==F){
    MF_Max <- ''
    n_C_max = ceiling((Q1ref.t-2-14)/14)  # Max nC if fully saturated with H and has either an O or N to ionize
    
    if(length(idx_mzPos1)>0 && s.unk.t[idx_mzPos1,'Int_N'] > 2){
      idx_mzPos1 <- idx_mzPos1[which(s.unk.t$Int_N[idx_mzPos1] == max(s.unk.t$Int_N[idx_mzPos1]))]
      nCf = max(s.unk.t[idx_mzPos1,'Int_N'])
      
      n_C_E = round(nCf*(0.9177) - 0.1546)  # Exact nC based on M+1
      #n_Si_hi = round(n_C_E/4.7/2)  # Check M+1 is too low to have n_C_max, 4.7 = ratio of C13 to Si29
      
      #if(n_Si_hi > 0){
      #  eList1 <- c(eList1,c('Si'))
      #  eList2 <- c(eList2,c('Si'))
      #  MF_Max <-  paste0(MF_Max,'Si',n_Si_hi)
      #}
      if(nCf < n_C_max){
        n_C_ff = 0.15 - Q1ref.t/8000  #fudge factor
        n_C_lo <- round(nCf*(0.9177 - n_C_ff) - 0.1546)
        n_C_hi <- round(nCf*(0.9177 + n_C_ff) - 0.1546)
      }else{
        idx_mzPos2 <- which(s.unk.t$mzIP==2)
        if(length(idx_mzPos2)>0){
          nCf = max(s.unk.t[idx_mzPos2,'Int_N'])
          n_C_lo <- round(nCf*(0.9177 - 3/50) - 0.1546)
        }else{n_C_lo = n_C_max - 10}
        n_C_hi = n_C_max
      }
    }else{  # Only when M+1 Int is <2 % is nC very small
      n_C_lo=0
      n_C_hi=5
    }
    if(!is.null(nC)){  # Set bounds based on nC input variable, allows direct specification of nC
      n_C_lo=nC-1
      n_C_hi=nC+1
    }
    # if low and hi are equal set hi to +/- 2 of that value
    if(n_C_hi == n_C_lo){
      n_C_hi = n_C_hi + 4
      n_C_lo = n_C_lo - 4
    } 
    # if low and hi are <= 5 apart set to +/- 2.5 of the mean of lo and hi
    if(n_C_hi - n_C_lo <= 5){
      n_C_mean = round(mean(n_C_hi,n_C_lo))
      n_C_hi = n_C_mean + 4
      n_C_lo = n_C_mean - 8
    } 
    # Need to correct in case M+1 peak is biased high (noise)
    if(n_C_lo > n_C_hi){
      n_C_lo <- n_C_hi - 5
    } 
    
    MF_Min <- paste0('C',max(n_C_lo,0),'H1')
    MF_Max <-  paste0('C',n_C_hi,MF_Max,'H200')
  }else{
    MF_Min <- 'H0'
    MF_Max <- 'H999'
    idx_C <- which(eList1=='C')
    if(length(idx_C) >0){eList1=eList1[-idx_C]}  #Remove if not in spec
    idx_C <- which(eList2=='C')
    if(length(idx_C) >0){eList2=eList2[-idx_C]}  #Remove if not in spec
    n_C_lo <- 0
    n_C_hi <- 1
  }
  
  if(n_C_lo<0){n_C_lo <- 0} #Lower bound must not be negative
  
  # Halogen/Special atom detection in spectrum
  Hal <- detectBrCl(h.unk=s.unk.t,Q1ref.t, z=1)   #Charge already handled
  n_Hal <- length(Hal)
  
  if(length(which(eList1=='Cl' | eList1=='Br' | eList1=='S')) !=0 ){
    for(h in 1:n_Hal){
      h.e <- colnames(Hal[h])
      idx_El <- which(eList1==h.e)
      if(length(idx_El) >0 && Hal[h] <1){eList1=eList1[-idx_El]}  #Remove if not in spec
      if(length(idx_El)==0 && Hal[h] >0){eList1=c(eList1,h.e)}  #Add if not in eList
      if(Hal[h] > 0){
        MF_Max <- paste0(MF_Max,h.e,Hal[h])}  #Set MFMax bounds 
    }
  }
  
  if(Hal$Fe==0){s.unk.t <- s.unk.t[which(s.unk.t$mzIP>=0),]}
  
  #Initialize Elements
  elements <- Rdisop::initializeElements(eList1)
  
  ##  Setup Filters 
  if(is.null(Filters$minElements) || !is.character(Filters$minElements)
     || is.na(Filters$minElements) || Filters$minElements == ""){
    Filters$minElements <- MF_Min
  }
  if(is.null(Filters$maxElements) || !is.character(Filters$maxElements)
     || is.na(Filters$maxElements) || Filters$maxElements == ""){
    Filters$maxElements <- MF_Max
    
  }
  
  if(is.null(elements) || ( Q1ref.t>800 )){
    
    print('Large MF --> Element List defaulted CHNOPSF')
    
    MF_Max <-  paste0('C',round(Q1ref.t/35))
    MF_Min <-  paste0('C',round(Q1ref.t/55))
    Filters$minElements <- MF_Min
    Filters$maxElements <- MF_Max
    
    elements = Rdisop::initializeElements(c('C','H','O','N','S','P','F'))
    
  }else if(is.character(elements)){
    #if character formula is given, e.g. "CHOPS"
    
    elements <- MassTools::makeMF(elements)
    
    elements <- Rdisop::initializeElements( names(elements)[elements>0])
    
    #print(paste0('MF range ',MF_Min, ' to ',MF_Max))
  }
  
  print(paste0('MF range ',MF_Min, ' to ',MF_Max))
  
  ##  Setup DBE Pre-Filter
  Filters$DBErange = c(-5,round(Q1ref.t/16))
  
  #how to handle it if there are no molecular formulas given filter conditions
  if(summarize){
    failReturn <- "No MF within bounds" 
  }else{
    failReturn <- NULL
  }
  
  
  #have to add an electron mass to mz for each charge to get uncharged mass 
  #[needed for correct ppm calculation in decomposeMass; maybe breaks nitrogen rule calculation]
  mm <- Rdisop::decomposeIsotopes(s.unk.t[,'mz'] + q*5.48579909070e-4 - q*1.007276, 
                                  s.unk.t$Int_N, 
                                  z = 0, 
                                  maxisotopes = 1, 
                                  ppm = ppm, 
                                  mzabs = 0,
                                  elements = elements,
                                  minElements = Filters$minElements,
                                  maxElements = Filters$maxElements)
  
  if(is.null(mm)){return(failReturn)}
  
  f1 <- if(!is.null(Filters$DBErange)){mm$DBE >= Filters$DBErange[1] & mm$DBE<= Filters$DBErange[2]}else{rep(TRUE,length(mm$DBE))}
  
  if(!is.null(Filters$parity) && Filters$parity %in% c('e','o')){
    f1 <- f1 & mm$parity == Filters$parity
  }
  
  if(!any(f1)){return(failReturn)}
  
  
  f1 <- which(f1)
  
  if(length(f1) == 0){return(failReturn)}
  
  #now reorder by mass difference
  f2 <- f1[order(abs( mm$exactmass[f1] - (Q1ref - q*1.00727647) ))]
  #f2 <- f1[order(abs(mm$exactmass[f1] - mzList[1]))]
  
  sfs <- MassTools::makeMF(mm$formula[f2], forcelist = TRUE)
  
  
  res <- data.frame(mz = (mm$exactmass[f2] - q*5.48579909070e-4),
                    MF = mm$formula[f2],
                    charge = 0,
                    RdisopScore = mm$score[f2],
                    unsat = mm$DBE[f2],
                    parity = mm$parity[f2],
                    error = (mm$exactmass[f2] - q*5.48579909070e-4) - (Q1ref.t - q*1.00727647),
                    nrule = mm$valid[f2],
                    stringsAsFactors = FALSE)
  
  res$ppm <- res$error*1e6/(Q1ref - q*1.00727647)
  
  if(nrow(res)<5){maxCounts=F}
  
  #Golden rule #1: Restriction for element numbers
  if(maxCounts){
    if(q==0){Da=Q1ref.t}else{Da <- Q1ref.t}
    
    if(Da < 500){
      maxCountLimit <- MassTools::consolidateMF(c(C = 29,H = 3*n_C_hi+1, N= round(n_C_hi*0.33)+1, O = max(round(n_C_hi*0.6),8)+1, P= 1, 
                                                  S= max(round(n_C_hi/10+2),5), F=3*n_C_hi, Cl = 6, Br = 4))
    }else if(Da <1500){
      maxCountLimit <- MassTools::consolidateMF(c(C = 66,H = 126, N= 25, O = 27, P= 2, S= 8, F =40, Cl = 11, Br = 8))
    }else if(Da < 2300){
      maxCountLimit <- MassTools::consolidateMF(c(C = 115,H = 236, N= 32, O = 63, P= 6, S= 8, F =16, Cl = 11, Br = 8))
    }else if(Da < 3000){
      maxCountLimit <- MassTools::consolidateMF(c(C = 162,H = 208, N= 48, O = 78, P= 6, S= 9, F =16, Cl = 11, Br = 8))
    }else{
      maxCountLimit <- NULL
    }
    
    if(!is.null(maxCountLimit)){
      
      maxCountLimit[maxCountLimit == 0] <- Inf
      
      res$maxCounts <- sapply(sfs, "<=", maxCountLimit)
      
      
    }else{
      res$maxCounts <- TRUE
    }
    
    if(!is.null(Filters$maxCounts) && Filters$maxCounts){
      sel <- res$maxCounts
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  #Golden rule #2: (Senior's third theorem): 
  if(SeniorRule){
    
    res$SENIOR3 <- sapply(sfs,getSenior3.MFobject2)
    
    if(!is.null(Filters$SENIOR3) && is.numeric(Filters$SENIOR3)){
      
      sel <- res$SENIOR3 >= Filters$SENIOR3
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  #Golden Rule #4: Hydrogen/Carbon element ratio check
  if(HCratio){
    
    res$HtoC <-  sapply(sfs,function(MF){
      
      MF["H"]/MF["C"]
      
    })
    
    if(!is.null(Filters$HCratio) && Filters$HCratio){
      
      sel <- is.na(res$HtoC) | (res$HtoC > 0.2 & res$HtoC < 3.1)
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
  }
  
  #Golden Rule #5: Heteroatom check
  if(moreRatios){
    res$NtoC <-  sapply(sfs,function(MF){
      MF["N"]/MF["C"]
    })
    res$OtoC <-  sapply(sfs,function(MF){
      MF["O"]/MF["C"]
    })
    res$PtoC <-  sapply(sfs,function(MF){
      MF["P"]/MF["C"]
    })
    res$StoC <-  sapply(sfs,function(MF){
      MF["S"]/MF["C"]
    })
    res$OStoP <-  sapply(sfs,function(MF){
      MF["P"]/(MF["O"]+ MF["S"])
    })
    res$FtoC <-  sapply(sfs,function(MF){
      MF["F"]/MF["C"]
    })
    res$CltoC <-  sapply(sfs,function(MF){
      MF["Cl"]/MF["C"]
    })
    res$BrtoC <-  sapply(sfs,function(MF){
      MF["Br"]/MF["C"]
    })
    
    res$OCS <-  sapply(sfs,function(MF){
      (MF["O"])/((MF["C"])+1+((MF["S"]/MF["S"])*4))
    })    
    
    res$HFClBrItoC <-  sapply(sfs,function(MF){
      (MF["H"]+MF["F"]+MF["Cl"]+MF["Br"]+MF["I"])/MF["C"]  #Typically, halogens replace H so this ratio is useful (ie: for PFAS)
    })
    
    
    
    if(!is.null(Filters$moreRatios) && Filters$moreRatios){
      
      sel <- ((is.na(res$NtoC) | res$NtoC < 1.3)
              & (is.na(res$OtoC) | res$OtoC < 1.2)
              & (is.na(res$PtoC) | res$PtoC < 0.3)
              & (is.na(res$StoC) | res$StoC < 0.8))
      
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  # Golden rule #6: element probability check
  if(elementHeuristic){
    
    res$elementHeuristic = sapply(sfs,function(MF){
      
      r1 <- na.omit(MF[c("N","O","P","S")])
      
      r2 <- r1[r1>1]
      
      #if two of these elements are not in the molecule, don't apply restrictions
      if(length(r2) <= 2){return(TRUE)}
      
      if(length(r2) == 4 
         && (r2["N"] >=10
             || r2["O"] >=20
             ||r2["P"] >=4
             ||r2["S"] >=3)){return(FALSE)}
      
      switch(paste(names(r2), collapse = ""),
             NOP = {if(min(r2) > 3 && (r2["N"] >=11 || r2["O"] >=22 || r2["P"] >=7)){return(FALSE)}},
             OPS = {if(min(r2) > 1 && (r2["S"] >=3 || r2["O"] >=14 || r2["P"] >=3)){return(FALSE)}},                                                                            NPS = {if(min(r2) > 1 && (r2["N"] >=4 || r2["S"] >=3 || r2["P"] >=3)){return(F)}},
             NPS = {if(min(r2) > 1 && (r2["N"] >= 4 || r2["S"] >= 3 || r2["P"] >= 3)) {return(FALSE)}}, 
             NOS = {if(min(r2) > 6 && (r2["S"] >=8 || r2["O"] >=14 || r2["N"] >=19)){return(FALSE)}})
      
      return(TRUE)
      
      
    })
    
    
    if(!is.null(Filters$elementHeuristic) && Filters$elementHeuristic){
      
      sel <- res$elementHeuristic
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  res <- res[order(-res$RdisopScore),]
  
  # If sulfur detected, require sulfur in MF
  if(Hal$S>0){res <- res[grep('S',res$MF),]}
  
  
  #final selection:
  
  if(!is.null(top) && is.numeric(top)){ 
    
    res <- res[seq(min(top,nrow(res))),]
    
  }
  
  if(summarize){
    
    return(paste(paste0(res$MF,
                        if(res$charge[1] > 0){paste0("(+",res$charge,")")}
                        else if(res$charge[1] < 0){paste0("(",res$charge,")")}
                        else{""}),
                 collapse = "|"))
  }
  
  return(res)
}


calcMF2 <- function (mz = 200.000659, z = 1, ppm = 5, top = NULL, elements = Rdisop::initializeCHNOPS(), 
                     maxCounts = TRUE, SeniorRule = TRUE, HCratio = TRUE, moreRatios = TRUE, 
                     elementHeuristic = TRUE, Filters = list(DBErange = c(-5, 
                                                                          40), minElements = "C0", maxElements = "C99999", parity = "e", 
                                                             maxCounts = TRUE, SENIOR3 = 0, HCratio = TRUE, moreRatios = TRUE, 
                                                             elementHeuristic = TRUE), summarize = FALSE, BPPARAM = NULL) 
{
  if (is.null(Filters$minElements) || !is.character(Filters$minElements) || 
      is.na(Filters$minElements) || Filters$minElements == 
      "") {
    Filters$minElements <- "C0"
  }
  if (is.null(Filters$maxElements) || !is.character(Filters$maxElements) || 
      is.na(Filters$maxElements) || Filters$maxElements == 
      "") {
    Filters$maxElements <- "C100"
  }
  if (is.null(elements)) {
    elements = initializeCHNOPS()
  } else if (is.character(elements)) {
    elements <- makeMF(elements)
    elements <- initializeElements(names(elements)[elements > 0])
  }
  if (length(mz) > 1) {
    rl <- bplapply(mz, MassTools::calcMF, z = z, ppm = ppm, top = top, 
                   elements = elements, maxCounts = maxCounts, SeniorRule = SeniorRule, 
                   HCratio = HCratio, moreRatios = moreRatios, elementHeuristic = elementHeuristic, 
                   Filters = Filters, summarize = summarize, BPPARAM = if (is.null(BPPARAM)) {
                     SerialParam()
                   }
                   else {
                     BPPARAM
                   })
    if (summarize) {
      return(unlist(rl))
    }
    else {
      return(rl)
    }
  }
  if (summarize) {
    failReturn <- ""
  }
  else {
    failReturn <- NULL
  }
  mm <- Rdisop::decomposeMass(mz + z * 0.00054857990907, z = z, maxisotopes = 1, 
                              ppm = ppm, mzabs = 0, elements = elements, minElements = Filters$minElements, 
                              maxElements = Filters$maxElements)
  if (is.null(mm)) {
    return(failReturn)
  }
  f1 <- if (!is.null(Filters$DBErange)) {
    mm$DBE >= Filters$DBErange[1] & mm$DBE <= Filters$DBErange[2]
  }
  else {
    rep(TRUE, length(mm$DBE))
  }
  if (!is.null(Filters$parity) && Filters$parity %in% c("e", 
                                                        "o")) {
    f1 <- f1 & mm$parity == Filters$parity
  }
  if (!any(f1)) {
    return(failReturn)
  }
  f1 <- which(f1)
  if (length(f1) == 0) {
    return(failReturn)
  }
  f2 <- f1[order(abs(mm$exactmass[f1] - z * 0.00054857990907 - 
                       mz))]
  sfs <- makeMF(mm$formula[f2], forcelist = TRUE)
  res <- data.frame(mz = mm$exactmass[f2] - z * 0.00054857990907, 
                    MF = mm$formula[f2], charge = z, RdisopScore = mm$score[f2], 
                    unsat = mm$DBE[f2], parity = mm$parity[f2], error = mm$exactmass[f2] - 
                      z * 0.00054857990907 - mz, nrule = mm$valid[f2], 
                    stringsAsFactors = FALSE)
  res$ppm <- res$error/mz * 1e+06
  
  #Golden rule #1: Restriction for element numbers
  if (maxCounts) {
    Da <- abs(mz/z)
    if (Da < 500) {
      maxCountLimit <- consolidateMF(c(C = 29, H = 72, 
                                       N = 10, O = 18, P = 4, S = 7, F = 15, Cl = 8, 
                                       Br = 5))
    }
    else if (Da < 1000) {
      maxCountLimit <- consolidateMF(c(C = 66, H = 126, 
                                       N = 25, O = 27, P = 6, S = 8, F = 16, Cl = 11, 
                                       Br = 8))
    }
    else if (Da < 2000) {
      maxCountLimit <- consolidateMF(c(C = 115, H = 236, 
                                       N = 32, O = 63, P = 6, S = 8, F = 16, Cl = 11, 
                                       Br = 8))
    }
    else if (Da < 3000) {
      maxCountLimit <- consolidateMF(c(C = 162, H = 208, 
                                       N = 48, O = 78, P = 6, S = 9, F = 16, Cl = 11, 
                                       Br = 8))
    }
    else {
      maxCountLimit <- NULL
    }
    if (!is.null(maxCountLimit)) {
      maxCountLimit[maxCountLimit == 0] <- Inf
      res$maxCounts <- sapply(sfs, "<=", maxCountLimit)
    }
    else {
      res$maxCounts <- TRUE
    }
    if (!is.null(Filters$maxCounts) && Filters$maxCounts) {
      sel <- res$maxCounts
      if (!any(sel)) {
        return(failReturn)
      }
      res <- res[sel, ]
      sfs <- sfs[sel]
    }
  }
  
  #Golden rule #2: (Senior's third theorem): 
  if(SeniorRule){
    
    res$SENIOR3 <- sapply(sfs,getSenior3.MFobject2)
    
    if(!is.null(Filters$SENIOR3) && is.numeric(Filters$SENIOR3)){
      
      sel <- res$SENIOR3 >= Filters$SENIOR3
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  #Golden Rule #4: Hydrogen/Carbon element ratio check
  if(HCratio){
    
    res$HtoC <-  sapply(sfs,function(MF){
      
      MF["H"]/MF["C"]
      
    })
    
    if(!is.null(Filters$HCratio) && Filters$HCratio){
      
      sel <- is.na(res$HtoC) | (res$HtoC > 0.2 & res$HtoC < 3.1)
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
  }
  
  
  #Golden Rule #5: Heteroatom check
  if(moreRatios){
    res$NtoC <-  sapply(sfs,function(MF){
      MF["N"]/MF["C"]
    })
    res$OtoC <-  sapply(sfs,function(MF){
      MF["O"]/MF["C"]
    })
    res$PtoC <-  sapply(sfs,function(MF){
      MF["P"]/MF["C"]
    })
    res$StoC <-  sapply(sfs,function(MF){
      MF["S"]/MF["C"]
    })
    res$OStoP <-  sapply(sfs,function(MF){
      MF["P"]/(MF["O"]+ MF["S"])
    })
    res$FtoC <-  sapply(sfs,function(MF){
      MF["F"]/MF["C"]
    })
    res$CltoC <-  sapply(sfs,function(MF){
      MF["Cl"]/MF["C"]
    })
    res$BrtoC <-  sapply(sfs,function(MF){
      MF["Br"]/MF["C"]
    })
    
    res$HFClBrItoC <-  sapply(sfs,function(MF){
      (MF["H"]+MF["F"]+MF["Cl"]+MF["Br"]+MF["I"])/MF["C"]  #Typically, halogens replace H so this ratio is useful (ie: for PFAS)
    })
    
    if(!is.null(Filters$moreRatios) && Filters$moreRatios){
      
      sel <- ((is.na(res$NtoC) | res$NtoC < 1.3)
              & (is.na(res$OtoC) | res$OtoC < 1.2)
              & (is.na(res$PtoC) | res$PtoC < 0.3)
              & (is.na(res$StoC) | res$StoC < 0.8))
      
      if(!any(sel)){return(failReturn)}
      
      res <- res[sel,]
      sfs <- sfs[sel]
      
    }
    
  }
  
  
  # Golden rule #6: element probability check
  if (elementHeuristic) {
    res$elementHeuristic = sapply(sfs, function(MF) {
      r1 <- na.omit(MF[c("N", "O", "P", "S")])
      r2 <- r1[r1 > 1]
      if (length(r2) <= 2) {return(TRUE)}
      if (length(r2) == 4 && (r2["N"] >= 10 || r2["O"] >= 20 || r2["P"] >= 4 || r2["S"] >= 3)) { return(FALSE)}
      switch(paste(names(r2), collapse = ""), 
             NOP = {if (min(r2) > 3 && (r2["N"] >= 11 || r2["O"] >= 22 || r2["P"] >= 7)) {return(FALSE)}}, 
             OPS = {if (min(r2) > 1 && (r2["S"] >= 3 || r2["O"] >= 14 || r2["P"] >= 3)) {return(FALSE)}},
             NPS = {if (min(r2) > 1 && (r2["N"] >= 4 || r2["S"] >= 3 || r2["P"] >= 3)) {return(FALSE)}}, 
             NOS = {if (min(r2) > 6 && (r2["S"] >= 8 || r2["O"] >= 14 || r2["N"] >= 19)) {return(FALSE)}})
      return(TRUE)
    })
    
    if (!is.null(Filters$elementHeuristic) && Filters$elementHeuristic) {
      sel <- res$elementHeuristic
      if (!any(sel)) {
        return(failReturn)
      }
      res <- res[sel, ]
      sfs <- sfs[sel]
    }
  }
  if (!is.null(top) && is.numeric(top)) {
    res <- res[seq(min(top, nrow(res))), ]
  }
  if (summarize) {
    return(paste(paste0(res$MF, if (res$charge[1] > 0) {
      paste0("(+", res$charge, ")")
    } else if (res$charge[1] < 0) {
      paste0("(", res$charge, ")")
    } else {
      ""
    }), collapse = "|"))
  }
  return(res)
}



detectBrCl <- function(h.unk, Q1ref, z=1){
  
  
  Q1_dist <- round((h.unk[,'mz']-Q1ref)*z)  # Setup M diff
  
  idx_mzPos2 <- which(Q1_dist==2)
  if(length(idx_mzPos2)>0){
    
    idx_mzPos0 <- min(which(h.unk$Int_N>40)) # Check later, may not work for unrefined spectra
    
    idx_mzPos2 <- idx_mzPos2[which(h.unk$Int_N[idx_mzPos2] == max(h.unk$Int_N[idx_mzPos2]))] #Get highest mzPos2
    
    n_Cl <- round(h.unk$Int_N[idx_mzPos2]*0.9/(h.unk$Int_N[idx_mzPos0]*0.3199))
    n_Br <- round(h.unk$Int_N[idx_mzPos2]*0.9/(h.unk$Int_N[idx_mzPos0]*0.9727))
    
    # Check for sulfur if Cl and Br absent
    if(n_Cl<1 && n_Br<1){
      n_S <- round(h.unk$Int_N[idx_mzPos2]/(h.unk$Int_N[idx_mzPos0]*0.04250))  
    }else{
      Int_N2_adj <- h.unk$Int_N[idx_mzPos2] - 100*0.8*(n_Cl*0.3199 - n_Br*0.9727)
      n_S <- round(Int_N2_adj/(h.unk$Int_N[idx_mzPos0]*0.04250))
    }
    
    Hal <- data.frame(Cl=min(n_Cl,7), Br=min(n_Br,5), S=min(n_S,3))  #Incorporate Max of special elements
    
  }else{Hal <- data.frame(Cl=0, Br=0, S=0)}
  
  
  #MFList <- MFList[which(MetaboCoreUtils::containsElements(MFList$MF, 'Br')==F),]   ## Remove MFs with Br
  #MFList <- MFList[which(MetaboCoreUtils::containsElements(MFList$MF, 'Cl')==F),]   ## Remove MFs with Cl
  
  if(h.unk[,'mz'][1] < 200){  #Only detect Boron and Iron when mz <200
    # Detect Boron
    idx_mzPosB <- mzIPidx_Max(spec=h.unk, Q1_dist, mzIP=-1, minInt=12.4) 
    if(length(idx_mzPosB)>0){
      n_B <- round(h.unk$Int_N[idx_mzPosB]/24.84)
      Hal <- cbind(Hal, B=n_B)
    }else{Hal <- cbind(Hal, B=0)}
    
    # Detect Iron
    idx_mzPosFe <- mzIPidx_Max(spec=h.unk, Q1_dist, mzIP=-2, minInt=3.15) 
    if(length(idx_mzPosFe)>0){
      n_Fe <- round(h.unk$Int_N[idx_mzPosFe]/6.37)
      if(length(which(Q1_dist==-1))!=0){n_Fe=0}
    }else{n_Fe=0}
    #Double check
  }else{
    Hal <- cbind(Hal, B=0)
    n_Fe=0
  }
  
  Hal <- cbind(Hal, Fe=n_Fe)
  
  return(Hal)
}

## Collect Molecular Formula with adducts and deducts
ionizeMF <- function(m, a=NA, d=NA){
  
  M <- MetaboCoreUtils::standardizeFormula(m)
  if(is.na(a)==F){
    a <- MetaboCoreUtils::standardizeFormula(a)
    M <- MetaboCoreUtils::addElements(M, a)}
  if(is.na(d)==F){
    d <- MetaboCoreUtils::standardizeFormula(d)
    if(MetaboCoreUtils::containsElements(M, d)){
      M <- MetaboCoreUtils::subtractElements(M, d)
    }else{
      print(paste('Cannot deduct',d,'from',M))
      M <- NA
    }
  }
  if(is.na(M)==F){
    M <- MassTools::makeMF(M)
    M <- MetaboCoreUtils::pasteElements(M[which(M!=0)])
    #M <- print(M)
  }else{M <- NA}
  #M <- countElements(M)
  #M <- pasteElements(M[[1]][-which(M[[1]]==0)])  #Drop 0 count elements
  
  return(list(MFinput=m, Mion=M))
} 


makeMF2 <- function(s, forcelist = F){
  
  s <- MassTools::reformatMFstring(s)
  elnums <- strsplit(s, "[A-Z]|[A-Z][a-z]")
  elnums <- lapply(elnums, function(el) {
    as.integer(el[-1])
  })
  elnames <- strsplit(s, "-[0-9]+|[0-9]+")
  elnums <- mapply(function(enu, ena) {
    names(enu) <- ena
    return(enu)
  }, enu = elnums, ena = elnames, SIMPLIFY = F)
  consolidated <- lapply(elnums, consolidateMF2)
  if (!forcelist && length(consolidated) == 1) {
    return(consolidated[[1]])
  }
  else {
    return(consolidated)
  }
}


consolidateMF2 <- function (x) {
  target <- getOption("MassTools.elements")
  target[unique(isotopes$element)] <- rep(0,length(unique(isotopes$element)))  ## VERY SLOW!
  
  #E <- c('Ga','Te','Ge','Sn')
  #target[E] <- rep(0,length(E))
  matched <- match(names(x), names(target))
  matched<-matched[!is.na(matched)]
  if (!any(duplicated(names(x)))&(sum(is.na(matched))<1)) {
    target[matched] <- x
  }
  else {
    for (i in unique(matched)) {
      target[i] <- sum(x[matched == i])
    }
  }
  class(target) <- "MFobject"
  return(target)
}


## Use EnviPat to generate centroided MS
MFtoSpectrum <- function(formulas, q=F, Res=F, ppmTol=F, patOnly=F){
  
  ## Formula modified to accept ppmTol or Resolving power @ m/z
  
  # data(isotopes) should already be loaded
  #formulas <- check_chemform(isotopes,formulas) crashes R
  
  pattern<-enviPat::isopattern(
    isotopes,
    formulas,
    threshold=1,
    charge=q,
    emass=0.00054858,
    plotit=F,
    algo=1,
    rel_to=1,
    verbose=F,
    return_iso_calc_amount=F)
  
  if(patOnly==T){
    np = length(pattern)
    if(np>1){
      centro <- list()
      for(p in 1:np){
        centro[[p]] <- data.frame(pattern[[p]][,1:2])
      }
      names(centro) <- names(pattern)
    }else{
      centro <- data.frame(pattern[[1]][,1:2])
    }
    
  }else{
    
    if(Res==F){
      
      profiles<-enviPat::envelope(
        pattern,
        ppm=T,
        dmz=ppmTol*200/1e6,
        env="Gaussian",
        verbose=F,
        plotit=F)
    }
    
    if(ppmTol==F){
      print(pattern)
      profiles<-enviPat::envelope(
        pattern,
        ppm=F,
        dmz='get',
        frac=1/4,
        env="Gaussian",
        resolution=Res,
        verbose=F,
        plotit=F)
    }
    
    
    centro<-vdetect(profiles,detect="centroid", plotit=F, verbose=F)
    
    if(length(centro)==1){
      centro <- data.frame(centro)
      colnames(centro) <- c('m/z','abundance')
    }
  }# end if PatOnly
  
  return(centro)
  
}


## --------------------------------------------------------------------------------
## Customized version of OrgMassSpecR::SpectrumSimilarity to add intensity precision

SpectrumSimilarity_custom <- function (spec.top, spec.bottom, dppm = 10, b = 10, top.label = NULL, 
                                       bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, int.prec=0,
                                       print.alignment = FALSE, print.graphic = TRUE, output.list = FALSE) 
{
  
  t = spec.top[1,1]*dppm/1e6  #build alignment tolerance from ppm tolerance (in amu)
  
  top_tmp <- data.frame(mz = spec.top[, 1], intensity = spec.top[,2])
  top_tmp$normalized <- round((top_tmp$intensity/max(top_tmp$intensity)) * 100, int.prec)
  top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= xlim[2])
  top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)
  top <- subset(top_plot, top_plot$intensity >= b)
  bottom_tmp <- data.frame(mz = spec.bottom[, 1], intensity = spec.bottom[,2])
  bottom_tmp$normalized <- round((bottom_tmp$intensity/max(bottom_tmp$intensity)) * 100, int.prec)
  bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & bottom_tmp$mz <= xlim[2])
  bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)
  bottom <- subset(bottom_plot, bottom_plot$intensity >= b)
  for (i in 1:nrow(bottom)) top[, 1][bottom[, 1][i] >= top[,1] - t & bottom[, 1][i] <= top[,1] + t] <- bottom[,1][i]
  alignment <- merge(top, bottom, by = 1, all = TRUE)
  if (length(unique(alignment[, 1])) != length(alignment[, 1])) 
    warning("the m/z tolerance is set too high")
  alignment[, c(2, 3)][is.na(alignment[, c(2, 3)])] <- 0
  names(alignment) <- c("mz", "int.ref", "int.unk")
  if (print.alignment == TRUE) {
    print(alignment)
  }
  if (x.threshold < 0) 
    stop("x.threshold argument must be zero or a positive number")
  alignment <- alignment[alignment[, 1] >= x.threshold, ]
  u <- alignment[, 2]
  v <- alignment[, 3]
  if(length(which(v>2)) > 1 ){
    similarity_score1 <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))
    
  }else{similarity_score1 <- as.vector((u %*% v)/(sqrt(sum(u^2)) * sqrt(sum(v^2))))*length(which(v>2))/length(which(u>2)) }  #major penalty for only one iso
  if (print.graphic == TRUE) {
    plot.new()
    plot.window(xlim = xlim, ylim = c(-125, 125))
    ticks <- c(-100, -50, 0, 50, 100)
    for (i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i],2), c(0, top_plot$intensity[i]), col = "blue", lwd=2, lend=1)
    for (i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i],2), c(0, -bottom_plot$intensity[i]), col = "red", lwd=2,lend=1)
    axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "intensity")
    axis(1, pos = -125)
    lines(xlim, c(0, 0), lwd=3)
    rect(xlim[1], -125, xlim[2], 125)
    mtext("m/z", side = 1, line = 2)
    mtext("intensity (%)", side = 2, line = 2)
    plot.window(xlim = c(0, 20), ylim = c(-10, 10))
    text(10, 9, top.label)
    text(10, -9, paste0(bottom.label,' MFscore=',round(similarity_score1,4)))
    
  }
  if (output.list == TRUE && print.graphic == TRUE) {
    headTailPlot <- function() {
      grid::pushViewport(grid::plotViewport(c(5, 5, 2, 2)))
      grid::pushViewport(grid::dataViewport(xscale = xlim, yscale = c(-125, 
                                                                      125)))
      grid::grid.rect()
      tmp <- pretty(xlim)
      xlabels <- tmp[tmp >= xlim[1] & tmp <= xlim[2]]
      grid::grid.xaxis(at = xlabels)
      grid::grid.yaxis(at = c(-100, -50, 0, 50, 100))
      grid::grid.segments(top_plot$mz, top_plot$intensity, top_plot$mz, 
                          rep(0, length(top_plot$intensity)), default.units = "native", 
                          gp = grid::gpar(lwd = 0.75, col = "blue"))
      grid::grid.segments(bottom_plot$mz, -bottom_plot$intensity, 
                          bottom_plot$mz, rep(0, length(bottom_plot$intensity)), 
                          default.units = "native", gp = grid::gpar(lwd = 0.75, 
                                                                    col = "red"))
      grid::grid.abline(intercept = 0, slope = 0)
      grid::grid.text("intensity (%)", x = grid::unit(-3.5, "lines"), 
                      rot = 90)
      grid::grid.text("m/z", y = grid::unit(-3.5, "lines"))
      grid::popViewport(1)
      grid::pushViewport(grid::dataViewport(xscale = c(0, 20), yscale = c(-10, 
                                                                          10)))
      grid::grid.text(top.label, grid::unit(10, "native"), grid::unit(9, 
                                                                      "native"))
      grid::grid.text(bottom.label, grid::unit(10, "native"), grid::unit(-9, 
                                                                         "native"))
      grid::popViewport(2)
    }
    p <- grid::grid.grabExpr(headTailPlot())
  }
  if (output.list == F && print.graphic==F) {return(similarity_score1)}
  if (output.list == T && print.graphic==T) {return(list(similarity.score = similarity_score1, alignment = alignment,plot = p))}
  if (output.list == T && print.graphic==F) {return(list(similarity.score = similarity_score1, alignment = alignment))}
}


mzIPidx_Max <- function(spec, Q1_dist, mzIP, minInt=0){  #Gets max Int pos for given mzIP
  
  idx_mzIP <- which(Q1_dist==mzIP & spec$Int_N >minInt)
  if(length(idx_mzIP)>0){
    idx_mzIP <- idx_mzIP[which(spec$Int_N[idx_mzIP] == max(spec$Int_N[idx_mzIP]))]
  }else{idx_mzIP <- NULL}
  
  return(idx_mzIP)
}


SpectrumSimilarity_custom2 <- function (obs = NULL, the = NULL, dppm = 2, int_prec = 0.3, 
                                        limit = 0, rnd_prec = 5) 
{
  
  # Normalize spectra, Intensity is column 2
  obs[,2] <- obs[,2]/max(obs[,2])
  the[,2] <- the[,2]/max(the[,2])
  
  # Align spectra on m/z. Choose 1 m/z from Obs that is within the ppmTol.  
  # Assumes observed spectrum may have superfluous peaks
  n_m <- nrow(the)
  b_obs <- rep(0,n_m)
  spec_align <- data.frame(m_ref=the[,1],i_ref=the[,2], m_obs=b_obs, i_obs=b_obs, m_ppm=b_obs)
  for(m in 1:n_m){
    dist_m <- abs(the[m,1] - obs[,1])*1e6/the[m,1]
    idx_m_min <- which(dist_m==min(dist_m) )
    if(length(idx_m_min)==1){
      spec_align[m,3] <- obs[idx_m_min,1]
      spec_align[m,4] <- obs[idx_m_min,2]
      spec_align[m,5] <- dist_m[idx_m_min]
    }
  }
  
  ## Mass residuals relative to mass error tolerance
  max_err_mz <- dppm * spec_align[,1]/10^6
  dmz <- 1 - abs(spec_align[,3] - spec_align[,1])/max_err_mz
  dmz[dmz > 1] <- 1
  
  ## Intensity residuals relative to Intensity tolerance
  max_err_int <- int_prec * spec_align[,2]
  dint <- 1 -  (abs(spec_align[,4] - spec_align[,2]))/max_err_int
  dint[dint > 1] <- 1
  
  ## Dot-product calculations
  u.m <- spec_align[,1]/max_err_mz
  v.m <- spec_align[,3]/max_err_mz
  u.i <- spec_align[,2]/max_err_int
  v.i <- spec_align[,4]/max_err_int
  
  out.m <- 1 - (1 - as.vector((u.m %*% v.m)/(sqrt(sum(u.m^2)) * sqrt(sum(v.m^2)))))*1e6
  out.i <- as.vector((u.i %*% v.i)/(sqrt(sum(u.i^2)) * sqrt(sum(v.i^2))))
  
  score1 <- weighted.mean(c(out.m,out.i),w=c(1,0)) #equal weighted
  score2 <- weighted.mean(c(out.m,out.i),w=c(0,1)) #int only weighted
  score3 <- ( max(0,weighted.mean(dmz,spec_align[,2])) + max(0,weighted.mean(dint,spec_align[,2]))  )/2
  #old method
  #score3<- round(1.01 - mean(sqrt(dmz * dint)), rnd_prec)
  
  
  return(list(Score1=round(out.i,rnd_prec),Score2=round(score3,rnd_prec), Alignment=spec_align))
}

MFnRuleValidate <- function(MFList, Type=1){
  pcMFs_N <- unlist(lapply(MFList$MF, nElemInMF, E='N')) # Get # of N for each MF
  MFList$nrule <- 'Invalid' #Default invalid
  if(Type==1){idx_NValids <- which(pcMFs_N%%2 == round(MFList$mz-MFList$charge)%%2)}
  else{idx_NValids <- which(pcMFs_N%%2 == round(MFList$mz)%%2)}
  MFList$nrule[idx_NValids] <- 'Valid'  
  MFList$nrule[which(pcMFs_N == 0)] <- 'Valid' #No N needs not the n-rule
  #pcMFs
  
  return(MFList)
}

## Count Element in MF
nElemInMF <- function(MF, E='H'){
  nE <- as.data.frame(MetaboCoreUtils::countElements(MF))[E,]
  if(is.na(nE)){nE <-0}
  return(nE)
}

## Compare MFs
compareMFs <- function(MFs){
  mod <- RMassBank::add.formula(MFs[1], RMassBank::multiply.formula(MFs[2], -1))
  if(mod==''){m <- T}else{m <- F}
  return(m)
}

