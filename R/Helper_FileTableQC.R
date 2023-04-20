#
#
#
#
#
#
#------ Establish function -------#
#' fileTableQC
#'
#' Check if the loaded output files (e.g. from MaxQuant), include all the columns required for the downstream assembly and include measured values (not NA or 0 for all features). This is a Quality Control (QC) check of the input files before proceeding to SmartPhos::readExperiment.
#'
#' @param fileTable Loaded fileTable using SmartPhos::generateInputTable
#'  @export
#'  
fileTableQC <- function(fileTable){
  fileTable_FP <- fileTable[fileTable$type == "proteome",]
  fileTable_PP <- fileTable[fileTable$type == "phosphoproteome",]
  
  #--- Full Proteome (FP) ----- 
  list_FP <- unique(fileTable_FP$fileName)
  
  for(i in list_FP){
    inputTab_FP <- data.table::fread(i, sep = "\t", check.names = TRUE)#Load data
    list_sampleName <- unique(fileTable_FP$sample)
    
    mylist <- list() #create an empty list
    
    for(k in list_sampleName){
      #Check if the required columns are in the DF
      if (((paste0("LFQ.intensity.", k) %in% colnames(inputTab_FP))==TRUE) & ((paste0("Intensity.", k) %in% colnames(inputTab_FP))==TRUE) & ((paste0("Razor...unique.peptides.", k) %in% colnames(inputTab_FP))==TRUE)){
        RequiredColumn <- TRUE
        
        # Check if the required columns have values, 0s or NA
        LFQ <- inputTab_FP%>%
          select(paste0("LFQ.intensity.", k))
        Intensity <- inputTab_FP%>%
          select(paste0("Intensity.", k))
        UniquePeptides <- inputTab_FP%>%
          select(paste0("Razor...unique.peptides.", k))
        
        if((all(UniquePeptides[,1] == 0))==TRUE |(all(Intensity[,1] == 0))==TRUE |(all(LFQ[,1] == 0))==TRUE | all(is.na(UniquePeptides[,1])) | all(is.na(Intensity[,1])) | all(is.na(LFQ[,1]))){
          AllValuesMissing <- TRUE
          #Add results to mylist:
          mylist[[k]] <- c(RequiredColumn, AllValuesMissing)
        } else{
          AllValuesMissing <- FALSE
          #Add results to mylist:
          mylist[[k]] <- c(RequiredColumn, AllValuesMissing)
        }
        
      }else{
        RequiredColumn <- FALSE
        mylist[[k]] <- RequiredColumn
      }
    }
    
    # Combine all vectors into a matrix
    RequiredColumns <- as.data.frame(do.call("rbind",mylist))%>% 
      rename("Required Columns Present?"=1,
             "All features NA/0?"=2)%>%
      rownames_to_column("sample")
    # Merge with input:
    fileTable_FP <- merge(fileTable_FP, RequiredColumns, by="sample", all.x=TRUE)
  }
  
  #---  Phosphoproteome (PP)--- 
  list_PP <- unique(fileTable_PP$fileName)
  
  for(i in list_PP){
    inputTab_PP <- data.table::fread(i, sep = "\t", check.names = TRUE)#Load data
    list_sampleName <- unique(fileTable_PP$sample)
    
    mylist <- list() #create an empty list
    
    for(k in list_sampleName){
      #Check if the required columns are in the DF
      if (((paste0("Localization.prob.", k) %in% colnames(inputTab_PP))==TRUE) & ((paste0("Score.diff.", k) %in% colnames(inputTab_PP))==TRUE) & ((paste0("Intensity.", k) %in% colnames(inputTab_PP))==TRUE)){
        RequiredColumn <- TRUE
        
        # Check if the required columns have values, 0s or NA
        Localization <- inputTab_PP%>%
          select(paste0("Localization.prob.", k))
        ScoreDiff <- inputTab_PP%>%
          select(paste0("Score.diff.", k))
        Intensity <- inputTab_PP%>%
          select(paste0("Intensity.", k))
        
        if((all(Intensity[,1] == 0))==TRUE |(all(ScoreDiff[,1] == 0))==TRUE |(all(Localization[,1] == 0))==TRUE | all(is.na(Intensity[,1])) | all(is.na(ScoreDiff[,1])) | all(is.na(Localization[,1]))){
          AllValuesMissing <- TRUE
          #Add results to mylist:
          mylist[[k]] <- c(RequiredColumn, AllValuesMissing)
        } else{
          AllValuesMissing <- FALSE
          #Add results to mylist:
          mylist[[k]] <- c(RequiredColumn, AllValuesMissing)
        }
        
      }else{
        RequiredColumn <- FALSE
        mylist[[k]] <- RequiredColumn
      }
    }
    
    # Combine all vectors into a matrix
    RequiredColumns <- as.data.frame(do.call("rbind",mylist))%>% 
      rename("Required Columns Present?"=1,
             "All features NA/0?"=2)%>%
      rownames_to_column("sample")
    # Merge with input:
    fileTable_PP <- merge(fileTable_PP, RequiredColumns, by="sample", all.x=TRUE)
  }
  
  #--- Messgae & Return ----- 
  #Filter:
  fileTable_FP_filt <-fileTable_FP[,c(2,1,3:7)]%>%
    subset(!`Required Columns Present?`==FALSE)%>%#Remove samples where not all required columns are available
    subset(!`All features NA/0?`==TRUE)#Remove samples in at least one of the required for all features either 0 or NA was detected
  
  fileTable_PP_filt<- fileTable_PP[,c(2,1,3:7)]%>%
    subset(!`Required Columns Present?`==FALSE)%>%#Remove samples where not all required columns are available
    subset(!`All features NA/0?`==TRUE)#Remove samples in at least one of the required for all features either 0 or NA was detected
  
  #Message:
  if((nrow(fileTable_FP)==nrow(fileTable_FP_filt))==FALSE){
    warning("fileTable$type == proteome included ", nrow(fileTable_FP), ". After filtering due to required columns missing and/or all features being NA/0 only ", nrow(fileTable_FP_filt), " samples are kept. Please check that this is due to the settings used in your analysis software that led to samples enriched for phospho sites not to be quantified for the full protein file."  )
  }else{
    message("All samples in fileTable$type == proteome have the required columns.")
  }
  
  if((nrow(fileTable_PP)==nrow(fileTable_PP_filt))==FALSE){
    warning("fileTable$type == phosphoproteome included ", nrow(fileTable_PP), ". After filtering due to required columns missing and/or all features being NA/0 only ", nrow(fileTable_PP_filt), " samples are kept. Please check that this is due to the settings used in your analysis software that led to samples not enriched for phospho sites not to be quantified for the phosphoproteome file."  )
  }else{
    message("All samples in fileTable$type == phosphoproteome have the required columns.")
  }
  
  #Return
  if((nrow(fileTable_FP)==nrow(fileTable_FP_filt))==FALSE | (nrow(fileTable_PP)==nrow(fileTable_PP_filt))==FALSE){
    message("fileTable as filtered and unfiltered version is returned into your environment. Please check your warnings!")
    fileTable_NotFilt <- rbind(fileTable_PP, fileTable_FP)
    fileTable_NotFilt <- fileTable_NotFilt[,c(2,1,3:7)]%>%
      mutate(`Removed?` = case_when((`Required Columns Present?` =="FALSE") | (`All features NA/0?`==TRUE)~ 'TRUE',
                                    TRUE ~ 'FALSE'))
    assign("fileTable_FiltParameters", fileTable_NotFilt, envir=.GlobalEnv)
    
    
    fileTable_Filt <-rbind(fileTable_PP_filt, fileTable_FP_filt)
    fileTable_Filt <-fileTable_Filt[,1:5]
    assign("fileTable_Filtered", fileTable_Filt, envir=.GlobalEnv)
    
  }else{
    message("fileTable has not been changed.")
  }
  
}
