#Perform blank correction, filter corrections, etc.
rm(list = ls ()) # clears R memory

library("staRdom")
library(tidyverse)
cores <- detectCores(logical = FALSE)
setwd("~")

#Run these functions every time
###############################################################################
removed_fun<-function (eem, sample, keep = FALSE, ignore_case = FALSE, verbose = TRUE) 
{
  stopifnot(class(eem) == "eemlist", is.character(sample) | 
              is.numeric(sample))
  sample_names <- unlist(lapply(eem, function(x) {
    x$sample
  }))
  if (is.numeric(sample)) {
    stopifnot(all(is_between(sample, 1, length(eem))))
    to_remove <- ifelse(rep(keep, length(sample)), setdiff(1:length(eem), 
                                                           sample), sample)
    eem[to_remove] <- NULL
    if (verbose) {
      cat(ifelse(keep, "Extracted sample(s):", "Removed sample(s):"), 
          sample_names[sample], "\n")
    }
  }
  if (is.character(sample)) {
    to_remove <- grepl(paste(sample, collapse = "|"), sample_names, 
                       ignore.case = ignore_case)
    eem[xor(to_remove, keep)] <- NULL
    if (verbose) {
      if (all(to_remove == FALSE)) {
        cat("Nothing to remove.")
      }
      else {
        cat(ifelse(keep, "Extracted sample(s):", "Removed sample(s):"), 
            sample_names[to_remove], "\n")
      }
    }
  }
  return(sample_names[to_remove])
}


#cl <- snow::makeCluster(2, setup_strategy = "sequential")
#clusterEvalQ(cl, .libPaths("C:/Users/primers/Documents/R/eems_library"))
#cl <- parallel::makeCluster(6, setup_timeout = 0.5)


eem_read_csv2 <- function(path, col = "ex", recursive = TRUE, is_blank_corrected = FALSE, is_scatter_corrected = FALSE, is_ife_corrected = FALSE, is_raman_normalized = FALSE, manufacturer = "unknown", ...){
  #path <- "./inst/extdata/EEMs"
  if(col == "em"){
    csv_func <- eem_csv2
  } else {
    csv_func <- eem_csv
  }
  eems <- eem_read(path, recursive = TRUE, import_function = csv_func)
  
  lapply(eems, function(eem){
    class(eem) <- "eem"
    attr(eem, "is_blank_corrected") <- is_blank_corrected
    attr(eem, "is_scatter_corrected") <- is_scatter_corrected
    attr(eem, "is_ife_corrected") <- is_ife_corrected
    attr(eem, "is_raman_normalized") <- is_raman_normalized
    attr(eem, "manufacturer") <- manufacturer
    eem
  }) %>%
    `class<-`("eemlist")
}



#check this tutorial for reference: https://cran.r-project.org/web/packages/staRdom/vignettes/PARAFAC_analysis_of_EEM.html







#data_folder<-choose.dir()
data_folder<-"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/spectroscopy/Sepulveda/corrected_data"
#data_folder<-"/Volumes/GoogleDrive/My Drive/USDA Wetlands/Aqualog_Data/Corrected Data"
#folder <- system.file("extdata/EEMs/", package = "staRdom") # folder containing example EEMs
eem_folder<-paste(data_folder,"/eem",sep="")
eem_list <- eem_read_csv2(eem_folder) # in case you use your own data, just replace folder by a path. e.g. "C:/folder/another folder" and use eem_read if you do not have plain csv tables to import



absorbance_path<-paste(data_folder,"/abs",sep="")
absorbance0 <- absorbance_read(absorbance_path,recursive=TRUE) # load csv or txt tables in folder
absorbance <- abs_blcor(absorbance0,wlrange = c(680,700))

#metatable <- system.file("extdata/metatable_dreem.csv",package = "staRdom") # path to example data, can be replaced by a path to your own data
#meta <- read.table(paste(data_folder,"/metatable.csv",sep=""), header = TRUE, sep = ",", dec = ".", row.names = 1) # load data
#dilution = "meta" # e.g. 1 for undiluted samples
#dil_sample_name_correction = FALSEc

#eem_listtemplate(eem_list, absorbance) %>%
#  write.csv(file="metatable.csv", row.names = FALSE)

abst<-as.data.frame(t(absorbance))
names(abst)<-abst[1,]
abst<-abst[2:nrow(abst),]



# adjust range of EEMs to cover correction vectors
#eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1]))
#eem_list <- eem_range(eem_list,ex = range(Excor[,1]), em = range(Emcor[,1][1],780))
eem_lista <- eem_range(eem_list,ex = c(250,600), em = range(Emcor[,1][1],600))

#eem_list2b <- eem_spectral_cor(eem_lista,Excor,Emcor)

eem_list2c <- eem_extend2largest(eem_lista, interpolation = 1, extend = FALSE, cores = cores)

eem_list3 <- eem_ife_correction(eem_lista,absorbance, cuvl = 1)

# blank subtraction
#eem_list4 <- eem_remove_blank(eem_list3)
#eem_overview_plot(eem_list2, spp=8)


#eem_list4b <- eem_raman_normalisation(eem_list3)
eem_list4b <- eem_raman_normalisation2(eem_list3, blank = "blank")
#eem_list <- eem_raman_normalisation(eem_list)


#SL<-read.csv("/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/USDA Wetlands/Stats/usda_wetlands_jake/data/basic_fl_stats_20230316_env_water copy.csv")



#samp_names<-eem_humification_index(eem_list4b)
samp_names<-eem_humification_index(eem_list4b)
to_rm<-samp_names$sample[!(samp_names$sample %in% SL$sample)]

eem_list5 <- eem_extract(eem_list3, to_rm, ignore_case = TRUE) #remove unwanted samples
abs_rm <- removed_fun(eem_list3, to_rm, ignore_case =TRUE)
absorbance2 <- dplyr::select(absorbance, -dplyr::matches(paste(abs_rm, collapse="|"), ignore.case = T))

list5 <- sapply(eem_list5,"[[","sample")
abs <- names(absorbance2)
abs[!(abs %in% list5)] #should just display "wavelength", other samples listed should be excluded using 'absadd' above


                            
#set scatter removal dimensions
remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
#remove_scatter_width <- c(30,50,85,50) #(nm above and below, ...)
remove_scatter_width <- c(50,35,45,35)

eem_list_rem <- eem_rem_scat(eem_list3, remove_scatter = remove_scatter, remove_scatter_width = remove_scatter_width)

#look at data
#eem_overview_plot(eem_list5[197:202], spp=6) # #before scatter is removed
#eem_overview_plot(eem_list_rem, spp=6) # band


#only plot samples/rows 50 to 55. plotting subsamples of the complete list
#eem_overview_plot(eem_list_rem[50:55], spp=6)
#eem_exlude()
#eem_overview_plot(eem_list5, spp=6)


#==========================================
#this section of code is for identifying and excluding samples in stardom that have already been reviewed
#slsrfl$Aqualog_Date<-as.Date(as.character(slsrfl$Date.run.on.aqualog),format="%Y%m%d")
#already_reviewed<-subset(slsrfl,Aqualog_Date!=as.Date("2021-10-03"))$Aqualog_Project_ID

#exclude <- list("ex" = c(),
#                "em" = c(),
#               "sample" = already_reviewed
#)
#eem_list_new_only<-eem_eclude(eem_list,exclude)
#eem_overview_plot(eem_list_new_only)
#end code subset
#==========================================

#linear interpolation to fill in the gaps
eem_list_rem_int <- eem_interp(eem_list_rem, cores = cores, type = 1, extend = FALSE)
##use if getting noise in models
#eem_overview_plot(eem_list_rem_int[171], spp=6)
#stop 171, go through rest 04/18

#eem_overview_plot(eem_list_rem_int[197:202], spp=6)
#eem_overview_plot(eem_list_rem_int[12:13], spp=6)

#282 flagged
#watch for 198, rnadom high value pixel
##This is what checking? they were blue

#dil_data <- meta["dilution"]

#eem_list_dil <- eem_dilution(eem_list_rem_int,dil_data)
eem_list_dil<-eem_list_rem_int #eems list, blanks removed, interpolated - no dilutions


##save eem list for easy re-loading: ======================================== ##
saveRDS(eem_list_dil, paste(gsub("/data/","",data_folder,fixed=TRUE),"/eem_list_dil_sepulveda.rds",sep=""))
#reload: 
#eem_list_dil <- read_rds("~/MRitter/Modeling/StaRdom/Corrected_Data/eem_list_dil_20221207_LW.rds")
##=========================================================================== ##

library(dplyr)
library(reshape2)
library(viridis)
comp_scans<-data.frame()
comp_scans[1,]<-NA
for(i in 1:length(eem_list_dil)){
eem_list_temp<-as.data.frame(eem_list_dil[[i]]$x)
row.names(eem_list_temp)<-(eem_list_dil[[i]]$em)
names(eem_list_temp)<-(eem_list_dil[[i]]$ex)
#comp_scans$em<-(eem_list_dil[[i]]$em)
keepers<-as.numeric(names(eem_list_temp))[as.numeric(names(eem_list_temp))>=240 & as.numeric(names(eem_list_temp))<=340]
eem_list_temp_plot<-eem_list_temp[,as.character(keepers)]
eem_list_temp_plot$rowSums<-rowSums(eem_list_temp_plot)
eem_list_temp_plot$rowSums
comp_scans<-bind_cols(comp_scans,data.frame(eem_list_temp_plot$rowSums))
names(comp_scans)[i]<-eem_list_dil[[i]]$sample
}
comp_scans$em<-(eem_list_dil[[i]]$em)

comp_scans_melt<-melt(comp_scans,id.vars=c("em"))
comp_scans_melt<-comp_scans_melt[!grepl("blank",comp_scans_melt$variable,fixed=TRUE),]
comp_scans_melt<-comp_scans_melt[!grepl("PFOS0.25uM (01)",comp_scans_melt$variable,fixed=TRUE),]
comp_scans_melt<-comp_scans_melt[!grepl("PFOS2.5uM (01)",comp_scans_melt$variable,fixed=TRUE),]
comp_scans_melt<-comp_scans_melt[!grepl("PFOS0.25uM (02)",comp_scans_melt$variable,fixed=TRUE),]

ggplot(comp_scans_melt,aes(em,value,color=variable))+
geom_point()+
scale_colour_viridis_d(option="A")


abst<-as.data.frame(t(absorbance2)) #extract abs stats #***changed from absorbance to absorbance2***#
names(abst)<-abst[1,]
abst<-abst[2:nrow(abst),]

abs_stats<-data.frame()
for(i in 1:nrow(abst)){
  data_subset<-as.numeric(abst[i,])
  abs_base<-data.frame(V1=as.numeric(names(abst)),V3=data_subset)
  
  abs_base$am1<-abs_base$V3*100*2.3025851
  abs_base$lnam1<-log(abs_base$am1)
  
  abs_275_295<-subset(abs_base,V1>=275 & V1<=295)
  if(sum(!is.na(abs_275_295$lnam1))>5){
    s275295<-lm(abs_275_295$lnam1~abs_275_295$V1)$coefficients[2]*-1
  }else{s275295<-NA}
  
  abs_350_400<-subset(abs_base,V1>=350 & V1<=400)
  abs_350_400$lnam1[is.na(abs_350_400$lnam1)]<-0
  if(sum(!is.na(abs_350_400$lnam1))>2){
    s350400<-lm(abs_350_400$lnam1~abs_350_400$V1)$coefficients[2]*-1
  }else{s350400<-NA}
  
  sratio<-s275295/s350400
  
  a254<-approx(abs_base[,1],abs_base$am1,xout=254)$y[1]
  a254_dec<-a254/2.3025851
  
  abs_280_450<-subset(abs_base,V1>=280 & V1<=450)
  if(sum(!is.na(abs_280_450$lnam1))>2){
    s280450<-lm(abs_280_450$lnam1~abs_280_450$V1)$coefficients[2]*-1
  }else{s280450<-NA}
  
  
  totala250450<-sum(approx(abs_base$V1,abs_base$V3, xout=c(250:450))[[2]])
  
  abs_stats_temp<-data.frame(sample=row.names(abst)[i],s275295=s275295,s350400=s350400,sratio=sratio,s280450=s280450,a254=a254,a254_dec=a254_dec,totala250450=totala250450)
  row.names(abs_stats_temp)<-"R"
  
  abs_stats<-bind_rows(abs_stats,abs_stats_temp)
}


#this generates the table of summary statistics and outputs.
fl_data<-cbind.data.frame(eem_biological_index(eem_list_dil),eem_coble_peaks(eem_list_dil),eem_fluorescence_index(eem_list_dil),eem_humification_index(eem_list_dil))
fl_data2<-fl_data[,c(1,2,4,5,6,7,8,10,12)]
fl_abs_data<-merge(fl_data2,abs_stats,by="sample")


#saves corrected fl and abs stats to file:

write.csv(fl_abs_data,paste(gsub("/data/","",data_folder,fixed=TRUE),"/basic_fl_stats_20240125_W.csv",sep=""),row.names=FALSE)

# saveRDS(fl_abs_data,paste(gsub("/data/","",data_folder,fixed=TRUE),"/basic_fl_stats_20221207_LW.RDS",sep=""))
FL_corrections <- summary(eem_list_dil)
# write.csv(FL_corrections,paste(gsub("/data/","",data_folder,fixed=TRUE),"/eems_list_20221207_LW.csv",sep=""),row.names=FALSE)


#------------------------------> Up to this point 10/11/2022 mr - all data (lake, wetland, soil, other): only shortest list of samples removed
#------------------------------> Up to this point 10/19/2022 mr - data (lake, wetland): list of samples from AQ scan log (Ritter_rm) removed

###here I am saving the corrected interpolated files ready for parafac analysis###

#I have commented the line that saves the file (as a RDS file which is a format for R). You can uncomment this if you correct additional scans and want to save a new dataset. Otherwise, you can just use the readRDS line to load the data and do PARAFAC analysis.

#saveRDS(eem_list_dil, paste(gsub("/data/","",data_folder,fixed=TRUE),"/Corrected_EEMs_20221019_LW.RDS",sep=""))



#saveRDS(eem_list_dil,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/USDA Wetlands/Stats/usda_wetlands_jake/data/usda_eems_w_cor.rds")


