#We are going to store our R packages for eem processing in a separate folder to avoid errors. Make a directory in a place where it will be secure (people sometimes use documents or their home account folder). Enter that directory location below. If you are using macos or linux and you want to keep the directory in your home folder, don't change anything.
#.libPaths("C:/Users/primers/Documents/R/eems_library")

#this will install staRdom. You only need to run this once when you are first installing packages. Make sure the directory indicated by lib is the same as the libPaths directory above
#install.packages("staRdom",lib="C:/Users/primers/Documents/R/eems_library")

library("staRdom")
cores <- detectCores(logical = FALSE)
cores<-6

#cl <- snow::makeCluster(2, setup_strategy = "sequential")
#clusterEvalQ(cl, .libPaths("C:/Users/primers/Documents/R/eems_library"))
#cl <- parallel::makeCluster(6, setup_timeout = 0.5)


eem_read_csv <- function(path, col = "ex", recursive = TRUE, is_blank_corrected = FALSE, is_scatter_corrected = FALSE, is_ife_corrected = FALSE, is_raman_normalized = FALSE, manufacturer = "unknown", ...){
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


#Load in Dataset
eem_list_dil_k0<-readRDS("/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/fl_dataset_20240609.rds")

eem_list_dil_k <- eem_extract(eem_list_dil_k0, c("nano", "miliq", "milliq", "mq", "blank","B1S1S1SEM","B1S6FawleySEM","B1S7SwitcesSEM","B1S8AdamsSEM","B1S9McBeeSEM","25X","15X","30X","50X","35X","120X","75X"),ignore_case = TRUE)



###here I am saving the corrected interpolated files ready for parafac analysis###

#I have commented the line that saves the file (as a RDS file which is a format for R). You can uncomment this if you correct additional scans and want to save a new dataset. Otherwise, you can just use the readRDS line to load the data and do PARAFAC analysis.
#saveRDS(eem_list_dil, paste(gsub("/data/","",data_folder,fixed=TRUE),"/Corrected_EEMs_20202702.rds",sep=""))

#this line loads the corrected dataset that I saved and is ready for PARAFAC analysis.
#eem_list_dil<-readRDS(paste(gsub("/data/","",data_folder,fixed=TRUE),"/Corrected_EEMs_20202702.rds",sep=""))

# minimum and maximum of numbers of components
dim_min <- 5
dim_max <- 5

nstart <- 20 # number of similar models from which best is chosen
maxit = 8000 # maximum number of iterations in PARAFAC analysis
ctol <- 10^-5 # tolerance in PARAFAC analysis

###  CHANGE FOR SEDIMENT PARAFAC ###
##if you have bad samples you need to exclude, you can do that here.##
#exclude1<-list(sample=c("B1S15S15","B1S1S1","B1S2S2","B1S3S3","B1S4S4","B1S5S5"))
#eem_list_dil2<-eem_exclude(eem_list_dil, exclude1)
#sample_remove<-c("name1","name2")
#exclude<-list("ex"=c(), "em"=c(), "sample"=sample_remove

#Running the initial model with 3 through 6 components. When we're done with the first model we will view the resulting components and choose based on visual assessment how many components we expect to fit in the final model.
#TRUE to scale to one (if have large range ie. soils + water)
#FALSE does not scale


eem_list_dil_k <- eem_extract(eem_list_dil_k, c(354:360), ignore_case = TRUE) #remove unwanted samples



# Run initial model manually, and verify:
pf1n <- eem_parafac(eem_list_dil_k, comps = seq(dim_min,dim_max), normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)

#eempf_plot_comps(pf1n, contour = TRUE, type = 1)


#Troubleshooting function:
paramod<-pf1n[[1]] #select one model you want it to run. based on dim_min/dim_max above. If min=3, pf1n[[3]]=5 component model
dataset<-eem_list_dil_k
normit<-TRUE #T/F if data should be normalized

parafac_autofit<-function(paramod,dataset,normit){
  cores <- detectCores(logical = FALSE)
  nstart <- 20 # number of similar models from which best is chosen
  maxit = 5000 # maximum number of iterations in PARAFAC analysis
  ctol <- 10^-5 # tolerance in PARAFAC analysis
  comps<-length(names(as.data.frame(paramod$A)))
  complete<-FALSE
  cpl <- eempf_leverage(paramod)	
  while(complete!=TRUE){
    
    scaled_leverages<-scale(log10(cpl$A),center=TRUE,scale=TRUE)
    auto_exclude<-list(sample=row.names(scaled_leverages)[scaled_leverages>1.96],em=c(),ex=c())
    print(auto_exclude)
    if(length(auto_exclude[[1]])>0){
      dataset <- eem_exclude(dataset, auto_exclude)
      splithalf <- splithalf(dataset, comps=comps, normalise = normit, rand = FALSE, cores = cores, nstart = nstart, maxit = maxit, ctol = ctol)
      splithalf_tcc <- splithalf_tcc(splithalf)
      print(splithalf_tcc)
      if(sum(splithalf_tcc[,3:4]<0.95)==0){complete<-TRUE}else{paramod <- eem_parafac(dataset, comps = comps, normalise = normit, const = c("nonneg", "nonneg", "nonneg"), maxit = maxit, nstart = nstart, ctol = ctol, cores = cores)
      cpl <- eempf_leverage(paramod[[1]])}}else{
        print("no model fit.")
        complete<-TRUE 
        paramod<-NULL}
  }
  results<-list(model=paramod,data=dataset)
  return(results)
}


#eempf_compare(validated_7comp[[1]], contour = TRUE)
validated_7comp<-parafac_autofit(paramod=paramod,dataset=dataset,normit=TRUE)
#validated_5comp<-parafac_autofit(paramod=pf1n[[3]],dataset=dataset,normit=TRUE)
#validated_6comp<-parafac_autofit(paramod=pf1n[[4]],dataset=dataset,normit=TRUE)
#summary(validated_5comp)



#now we've validated a 5 component model with sediment and water column samples. We need to run a high tolerance model and then output the results from that model

#saveRDS(validated_4comp, "~/...")
final_dataset <- validated_7comp[[2]] #[[2]]= second list element, all of the data including outliers that were previously excluded

#validated the 3 component model based on eem_list_dil_ex1
#now running at higher tolerance
ctol <- 10^-9 # tolerance in PARAFAC analysis
pf1n_hightole <- eem_parafac(final_dataset, comps = 5, normalise = TRUE, const = c("nonneg", "nonneg", "nonneg"), maxit = 11000, nstart = nstart, ctol = ctol, cores = cores)
eempf_compare(pf1n_hightole)

#we're adding in the samples that were excluded during the outlier process and exporting all the sample scores.
#pf3n_hightole<-readRDS("GIW_DOC_20190529.rds")
pf2n_wOutliers <- A_missing(eem_list_dil_k, pfmodel = pf1n_hightole[[1]], cores = cores)
pf2n_wOutliers

dataout<-as.data.frame(pf2n_wOutliers$A)
dataout$sample<-row.names(dataout)

write.csv(dataout,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_scores_20240611.csv",row.names=FALSE)

saveRDS(pf2n_wOutliers,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_scores_20240611_fullobj.rds")

saveRDS(pf1n_hightole,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_fullmodel.rds")



pf1n_hightole<-readRDS("/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_fullmodel.rds")


para_scores<-as.data.frame(pf1n_hightole[[1]]$A)
para_scores$sample<-row.names(para_scores)
write.csv(para_scores,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_scores_20240611.csv",row.names=FALSE)


#save PARAFAC for Openfluor - DO NOT OPEN WITH EXCEL after original formatting
eempf_openfluor(pf1n_hightole[[1]], file = "/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_openfluor.txt")


#save PARAFAC scores
write.csv(pf2n_wOutliers$A,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/USDA Wetlands/Data/Instrument_Data/Aqualog/Aqualog_Data/Corrected Data/PARAFAC Scores Water and Soil 4 Component USDA Wetlands 20240125.csv")

#save PARAFAC for Openfluor - DO NOT OPEN WITH EXCEL after original formatting
eempf_openfluor(pf1n_hightole[[1]], file = "/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/USDA Wetlands/Data/Instrument_Data/Aqualog/Aqualog_Data/Corrected Data/PARAFAC Scores Water and Soil 4 Component USDA Wetlands 20240125.txt")

saveRDS(pf1n_hightole,"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/USDA Wetlands/Data/Instrument_Data/Aqualog/Aqualog_Data/Corrected Data/PARAFAC Scores Water and Soil 4 Component USDA Wetlands 20240125.rds") # After save, change decimal notation so there is no "e" -> instead 9 decimal places


mod<-readRDS("/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/EcosystemScience/USDA Wetlands/Stats/usda_wetlands_jake/data/all_in_parafac_20240611/all_in_parafac_scores_20240611_fullobj.rds")

mod %>%
  ggeem(contour = TRUE)

#eempf_report(pf2n_wOutliers$A, export = "SuperEEM_Report_3Comp.txt", eem_list = eem_list_dil, splithalf = TRUE, shmodel = sh_ex1_c3, performance = TRUE, spp=1, cores = cores )

#eempf_report(pf2n_hightole[[1]], export = "parafac_report.html", eem_list = eem_list, shmodel = sh, performance = TRUE)

#Peak identification/Association Notes:hea
#--------------------------------------#
# peaks within 10-15nm of another study are acceptable
# parentheses around the secondary peak, primary (largest) peak no parentheses


