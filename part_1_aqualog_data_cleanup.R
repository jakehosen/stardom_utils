options(scipen=999)


choose.dir<- function() {
	system("osascript -e 'tell app \"RStudio\" to POSIX path of (choose folder with prompt \"Choose Folder:\")' > /tmp/R_folder",
			intern = FALSE, ignore.stderr = TRUE)
	p <- system("cat /tmp/R_folder && rm -f /tmp/R_folder", intern = TRUE)
	return(ifelse(length(p), p, NA))
}


convert_aqualog<-function(){
library(tcltk)
library(utils)

#choose folder where corrected data will be stored.
cat("\nchoose destination folder")
#destination_folder<-choose.dir()
destination_folder<-"/Users/jdh/Library/CloudStorage/GoogleDrive-jakehosen@gmail.com/My Drive/spectroscopy/Sepulveda/corrected_data"

x<-"yes"
while(x=="yes"|x=="y"){
#choose folder where original aqualog data are stored.
cat("\nchoose folder where original aqualog data are stored.")
#main_folder<-choose.dir()
main_folder <- dirname(file.choose())
#main_folder<-"/Volumes/GoogleDrive/My Drive/USDA Wetlands/Aqualog_Data/Raw Data/"

cat("choose name for project")
#proj_name<-readLines(n=1)
proj_name<-gsub('.*/ ?(\\w+)', '\\1', main_folder)

abs_files<-list.files(main_folder,pattern="*Abs Spectra Graphs.dat")
#abs_file<-read.table("/Users/prime/gdrive/uf_aqualog/uf_aqualog/20190408/20190408_Group1/abs/B1S15S15ABS.csv")

#gsub(".*/(.*?)/$", "\\1", main_folder)
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),"/abs",sep=""))
#dir.create(paste(destination_folder,gsub(".*/(.*?)/$", "\\1", main_folder),"/eem",sep=""))
dir.create(paste(destination_folder,"/abs",sep=""))
dir.create(paste(destination_folder,"/eem",sep=""))
#dir.create(paste(destination_folder,"/abs/",gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
#dir.create(paste(destination_folder,"/eem/",gsub(".*/(.*?)/$", "\\1", main_folder),sep=""))
#abs_cor_folder<-paste(destination_folder,"/abs/",gsub(".*/(.*?)/$", "\\1", main_folder),sep="")
#eem_cor_folder<-paste(destination_folder,"/eem/",gsub(".*/(.*?)/$", "\\1", main_folder),sep="")

dir.create(paste(destination_folder,"/abs/",proj_name,sep=""))
dir.create(paste(destination_folder,"/eem/",proj_name,sep=""))
abs_cor_folder<-paste(destination_folder,"/abs/",proj_name,sep="")
eem_cor_folder<-paste(destination_folder,"/eem/",proj_name,sep="")



for(i in 1:length(abs_files)){
	abs_temp0<-read.table(paste(main_folder,"/",abs_files[i],sep=""),sep="\t")[1:191,]
	abs_temp<-abs_temp0[-c(1:3),c("V1","V10")]
	write.table(abs_temp,paste(abs_cor_folder,"/",gsub(" - Abs Spectra Graphs.dat",".csv",abs_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}


sem_files<-list.files(main_folder,pattern="*Blank Waterfall Plot.dat")

for(i in 1:length(sem_files)){
	abs_temp<-read.table(paste(main_folder,"/",sem_files[i],sep=""),sep="\t")
#	abs_col<-ncol(abs_temp)
	abs_temp2<-abs_temp[-c(2),]
#	names(abs_temp2)<-abs_temp[1,]
	abs_temp2[1,1]<-"Wavelength"
	write.table(abs_temp2,paste(eem_cor_folder,"/",gsub(" - Sample - Blank Waterfall Plot.dat",".csv",sem_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")
}


bem_files<-list.files(main_folder,pattern="*Waterfall Plot Blank.dat")
bem_temp<-read.table(paste(main_folder,"/",bem_files[1],sep=""),sep="\t")
bem_temp2<-bem_temp[-c(1:2),]
names(bem_temp2)<-bem_temp[1,]
bem_temp2[1,1]<-"Wavelength"
write.table(bem_temp2,paste(eem_cor_folder,"/",gsub(" - Sample - Blank Waterfall Plot.dat","_blank.csv",bem_files[i],fixed=TRUE),sep=""),row.names=FALSE,col.names=FALSE,sep=",")

print("\nCorrect Another Folder?")
x<-readLines(n=1)
print(x)
}
}


convert_aqualog()

#==========================================
#this section of code is for identifying and excluding samples in stardom that have already been reviewed
#slsrfl$Aqualog_Date<-as.Date(as.character(slsrfl$Date.run.on.aqualog),format="%Y%m%d")
#already_reviewed<-subset(slsrfl,Aqualog_Date!=as.Date("2021-10-03"))$Aqualog_Project_ID

#exclude <- list("ex" = c(),
#                "em" = c(),
#                "sample" = already_reviewed
#)
#eem_list_new_only<-eem_eclude(eem_list,exclude)
#eem_overview_plot(eem_list_new_only)
#end code subset
#==========================================
