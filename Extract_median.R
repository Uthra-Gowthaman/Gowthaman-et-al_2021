#install and load the flowcore package from www.bioconductor.org
library(flowCore)

#Load the plate data and specify the location
plate<-"xxx"
fcs_files_folder<-"xxx"

# Listing of the fcs files,
temp = list.files(paste("00_FCS-files",fcs_files_folder,sep="/"),pattern="*fcs",full.names=F)

# Reordering of the wells (A1 as first and P24 as last)
g<-c()
for(i in 1:length(temp)){
  g_i<-grep(sort(gsub(".*._","",temp))[i],temp)	
  g<-c(g,g_i)
}
temp<-temp[g]

# Read the fcs files from the plate
fcs<-lapply(paste("00_FCS-files",fcs_files_folder,temp,sep="/"),read.FCS)

# Creating a new list 'fcs_unfiltered'
fcs_unfiltered<-fcs

#Defining a minimum cut-off cell number to consider data in the well
cell_number<-500
less500<-c()

# Looping the fcs files for filtering
for(i in 1:length(fcs)){
  
  if(nrow(fcs[i][[1]]@exprs)>cell_number){ 
    less500_i<-"F"
    
    #quantiles
    q25_fsc<-quantile(fcs[i][[1]]@exprs[,"FSC-H"] , 1/4)
    q75_fsc<-quantile(fcs[i][[1]]@exprs[,"FSC-H"] , 3/4)
    q25_ssc<-quantile(fcs[i][[1]]@exprs[,"SSC-H"] , 1/4)
    q75_ssc<-quantile(fcs[i][[1]]@exprs[,"SSC-H"] , 3/4)
    
    if(length(which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc)) >1){
      
      fcs[i][[1]]@exprs<-fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),]  
      
    } else { fcs[i][[1]]@exprs<-rbind(fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),] ,fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),] ) }
    
  } else {less500_i<-"T"}
  less500<-c(less500,less500_i)
  
  if(i %in% seq(0,100,10)) { print(i) }
}

# Extracting the plate min, max and median values mCherry and YFP
yall<-c()
xall<-c()
for(i in 1:length(fcs)){
  if(nrow(fcs[i][[1]]@exprs)>cell_number){ 
    yall_i<-log10(fcs[i][[1]]@exprs[,"B535/345-H"])
    xall_i<-log10(fcs[i][[1]]@exprs[,"mCherry-H"])
    yall<-c(yall,yall_i)
    xall<-c(xall,xall_i)
  }
}
xmedian<-median(xall)
ymedian<-median(yall)
xmin<-min(xall)
ymin<-min(yall)
xmax<-max(xall)
ymax<-max(yall)


#Plots of the entire plate
filename<-paste("01_plate view/",plate,".png",sep="")
png(filename,width=12,height=8,res=200,units="in")

#vector
leg<-paste(rep(c("A","B","C","D","E","F","G","H"),each=12),seq(1,12,1),sep="")

#Defining the plot values
par(mfrow=c(8,12),mar=c(1.5,1.5,0.5,0.5),bty="n")

xmin=1
xmax=4
ymin=1
ymax=4
xlim<-c(floor(xmin),ceiling(xmax))
ylim<-c(floor(ymin),ceiling(ymax))

makeTransparent<-function(someColor, alpha=3)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1],green=curcoldata[2],blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
color1<-makeTransparent("black")
col<-rep(color1,96)

# Looping of the fcs files to plot log10 signal median
for(i in 1:length(fcs)){
  
  if(less500[i] == "F"){
    #removal of negative values to avoid "-Inf" after log transformation
    fcs[i][[1]]@exprs<-fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"mCherry-H"]>0 & fcs[i][[1]]@exprs[,"B535/345-H"]>0),]
    x<-log10(fcs[i][[1]]@exprs[,"mCherry-H"])
    y<-log10(fcs[i][[1]]@exprs[,"B535/345-H"])
    
    plot(x,y,cex=0.4,xaxt="n",yaxt="n",xlab="",ylab="",pch=19,xlim=xlim,ylim=ylim,col=col,cex.axis=0.3)
    legend("topleft",legend=leg[i],bty="n",cex=0.5,inset=c(-0.18,-0.1),text.col="grey")
    
    #median line
    mm<-median(x)
    ym<-median(y)
    abline(h=ym,col=rgb(0.15,0.15,0.15,maxColorValue=1),lwd=0.5,lty=3)
    abline(v=mm,col=rgb(0.15,0.15,0.15,maxColorValue=1),lwd=0.5,lty=3)
    
    #axis line
    axis(side=1,at=seq(floor(xmin),ceiling(xmax),1),cex.axis=0.4,tcl=-0.1,labels=rep("",length(seq(floor(xmin),ceiling(xmax),1))),col="grey",col.ticks="grey")
    axis(side=1,at=seq(floor(xmin),ceiling(xmax),1),cex.axis=0.4,line=-1.2,lwd=0,col.axis="grey")
    axis(side=2,at=seq(floor(ymin),ceiling(ymax),1),cex.axis=0.4,tcl=-0.1,labels=rep("",length(seq(floor(ymin),ceiling(ymax),1))),col="grey",col.ticks="grey")
    axis(side=2,at=seq(floor(ymin),ceiling(ymax),1),cex.axis=0.4,line=-0.9,lwd=0,col.axis="grey")
    mtext(side=1,line=0.4,"log10(mCherry-H)",cex=0.3,col="grey")
    mtext(side=2,line=0.6,"log10(YFP-H)",cex=0.3,col="grey")
    #rect(xleft=log10(q40_fsc),xright=log10(q60_fsc),ybottom=log10(q40_ssc),ytop=log10(q60_ssc),border="red2",lwd=2.5)
    
  }
    else{  
    plot(1,1,cex=0.4,xaxt="n",yaxt="n",xlab="",ylab="",pch=19,xlim=xlim,ylim=ylim,col="white",cex.axis=0.3)
    legend("topleft",legend=leg[i],bty="n",cex=0.5)
    text(x=(xlim[2]+xlim[1])/2,y=(ylim[2]+ylim[1])/2,"not enough cells",cex=0.25,col="grey")
  }
    if(i == 1){  
    date<-fcs[i][[1]]@description["$DATE"][[1]]
    time<-fcs[i][[1]]@description["$BTIM"][[1]]
    legend("bottomleft",legend=c(date,time),bty="n",cex=0.5,text.col="grey")
  }  
    legend("bottomright",legend=nrow(fcs[i][[1]]),bty="n",cex=0.5,text.col="grey")
  
  if(i %in% seq(0,100,10)) { print(i) }
}

# restoring the 'fcs unfiltered' files for extracting values
fcs<-fcs_unfiltered

ncell_before<-unlist(lapply(fcs,nrow))
ncell_before

file_name<-paste("02_cell numbers/cell_numbers_",plate,".pdf",sep="")
pdf(file_name)
plot(ncell_before)
abline(h=median(ncell_before),lwd=2.5,col="grey")
text(x=48,y=median(ncell_before),round(median(ncell_before)),font=2,cex=1.5)
dev.off()
sink(paste("02_cell numbers/cell_numbers_",plate,".txt",sep=""))
sum(ncell_before)
sink()

cell_number<-500
less500<-c()
for(i in 1:length(fcs)){
  if(nrow(fcs[i][[1]]@exprs)>cell_number){ 
    less500_i<-"F"
    q25_fsc<-quantile(fcs[i][[1]]@exprs[,"FSC-H"] , 1/4)
    q75_fsc<-quantile(fcs[i][[1]]@exprs[,"FSC-H"] , 3/4)
    q25_ssc<-quantile(fcs[i][[1]]@exprs[,"SSC-H"] , 1/4)
    q75_ssc<-quantile(fcs[i][[1]]@exprs[,"SSC-H"] , 3/4)
    
    if(length(which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc)) >1){
      
      fcs[i][[1]]@exprs<-fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),]  
      
    } else { fcs[i][[1]]@exprs<-rbind(fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),] ,fcs[i][[1]]@exprs[which(fcs[i][[1]]@exprs[,"SSC-H"]>q25_ssc & fcs[i][[1]]@exprs[,"SSC-H"]<q75_ssc & fcs[i][[1]]@exprs[,"FSC-H"]>q25_fsc & fcs[i][[1]]@exprs[,"FSC-H"]<q75_fsc),] ) }
    
  } else {less500_i<-"T"}
  less500<-c(less500,less500_i)
  
  if(i %in% seq(0,100,10)) { print(i) }
}

#Extraction of cell numbers after filtering
ncell_after_filtering<-unlist(lapply(fcs,nrow))
write.csv(cbind(ncell_before,ncell_after_filtering),paste("02_cell numbers/cell_numbers_every_well_",plate,".csv",sep=""))

#Defining cell number cut-off after filtering
cell_number2<-100
mCherry_medians<-c()
YFP_medians<-c()
for(i in 1:length(fcs)){
  if(nrow(fcs[i][[1]]@exprs)>cell_number2 & less500[i]=="F"){  
    
    mCherry_medians_i<-log10(median(fcs[i][[1]]@exprs[,"mCherry-H"]))
    mCherry_medians<-c(mCherry_medians,mCherry_medians_i)
    YFP_medians_i<-log10(median(fcs[i][[1]]@exprs[,"B535/345-H"]))
    YFP_medians<-c(YFP_medians,YFP_medians_i)
  }
  else{
    mCherry_medians_i<-NaN
    mCherry_medians<-c(mCherry_medians,mCherry_medians_i)
    YFP_medians_i<-NaN
    YFP_medians<-c(YFP_medians,YFP_medians_i)
  }
}

#Create data.frame with medians
medians<-data.frame(mCherry_medians,YFP_medians)
names(medians)<-c("mCherry","YFP")

#Adding columns
well<-paste(rep(c("A","B","C","D","E","F","G","H"),each=12),1:12)
medians$well<-well[1:nrow(medians)]
medians$cell_number_before<-ncell_before
medians$cell_number_after<-ncell_after_filtering
head(medians)

file_name<-paste("03_medians csv-files/medians_",plate,".csv",sep="")
write.csv(medians,file_name,row.names=F)

dev.off()

