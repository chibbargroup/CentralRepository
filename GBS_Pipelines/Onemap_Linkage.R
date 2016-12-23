#!/usr/bin/env Rscript
#Multicore Linkage Analysis script using R package onemap
#Onemap documention is located at ./Documents/Tutorial_Onemap_complete_version.pdf
#Written May 13, 2015
#Craig Irvine (cri646@mail.usask.ca)
#version 0.5

#Changelog
#28May15 : Ver 0.2 - Moved save results into groupProcessing so results are saved as they are completed.
#		   - Added safe and forced output files for order.seq method
#03June15 : Ver 0.3 - Added command line options to enable testing of LOD values before mapping groups
#		    - Added Options --lod, --maxrf and --test to set LOD value and run group size tests.
#02July15 : Ver 0.4 - Added options --oc and --f2 to allow for setting of input file from the command line
#		    - Added Input File to sink commands to allow for multiple runs in the same directory
#10Nov15 : Ver 0.5  - Added option for saving the forced or save version of Allgroup_maps.mapmaker
#                   - Fixed wrong input varable for compare method of groupProcessing.
#                   - Added LG.f2.all.safe assignment in compare method of groupProcessing
#                   - Fixed wrong placement of closing bracket in order.seq method of groupProcessing
#                   - Added writting of Allgroup_maps.maps and Allgroup_maps.png after all processing is completed
#                       for use with R/QTL.
#4Jan16 : Ver 0.6   - Fixed problem with ord.seq method usage, THRES options should have matched the LOD setting.
#                   - Added LOD and max.rf options to rf.2pts usage.
#                   - Fixed writting of Allgroup mapping and png files.
#7Jan16 : Ver 0.7   - Added global lists to hold both safe and forced maps for saving at the end of the process.

sink(file="Onemap_linkage_log.txt",split=TRUE)

cat('Loading Libraries')
library(onemap)
library(qtl)
library(parallel)


cat('Clearing Workspace\n')
rm(list=ls()) #clear workspace

#defaut options
opt_lod <- 3.0
opt_maxrf <- 0.5
opt_testrun <- FALSE
opt_safe <- FALSE
opt_f2_cross_input <- ""
opt_outcross_input <- ""
log_name <- ""


LG.f2.all.vec <- list()
LG.f2.vec <- list()
G.LG.f2.ord.all <- list()
G.LG.f2.ord.safe <- list()




#Group processing function is used by mclapply to fork the marker ordering
#processes and speed up completion time.  Uses global var LG.f2.vec and 
#passed index values to perform operations.
groupProcessing <- function(inputIndexs)
{

  sink(file=paste0(log_name,"_MG_",inputIndexs,"_log.txt"))

  cat(paste0("**Marker Group: ",inputIndexs," - groupProcessing Started\n"))

  if (as.integer(length(LG.f2.vec[[inputIndexs]][[1]])) < 10) #number of markers in a group.
  {
    cat(paste0("Marker Group: ",inputIndexs," - Marker Number less than 10, compare method\n"))
    comp <- compare(LG.f2.vec[[inputIndexs]])
    
    if (as.integer(length(LG.f2.vec[[inputIndexs]][[1]])) < 3)
    {
      LG.f2.all <- comp
      LG.f2.all.safe <- comp

    } else
    {
      LG.f2.all <- make.seq(comp,1,1)
      LG.f2.all.safe <- make.seq(comp,1,1)
    }
    
    write.map(LG.f2.all, paste0("MarkerGroup_",inputIndexs,".map"))
    png(filename=paste0("MarkerGroup_",inputIndexs,".png"),height=960,width=960,pointsize=14)
    draw.map(LG.f2.all,names=TRUE,grid=TRUE,cex.mrk=1.4)
    dev.off()
    
  } else
  { 
    cat(paste0("Marker Group: ",inputIndexs," - Marker Number greater than 10, order.seq method\n"))
    LG.f2.ord <- order.seq(input.seq=LG.f2.vec[[inputIndexs]], n.init = 5, subset.search = "twopt", twopt.alg = "rcd", THRES = opt_lod, draw.try = FALSE, wait = 0) #draw.try and wait should not be changed due to the multiprocessing nature of the function.  If require set mc.core in 
		     #mclapply equal to 1
    LG.f2.all <- make.seq(LG.f2.ord,"force") # force or safe
    LG.f2.all.safe <- make.seq(LG.f2.ord,"safe") # force or safe
    
    
    #ripple.seq(LG.f2.all, ws=5, LOD=3)#checks for a better marker ordering
  

    write.map(LG.f2.all, paste0("MarkerGroup_",inputIndexs,"_forced.map"))
    png(filename=paste0("MarkerGroup_",inputIndexs,"_forced.png"),height=960,width=960,pointsize=14)
    draw.map(LG.f2.all,names=TRUE,grid=TRUE,cex.mrk=1.4)
    dev.off()
    
    write.map(LG.f2.all.safe, paste0("MarkerGroup_",inputIndexs,"_safe.map"))
    png(filename=paste0("MarkerGroup_",inputIndexs,"_safe.png"),height=960,width=960,pointsize=14)
    draw.map(LG.f2.all.safe,names=TRUE,grid=TRUE,cex.mrk=1.4)
    dev.off()
  }
  
  cat(paste0("**Marker Group: ",inputIndexs," - groupProcessing Complete\n"))
  
  sink()
  
  
  G.LG.f2.ord.all <- list(LG.f2.ord.all,LG.f2.all)
  G.LG.f2.ord.safe <- list(LG.f2.ord.safe,LG.f2.all.safe)
  
  
  if (opt_safe == TRUE)
    return(LG.f2.all.safe)
  else
    return(LG.f2.all)
}


#Script starts here.
cat('Options\n')
args <- commandArgs(trailingOnly=TRUE)
print(args)

set_lod <- FALSE
set_maxrf <- FALSE
set_f2 <- FALSE
set_outcross <- FALSE

for (opts in args)
{
  
  #Set options
  if (set_lod == TRUE)
  {
    opt_lod <- as.numeric(opts)
    set_lod <- FALSE    
  } else if (set_maxrf == TRUE)
  {
    opt_maxrf <- as.numeric(opts)
    set_maxrf <- FALSE
  } else if (set_outcross == TRUE)
  {
    opt_outcross_input <- opts
    log_name <- opts
    set_outcross <- FALSE
    
  } else if (set_f2 == TRUE)
  {
    opt_f2_cross_input <- opts
    log_name <- opts
    set_f2 <- FALSE
  

  
  #check options
  } else if (opts == '--lod')
  {
    set_lod <- TRUE
  } else if (opts == '--maxrf')
  {
    set_maxrf <- TRUE
  } else if (opts == '--test')
  {
    opt_testrun <- TRUE
  }else if (opts == '--f2')
  {
    set_f2 <- TRUE
  }else if (opts == '--oc')
  {
    set_outcross <- TRUE
  }
  else if (opts == '--safe')
  {
    opt_safe <- TRUE
  }
  
  

}




cat('Reading Data File\n')
Input_data <- ""

if (opt_f2_cross_input != "")
{
  Input_data <- read.mapmaker(file=opt_f2_cross_input) #f2 crosses
}else if (opt_outcross_input != "")
{
  Input_data <- read.outcross(file=opt_outcross_input) #Outcross
}else
{
  cat('No input file specified\n')
  exit(1)
}


cat('Two Point Analysis\n')
twopts <- rf.2pts(Input_data, LOD=opt_lod, max.rf=opt_maxrf)

cat('Marker Grouping\n')
mark.all.f2 <- make.seq(twopts,"all")
(LGs.f2 <- group(mark.all.f2, LOD=opt_lod, max.rf=opt_maxrf)) #display it for the record


if (opt_testrun == FALSE)
{
  type="kosambi" # can be haldane or kosambi
  cat(paste0("Mapping Function: ",type,"\n"))
  set.map.fun(type)

  NumGroup <- LGs.f2["n.groups"] #Accesses number of groups in LGs.f2 class
  cat(paste0("Number of groups: ",NumGroup,"\n"))
  cat(paste0("Vectoring Groups\n"))
  for (i in seq(1,as.integer(NumGroup))) 
  { #create the vector for multiprocessing.
    (LG.f2.vec[[i]] <- make.seq(LGs.f2,i));
  }

  vecIndexs <- 1:as.integer(NumGroup)
  #use multiprocessing to process all groups
  cat(paste0("Starting groupProcessing()\n"))
  (LG.f2.all.vec = mclapply(vecIndexs,groupProcessing,mc.preschedule=FALSE,mc.cores=8))
  
  
  G.LG.f2.ord.all
  
  write.map(G.LG.f2.ord.all,"AllGroups_all.map")
  png(filename="AllGroups_all.png",height=960,width=960,pointsize=14)
  draw.map(LG.f2.all.vec,names=TRUE,grid=TRUE,cex.mrk=1.4)
  dev.off()
  
  
  write.map(G.LG.f2.ord.safe,"AllGroups_safe.map")
  png(filename="AllGroups_safe.png",height=960,width=960,pointsize=14)
  draw.map(LG.f2.all.vec,names=TRUE,grid=TRUE,cex.mrk=1.4)
  dev.off()
  
} else
{
  cat('Testrun Finished')
}

cat("Analysis Complete\n")

sink()



