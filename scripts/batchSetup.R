##########################
#parameters

exptNameIn = "mapSetA_exptMinExample"
nameAddBit = "Test" #identifier for this run

##############################
#Load packages and scripts
#Note: for parallel processing need to use package installed from github rather than local dev version.
#Ensure any recent changes to the package have been committed before proceeding
if(0){# reinstall package if necessary.
  Sys.setenv(TZ = "UTC")
  #deps <- tools::package_dependencies("LSTDConnect", recursive = TRUE)
  #update.packages(oldPkgs = unlist(deps), ask = FALSE)
  devtools::install_github("LandSciTech/LSTDConnect")
}
library(LSTDConnect) ###
library(DBI); library(RSQLite)
library(tidyverse)
library(dbplyr);library(stringr)
library(raster)
library(samc)
source("scripts/batchHelperFns.R")

##############################
#Create empty SQLite database -  done once prior to starting the batch run.
exptName= paste0(exptNameIn,nameAddBit)
exptDir = paste0("output/",exptName)
exptTab = read.csv(paste0("./input/",exptNameIn,".csv"))

exptTab = exptTab[order(exptTab$scnID),]
exptTab$status = 'not done'
exptTab$sampleRowID = seq(1:nrow(exptTab))

#create blank output table
outputTab = data.frame(scnID=exptTab$scnID,sampleRowID=exptTab$sampleRowID,outVal = NA, resultsPath = NA,elapsedTime=NA)
outputTab = subset(outputTab,scnID<0)

myDBName = initiateExpt(exptName,exptTab,outputTab)
#Now we have a database with an experimental design table (one row per map processing event) and an empty output table.

expectFile = paste0(exptDir,"/",exptName,"DB.sqlite")
if(!file.exists(expectFile)){
  stop("Database creation failed.")
}
#proceed to batchRun.R