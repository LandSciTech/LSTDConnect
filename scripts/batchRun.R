### Start or restart parallel process to populate SQL Database ##
#Note: this script will resume processing the status = 'not done' rows in the 'expt' table. Rows that were 'in process' at the time of the interruption will not be handled by this script. Use 'tidyInProcessRows.R' to collect and finish processing these layers at the end of the run.

# Requires the user to input the 'experiment name' to access the appropriate SQL database initiated at the start of the run.
#   Run batchSetup.R to create the empty database.
#   LSTDConnect/output/'experimentname'/'experimentname'DB.sqlite

# Parameters
###################################################
#set the name of the experiment to resume.
exptName = "mapSetA_exptMinExampleTest"

#set the database path.
myDBName= paste0("output/", exptName,"/",exptName,"DB.sqlite")

#set the batch size
Nbatch=2

#########################
#load library and source files into the environment.
#Note: for parallel processing need to use package installed from github rather than local dev version.
#Ensure any recent changes to the package have been committed before proceeding
if(0){# reinstall package if necessary.
  Sys.setenv(TZ = "UTC")
  #deps <- tools::package_dependencies("LSTDConnect", recursive = TRUE)
  #update.packages(oldPkgs = unlist(deps), ask = FALSE)
  devtools::install_github("LandSciTech/LSTDConnect")
}
library(LSTDConnect)
library(DBI); library(RSQLite);library(tidyverse);library(dbplyr);
library(stringr);library(raster)
library(parallel)
library(samc)
source("scripts/batchHelperFns.R")

################################################################
#run one to check for obvious errors
dbL <- DBI::dbConnect(RSQLite::SQLite(), myDBName)
exptSel = tbl(dbL,sql("SELECT * FROM expt")) %>%
                arrange(sampleRowID) %>% head(n=100) %>% collect()
DBI::dbDisconnect(dbL)
lcdir="./nonVersionedInput/LayerCatalog/"
omaps<- applyBufferedSum(indir = lcdir, expt.row = exptSel[2,])
plot(omaps$bufsum.Patch)#;plot(omaps$fPatch,add=T)

check = processingFunction(Nbatch=2,myDBName=myDBName,exptName=exptName,
                           innerFn=applyBufferedSum,initWait=0,continue=F,printError=T)
#see example output
example= raster(paste0("output/", exptName,"/maps/2.tif"))
plot(example)

#####################
#run in parallel
no_cores <- detectCores() - 1
results=c(1:no_cores)

# Initiate cluster
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(DBI))
clusterEvalQ(cl, library(dbplyr))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl,library(samc))
clusterEvalQ(cl, library(LSTDConnect))
clusterEvalQ(cl, library(raster))

#Apply processingFunction in parallel.
#See batchHelperFns.R for processingFunction and applyBufferedSum code.
#If database access conflicts try increasing initWait parameter.
testResPar = parLapply(cl,results,processingFunction,Nbatch=Nbatch,myDBName=myDBName,exptName=exptName,
                       innerFn=applyBufferedSum,initWait=10,continue=T)
stopCluster(cl)
#Output
#####################################
# 1. writes map to /maps folder under the '/experiment name' directory. Map filename is the scnID for that particular scenario.
# 2. updates status field in 'expt' table in SQL DB to 'done'
# 3. writes row to 'output' table in SQL DB with scnID, mean map value, path to map, and ############

########################
#Notes:

#The problem with this implementation is SQLite database access conflicts. 
#For now we have done hokey things to force delays and stagger timing - see processingFunction in batchHelperFns.R for details.
#It doesn't work well.

#It would be helpful to have a general solution to the problem of batch processing scenarios. 
#The need arises often and is not specific to connectivity analysis. 
#Note a key feature of a good batch processsing solution is no need to manually assign sets of scenarios to threads.
#The algorithm should be able to figure out what is next on the cue, 
#and each thread should continue doing the next thing until all jobs are done.

