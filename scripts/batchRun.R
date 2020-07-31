### Start or restart parallel process to populate SQL Database ##
#Note: this script will resume processing the status = 'not done' rows in the 'expt' table. Rows that were 'in process' at the time of the interruption will not be handled by this script. Use 'tidyInProcessRows.R' to collect and finish processing these layers at the end of the run.

# Requires the user to input the 'experiment name' to access the appropriate SQL database initiated at the start of the run.
#   Run batchSetup.R to create the empty database.
#   LSTDConnect/output/'experimentname'/'experimentname'DB.sqlite

#load library and source files into the environment.
library(DBI); library(RSQLite);library(tidyverse);library(dbplyr);
library(stringr);library(raster)
#Use local dev version of LSTDConnect package.
devtools::document();devtools::load_all()
#library(LSTDConnect) 
source("scripts/batchHelperFns.R")
library(parallel)

# Parameters
###################################################
#set the name of the experiment to resume.
exptName = "mapSetA_exptMinExampleTest"

#set the database path.
myDBName= paste0("output/", exptName,"/",exptName,"DB.sqlite")

#set the batch size
Nbatch=2

################################################################
#connect to database.
myDB<- dbConnect(RSQLite::SQLite(), myDBName)

no_cores <- detectCores() - 1
results=c(1:no_cores)

# Initiate cluster
cl <- makeCluster(no_cores)
clusterEvalQ(cl, library(DBI))
clusterEvalQ(cl, library(dbplyr))
clusterEvalQ(cl, library(dplyr))
clusterEvalQ(cl, source("scripts/batchHelperFns.R"))
clusterEvalQ(cl,library(samc))
clusterEvalQ(cl, library(LSTDConnect))
clusterEvalQ(cl, library(raster))

#Apply processingFunction in parallel.
#See batchHelperFns.R for processingFunction and applyBufferedSum code.
testResPar = parLapply(cl,results,processingFunction,Nbatch=Nbatch,myDBName=myDBName,exptName=exptName,
                       innerFn=applyBufferedSum,continue=T)
stopCluster(cl)
dbDisconnect(myDB)
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

