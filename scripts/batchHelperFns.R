
processingFunction<-function(threadId,Nbatch,myDBName,exptName,innerFn,continue,initWait=240,printError=F,lcdir="./nonVersionedInput/LayerCatalog/"){
  #threadID=1;innerFn = applyBufferedSum;Nbatch=1;continue = T;lcdir="./nonVersionedInput/LayerCatalog/"
  batchDir =paste0("output/",exptName,"/tabs")
  mapDir = paste0("output/",exptName,"/maps")

  dir.create(batchDir,recursive = T)
  dir.create(paste0(batchDir,"/errorReport"),recursive = T)

  dir.create(mapDir, recursive = T)

  first=T
  batchNames=c()
  waitTime = 20*runif(1)+10

  Sys.sleep(runif(1)*initWait)
  while(continue|first){
    first=F
    dbL <- try(DBI::dbConnect(RSQLite::SQLite(), myDBName),silent=T)
    if(inherits(dbL,"try-error")){
      Sys.sleep(waitTime)
      waitTime=waitTime+10
      next
    }
    #pull & modify exptl table
    #note there is a split second between collecting & dbSendQuery when database could be modified by another process.

    exptSel = try(tbl(dbL,sql("SELECT * FROM expt")) %>% filter(status=="not done") %>%
      arrange(sampleRowID) %>% head(n=Nbatch) %>% collect(),silent=T)
    #exptSel = tbl(dbL,sql("SELECT * FROM expt")) %>% filter(status=="in progress") %>%  arrange(sampleRowID) %>% head(n=Nbatch) %>% collect()

      #exptSel = tbl(dbL,sql("SELECT * FROM expt")) %>% filter(status=="in progress") %>%  arrange(sampleRowID) %>% head(n=Nbatch) %>% collect()
    #exptSel = tbl(dbL,sql("SELECT * FROM expt"))%>% filter(status=="not done") %>% arrange(sampleRowID) %>% collect()

    #subset(exptSel,sampleRowID==1214)
    if(nrow(exptSel)==0){continue=F;break}

    #exptSel=subset(exptSel,paArea==max(exptSel$paArea))

    resp = try(dbSendQuery(dbL, paste("UPDATE expt SET status = 'in progress' WHERE scnID IN (",
                           paste0(shQuote(exptSel$scnID), collapse=","), ")" )),silent=T)

    numTrys = 0
    while(numTrys<10){
      try(DBI::dbDisconnect(dbL),silent=T)
      if(!inherits(resp,"try-error")){break}
      numTrys=numTrys+1
      Sys.sleep(4)
    }  
    
    if(inherits(exptSel,"try-error")){
      Sys.sleep(waitTime)
      waitTime=waitTime+10
      #waitTime=waitTime*(1+runif(1))
      next
    }


    #Output table for this batch of scenarios.
    scn.outputTab = data.frame(scnID=exptSel$scnID,sampleRowID = exptSel$sampleRowID,outVal = NA, resultsPath = NA,elapsedTime=NA)
    #outputBit<- scn.outputTab

    batchName = paste0("rows",min(exptSel$sampleRowID),"to",max(exptSel$sampleRowID))

    for(i in 1:nrow(exptSel)){
      #exptSel$rowID=seq(1:nrow(exptSel));exptSel$rowID[exptSel$paArea==min(exptSel$paArea)]
      #exptSel$rowID=seq(1:nrow(exptSel));exptSel$rowID[exptSel$sampleRowID==7839]

      #i =2;exptSel[i,]$paArea
      start.time<-Sys.time()
      #source("scripts/batchHelperFns.R")
      #exptSel[i,"kernelScale"]=80

      omaps<- try(innerFn(indir = lcdir, expt.row = exptSel[i,]),silent=T)
      if(inherits(omaps,"try-error")){
        #TO DO: handle memory errors differently - crash thread rather than recovering.
        scn.outputTab$resultsPath[scn.outputTab$scnID==exptSel$scnID[i]]=omaps[1]


        #write error reports to table directory.
        write.csv(x= data.frame("scenario" = c(exptSel$scnID[i]), "error" = omaps[1]), file = paste0(batchDir,"/errorReport/", exptSel$scnID[i],".csv"), row.names = F)

        if(printError){
          stop(omaps[1])
        }
        if(grepl("cannot allocate vector of size",scn.outputTab$resultsPath[scn.outputTab$scnID==exptSel$scnID[i]],fixed=T)){
          next
        }else{next}

      }
      omaps$fPatch[is.na(omaps$fPatch)]=0
      map = omaps$bufsum.Patch
      map[omaps$fPatch!=1]=NA
      meanValue <-raster::cellStats(map,stat="mean")

      #print(paste("Mean value of PA ", exptSel$paId[i], "is ", meanValue))
      try(writeRaster(omaps$bufsum.Patch, filename = paste0(mapDir,"/", exptSel$scnID[i], ".tif"), overwrite=TRUE),silent=T)

      scn.outputTab$elapsedTime[scn.outputTab$scnID == exptSel$scnID[i]]=Sys.time()-start.time
      scn.outputTab$outVal[scn.outputTab$scnID == exptSel$scnID[i]]=meanValue
      scn.outputTab$resultsPath[scn.outputTab$scnID==exptSel$scnID[i]]=paste0(mapDir,"/", exptSel$scnID[i],".tif")

      #resultsPath<- paste0(mapDir,"/", exptSel$scnID[i])
    }

    try(write.csv(scn.outputTab, file = paste0(batchDir, "/",batchName,"_out.csv"),row.names=F),silent=F)

    #modify exptl table
    batchNames = c(batchNames,batchName)

    numTrys = 0
    while(numTrys<2){
      dbL <- try(DBI::dbConnect(RSQLite::SQLite(), myDBName),silent=T)

      #only set done for runs without errors
      noErrorRuns = subset(scn.outputTab,!is.na(scn.outputTab$resultsPath)&&!grepl("Error",scn.outputTab$resultsPath))$scnID
      if(length(noErrorRuns)>0){
        resp = try(dbSendQuery(dbL, paste("UPDATE expt SET status = 'done' WHERE scnID IN (",
                                          paste0(shQuote(noErrorRuns), collapse=","), ")" )),silent=T)
      }
      resp = try(dbWriteTable(dbL,"output",scn.outputTab,append=T),silent=T)

      if(!inherits(resp,"try-error")){break}
      numTrys=numTrys+1
      Sys.sleep(waitTime)
      waitTime=waitTime+10
      #waitTime=waitTime*(1+runif(1))

    }
    numTrys = 0
    while(numTrys<10){
      try(DBI::dbDisconnect(dbL),silent=T)
      if(!inherits(resp,"try-error")){break}
      numTrys=numTrys+1
      Sys.sleep(4)
    }  
    #if memory size error, crash thread
    memSizeErrors = subset(scn.outputTab,grepl("cannot allocate vector of size",resultsPath))
    if(nrow(memSizeErrors)>0){
      stop(memSizeErrors$resultsPath[1])
    }

    if(inherits(resp,"try-error")){next}

    exptSel=NULL
    ##end continue for single thread testing (non-parallel)

  }

  #return(data.frame(threadId=threadId,batchNames=batchNames))
  return(batchNames)
}



applyBufferedSum<- function(indir,expt.row,maxFocalPatchCells=15000){


  # expt.row<-getScenario(dbName = "./output/mapSetA_exptPriorityLimit9RunParc4DB.sqlite", metricVariant = 'parc');
  #expt.row<- expt.row[1,]
  # indir='./nonVersionedInput/LayerCatalog/'

  #exptSel$rowID=seq(1:nrow(exptSel));exptSel$rowID[exptSel$paArea==min(exptSel$paArea)]
  #exptSel$rowID=seq(1:nrow(exptSel));exptSel$rowID[exptSel$paID==4354]

  #indir = lcdir; expt.row = exptSel[1,]

  #unique(subset(exptSel,paArea>3000000,select=c(paID,paArea,nameEco)))

  #data.frame(expt.row)
  #load the required layers into environment

  layer.dir<- paste0(indir,"pa",expt.row$paID)

  #layer.dir<-"./nonVersionedInput/LayerCatalog/pa4569"

  if(!dir.exists(layer.dir)){
    stop("Layer not included in catalogue:",expt.row$paID)
  }
  #patch raster
  #fPatch<- paste0(layer.dir, "/FocalPatch.tif") %>% raster::raster()
  #QUESTION: why is there a landscape without a focal patch??? Fix this. All sampled PAs should be real things...
  #NOTE: not necessary to save focal patch as a separate layer - save disk space in layer catalogue, and read/write time.

  #adjacent Protected Areas
  adjPAs<- paste0(layer.dir, "/patches.tif") %>% raster::raster()

  fPatch = adjPAs<0
  #fPatch[adjPAs==4569]=1
  fPatch[adjPAs==expt.row$paID]=1
  if(maxValue(fPatch)==0){
    warning("No focal patch in this landscape. Why? Fix this!")
  }

  adjPAs[is.na(adjPAs)]=0

  #quality raster
  qSurface<- grep(pattern = expt.row$qualitySurface, x = list.files(path = layer.dir, pattern = '.tif', recursive = T, full.names = TRUE), ignore.case = TRUE, value = TRUE ) %>% raster::raster()

  #plot(qSurface);plot(fPatch,add=T)
  #qSurface<-grep(pattern = 'kennedy', x = list.files(path = layer.dir, pattern = '.tif', recursive = T, full.names = TRUE), ignore.case = TRUE, value = TRUE ) %>% raster::raster()

  npPenalty<- expt.row$nonProtectedPenalty
  #npPenalty=0.75

  qSurface[adjPAs==0] = qSurface[adjPAs==0]*npPenalty
  #plot(qSurface)

  waterbuf<- expt.row$waterBuffer
  #waterbuf = 1

  if(waterbuf!=0){

    #NOTE: the 'Patch' map for each protected area was updated to include water (value=10) and terrestrial area outside 500km     buffered study area (value=5)
    # To add water cells, I used the Kennedy human disturbance raster as a template.
    # Please see './scripts/WorkflowSteps/Step3AddWaterFeatures.R' for the workflow I used for this process.

    qSurface[adjPAs==10]=1111 #placeholder value for water
    qSurface[adjPAs==5]=2222 #placeholder value for terrestrial area outside buffered study zone
    qSurface[is.na(qSurface)]=3333 #placeholder value for whitespace around the study area map; not included in buffer

    qSurface[qSurface==1111]=NA #set the water back to NA. these are the cells we want to buffer 'into'.

    wBuff= lsBuffer(!is.na(qSurface)&qSurface!=3333, waterbuf*res(qSurface)[1]) #run the buffer 'into' the water.

    wBuff[!is.na(qSurface)]=NA

    qSurface[wBuff==1]=100
    qSurface[qSurface==2222]=NA #revert whitespace back to NA
    qSurface[qSurface==3333]=NA #revert whitespace back to NA

  }

  landqualitymin<- expt.row$landQualityMin
  #landqualitymin=100
  if(landqualitymin>0){
    qSurface[qSurface<landqualitymin]=landqualitymin
  }

  #velox will report NA for any neighbourhood that includes an NA. So set water quality to 0, not NA
  qSurface[(adjPAs==10)&(is.na(qSurface))]=0
  qSurface[is.na(qSurface)]=0 #velox handles NA values very oddly - they seem to propogate through in unexpected ways. Better to remove them entirely.

  cellDim = res(qSurface)[1]
  dbar=expt.row$kernelScale*res(qSurface)[1]/cellDim
  k = exponentialKernel(dbar,negligible=10^-6)
  maxDistHold = (nrow(k)-1)/2
  #reduce size of problem by cutting off values > exp kernel size from park.
  cBuff = 500
  cutD = cBuff-maxDistHold
  qSurfaceC=qSurface
  qSurfaceC[1:cutD,]=NA
  qSurfaceC[(nrow(qSurfaceC)-cutD):nrow(qSurfaceC),]=NA
  qSurfaceC[,1:cutD]=NA
  qSurfaceC[,(ncol(qSurfaceC)-cutD):ncol(qSurfaceC)]=NA

  if(expt.row$metricVariant %in% c('uniform', 'Uniform')){
    bufsum.Patch=applyKernel(quality = qSurfaceC,
                                               kernel = 'Uniform',
                                               d=expt.row$kernelScale*res(qSurface)[1],useAveDist=T)
  }
  if(expt.row$metricVariant %in% c('exponential', 'Exponential')){
    bufsum.Patch= applyKernel(quality = qSurfaceC,
                                                kernel = 'Exponential',
                                                d= expt.row$kernelScale*res(qSurface)[1],negligible=10^-6)
  }
  if(grepl('samc',expt.row$metricVariant)){
    #velox will report NA for any neighbourhood that includes an NA. So set water quality to 0, not NA
    backgroundMortality =0.00001
    tSet= data.frame(t=c(4,23,91,361,1455),meanDisplacement=c(2,5,10,20,40))
    t =tSet$t[round(tSet$meanDisplacement)==expt.row$kernelScale]

    oneMap = qSurface;oneMap[oneMap!=1]=1
    if(!grepl("samcM",as.character(expt.row$metricVariant))){
      ds = expt.row$metricVariant;ds=gsub("B","",ds)
      ds=gsub("R","",ds)
      shapeMult = as.numeric(gsub("samc","",as.character(ds)))
      cCost =1+(shapeMult-1)*(100-qSurface)/100
      cMort <- oneMap;cMort[cMort==1]=backgroundMortality

    }
    if(grepl("samcM",expt.row$metricVariant)|(grepl("samcB",expt.row$metricVariant))){
      if(grepl("samcM",as.character(expt.row$metricVariant))){
        shapeMult = as.numeric(gsub("samcM","",as.character(expt.row$metricVariant)))/300
        cCost=oneMap
      }
      if(grepl("samcB",as.character(expt.row$metricVariant))){
        shapeMult = as.numeric(gsub("samcB","",as.character(expt.row$metricVariant)))/300
      }

      #shapeMult=0.01
      cMort <- oneMap;cMort[cMort==1]=backgroundMortality
      cMort = cMort+shapeMult*(1-cMort)*((100-qSurface)/100)
    }

    fPatch[(adjPAs==10)&(is.na(fPatch))]=0 #parc calculations include water within protected areas. Probably that should change...

    cCost[is.na(qSurfaceC)]=NA
    cMort[is.na(qSurfaceC)]=NA
    occ_prob_data=qSurfaceC
    occ_prob_data = as.matrix(occ_prob_data)

    samc_obj <- samc(as.matrix(cCost), as.matrix(cMort), tr_fun = function(x) 1/mean(x))

    short_disp = distribution(samc_obj,occ_prob_data,time=t)
    short_disp_map <- map(samc_obj, short_disp)
    bufsum.Patch = oneMap; bufsum.Patch@data@values = short_disp_map@data@values #+short_mort_map
    #plot(qSurface);plot(bufsum.Patch,add=T)
  }
  if(grepl('parc',expt.row$metricVariant)){
    #expt.row$kernelScale = 5
    #calculate maxDist from negligible value +alpha
    #qSurface=100-qSurface
    #temp  = matrix(1,nrow=1000000,ncol=10000000)

    shapeMult = as.numeric(gsub("parc","",as.character(expt.row$metricVariant)))
    #shapeMult = as.numeric(gsub("samcM","",as.character(expt.row$metricVariant)))

    exCost = 1+(shapeMult-1)*(100-qSurface)/100 #Drielsma 2007 example scales cost from 1 to 2.
    exCost[adjPAs==10]=shapeMult

    fPatch[(adjPAs==10)&(is.na(fPatch))]=0 #parc calculations include water within protected areas. Probably that should change...
    #qSurface[(adjPAs==10)&(is.na(qSurface))]=0

    focalCells<-ncell(fPatch[fPatch==1])

    if(focalCells>maxFocalPatchCells){
      bufsum.Patch<-splitParcScenarios(fPatch=fPatch, qSurface=qSurface, exCost=exCost, maxFocalCellsPerSegment=maxFocalPatchCells, segCount=4, cellDim=cellDim, dbar=dbar, k=k, maxDistHold=maxDistHold)
    }else{
      #maxDistHold=300
      testPC = parcConnectedness(x=stack(fPatch,exCost,qSurface),maxDist=maxDistHold*res(qSurface)[1],alpha=2/(dbar*res(qSurface)[1]),memoryLimit=0,stopOnMemoryLimit=T,neighbourhood="octagon")
      bufsum.Patch = testPC@M_i
      #plot(bufsum.Patch)
    }
  }

  #TO DO: save little output map to batchDir, appropriately named. And add relative path to output map in result table

  #bufsum.Mask = bufsum.Patch[(exPatch==0)|is.na(exPatch)]=NA
  #calculate mean value of cells within park boundary
  # mean.bufsum.patch<- getValues(bufsum.Patch)
  # mean.bufsum.patch<- mean.bufsum.patch[!is.na(mean.bufsum.patch)]
  # mean.bufsum.patch<- mean(mean.bufsum.patch)

  return(list(bufsum.Patch=bufsum.Patch,fPatch=fPatch))
}

splitParcScenarios<-function(fPatch, qSurface, exCost, maxFocalCellsPerSegment, segCount=4,cellDim, dbar, k, maxDistHold){
  #maxFocalCellsPerSegment=20000; segCount=4

  paCellCount<- ncell(fPatch[fPatch==1])

  #break up the patch layer until each segment has less than the max threshold of focal cells.
  while(paCellCount > maxFocalCellsPerSegment){

    patchSplit<- SpaDES.tools::splitRaster(fPatch, nx=0.5*segCount, ny=0.5*segCount)

    fpCells<-lapply(patchSplit, FUN = function(x){
      fpCellcount<- ncell(x[x==1])
      return(fpCellcount)}) %>% unlist()

    paCellCount<-max(fpCells)
    segCount=segCount+2
  }

  splitValues<- list()

  # discard segments that don't contain PA cells; make an "index" vector for subsetting the split layers to exclude chunks that contain no focal patch cells.
  for(i in 1:length(patchSplit)){
    p<- patchSplit[[i]]
    vals<-unique(getValues(p))
    if(1 %in% vals){
      splitValues<-append(splitValues,i)
    }
  }

  splitValues<-unlist(splitValues)
  patchSplit<-patchSplit[splitValues];splitValues=NULL

  #run parcConnectedness on each segment which has PA cells.
  patchSplitPC<-lapply(patchSplit, FUN= function(x){
    segmentLayer<-extend(x, qSurface)
    mask<-fPatch
    ### each segment should contain a portion of the focal patch & '0' values for the rest of the landscape.
    mask[!is.na(mask)]=0
    mask[segmentLayer==1]=1
    segmentLayer=mask;rm(mask)

    pc<-parcConnectedness(stack(segmentLayer,exCost,qSurface), maxDist=maxDistHold*res(qSurface)[1], alpha=2/(dbar*res(qSurface)[1]), memoryLimit=0, stopOnMemoryLimit=T)
    return(pc)
  })

  #list of M_i layers from the resulting parcConnectedness objects.
  pcRasters<- lapply(patchSplitPC, FUN = function(x){
    m_i<-x@M_i[[1]]
  })

  #merge into a single map
  bufsum.patch<- do.call(merge,pcRasters)

  #bsum.patch<-mosaic(pcRasters[[1]],pcRasters[[2]],pcRasters[[3]],pcRasters[[4]],fun='max')

  return(bufsum.patch)
}


initiateExpt<-function(exptName,exptTab,outputTab){
  #set up SQLite database
  #Here is a useful R database tutorial:
  #https://datacarpentry.org/R-ecology-lesson/05-r-and-databases.html
  #exptName =paste0(exptName,"Test")
  exptDir =  paste0("output/",exptName)
  if(dir.exists(exptDir)){stop("Directory ",exptDir,"already exists. Please remove the directory or change exptName to initiate a new experiment.")}
  dir.create(exptDir,recursive = T)
  myDBName <- paste0(exptDir,"/",exptName,"DB.sqlite")

  myDB <- dbConnect(RSQLite::SQLite(),myDBName)

  dbWriteTable(myDB, "expt", exptTab)
  dbWriteTable(myDB, "output", outputTab)
  dbListTables(myDB)
  dbDisconnect(myDB)#Don't want stray database connections accumulating.
  return(myDBName)
}

#Function for tidying up 'in progress' rows at the end of a processing run.

tidyInProgressResults<- function(threadId,Nbatch,myDBName,exptName,innerFn,continue,lcdir="./nonVersionedInput/LayerCatalog/"){

  batchDir =paste0("./output/",exptName,"/tabs")
  mapDir = paste0("./output/",exptName,"/maps")

  if(dir.exists(batchDir) &  dir.exists(mapDir)!=TRUE){
    stop("Check the experiment name -- cannot write to the results path.")
  }

  first=T
  batchNames=c()
  while(continue|first){
    first=F
    dbL <- try(DBI::dbConnect(RSQLite::SQLite(), myDBName),silent=T)
    if(inherits(dbL,"try-error")){
      Sys.sleep(0.2)
      next
    }
    #pull & modify exptl table
    #note there is a split second between collecting & dbSendQuery when database could be modified by another process.


    exptSel = tbl(dbL,sql("SELECT * FROM expt")) %>% filter(status!="done") %>% arrange(sampleRowID) %>% head(n=Nbatch) %>% collect()
    if(nrow(exptSel)==0){continue=F;break}

    if(inherits(exptSel,"try-error")){
      Sys.sleep(0.2)
      next
    }

    numTrys = 0
    while(numTrys<10){
      resp = try(dbSendQuery(dbL, paste("UPDATE expt SET status = 'tidying in progress' WHERE scnID IN (",
                                        paste0(shQuote(exptSel$scnID), collapse=","), ")" )),silent=T)

      if(!inherits(resp,"try-error")){break}
      numTrys=numTrys+1
      Sys.sleep(0.2)

    }
    try(dbDisconnect(dbL),silent=T)

    if(inherits(resp,"try-error")){next}

    #Output table for this batch of scenarios.
    scn.outputTab = data.frame(scnID=exptSel$scnID,sampleRowID = exptSel$sampleRowID,outVal = NA, resultsPath = NA,elapsedTime=NA)

    batchName = paste0("rows",min(exptSel$sampleRowID),"to",max(exptSel$sampleRowID))

    for(i in 1:nrow(exptSel)){
      #i =1
      start.time<-Sys.time()

      omaps<- try(innerFn(indir = lcdir, expt.row = exptSel[i,]),silent=T)
      if(inherits(omaps,"try-error")){
        scn.outputTab$resultsPath[scn.outputTab$scnID==exptSel$scnID[i]]=omaps[1]
        next
      }
      omaps$fPatch[is.na(omaps$fPatch)]=0

      map = omaps$bufsum.Patch
      map[omaps$fPatch!=1]=NA
      meanValue <-raster::cellStats(map,stat="mean")

      #print(paste("Mean value of PA ", exptSel$paId[i], "is ", meanValue))
      try(writeRaster(omaps$bufsum.Patch, filename = paste0(mapDir,"/", exptSel$scnID[i], ".tif"), overwrite=TRUE),silent=T)

      scn.outputTab$elapsedTime[scn.outputTab$scnID == exptSel$scnID[i]]=Sys.time()-start.time
      scn.outputTab$outVal[scn.outputTab$scnID == exptSel$scnID[i]]=meanValue
      scn.outputTab$resultsPath[scn.outputTab$scnID==exptSel$scnID[i]]=paste0(mapDir,"/", exptSel$scnID[i],".tif")

      #resultsPath<- paste0(mapDir,"/", exptSel$scnID[i])
    }
    try(write.csv(scn.outputTab, file = paste0(batchDir, "/",batchName,"_out.csv"),row.names=F),silent=F)

    #modify exptl table
    batchNames = c(batchNames,batchName)

    numTrys = 0
    while(numTrys<1000){
      dbL <- try(DBI::dbConnect(RSQLite::SQLite(), myDBName),silent=T)

      resp = try(dbSendQuery(dbL, paste("UPDATE expt SET status = 'done' WHERE scnID IN (",
                                        paste0(shQuote(exptSel$scnID), collapse=","), ")" )),silent=T)
      resp = try(dbWriteTable(dbL,"output",scn.outputTab,append=T),silent=T)

      if(!inherits(resp,"try-error")){break}
      numTrys=numTrys+1
      Sys.sleep(0.2)
    }
    try(DBI::dbDisconnect(dbL),silent=T)

    if(inherits(resp,"try-error")){next}

    exptSel=NULL
    ##end continue for single thread testing (non-parallel)

  }

  #return(data.frame(threadId=threadId,batchNames=batchNames))
  return(batchNames)
}

testInner<-function(indir, expt.row){
  layer.dir<- paste0(indir,"pa",expt.row$paID)

  #layer.dir<-"./nonVersionedInput/LayerCatalog/pa4569"

  if(!dir.exists(layer.dir)){
    stop("Layer not include in catalogue:",expt.row$paID)
  }
  #patch raster
  #fPatch<- paste0(layer.dir, "/FocalPatch.tif") %>% raster::raster()
  #QUESTION: why is there a landscape without a focal patch??? Fix this. All sampled PAs should be real things...
  #NOTE: not necessary to save focal patch as a separate layer - save disk space in layer catalogue, and read/write time.

  #adjacent Protected Areas
  adjPAs<- paste0(layer.dir, "/patches.tif") %>% raster::raster()

  fPatch = adjPAs<0
  #fPatch[adjPAs==4569]=1
  fPatch[adjPAs==expt.row$paID]=1
  return(list(bufsum.Patch=fPatch,fPatch=fPatch))
}
if(0){

  checkDir = paste0(mapDir,"/", exptSel$scnID[i], ".tif")
  if(file.exists(checkDir)){

    checkNA = raster(checkDir)

    layer.dir<- paste0(lcdir,"pa",exptSel$paID[i])

    adjPAs<- paste0(layer.dir, "/patches.tif") %>% raster::raster()
    fPatch = adjPAs<0
    #fPatch[adjPAs==4569]=1
    fPatch[adjPAs==exptSel$paID[i]]=1
    sumNA = cellStats((fPatch==1)&is.na(checkNA),"sum")
  }else{
    sumNA=1
  }
  if(sumNA==0){
    omaps = list(fPatch=fPatch,bufsum.Patch=checkNA)
  }else{
    omaps<- try(innerFn(indir = lcdir, expt.row = exptSel[i,]),silent=T)
  }
}

