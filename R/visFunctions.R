#' @export
plotComposites<-function(outDir,outName,panelVars,xVars,yVars,renameList=list(cropland="cropland",urban="urban",cutblock="cutblock"),borderWidth=3,borderValues=NULL,skipMaps=F,addName="",yLab=NULL,xLab=NULL,responseName=NULL,addMaps=NULL,landclassColours =NULL,rep=NULL,includeCombos=NULL){
  #renameList=list(cropland="cropland",urban="urban",cutblock="cutblock")
  #borderWidth=3; skipMaps=F;addName="costHumanFootprint";borderValues=NULL
  #borderValues=testSetN;skipMaps=T;addName=aAddName;yLab="Anthropogenic disturbance amount & type";xLab="Landscape type. All are 50% forest.";responseName = "Colonization potential"
  #renameList=list(cropland="cropland",urban="urban",cutblock="cutblock");borderWidth=3;borderValues=NULL;skipMaps=F;addName="";yLab=NULL;xLab=NULL;responseName=NULL;addMaps=NULL;landclassColours =NULL;rep=1

  #load maps - see CompositeMap.R for map creation example
  if(!skipMaps){
    #if (is.null(inMaps)){
      composites=NULL
      if(!is.null(rep)){
        load(paste0(outDir,"/",outName,"_rep",rep,".RData"))

      }else{
        load(paste0(outDir,"/",outName,".RData"))
      }
    #}else{
    #  composites=inMaps
    #}
  }
  combos = read.csv(paste0(outDir,"/",outName,"_combos.csv"))
  landclasses = read.csv(paste0(outDir,"/",outName,"_landclasses.csv"))
  missed=setdiff(names(combos),c(panelVars,xVars,yVars))
  combos=getName(combos)
  #str(combos)
  #tail(combos)
  if(!is.null(includeCombos)){
    combos=merge(combos,includeCombos)
    for (nn in names(includeCombos)){
      #nn="anthType"
      combos[[nn]]=as.character(combos[[nn]])
      combos[[nn]]=factor(combos[[nn]],levels = levels(includeCombos[[nn]]))
    }
  }

  if(!is.null(borderValues)){
    if(!is.element("getName",names(borderValues))){stop("Expecting getName column in border values")}
    if(!is.element("borderVal",names(borderValues))){stop("Expecting borderVal column in border values")}

    missingBorderVals = setdiff(combos$getName,borderValues$getName)
    if(length(missingBorderVals)>0){warning("Missing border values. Should have one for each row of combos")}
    combos$borderVal=NULL
    combos=merge(combos,subset(borderValues,select=c("getName","borderVal",intersect(panelVars,names(borderValues)))),all.x=T)
  }

  landclasses=setLandclassColours(landclasses,landclassColours)

  allPanels = unique(subset(combos,select=panelVars))
  if(skipMaps){
    allPanels$panelName=""
    if(length(panelVars)==1){
      allPanels$panelName=allPanels[[panelVars]]
    }else{
      for(p in panelVars){
        #p=panelVars[[1]]
        allPanels$panelName=paste(allPanels$panelName,allPanels[[p]],sep=".")
      }
    }

    cCombos=combos
    #cCombos$borderVal = runif(nrow(cCombos))*10000

    xLabs =getLabs(cCombos,xVars,renameList);names(xLabs)[names(xLabs)=="lab"]="xLab"
    yLabs=getLabs(cCombos,yVars,renameList);names(yLabs)[names(yLabs)=="lab"]="yLab"
    cCombos$xLab=NULL;cCombos$yLab=NULL;cCombos$panelName=NULL
    cCombos=merge(cCombos,xLabs)
    cCombos=merge(cCombos,yLabs)
    cCombos=merge(cCombos,allPanels)

    base = ggplot(cCombos,aes(x=xLab,y=yLab,fill=borderVal))+geom_tile()+facet_wrap(~panelName)
    if(!is.null(xLab)){base=base+xlab(xLab)}
    if(!is.null(yLab)){base=base+ylab(yLab)}
    if(!is.null(responseName)){base=base+ labs(fill=responseName) }
    base=base+theme(axis.text.x = element_text(angle = -45,hjust = 0))
    base=base+scale_fill_viridis_c(option="magma")#scale_fill_distiller(palette="YlOrRd")
    pdf(paste0(outDir,"/",outName,addName,".pdf"),width=10,height=6)
    print(base)
    dev.off()

  }else{
    if(nrow(allPanels)==0){
      allPanels=data.frame(panel="all")
    }
    for(i in 1:nrow(allPanels)){
      #i=3
      cPanel = allPanels[i,]
      cPanelName=paste(t(cPanel),collapse="_")

      print(cPanelName)
      cCombos = merge(combos,cPanel)

      outMapList = plotXYSet(composites,cCombos,xVars,yVars,borderWidth,renameList=renameList)
      fullMap=outMapList$fullMap;xLabs=outMapList$xLabs;yLabs=outMapList$yLabs;cCombos=outMapList$cCombos
      fullMap=setRat(fullMap,landclasses)

      if(!is.null(addMaps)){
        addMapList = plotXYSet(addMaps,cCombos,xVars,yVars,borderWidth,renameList=renameList)
        addMap=addMapList$fullMap
        addMap[addMap<1]=NA
        aP = rasterVis::levelplot(addMap,xlab="", ylab="",
                                  scales=list(x=list(at=(xLabs$startX+0.5*ncol(composites[[1]]))*res(fullMap)[1],labels=xLabs$xLab,rot=15),
                                              y=list(at=(yLabs$startY+0.5*nrow(composites[[1]]))*res(fullMap)[1],labels=yLabs$yLab,rot=55)),
                                  maxpixels=2e6,margin=F,colorkey=list(space="right"))
      }

      pdf(paste0(outDir,"/",outName,cPanelName,addName,".pdf"),width=10,height=6.5)
      bP = rasterVis::levelplot(fullMap, col.regions=levels(fullMap)[[1]]$color, xlab="", ylab="",
                      scales=list(x=list(at=(xLabs$startX+0.5*ncol(composites[[1]]))*res(fullMap)[1],labels=xLabs$xLab,rot=15),
                                  y=list(at=(yLabs$startY+0.5*nrow(composites[[1]]))*res(fullMap)[1],labels=yLabs$yLab,rot=55)),
                      maxpixels=2e6)
      if(!is.null(addMaps)){
        #print(plot(bP+aP))
        print(plot(aP + as.layer(bP, under = TRUE)))
      }else{
        print(plot(bP))
      }
      dev.off()
    }
  }
}
#' @export
setRat<-function(fullMap,landclasses,naCol="white"){
  fullMap = raster::as.factor(fullMap)
  rat <- raster::levels(fullMap)[[1]];rat$Landclass=NULL;rat$color=NULL
  lcs = unique(subset(landclasses,select=c(Landclass,LCID,color)))
  names(lcs)[names(lcs)=="LCID"]="ID"
  rat=merge(rat,lcs,all.x=T)
  rat$color=as.character(rat$color)
  rat$color[is.na(rat$Landclass)]=naCol
  levels(fullMap) <- rat
  return(fullMap)
}
plotXYSet<-function(composites,cCombos,xVars,yVars,borderWidth,renameList=list(),skipMaps=F){
  #cCombos = merge(combos,cPanel);skipMaps=T

  xLabs =getLabs(cCombos,xVars,renameList);names(xLabs)[names(xLabs)=="lab"]="xLab"
  yLabs=getLabs(cCombos,yVars,renameList);names(yLabs)[names(yLabs)=="lab"]="yLab"

  if(skipMaps){
    xDim=1;yDim=1;cRes=2
  }else{
    xDim = dim(composites)[1]
    yDim = dim(composites)[2]
    cRes=res(composites[[1]])[1]

  }

  fullWidth = nrow(xLabs)*(xDim+2*borderWidth)
  fullHeight = nrow(yLabs)*(yDim+2*borderWidth)

  fullMap=NULL
  fullMap = raster(nrow=fullHeight*cRes,ncol=fullWidth*cRes,xmn=1,xmx=fullWidth*cRes,ymn=1,ymx=fullHeight*cRes)
  res(fullMap)=c(cRes,cRes)

  xLabs$startX=1+seq(0,nrow(xLabs)-1)*(xDim+2*borderWidth)
  yLabs$startY=1+seq(0,nrow(yLabs)-1)*(yDim+2*borderWidth)

  cCombos$startX=NULL;cCombos$startY=NULL
  cCombos=merge(cCombos,xLabs)
  cCombos=merge(cCombos,yLabs)

  #TO DO - set border values for each combo
  if(!is.element("borderVal",names(cCombos))){
    cCombos$borderVal=-1
  }

  #for each row of cCombo, set map portion
  for(rr in 1:nrow(cCombos)){
    #rr=1
    r=cCombos[rr,]
    #r=subset(cCombos,(startX==1)&(startY==1))
    fullMap[dim(fullMap)[1]-(r$startY:(r$startY+yDim+2*borderWidth-1)),r$startX:(r$startX+xDim+2*borderWidth-1)]=r$borderVal

    if(!skipMaps){
      cm = names(composites)[grepl(r$getName,names(composites),fixed=T)]
      if(length(cm)>1){stop("something is wrong. Only one map of this name should be included in map set.")}
      if(length(cm)==0){next}
      cMap = composites[[cm]]#pull appropriate composite map
      fullMap[dim(fullMap)[1]-((r$startY+borderWidth):(r$startY+borderWidth+yDim-1)),(r$startX+borderWidth):(r$startX+borderWidth+xDim-1)]=cMap
    }
  }

  return(list(fullMap=fullMap,xLabs=xLabs,yLabs=yLabs,cCombos=cCombos))
  #figure out how to set x, y labels and plot.
}

getLabs<-function(cCombos,xVars,renameList){
  #xVars=yVars
  #assume combo columns are factors with order set by levels.
  xLabs = unique(subset(cCombos,select=xVars))
  xLabs$lab = "__"
  for(v in xVars){
    xLabs$lab=paste(xLabs$lab,xLabs[[v]],sep=".")
  }
  xLabs$lab=gsub("__.","",xLabs$lab,fixed=T)

  xLabs$order = 0
  cf=1
  for(i in length(xVars):1){
    #i=2
    xx=xVars[i]
    xLabs$order=xLabs$order+as.numeric(xLabs[[xx]])*cf
    cf = cf*10

  }
  xLabs=xLabs[order(xLabs$order),]
  head(xLabs)

  if(length(renameList)>0){
    for(l in 1:length(renameList)){
      xLabs$lab=gsub(names(renameList)[l],renameList[[l]],xLabs$lab,fixed=T)

    }
  }
  xLabs$lab=factor(xLabs$lab,levels=xLabs$lab)

  xLabs$order=NULL

  return(xLabs)
}
