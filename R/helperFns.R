#' @export
setLandclassColours<-function(landclasses,landclassColours=NULL){
  #http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
  if(is.null(landclassColours)){
    landclasses$color=""
    landclasses$color[landclasses$Landclass=="water"]= "cadetblue2"#blue
    landclasses$color[landclasses$Landclass=="road"]="black"
    landclasses$color[landclasses$Landclass=="road P"]="black"
    landclasses$color[landclasses$Landclass=="urban"]="black"# dark grey
    landclasses$color[landclasses$Landclass=="barren"]="gray90"# light grey
    landclasses$color[landclasses$Landclass=="barren P"]="gray100"# light grey
    landclasses$color[landclasses$Landclass=="cropland"]="gray50"
    landclasses$color[landclasses$Landclass=="cutblock"]="gray70"
    landclasses$color[landclasses$Landclass=="grass"]="lightgoldenrod2"
    landclasses$color[landclasses$Landclass=="grass P"]="lightgoldenrod1"
    landclasses$color[landclasses$Landclass=="young forest"]="darkorange2"
    landclasses$color[landclasses$Landclass=="young forest P"]="darkorange1"
    landclasses$color[landclasses$Landclass=="old forest"]="green3"
    landclasses$color[landclasses$Landclass=="old forest P"]="green2"
  }else{
    landclasses$color=NULL
    landclasses = merge(landclasses,landclassColours)
  }
  return(landclasses)
}

#' @export
stripX<-function(x,parts=1,sep="_"){
  #x=testSet$landscape[1]
  xparts=strsplit(x,sep)[[1]]
  xparts=xparts[(parts+1):length(xparts)]
  return(paste(xparts,collapse=sep))
}

mergeParcConnectedness<-function(cPC,outPC=NULL,memoryLimit,clumpName = "c"){
  #outPC=NULL;memoryLimit=1;clumpName=ic
  #outPC=cPC
  cPC@Mbar_p$clump=clumpName
  cPC@Mbar$clump=clumpName

  if(is.null(outPC)){
    return(cPC)
  }

  outPC@patches=raster::merge(cPC@patches,outPC@patches)
  if(as.numeric(object.size(outPC@d_ij))/10^9<memoryLimit){
    #names(cPC@d_ij)="clump2"
    outPC@d_ij = c(outPC@d_ij,cPC@d_ij)
    outPC@H_j = c(outPC@H_j,cPC@H_j)
  }else{
    addNames = names(cPC@d_ij)
    addBit = list("Memory limit exceeded.")
    for(adn in addNames){
      names(addBit)=adn
      outPC@d_ij=c(outPC@d_ij,adn)
      outPC@H_j=c(outPC@H_j,adn)
    }
  }

  addAlphas= names(cPC@M_i)
  for(aName in addAlphas){
    #aName = addAlphas[1]
    if(!is.element(aName,names(outPC@M_i))){
      outPC@M_i=raster::addLayer(cPC@M_i[[aName]])
    }else{
      outPC@M_i[[aName]]=raster::merge(outPC@M_i[[aName]],cPC@M_i[[aName]])
    }
  }

  outPC@Mbar_p=rbind(outPC@Mbar_p,cPC@Mbar_p)
  outPC@Mbar = rbind(outPC@Mbar,cPC@Mbar)
  return(outPC)
}

getVs<-function(inMap,id=T,locs=NULL,omit0=T){
  #TO DO: likely there is a more efficient straightforward way to do this. Fix if it matters.
  #inMap = x$patches; id = F; locs = names(R_iv);omit0=T

  ras.matrix<-raster::as.matrix(inMap)
  weight<-c(t(ras.matrix))
  weight<-data.frame(weight=weight)
  weight$id = seq_len(nrow(weight))
  subset(weight,weight>0)
  weight = subset(weight,!is.na(weight))
  if(omit0){
    weight = subset(weight,weight>0)
  }
  if(!is.null(locs)){
    #locs = c(3,7)
    locs=as.numeric(locs)
    missingLocs = setdiff(locs,weight$id)
    if(length(missingLocs)>0){
      stop("inMap does not include these selected locations:",paste(missingLocs,collapse=","))
    }
    weight= merge(data.frame(id=locs),weight)
  }
  if(id){
    return(weight$id)
  }else{
    return(weight$weight)
  }
}
