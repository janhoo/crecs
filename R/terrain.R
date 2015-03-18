wtf<-function(r)
{
	stopifnot(inherits(r, "RasterLayer"))
	return(as.matrix(r))
}


#' @title bathytopografic position index
#' 
#' @description Returns different tpi's.
#' 
#' @param dem the elevation model as \code{RasterLayer}, \code{SpatialGridDataFrame} or \code{SpatialPointsDataFrame} 
#' @param orad outer circle radius in map units for which the tpi is calculated
#' @param irad inner circle radius. If given, an annulus is used fot tpi caculation. If 0, the inner radius contains exacly one pixel. Default is \code{NULL}
#' @param slope threshold to delineate flat from slope. Unit is 'degrees'. In case of default \code{NULL}, single standart deviation of slopes is used
#' @param fac modifier to slope threshold. Applies only if \code{simple=TRUE}.
#' @param simple output only 3 categories: (1) valley; (2) slope/flat; (3) ridge/top. Overrides \code{simpler}
#' @param simpler simplify the non-simple schema to 6 categories:  (1) valley; (2) lower slope; (3) flat; (4) middle slope; (5) upper slope; (6) ridge/top
#' @param scale apply preemptive \code{\link{scale}}
#' @param inv switch sign (in case of bathemetry with positive depth values)
#' @param verbose be verbose and plot 
#' @details #FIXME
#' #FIXME
#' @examples #FIXME
#' 
bpi<-function(dem,orad,irad=NULL,slope=NULL,fac=1,simple=TRUE,simpler=TRUE,scale=FALSE,inv=FALSE,verbose=FALSE)
{
	stopifnot(class(slope)%in%c("NULL","numeric"))
	stopifnot(class(irad)%in%c("NULL","numeric"))
	stopifnot(class(orad)=="numeric")
	stopifnot(class(dem)%in%c("RasterLayer","SpatialGridDataFrame","SpatialPixelsDataFrame"))
	if(verbose){
		bpi.colors <-  colorRampPalette(c("steelblue4", "turquoise1", "yellowgreen", "yellow",   "darkorange", "red"))
	}
	if(class(dem)!="RasterLayer"){
		dem<-raster(dem)
	}
	if(inv){
		dem<-dem*(-1)
	}
	
	mo<-focalWeight(dem,d=orad,type="circle")
	if(is.null(irad)){
		M<-mo
	} else {
		# calc annulus matrix M
		mi<-focalWeight(dem,d=irad,type="circle")
		xoffset<-(ncol(mo)-ncol(mi))/2
		yoffset<-(nrow(mo)-nrow(mi))/2
		ro<-raster(mo) 
		ri<-raster(mi)
		extent(ro)<-extent(0,ncol(mo),0,nrow(ro))
		extent(ri)<-extent(xoffset,xoffset+ncol(mi),yoffset,yoffset+nrow(mi))
		
		m1<-as.matrix(extend(x=ri,y=extent(ro),value=0))
		m2<-m1/max(m1)
		mo2<-mo/max(mo)
		M1<-mo2-m2
		M<-M1/sum(M1)
	}
	# focal mean
	if(scale){
		dem<-scale(dem)
	}
	ra.tpi<-dem-focal(x=dem,w=M,na.rm=T,pad=F)
	# normalize
	if(scale){
		sdf<-1
	} else {
		sdf<-cellStats(ra.tpi, stat="sd",na.rm=TRUE, asSample=F)
	}
	
	if(simple){
		m<-matrix(ncol=3,byrow = TRUE, data=c(-Inf,(-sdf)*fac,1, (-sdf)*fac,sdf*fac,2,sdf*fac,Inf,4)) #(0) Valley; (1) slope/flat; (4) ridge
		BPI3<-reclassify(x=ra.tpi,rcl=m,right=TRUE)
		if(verbose){
			plot(BPI3,col=bpi.colors(3))
		}
		return(BPI3)
	} 
	m<-matrix(ncol=3,byrow = TRUE, data=c(-Inf,(-sdf),1, (-sdf),(0.5*(-sdf)),2,(0.5*(-sdf)),(0.5*sdf),3,(0.5*sdf),sdf,5,sdf,Inf,6)) # (1) valley; (2) lower slope; (3) flat [4-middle slope] ;(5) upper slope; (6) ridge
	ra.tpi.cl<-reclassify(x=ra.tpi,rcl=m,right=TRUE)               # (1) valley; (2) lower slope; (3) flat [4-middle slope] ;(5) upper slope; (6) ridge
	ra.slope<-terrain(dem,opt= "slope",neighbors = 8,unit = 'degrees')
	if(is.null(slope)){
		slope.threshold<-cellStats(ra.slope, stat="sd",na.rm=TRUE, asSample=F)
	} else {
		slope.threshold<-slope
	}
	ra.slope.cl<-reclassify(ra.slope,c(0,slope.threshold,0,slope.threshold,Inf,0.5))
	B<-ra.tpi.cl+ra.slope.cl
	BPI<-reclassify(B,c(3.1,3.9,4)) # middle slope
	if(simpler){
		BPI<-floor(BPI)
	} 
	if(verbose){
		plot(BPI,col=bpi.colors(length(table(BPI@data@values))))
	}
	return(BPI)
}




#' @title landform analysis indices
#' 
#' @description Returns landform analyses
#' 
#' @param dem the elevation model as \code{RasterLayer} 
#' @param large.bpi large scale \code{\link{bpi}} raster. Must use \code{simple=TRUE}
#' @param small.bpi small scale \code{\link{bpi}} raster. Must use \code{simple=TRUE}
#' @param slope threshold to delineate flat from slope. Unit is 'degrees'. In case of default \code{NULL}, single standart deviation of slopes is used
#' @param factors make factors (also include landform names)
#' @param verbose be verbose and plot 
#' @param trim trim the output raster to non-NA extent 
#' @details #FIXME
#' #FIXME
#' @examples #FIXME
#' 
lfi<-function(dem,large.bpi,small.bpi,slope.threshold=NULL,factors=TRUE,verbose=FALSE,trim=F)
{
	rll<-large.bpi*4
	sl<-terrain(dem,opt= "slope",neighbors = 8,unit = 'degrees')
	if(is.null(slope.threshold)){
		slth<-cellStats(sl, stat="sd",na.rm=TRUE, asSample=F)
	} else {
		slth<-slope.threshold
	}
	slc<-reclassify(sl,c(0,slth,0,slth,Inf,0.5))
	R<-small.bpi+rll+slc
	df<-data.frame(from=c(-Inf,5.9,7.9,8.9,9.3,9.9,10.3,11.9,12.3,16.9,17.9,19.9),
				   to=c(5.9,7.9,8.9,9.3,9.9,10.3,11.9,12.3,16.9,17.9,19.9,Inf),
				   becomes=c(1,5,9,2,3,6,7,10,11,4,8,12))
	m<-as.matrix(df)
	r<-reclassify(R,m,right=T)
	if(trim){
		r<-trim(r)
	}
	if(verbose){
		slope <- terrain(dem, opt='slope')
		aspect <- terrain(dem, opt='aspect')
		hill <- hillShade(slope, aspect, angle=40, direction=90)
		if(trim){
			hilli<-crop(hill,extent(r))
		} else {
			hilli<-crop(hill,extent(trim(r)))
		}
		lfi.colors <-  colorRampPalette(c("darkviolet", "orchid3", "orchid", "royalblue", "turquoise1", "mediumseagreen", "olivedrab1", "yellow", "wheat3", "sienna1", "sienna2", "sienna4"))
		plot(hilli, col=grey(0:100/100), legend=FALSE)
		plot(r,col=lfi.colors(12), alpha=0.35, legend=F, add=T)
		legend(x="topleft",title="Landforms", inset=.05, bty = "n", horiz=F,cex=0.7,fill=(lfi.colors(12)),c("canyons, deeply incised streams","local valleys in plains","lateral midslope incised drainages","upland incised drainages, headwaters","U-shape valleys","broad flat plains","open slopes","upper slopes, mesas","Local ridges in valleys","Local ridges in plains","lateral midslope drainage divides","mt tops, high ridges"))
		#plot(r,col=(heat.colors(12)), alpha=0.50, legend=F, add=T)
	}
	# factoring
	if(factors){
		landforms<-c("canyons, deeply incised streams","local valleys in plains","lateral midslope incised drainages","upland incised drainages, headwaters","U-shape valleys","broad flat plains","open slopes","upper slopes, mesas","Local ridges in valleys","Local ridges in plains","lateral midslope drainage divides","mt tops, high ridges")
		f<-ratify(r)
		rat <- levels(f)[[1]]
		rat$landcover<-landforms[levels(f)[[1]]$ID]
		levels(f) <- rat
		r<-f
	}
	return(r)
}








