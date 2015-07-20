# corresponds to func.dataset.prototyptyping.R

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#         preppredictors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preppredictors<-function(ecospec=wspe1, predictor=predictor, rel=relevant.cols,KILL=TRUE,ISOLATE=FALSE){
	Ap<-predictor
	COL<-match(rel,names(ecospec))
	for(i in COL){
		par<-names(ecospec)[i]
		cat(paste(par,"\n"))
		maximum<-repmax<-max(ecospec[[par]])
		minimum<-repmin<-min(ecospec[[par]])
		if(KILL){repmax<-repmin<-NA} #replace values out of paramater space with NA or min/max
		try(Ap@data[which(Ap@data[,par]> maximum),par]<-repmax)
		try(Ap@data[which(Ap@data[,par]< minimum),par]<-repmin)			
	}
	if(ISOLATE){
		return(Ap[seq(1,ncol(Ap),1)%in%match(rel,names(Ap))])
	}else{
		return(Ap)
	}
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# eliminate ambiguous stations 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elim.ambi.stations<-function(data,species,species.genus.list,verbose=FALSE,invert=FALSE)
{
	#e.g. species.genus.list=list(c("AbraN","Abra"),c("NereisD","Nereis"),c("AmphiuraF","Amphiura"))
	index<- seq_along(species.genus.list)[sapply(species.genus.list, FUN=function(X) species %in% X[1])]
	if(length(index)<1)
	{
		if (verbose) cat("elim.ambi.stations: \t","nothing to do \n")
		return(data)
	}
	genus <-species.genus.list[[index]][2]
	r<-which(!is.na(data[[genus]]) & is.na(data[[species]]))
	if (verbose) cat("elim.ambi.stations: \t",length(r),"ambiguous occurrences ")
	if(length(r)<1)
	{
		data
	} else {
		if(invert){
			if (verbose) cat("found\n")
			data[r,]
		} else {
			if (verbose) cat("removed\n")
			data[-r,]
		}
	}
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# aggregate.stations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aggregate.stations <- function(data,raster,species,xtracols=NULL)
{
	stopifnot(class(data)=="SpatialPointsDataFrame")
	stopifnot(class(raster)=="RasterLayer")
	stopifnot(class(species)=="character")
	stopifnot( class(xtracols) %in% c("NULL","character") )
	rid<-cellFromXY(raster,data)
	xynames <- colnames(slot(data,"coords"))
	df<-data.frame(data)
	names <- c(xynames, xtracols, species)
	tmp <- by(df[,names], rid,function(x) colMeans(x,na.rm=TRUE))
	df<-data.frame(do.call(rbind,tmp)) 
	SpatialPointsDataFrame(df[,xynames],data=df[,c(xtracols, species),drop=FALSE], proj4string=CRS(proj4string(data)))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# returns (sp)df of response variables 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
species.response <- function(species,response,ws=ws,ds=ds,species.cols,name.species=TRUE, verbose=FALSE)
{
	#uses isolate.species(),na2zero()
	
	# ws/ds (sp)df with mass/counts
	stopifnot( any(response %in% c("mass" , "count" , "prevalence", "presence")) ) 
	stopifnot( any(class(species.cols) %in% c("integer" , "character")) )
	# stopifnot(nrow(ws)==nrow(ds))
	
	# declarations:
	ZERO <- 1e-7
	
	tmp <- ds[,-species.cols,drop=FALSE]
	tmpcols <- names(tmp)
	if( any(response %in% c("count" , "prevalence", "presence")) )
	{
		d <- isolate.species(species=species,data=ds)[[species]]
		if("count"%in%response) 
		{
			tmp$count <- na2zero(d)
		}
		if("prevalence"%in%response)
		{
			tmp$prevalence <- sapply(d,function(x){ifelse(is.na(x),0,1)})
		}
		if("presence"%in%response)
		{
			tmp$presence <- sapply(d,function(x){ifelse(is.na(x),0,1)})
		}
	}	
	if("mass"%in%response)
	{
		tmp$mass <- isolate.species(species=species,data=ws)[[species]]
		# do the correction for 0 meaning unknown for mass
		if (verbose) cat("species.response: \t\t removing",length(which(tmp$mass<ZERO)),"occurrencees due to zeroness\n")
		tmp <- tmp[which( is.na(tmp$mass) | tmp$mass > ZERO ),]
		tmp$mass <- na2zero(tmp$mass)
		
		
	}
	if(name.species && length(response)==1) names(tmp)[ncol(tmp)]<-species
	
	return(tmp[,c(tmpcols, response)])
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# sample background data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title sample background aka pseudo-absence data
#' 
#' @description takes \code{SpatialPoints} or \code{SpatialPointsDataFrame} and add  values from  
#' FIXME 
#' @param presence \code{SpatialPoints} of presence data 
#' @param raster \code{RasterLayer} is required for sampling. If \code{buffer=NULL}, the background is sampled without the pixels with presence data. 
#' @param area \code{SpatialPolygons} or  \code{SpatialPolygonsDataFrame} with an area within \code{raster} to be sampled (cmp. \code{extent})
#' @param sample.size The sample size
#' @param extent Eiter \code{raster(default), extent, or area}. Sample the whole area (extent), only where raster vals are !=NA (rater) or confined to an area (area).
#' @param buffer Perimeter around presence data where background is not sampled. Distance unit of buffer: see \code{\link{raster::buffer}}.
#' @param verbose Be verbose
#' @param seed Seed for random background sampling 
#' @param xtracols transfer given features from \code{presence}--input data \code{SpatialPointsDataFrame} to returned \code{SpatialPointsDataFrame}.
#' @param out If \code{spdf}(default), a SpatialPointDataFrame is returned with a column "presence" with 1:presence 0:absence. If \code{sp}, only the \code{SpatialPoints} of the background are returned.  Overrides \code{xtracols}
#' @param colname New colnames for \code{xtracols}. \code{length(colname)==length(xtracols)}.
#' @return \code{SpatialPointsDataFrame} or \code{SpatialPoints}. See \code{out}.
#' @details If \code{FALSE}, preexisting columns in \code{points} will be overwritten by column creates by looking up values in raster layers with the same name.
#' #FIXME
#' @examples #FIXME
#' 
pseudo.absence<-function(presence,raster,area=NULL,sample.size=1000,extent="raster",buffer=NULL,chunks=NULL,verbose=TRUE,seed=NULL, xtracols=NULL,out="spdf",colname=NULL)
{
	# uses ordinal()
	library(raster)
	# distance unit of buffer: see ?raster::buffer()
	# presence, raster, and area must have same crs!
	# out="spdf" -> SpatialPointsDataFrame incl. presence (field "presence" with 1:presence 0:absence)
	# out="sp" -> SpatialPoints of pseudo.absence data
	if(!is.null(presence))
	{
		stopifnot( any( class(presence)=="SpatialPoints" , class(presence)=="SpatialPointsDataFrame" ))
		stopifnot( proj4string(raster)==proj4string(presence) )
		if(extent!="raster") {stopifnot(  proj4string(area)==proj4string(presence) ) }
	}
	stopifnot( any( class(sample.size)=="numeric" , class(sample.size)=="integer" ))
	stopifnot( any( extent=="raster" , extent=="extent",extent=="area" ))
	stopifnot( any( extent=="raster" , class(area)=="SpatialPolygons" , class(area)=="SpatialPolygonsDataFrame" ))
	stopifnot( any( is.null(buffer) , class(buffer)=="numeric" ))
	stopifnot( any( is.null(chunks) , class(chunks)=="numeric" ))
	stopifnot( any( is.null(colname) , class(colname)=="character" ))
	stopifnot( any( is.null(xtracols) , all(xtracols%in%names(presence)) ))
	sample.size <- floor(sample.size)
	
	raster.x <- raster
	raster.x[!is.na(raster.x)] <- 1   # H function
	raster.na <- setValues(raster,rep(NA,length(raster)))
	
	# mark areas to exclude from sampling and maybe buffering
	rid <- cellFromXY(raster.na, presence)
	rid <- rid[!duplicated(rid)]
	
	if(is.null(presence))
	{
		rbuff <- raster.na
	} else	{
		if(is.null(buffer))
		{
			raster.na[rid] <- 1
			rbuff <- raster.na
		} else {
			if(verbose) cat(paste("buffering",length(rid),"points"))
			if(is.null(chunks))
			{
				cat("\n")
				raster.na[rid] <- 1
				rbuff <- buffer(raster.na,width=buffer)	
			} 
			else 
			{
				rbuff<-raster.na
				s<-seq(1,length(rid),chunks)
				if(s[length(s)]!=length(rid)){s<-c(s,length(rid))}
				s[length(s)]<-s[length(s)]+1
				cat(" in chunks of",chunks,"\n")
				for( i in 1:(length(s)-1))
				{
					cat(paste0(" ... ",ordinal(i)," round ",s[i],"-",s[i+1]-1,"\n"))
					raster.na <- setValues(raster,rep(NA,length(raster)))
					raster.na[rid[s[i]:(s[i+1]-1)]] <- 1
					raster.na <- buffer(raster.na,width=buffer)
					rbuff<-cover(rbuff,raster.na)
				}
			}
		}
	}
	rbuff[is.na(rbuff)] <- 2
	
	# crop to extent
	switch(extent,
		   raster={ 	bucket <- rbuff
		   		  bucket[is.na(raster)] <- NA
		   },
		   extent={	ref <- crop(raster,extent(area))
		   		 bucket <- crop(rbuff,extent(area))
		   		 bucket[is.na(ref)] <- NA
		   },
		   area={		ma <- mask(rbuff,area)
		   		ma[is.na(raster)] <- NA
		   		bucket <- crop(ma,extent(area))
		   })
	bucket[bucket<2] <- NA
	
	#sampling
	absence<-as(bucket,"SpatialPoints")
	if(!is.null(seed)) { set.seed(seed) }
	abs<-absence[sample(length(absence),sample.size),]
	if(verbose)	{ plot(bucket,legend=FALSE) 
				  plot(abs,cex=0.5,pch=4,add=T) 
				  if(!is.null(presence)){
				  	plot(presence,cex=0.3,pch=20,col="red",add=T)
				  }
	}
	
	# output format
	if(out=="spdf")
	{
		q<-rbind( as(presence,"SpatialPoints"),abs)
		row.names(q)<-1:length(row.names(q))
		if(!is.null(xtracols))
		{	
			P<-data.frame(presence)[,xtracols,drop=FALSE]
			A<-data.frame(matrix(data=rep(0,length(abs)*ncol(P)),ncol=ncol(P)))
			names(A)<-names(P)
			abs<-SpatialPointsDataFrame( q,data=rbind(P, A))
		} else {
			abs<-SpatialPointsDataFrame( q,data=data.frame(presence=c(rep(0,length(abs)),rep(1,nrow(presence)))))
			if(!is.null(colname)) { names(abs)<-colname }
		}
	} 
	invisible(abs)
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# add predictors to your spatial points
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' @title read predictors and add respective values to your spatial points
#' 
#' @description takes your response as \code{SpatialPoints} or \code{SpatialPointsDataFrame} and add  values from  
#' a \code{RasterLayer}, \code{RasterStack}, or \code{RasterBrick} at those positions. New columns will have the names of the respective raster layers.
#' @param points \code{SpatialPoints} or \code{SpatialPointsDataFrame} 
#' @param environ \code{RasterLayer}, \code{RasterStack}, or \code{RasterBrick}
#' @param dropcols. if \code{TRUE}, remove preexisting columns in \code{points} that also exist in \code{environ} under the same name.
#' @return \code{SpatialPointsDataFrame}
#' @details If \code{FALSE}, preexisting columns in \code{points} will be overwritten by column creates by looking up values in raster layers with the same name.
#' #FIXME
#' @examples #FIXME
#' 
get.environ<-function(points, environ,dropcols=TRUE)
{
	library(raster)
	stopifnot( class(environ) %in% c("RasterLayer","RasterStack","RasterBrick") )
	stopifnot(class(points) %in% c("SpatialPoints","SpatialPointsDataFrame"))
	if(class(dropcols)=="logical") 
	{
		if(dropcols && class(points)=="SpatialPointsDataFrame")
		{
			dropcols<-names(points)
		} 
	} 
	
	rip<-cellFromXY(environ , points)
	df<-data.frame(environ[rip])[!(names(environ)%in%dropcols)] # cols in dropcols 
	if(class(points)=="SpatialPoints"){
		SpatialPointsDataFrame(points,data=df)	
	} else {
		#SpatialPointsDataFrame(points,data=cbind(data.frame(points)[1],data.frame(environ[rip])))
		SpatialPointsDataFrame(points,data=cbind(data.frame(points,drop=FALSE)[names(points)],df))
	}	
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# outlier elimination (only at roof)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elim.outliers<-function(data,species=NULL,na.rm=FALSE,threshold=NULL,verbose=TRUE)
{ # try threshold=1e-7
	library(sp)
	stopifnot( class(data) %in% c("data.frame","SpatialPointsDataFrame", "SpatialGridDataFrame") ) #sgdf makes not much sense anyway‚
	stopifnot( class(species) %in% c("NULL","character") )
	stopifnot( class(threshold) %in% c("NULL","numeric") )
	stopifnot( class(na.rm)=="logical" )
	
	
	if( is.null(species) ) # last column is species
	{	
		species <- names(data)[ncol(data)]
	} else stopifnot(species%in%names(data))
	
	string<-NULL
	if( NA %in% data[[species]] )
	{
		if( na.rm ) 
		{
			string<-"removing NA's."
			data<-data[!is.na(data[[species]]),]
		} else {
			string<-"converting NA to 0."
			data[is.na(data[[species]]),species]<-0
		}
	}
	
	if( is.null(threshold) ) base<-data[[species]] 
	else base <- data[[species]][which( data[[species]]>threshold)]
	sy<-summary(base)
	iqr<-sy[5]-sy[2]; 		# inter-quartile range
	ant<-sy[5]+1.5*iqr; 	# antenna
	r<-which( data[[species]]<ant)
	if(verbose) cat("elim.outliers: \t\t\t",string,nrow(data)-length(r),"outliers\n")
	data[ r,]
}


