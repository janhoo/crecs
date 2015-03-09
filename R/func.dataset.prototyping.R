#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Contents
#
# prototyping functions :
# utilized by prepare.dataset() 
# plus for
# elim.ambi.stations() & preppredictors() & fn()
#
# 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# corresponds to func.dataset.R


makefn.fn<-function(species,dataset,model,predictors)
{
	fn<-function(string)
	{
		return(paste(species,string,dataset,model,paste0(predictors,collapse="."),"bin",sep="."))
	}
}

#         preppredictors
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.preppredictors<-function(spdf, predictors,KILL=FALSE,ISOLATE=FALSE,verbose=FALSE)
{
	preppredictors<-function(data){
		library(sp)
		COL<-match(predictors,names(data))
		for(i in COL){
			par<-names(data)[i]
			if (verbose) cat(paste(par,"\n"))
			maximum<-repmax<-max(data[[par]])
			minimum<-repmin<-min(data[[par]])
			if(KILL){repmax<-repmin<-NA} #replace values out of paramater space with NA or min/max
			try(spdf@data[which(spdf@data[,par]> maximum),par]<-repmax)
			try(spdf@data[which(spdf@data[,par]< minimum),par]<-repmin)      
		}
		if(ISOLATE){
			return(spdf[seq(1,ncol(spdf),1)%in%match(predictors,names(spdf))])
		}else{
			return(spdf)
		}
	}
	return(preppredictors)
}

# eliminate ambiguous stations 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.elim.ambi.stations<-function(species.genus.list,verbose=FALSE,invert=FALSE)
{
	elim.ambi.stations<-function(data,species)
	{
		library(sp)
		#e.g. species.genus.list=list(c("AbraN","Abra"),c("NereisD","Nereis"),c("AmphiuraF","Amphiura"))
		index<- seq_along(species.genus.list)[sapply(species.genus.list, FUN=function(X) species %in% X[1])]
		if (verbose) cat("elim.ambi.stations:     ")
		if(length(index)<1)
		{
			if (verbose) cat("nothing to do \n")
			return(data)
		}
		genus <-species.genus.list[[index]][2]
		r<-which(!is.na(data[[genus]]) & is.na(data[[species]]))
		if (verbose) cat(length(r),"ambiguous occurrences ")
		if(length(r)<1)
		{
			if (verbose) cat("\n")
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
	return(elim.ambi.stations)
}


# aggregate.stations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.aggregate.stations <- function(raster,species,xtracols=NULL)
{
	aggregate.stations <- function(data)
	{
		library(sp)
		library(raster)
		stopifnot(class(data)=="SpatialPointsDataFrame")
		stopifnot(class(raster)=="RasterLayer")
		stopifnot(class(species)=="character")
		stopifnot( any(class(xtracols) %in% c("NULL","character")) )
		
		rid<-cellFromXY(raster,data)
		xynames <- colnames(slot(data,"coords"))
		df<-data.frame(data)
		names <- c(xynames, xtracols, species)
		tmp <- by(df[,names], rid,function(x) colMeans(x,na.rm=TRUE))
		df<-data.frame(do.call(rbind,tmp)) 
		SpatialPointsDataFrame(df[,xynames],data=df[,c(xtracols, species),drop=FALSE], proj4string=CRS(proj4string(data)))
	}
	return(aggregate.stations)
}



# returns (sp)df of response variables 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.species.response<-function(response, species.cols,name.species,verbose)
{
	species.response <- function(species, WS, DS)
	{
		
		#uses isolate.species()
		# WS/DS (sp)df with mass/counts
		library(sp)
		stopifnot( any(response %in% c("mass" , "count" , "prevalence", "presence")) ) 
		stopifnot( any(class(species.cols) %in% c("integer" , "character")) )
		stopifnot(nrow(WS)==nrow(DS))
		
		# declarations:
		ZERO <- 1e-7
		
		tmp <- DS[,-species.cols,drop=FALSE]
		tmpcols <- names(tmp)
		if( any(response %in% c("count" , "prevalence")) )
		{
			d <- isolate.species(species=species,data=DS)[[species]]
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
			tmp$mass <- isolate.species(species=species,data=WS)[[species]]
			# do the correction for 0 meaning unknown for mass
			if (verbose) cat("species.response:       removing",length(which(tmp$mass<ZERO)),"occurrencees due to zeroness\n")
			tmp <- tmp[which( is.na(tmp$mass) | tmp$mass > ZERO ),]
			tmp$mass <- na2zero(tmp$mass)
			
			
		}
		if(name.species && length(response)==1) names(tmp)[ncol(tmp)]<-species
		
		return(tmp[,c(tmpcols, response)])
	}
	return(species.response)
}




# sample background data
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.pseudo.absence<-function(raster,area,extent,buffer,verbose,chunks=NULL,xtracols=NULL,out="sp",colname=NULL)
{
	pseudo.absence<-function(presence,sample.size,seed)
	{
		# uses ordinal()
		library(sp)
		library(raster)
		library(gdata)
		# distance unit of buffer: see raster::buffer()
		# presence, raster, and area must have same crs!
		# out="spdf" -> SpatialPointsDataFrame incl. presence (field "presence" with 1:presence 0:absence)
		# out="sp" -> SpatialPoints of pseudo.absence data
		stopifnot( any( class(presence)=="SpatialPoints" , class(presence)=="SpatialPointsDataFrame" ))
		stopifnot( any( class(sample.size)=="numeric" , class(sample.size)=="integer" ))
		stopifnot( any( extent=="raster" , extent=="extent",extent=="area" ))
		stopifnot( any( extent=="raster" , class(area)=="SpatialPolygons" , class(area)=="SpatialPolygonsDataFrame" ))
		stopifnot( proj4string(raster)==proj4string(presence) )
		if(extent!="raster") stopifnot(  proj4string(area)==proj4string(presence) ) 
		stopifnot( any( is.null(buffer) , class(buffer)=="numeric" ))
		stopifnot( any( is.null(chunks) , class(chunks)=="numeric" ))
		stopifnot( any( is.null(colname) , class(colname)=="character" ))
		stopifnot( any( is.null(xtracols) , all(xtracols%in%names(presence)) ))
		sample.size <- floor(sample.size)
		
		if (verbose) cat("pseudo.absence:        ","drawing",sample.size,"samples\n")
		raster.x <- raster
		raster.x[!is.na(raster.x)] <- 1   # H function
		raster.na <- setValues(raster,rep(NA,length(raster)))
		
		# mark areas to exclude from sampling and maybe buffering
		rid <- cellFromXY(raster.na, presence)
		rid <- rid[!duplicated(rid)]
		if(is.null(buffer))
		{
			raster.na[rid] <- 1
			rbuff <- raster.na
		} else {
			if(verbose) cat(paste0("                        ","buffering ",length(rid)," points"))
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
					if(verbose) cat(paste0("\r                        ",ordinal(i)," chunk: ",s[i],"-",s[i+1]-1,"      "))
					
					raster.na <- setValues(raster,rep(NA,length(raster)))
					raster.na[rid[s[i]:(s[i+1]-1)]] <- 1
					raster.na <- buffer(raster.na,width=buffer)
					rbuff<-cover(rbuff,raster.na)
				}
				if(verbose) cat("\n")
			}
		}
		rbuff[is.na(rbuff)] <- 2
		
		# crop to extent
		switch(extent,
			   raster={   bucket <- rbuff
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
		if(verbose)	{ plot(bucket,legend=FALSE); plot(abs,cex=0.5,pch=4,add=T); plot(presence,cex=0.3,pch=20,col="red",add=T)}
		
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
	return(pseudo.absence)
}

# add predictors to your spatial points
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.get.environ<-function(environ,dropcols=TRUE)
{
	get.environ<-function(points)
	{
		library(sp)
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
		df<-data.frame(environ[rip])[!(names(environ)%in%dropcols)]
		if(class(points)=="SpatialPoints"){
			SpatialPointsDataFrame(points,data=df)  
		} else {
			#SpatialPointsDataFrame(points,data=cbind(data.frame(points)[1],data.frame(environ[rip])))
			SpatialPointsDataFrame(points,data=cbind(data.frame(points,drop=FALSE)[names(points)],df))
		}	
	}
	return(get.environ)
}

# outlier elimination (only at roof)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makefn.elim.outliers<-function(species=NULL,na.rm=FALSE,threshold=NULL,verbose=TRUE)
{
	elim.outliers<-function(data)
	{ 
		# try threshold=1e-7
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
		if(verbose) cat("elim.outliers:         ",string,nrow(data)-length(r),"outliers\n")
		data[ r,]
	}
	return(elim.outliers)
}

