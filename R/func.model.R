#' @title Make a (species distribution) regression or classification model
#' 
#' @description Returns the model, the performance metrics, and/or distribution maps depending on arguments.
#' 
#' @param data \code{SpatialPointsDataFrame} containing response and predictors
#' @param method SDM methos used: "gam", "rf", "gbm", "max", or "gbm.step". See details for details
#' @param responsetype Type of the response: "count", "continuous", or "presence". Presence denotes bimodal responses [0|1]
#' @param response Column name of response in data argument
#' @param predictors Column names of predictors to use in data argument
#' @param secondary Column name in data. Use this to calculate model performance metrics from instead of response. Default is NULL.
#' @param enviStack \code{RasterStack} of predictors. Used to calculate SD map. See details.
#' @param enviPix \code{SpatialPixelsDataFrame} of predictors. enviPix<-as(enviStack,"SpatialPixelsDataFrame"). Only for performance. See details.
#' @param seed Integer. For reproduceability.
#' @param aggregated logical or \code{NULL} Default is  \code{NULL}. This is for ensemble calculations and added to metrics in case of not being \code{NULL}
#' @param pseudoabsence logical or \code{NULL} Default is \code{NULL}. Same as above
#' @param gbm.trees gbm.trees param for dismo::gbm
#' @param maxargs argument to pass tp maxent
#' @param model \code{logical} return ONLY the model. Default is \code{FALSE}. Overrides args below.
#' @param prediction return ONLY the prediction. Default is \code{FALSE}. Overrides args below.
#' @param moran add spatial autocorrelation metric to model metric. Default is \code{FALSE}. See details.
#' @param flat return model performance metrics only (as \code{data.frame})
#' @param rast return \code{list} of metric and raster
#' @param rast.occ additionally adds bimodal occurence map. see details for threshold calculation.
#' @param ... ellipsis is used to pass arguments to subsequent functions like  \code{threshold.def} and \code{moRange}. See \code{\link{metrics}} \code{\link{moranii}} for details
#' @return A \code{model}, or a \code{data.frame} or a \code{list} depending on arguments
#' @details SDM methos used are "gam", "rf", "gbm", "max", or "gbm.step".
#' If you want spatial autocorrelation metrics you probably need to pass additional arguments to \code{\link{moranii}}. 
#' Note that calculations may take very long depending on number of points and parametrization.
#'
#' \code{enviStack} and \code{enviPix} are the same in different data types. 
#' It is sufficient if you supply only one--the other will be generated. 
#' For big data sets and repetitive tasks it may be worthwhile to pass both to increase performance
#' @examples 
#' data<-get.environ(species,deutschebucht)
#' mo<-model(data=data,method="rf",responsetype = "continuous", response = "species1", predictors = c("mgs","mud","depth"),enviStack = deutschebucht, model=TRUE)
#' qqplot(data$species1,mo$predicted)
#' 
#' mo<-model(data=data,method="rf",responsetype = "continuous", response = "species1", predictors = c("mgs","mud","depth"),enviStack = deutschebucht, rast=TRUE)
#' par(mfrow=c(1,2))
#' plot(mo$raster$full,main="prediction")
#' plot(raster(deutschebucht,layer=match("species1",names(deutschebucht))),main="true distribution")
#' 
model <- function(data, 
                  method,             
                  responsetype,      
                  response, 
                  predictors, 
                  secondary=NULL, 
                  enviStack=NULL, 
                  enviPix=NULL, 
                  seed=NULL, 
				  aggregated=NULL, 
                  pseudoabsence=NULL, 
                  gbm.trees=2000, 
				  maxargs=c("outputformat=logistic", "defaultprevalence=0.5"),
				  model=FALSE,
				  prediction=FALSE,
				  flat=FALSE,
				  moran=FALSE,
                  rast=TRUE,
				  rast.occ=TRUE,
                  ...)
{
	# Metrics and prediction map of Species Distribution Model specified as list
	#
	# Args:
	#  data             [<SpatialPointsDataFrame>] 
	#  method           [<character>] modeling method
	#  responsetype     [<character>] 
	#  response         [<character>] column name of data 
	#  predictors       [<character>] column names of data
	#  secondary        [<character>] column name of data; compare with prediction for cor() in metrics()
	#  enviStack        [<RasterStack>]
	#  enviPix          [<SpatialPixelsDataFrame>]
	#  seed             [<numeric>] seed
	#  aggregated       [<logical>] is aggregated
	#  pseudoabsence    [<logical>] absences are pseudoabsences
	#  gbm.trees        [<numeric>] 
	#  OUTPUT OPTIONS   
	#  model            [<logical>] output model only (overides options below)
	#  prediction		[<logical>] output prediction only (overides options below)
	#  moran            [<logical>] add spatial autocorrelation metric to model metric
	#  flat             [<logical>] output model metric as data.frame (overides options below)
	#  rast             [<logical>] output list of metric and raster
	#  rast.occ			[<logical>] raster becomes list of of rasters of model and absence/presence prediction
	#  ...              arguments for moranii() & makefn.free.model (moParams=c("response","residuals"),moRange=list(c(0,0.5),c(0.5,1)),check.names=TRUE)
	
	#dots <- list(...)
	stopifnot( !all(    c(is.null(enviStack), is.null(enviPix)) ))
	if(is.null(enviStack)){
		enviStack <- stack(enviPix)
	}
	if(is.null(enviPix)){
		enviPix <- as(enviStack,"SpatialPixelsDataFrame")
	}	
	if(is.null(seed)){
		seed <- as.integer(runif(1)*100000)
	}
	if(!is.null(secondary)){
		secondary <- data[[secondary]]
	}
	Model <- makefn.free.model(method=method, response=response, responsetype=responsetype, predictors=predictors, Predictor=enviPix, data.names=names(data), gbm.trees=gbm.trees, ...)
	Predict <-makefn.free.predict(responsetype=responsetype, method=method, environ=enviStack, maxargs=maxargs)
	set.seed(seed)
	fit <- Model(data=data)
	if(model){
		return(fit)
	}
	P <- Predict(fit,data)
	if(prediction){
		data$prediction<-P
		data$residuals<-data[[response]]-P
		return(data)
	}
	M <- metrics(response=data[[response]], prediction=P, secondary=secondary, ...)
	
	if(moran){
		dat <- data[,response]
		dat$response <- data[[response]]
		dat$residuals <- data[[response]]-P
		M <- c(M,moranii(dat, ...))
	}
	metric <- data.frame(list(
		as.list(c(response=response, method=method, responsetype=responsetype, predictors=paste0(predictors,collapse="."))),
		as.list(c(aggregated=ifelse(class(aggregated)=="logical",aggregated,FALSE),pseudoabsence=ifelse(class(pseudoabsence)=="logical",pseudoabsence,FALSE))),
		as.list(M)),stringsAsFactors=FALSE)
	if(flat){
		return(metric)
	}
	if(rast){
		enviPix$full <- Predict(fit,enviPix)
		full <- raster(enviPix,match("full",names(enviPix)))
		if(rast.occ){
			t<-M[["threshold"]]
			ap <- full
			ap@data@values<-presence.absence(ap@data@values,t=t)
			full <- list(full=full,ap=ap)
		}
	} else {
		full <- NULL
	}
	return(list(metric=metric,raster=full))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                     
# crossvalid :: cross validation routine 
#
#' @title Cross validation routine for (species distribution) regression or classification models
#' 
#' @description repeated k-fold cross-validation. Calculate model preformance metrics, disribution, standard deviation and occurence maps.
#' 
#' @param Ncv Number of cv repititions.
#' @param Kfold Number of folds. Usually take 10 or 5.
#' @param data \code{SpatialPointsDataFrame} containing response and predictors
#' @param method SDM methos used: "gam", "rf", "gbm", "max", or "gbm.step". See details for details
#' @param responsetype Type of the response: "count", "continuous", or "presence". Presence denotes bimodal responses [0|1]
#' @param response Column name of response in data argument
#' @param predictors Column names of predictors to use in data argument
#' @param secondary Column name in data. Use this to calculate model performance metrics from instead of response. Default is NULL.
#' @param strata criterion for fold generation. #FIXME
#' @param enviStack \code{RasterStack} of predictors. Used to calculate SD map
#' @param enviPix \code{SpatialPixelsDataFrame} of predictors. enviPix<-as(enviStack,"SpatialPixelsDataFrame"). Only for performance.
#' @param seed Integer. You probably want reproduceability. Note that Maxent's pseudoabsence generation can't be seeded this way--so expect those results to vary
#' @param aggregated logical or \code{NULL} Default is  \code{NULL}. This is for ensemble calculations and added to metrics in case of not being \code{NULL}
#' @param pseudoabsence logical or \code{NULL} Default is \code{NULL}. Same as above
#' @param gbm.trees gbm.trees param for dismo::gbm
#' @param maxargs argument to pass tp maxent
#' @param flat return model performance metrics only (as \code{data.frame})
#' @param rast return \code{list} of metric and raster
#' @param ... ellipsis is used to pass arguments to subsequent functions like  \code{threshold.def} (see \code{\link{metrics}}) or \code{moranii()} & \code{makefn.free.model()} (moParams=c("response","residuals"),moRange=list(c(0,0.5),c(0.5,1)),check.names=TRUE)
#' 
#' @examples 
#' data<-get.environ(species,deutschebucht)
#' cv<-crossvalid(Ncv=1,Kfold=5,data=data,method="rf",responsetype = "continuous", response = "species1", predictors = c("mgs","mud","depth"),enviStack = deutschebucht,seed=23, check.names = FALSE)
#'  
#' # aggregated results
#' cv$metric.agg
#' 
#' # plot result maps
#' par(mfrow=c(2,2))
#' plot(cv$rmean,main="cv prediction")
#' plot(cv$sd,main="standart deviation")
#' plot(raster(deutschebucht,layer=match("species1",names(deutschebucht))),main="true distribution")
#' plot(cv$rbin,main="prob. of occurrence")
#' 
crossvalid <- function(Ncv, 
                       Kfold, 
                       data, 
                       method, 
                       responsetype, 
                       response, 
                       predictors, 
                       secondary=NULL, 
                       strata=NULL, 
                       enviStack=NULL,
                       enviPix=NULL,
                       seed,
					   aggregated=NULL, 
					   pseudoabsence=NULL,
                       gbm.trees=2000,
					   maxargs=c("outputformat=logistic", "defaultprevalence=0.5"),
                       rast=TRUE,
                       flat=FALSE,
					   ...)
{
	# Metrics and prediction map of Species Distribution Model specified as list
	#
	# Args:
	#  Ncv              [<numeric>] Number of cross-validation runs
	#  Kfold            [<numeric>] Number of folds
	#  data             [<SpatialPointsDataFrame>] 
	#  method           [<character>] modeling method
	#  responsetype     [<character>] 
	#  response         [<character>] column name of data 
	#  predictors       [<character>] column names of data
	#  secondary        [<character>] column name of data; compare with prediction for cor() in metrics()
	#  strata           [<character>] column name of data; stratum used to create folds: a vector or factor with sub-groups
	#  enviStack        [<RasterStack>]
	#  enviPix          [<SpatialPixelsDataFrame>]
	#  seed             [<numeric>] seed
	#  aggregated       [<logical>] is aggregated
	#  pseudoabsence    [<logical>] absences are pseudoabsences
	#  gbm.trees        [<numeric>] 
	#  rast             [<logical>] output raster
	#  flat             [<logical>] no list, just metric as data.frame
#  ...              arguments for makefn.free.model
	library(dismo)
	stopifnot( !all(    c(is.null(enviStack), is.null(enviPix)) ))
	if(is.null(enviStack)){
		enviStack <- stack(enviPix)
	}
	if(is.null(enviPix)){
		enviPix <- as(enviStack,"SpatialPixelsDataFrame")
	}	
	sec<-NULL
	if(flat){
		rast<-FALSE
	}
	if(is.null(strata)){
		strata<-"presence"
		if("presence" %in% names(data)){
			stop("presence exists in data - use strata=\"presence\"")
		} else {
			data$presence<-presence.absence(data[[response]])
		}
	}
	if(rast){
		rbin<-raster(enviStack,layer=1)*0
	}
	tmpfiles<-vector(mode = "character",length = Ncv)
	
	set.seed(seed)
	seed<-as.integer(runif(Ncv,min=10000,max=99999)*1000)
	Model<-makefn.free.model(method=method, response=response, responsetype=responsetype, predictors=predictors, Predictor=enviPix, data.names=names(data), gbm.trees=gbm.trees, ...)
	Predict<-makefn.free.predict(responsetype=responsetype, method=method, environ=enviStack, maxargs=maxargs)
	for(n in 1:Ncv){
		set.seed(seed[n])
		group <- kfold(x=data,k=Kfold,by=data[[strata]]) 
		for(k in 1:Kfold){
			#cat(n," ",k,"\n")
			test  <- data[group == k,]
			train <- data[group != k,]
			
			set.seed(seed[n])
			fit<-Model(data=train)
			
			P<-Predict(fit,test)
			if(!is.null(secondary)){
				sec <- data[[secondary]][group == k]
			}
			M<-metrics(response=test[[response]],prediction=P,secondary=sec, ...)
			metric<-data.frame(list(
				as.list(c(response=response,method=method,responsetype=responsetype,predictors=paste0(predictors,collapse="."))),
				as.list(c(aggregated=ifelse(class(aggregated)=="logical",aggregated,FALSE),pseudoabsence=ifelse(class(pseudoabsence)=="logical",pseudoabsence,FALSE))),
				as.list(M)),stringsAsFactors=FALSE)
			
			
			if (n+k<3){ # init
				result<-allocateDf(metric,Ncv*Kfold)
				if(rast){
					Ra<-stack(replicate(Kfold,rbin))
				}
			}
			result[(n-1)*Kfold+k,]<-metric[1,]
			if(rast){
				t<-M[["threshold"]]
				enviPix$fit<-Predict(fit,enviPix)
				ro<-ra<-raster(enviPix,layer=match("fit",names(enviPix)))
				ro@data@values<-presence.absence(ro@data@values,t=t)
				rbin<-rbin+ro   #certainty
				Ra[[k]]<-ra
			}
		}
		if(rast){
			tmpfiles[n]<-tempfile(pattern="tmp",fileext = ".bin")
			save(Ra,file=tmpfiles[n])
			if(n==1){
				rmean<-sum(Ra)/(Ncv*Kfold)
			} else {
				rmean<-rmean+sum(Ra)/(Ncv*Kfold)
			}
		}
	}
	
	
	if(rast){
		# raster
		rbin<-(rbin/(Ncv*Kfold/2))-1
		sd<-sd.raster(N=Ncv,K=Kfold,tmpfiles=tmpfiles,rmean=rmean)
		names(rbin)<-"rbin"
		names(rmean)<-"rmean"	
		names(sd)<-"sd"
	} else {
		rbin <- rmean <- sd <- NULL
	}
	if(flat){
		return(result)
	} else {
		return(list(metric.agg=aggregate.results(result),metric=result,rbin=rbin,rmean=rmean,sd=sd))
	}
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                     
# nill models 
#
# Metrics of null model cross validation 
#
# Args:
#  data             [<SpatialPointsDataFrame>] 
#  method           [<character>] modeling method
#  predictors       [<character>] column names of data
#  enviStack        [<RasterStack>]
#  enviPix          [<SpatialPixelsDataFrame>]
#  Ncv              [<numeric>] Number of cross-validation runs
#  Ncv              [<numeric>] Number of cross-validation runs
#  Kfold            [<numeric>] Number of folds
#  seed             [<numeric>] seed
#  prob             [<numeric> or <character>] if <char> the prevalence of data[[prob]] is used as prob

nullModel <- function (data, 
					   method="rf", 
					   predictors, 
					   enviStack, 
					   enviPix, 
					   N0=100, 
					   Ncv=10, 
					   Kfold=5, 
					   seed,  
					   prob=0.5) {
	set.seed(seed)
	Seed<-as.integer(runif(N0,min=10000,max=99999)*1000)
	if (class(prob)=="character"){
		prob<-sum(presence.absence(data[[prob]]))/nrow(data)
	}
	
	for(i in 1:N0){
		set.seed(Seed[i])
		data$presence<-rbinom(nrow(data),1,prob)
		m <-crossvalid(Ncv=Ncv,Kfold=Kfold,data=data,method=method,responsetype="presence",response="presence",predictors=predictors,strata="presence",secondary="presence",enviStack=enviStack,enviPix=enviPix,seed=Seed[i],maxargs=c("outputformat=logistic", "defaultprevalence=0.5"),gbm.trees=2000,rast=FALSE,flat=TRUE)
		cat <- c(rep(FALSE,(i-1)*Ncv*Kfold),rep(TRUE,Ncv*Kfold),rep(FALSE,Ncv*Kfold*(N0-i)))
		if (i<2){ # init
			nullmodelmetric<-allocateDf(m,Ncv*Kfold*N0)
		}
		nullmodelmetric[cat,]<-m
	}
	return(aggregate.results(nullmodelmetric))
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                                     


moranii <- function(sp,moParams=c("response","residuals"),moRange=list(c(0,0.5),c(0.5,1)), ...)
{
	# moran() calculates spatial autocorrelation (MoranII) of response and residuals for given ranges
	# only for lonlat projs. FIXME
	# Args:
	#    sp: spatialPoints or SPDF 
	#    response,residuals [<numeric>] input for spdep::sp.correlogram
	#    range [<list>] range as for spdep::dnearneigh
	#
	# Uses:
	#   clean.null.list()
	#browser()
	library(spdep)
	library(rgdal)
	stopifnot(class(sp)=="SpatialPointsDataFrame")
	result<-NULL
	p.adj.method="holm"
	spcor<-function(){
		m <- sp.correlogram(nb, as.vector(sp[[P]]), order=2, method="I", zero.policy=TRUE)
		res <- as.matrix(m$res)
		ZI <- (res[, 1] - res[, 2])/sqrt(res[, 3])
		pv <- p.adjust(2 * pnorm(abs(ZI), lower.tail = FALSE), method = p.adj.method)			
		mo<-res[1,1]
		p<-pv[1]/2
		return(c(mo,p))
	}
	
  
  
  for(R in moRange){
    nb <- dnearneigh(sp, R[1], R[2])
    
    for(P in moParams){
      m <- tryCatch(spcor(),error=function(x){return(c(NA,NA))})

      
      result <- c(result, m[1])
      names(result)[length(result)] <- paste("M",P,paste(R,collapse="_"),sep=".")
      result <- c(result, m[2])
      names(result)[length(result)] <- paste("p",P,paste(R,collapse="_"),sep=".")
    }
  }
  # if(any(sapply(result,is.na))){stop("NA in metrics")}
  return(result)
}



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#   
# code testing 
if(FALSE){
  nb<-dnearneigh(data, 0, 50)
  sp.correlogram(nb, as.vector(data$density), order=2, method="I")
  
  
  print(subs.correlog, p.adj.method="holm")
  #quartz(title="Correlogram of substrate density")
  plot(subs.correlog)
  
  # check catchment area of nb
  nn<-sapply(nb,FUN=function(x){length(x)}) # number of neighbors
  n<-match(max(nn),nn) # point with most neighbors
  load(file=file.path(WD,"data.geo/v.l.europeS.shape.bin"))
  plot(data,pch=20,cex=0.2)
  plot(data[n,],pch=20,cex=1,col="green",add=T)
  plot(data[nb[n][[1]],],pch=20,cex=0.2,col="red",add=T)
  plot(europeS.shape,add=T)
  
}

