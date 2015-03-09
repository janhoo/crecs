moranii <- function(sp,moParams=c("response","residuals"),moRange=list(c(0,0.5),c(0.5,1)))
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
	
	library(spdep)
	library(rgdal)
	stopifnot(class(sp)=="SpatialPointsDataFrame")
	result<-NULL
	p.adj.method="holm"
	

	
	
	for(a in moRange){
		nb <- dnearneigh(sp, a[1], a[2])
		
		for(b in moParams){
			m <- sp.correlogram(nb, as.vector(sp[[b]]), order=2, method="I", zero.policy=TRUE)
			
			res <- as.matrix(m$res)
			ZI <- (res[, 1] - res[, 2])/sqrt(res[, 3])
			pv <- p.adjust(2 * pnorm(abs(ZI), lower.tail = FALSE), method = p.adj.method)			
			mo<-res[1,1]
			p<-pv[1]/2
			
			result <- c(result, mo)
			names(result)[length(result)] <- paste("M",b,paste(a,collapse="-"),sep="*")
			result <- c(result, p)
			names(result)[length(result)] <- paste("p",b,paste(a,collapse="-"),sep="*")
		}
	}
    # if(any(sapply(result,is.na))){stop("NA in metrics")}
    return(result)
}


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
