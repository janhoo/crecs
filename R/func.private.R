
sd.raster<-function(N,K,tmpfiles,rmean)
{
	for(i in 1:N){
		load(file=tmpfiles[i])
		if(i==1)
		{
			temp<-sum((Ra-rmean)^2)
		} else {
			temp<-temp+sum((Ra-rmean)^2)
		}
	}
	return(sqrt(temp/((N*K)-1)))
}


aggregate.raster<-function(N,K,tmpdir,tmpstack,verbose=TRUE)
{
	load(file=file.path(mapdir,fn("rmean")))
	for(i in 1:N){
		load(file=file.path(tmpdir,paste0(tmpstack,i,".bin")))
		if(i==1)
		{
			temp<-sum((Ra-rmean)^2)
		} else {
			temp<-temp+sum((Ra-rmean)^2)
		}
		if(verbose) catline(string=string,add=paste0("std.dev ",i,"/",N),sep=" ")
	}
	rsd<-sqrt(temp/((N*K)-1))
	save(rsd,file=file.path(mapdir,fn("rsd")))
	return("|")
}



#aggregate.results<-function(res)
#{
#	cat<-sapply(res, is.numeric)
#	num<-(1:ncol(res))[cat]
#	fac<-(1:ncol(res))[!cat]
#	fac.names<-names(res)[fac]
#	
#	MEAN<-by(res[,num],INDICES=res[,fac.names],function(x) colMeans(x,na.rm=TRUE))
#	mat1<-do.call(rbind,MEAN) 
#	SD<-by(res[,num],INDICES=res[,fac.names],function(x) sapply(x,FUN=function(y){sd(y,na.rm=T)}))
#	mat2<-do.call(rbind,SD) 
#	colnames(mat2)<-paste0(colnames(mat2),".sd")
#	ex<-expand.grid(attributes(MEAN)$dimnames)
#	#names(ex)<-fac.names
#	cbind(ex,mat1,mat2)
#}

aggregate.results<-function(res)
{
	cat<-sapply(res, is.numeric)
	MEAN<-by(res[cat],INDICES=res[!cat],function(x) colMeans(x,na.rm=TRUE))
	mat1<-do.call(rbind,MEAN) 
	SD<-by(res[cat],INDICES=res[!cat],function(x) sapply(x,FUN=function(y){sd(y,na.rm=T)}))
	mat2<-do.call(rbind,SD) 
	colnames(mat2)<-paste0(colnames(mat2),".sd")
	ex<-expand.grid(attributes(MEAN)$dimnames,stringsAsFactors = FALSE)
	N<-nrow(res)
	cbind(N,ex,mat1,mat2)
}



