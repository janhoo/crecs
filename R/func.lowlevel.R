#
# Contents
#
# low level funcs:
#
#
allocateDf<- function(sourceDF, length) {
	return(data.frame(lapply(sapply(sourceDF, class), vector, length = length),stringsAsFactors=FALSE))
}



hret<-function(since) #human readable elapsed time
{
	elapsed<-proc.time()[3]-since
	if(elapsed<1) {return(paste(round(elapsed,digits=4),"sec.",collapse=""))}
	if(elapsed<100) {return(paste(round(elapsed,digits=1),"sec.",collapse=""))}
	if(elapsed<60*300) {return(paste(round(elapsed/60,digits=1),"min."))}
	if(elapsed<60*60*48) {return(paste(round(elapsed/3600,digits=1),"hours"))}
	else {return(paste(round(elapsed/86400,digits=1),"days"))}
}


isolate.species<-function(species,data,species.cols=NULL,keep.cols=NULL)
{  
	# in output, species is always in last column 
	# species.cols are removed
	# keep.cols are kept
	m<-match(species,names(data))
	if(!is.null(species.cols)) { 
		keep.cols <- c((1:ncol(data))[-(species.cols)] , keep.cols) } 
	if(!is.null(keep.cols)) { 
		keep.cols<-keep.cols[m!=keep.cols]}
	keep.cols <- c( keep.cols , m )
	data[,keep.cols[!duplicated(keep.cols)]]
}

gear.ssep<-function(data,species,trawlspecies,verbose=TRUE){
	if(trawlspecies)
	{
		
		r <- which(data$gear=="trawl" | data[[species]]>0)
	}else{
		r <- which(data$gear=="grab")
	}
	if(verbose) cat("gear.ssep:             ",length(r),"out of",nrow(data),"\n")
	data[r,]
}

ordinal<-function(num){
	library(gdata)
	n<-abs(as.integer(num))
	c<-as.character(n)
	j<-nchar(c)
	i<-ifelse((j>1 && n<14),j-1,j)
	paste0(num,case(substr(c, i, j),"st"="1","nd"="2","rd"="3",default="th"))
}

na2zero <- function(x)
{
	ifelse(is.na(x),0,x)
}

presence.absence <- function(x,t=1e-09,logical=FALSE,inverse=FALSE)
{
	A<-ifelse(inverse,0,1)
	B<-ifelse(inverse,1,0)
	if(logical){
		A <- A==TRUE  # T if A==1 else F
		B <- B==TRUE
	}
	ifelse(x>t,A,B)	
}

presence<-function(a,zero=1e-09){ # accepts vectors,dataframes, lists
	sapply(a,function(x) {x[which(x>zero)]<-1;x[which(!(x>zero))]<-NA;return(x)})
}

null.device<-function()
{
	switch(.Platform$OS.type,
		   windows={"NUL"},
		   unix={"/dev/null"}
	)
}

catline<-function(string=NULL,add,increment=NULL,sep = "",cr="\r")
{
	s<-paste(string,add,sep=sep)
	cat(paste0(cr,increment,s,"       "))
	invisible(s)
}

normalize<-function(x)
{
	(x - min(x))/(max(x)-min(x))
}

clean.null.list<-function(List)
{
	for(i in length(List):1){
		if(is.null(List[[i]])){
			List[[i]]<-NULL
		}
	}
	return(List)
}
