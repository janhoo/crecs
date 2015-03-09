#
# Contents
#
# generic helper funcs:
#

printf <- function(...) invisible(print(sprintf(...)))
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## add trancparency to color vector of either col=c("red","blue") or col=c("#FF000080" ,"#0000FF80")
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
transluc<-function(col="black",alpha=0.5){
	transparency<-round(255*alpha)
	if(class(col)=="character"){
		if(grepl("^#",col[1])){
			paste0(substr(col, 1, 7),as.hexmode(transparency))
		} else {
			v<-col2rgb(col)
			s<-sapply(data.frame(v),function(x){rgb(x[1],x[2],x[3],transparency, maxColorValue=255)})
			s
		}
	} else {
		cat("Use colors in the form of either col=c(\"red\",\"blue\") or col=c(\"#FF000080\" ,\"#0000FF80\")","\n","Done nothing.\n")
	}
}
# testing
if(FALSE){
	col<-c(red="red",  blue="blue",green="green",red2="red")
	(c2<-transluc(col))
	(c3<-transluc(col=c2))
	plot(c(1,2,3,4),col=c3)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# FIXME combine 2 rows in df with "\pm"  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ampersand <- function(df,cols,digits=NULL){
	stopifnot(all(length(cols)==2),class(df)=="data.frame")
	if(class(cols)=="character") cols<-match(cols,names(df))
	if(is.null(digits)){
		df[,cols[1]]<-paste(df[,cols[1]],"$\\pm$",df[,cols[2]],sep=" ")
	} else {
		fmt<-paste0("%.",digits,"f")
		df[,cols[1]]<-paste(sprintf(fmt,df[,cols[1]]),"$\\pm$",sprintf(fmt,df[,cols[2]]),sep=" ")
	}
	df[,-cols[2]]	
}
# testing
if(FALSE){
	(df<-data.frame(a=runif(10,min = 5, max = 10),b=runif(10),c=runif(10,min = 50, max = 100),d=runif(10,min = 2, max = 3)))
	cols<-c("a","b")
	(   MAT<-ampersand(df=df,cols=cols,digits=1)	)
	library(xtable)
	xtable(MAT)
	print(xtable(MAT),include.rownames=FALSE,sanitize.text.function = identity)
	#FIXME ther are sanitize.rownames.function, sanitize.colnames.function, sanitize.text.function
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%