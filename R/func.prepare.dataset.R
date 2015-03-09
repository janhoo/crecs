
prepare.dataset<-function(WS,DS,species,response,aggregated,pseudoabsence,seed=NULL,physio=FALSE)
{
    stopifnot( class(aggregated)=="logical" )
    stopifnot( class(pseudoabsence)=="logical" )
    stopifnot( any(response %in% c("mass" , "presence")) )
    data<-species.response(species=species,WS=WS,DS=DS)
    if(aggregated) # agg
    {
        agg<-aggregate.stations(data=data) 
        agg$presence <-presence.absence(agg[["prevalence"]])
        data <- get.environ(points= agg)
    }
    if(pseudoabsence) 
    {   
        presence<-data[which(data[["presence"]]>0),]
        sample.size<-length(which(data[["presence"]]<1))
        ps.abs<-pseudo.absence(presence=presence,sample.size=sample.size,seed=seed)
        ps.abs<-get.environ(points=ps.abs)
        ps.abs$mass<-ps.abs$presence<-0
        ps.pres<-presence[names(ps.abs)]
        data<-rbind(ps.pres,ps.abs)
        row.names(data)<-1:length(row.names(data))
    }
    if(response=="mass") 
    {
    	if(physio){
    		return(elim.outliers(data=massage(species=species,data=data)))
    	} else {
    		return(elim.outliers(data=data))
    	}
    } else {
        return(data)
    }
}

