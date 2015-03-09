
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

massage<-function(species,data)
{
	# some very ungeneric data cleaning routine
	lines<-switch(species,
				  PholoeB={
				  	c(which(log(data$mgs)>6.4 & data$mass>30) , 
				  	  which(log(data$mgs)>5.5 & data$mass>100) , 
				  	  which(((data$mud^.25)>2.3) & data$mass>110) )
				  },
				  CorbulaG={		
				  	c( which(log(data$mgs)>5.9 & data$mass>850) ,
				  	   which(((data$mud^.25)>2.25) & data$mass>2500))
				  },
				  AonidesP={
				  	which(((data$mud)^.25<(1.9)) & data$mass>120) 
				  },
				  OpheliaL={
				  	which(((data$mud)^.25<(0.1)) & data$mass>30)
				  },
				  OphiodromusF={
				  	which((log(data$mgs)>(5.6)) & data$mass>65)
				  },
				  HarpiniaA={
				  	c(which((log(data$mgs)>(5.3)) & data$mass>2)  , 
				  	  which((log(data$mgs)<(4.15)) & data$mass>44) )  	   	
				  },
				  NephtysCi={
				  	c(which((log(data$mgs)<(4.3)) & data$mass>3000) , 
				  	  which((log(data$mgs)>(7)) & data$mass>200))
				  },
				  TellinaFa ={	
				  	c(which(((data$mud)^.25>(1.8)) & data$mass>200),
				  	  which(((data$mud)^.25<(0.3)) & data$mass>200),
				  	  which((log(data$mgs)>(5.9)) & data$mass>1000),
				  	  which((log(data$mgs)>(6.2)) & data$mass>1))		
				  },
				  NuculaNi ={	
				  	c(which(((data$depth)<(24)) & data$mass>2000)	,
				  	  which(((data$depth)>(47)) & data$mass<200)	, 
				  	  which((log(data$mgs)>(5.6)) & data$mass>1000), 
				  	  which((log(data$mgs)<(4.2))), 
				  	  which(((data$mud^.25)>(2.2))), 
				  	  which(((data$depth)>(41)) ))
				  },
				  NephropsN ={	
				  	c(which(((data$gsd)>(1)) & data$mass>20),
				  	  which(((data$mud)>(20)) & data$mass>20 ),
				  	  which(((data$depth)>(43))  ))
				  	
				  },
				  Upogebia ={		
				  	c(which(((data$depth)>(43)) ),
				  	  which(((data$depth)<(35)) & data$mass>150 ),
				  	  which((log(data$mgs)>(5.2)) & data$mass>100),
				  	  which((log(data$ssWav)<(-0.2))))
				  },
				  Callianassa ={		
				  	c(which((log(data$mgs)>(5.23)) & data$mass>450),
				  	  which((log(data$mgs)<(4.5)) & data$mass>1500),
				  	  which((log(data$mgs)<(4.2)) & data$mass>150),
				  	  which(((data$mud^.25)<(1.29)) & data$mass>300),
				  	  which(((data$mud^.25)>(2.15)) & data$mass>10))
				  },
				  AmphiuraF ={
				  	c(which((log(data$mgs)>(5.7)) & data$mass>25000),
				  	  which((log(data$mgs)<(4.2)) & data$mass>17000),
				  	  which(((data$mud^.25)>(2.15))),
				  	  which(((data$mud^.25)<(0.35)) & data$mass>10000))
				  },
				  EchinoC ={
				  	c(which((log(data$mgs)>(6.4)) & data$mass>45000),
				  	  which((log(data$mgs)>(6)) & data$mass>200000),
				  	  which((log(data$mgs)<(4.15)) & data$mass>30000),
				  	  which(((data$mud^.25)>(2.15)) & data$mass>20000),
				  	  which(((data$mud^.25)<(0.35)) & data$mass>10000))
				  },
				  AbraN={
				  	c(which((log(data$mgs)<(4.5)) & data$mass>20),
				  	  which((log(data$mgs)>(5.5)) & data$mass>400),
				  	  which(((data$mud^.25)<(1.2)) & data$mass>170),
				  	  which(((data$mud^.25)>(1.93)) & data$mass>420))
				  },
				  {
				  	integer()
				  }
	)
	if(length(lines)==0L)
	{
		return(data)
	} else {
		return(data[-lines,])
	}
}
