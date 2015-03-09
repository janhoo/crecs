#
# Closures for model cv
# Creates generic method for cv procedure
#
# Usage:
# Model<-makefn.model(model,response,predictors,Predictor,data.names)
# fit<-Model(data)
#
#
#

makefn.model<-function(model,response,predictors,Predictor,data.names,gbm.trees)
{
    #uses: null.device()
    opts.model<-c("gam" , "rf" , "gbm", "max", "gbm.step","glm")
    opts.response<-c("mass" , "presence" )
    opts.predictors<-c("depth","mgs","gsd","mud","ssCav","ssCmax","ssWav","ssWmax","wtemp","atemp","sal","pp","trawlAgg","trawlAv","dist","ship")
    stopifnot( model %in% opts.model ) 
    stopifnot( response %in% opts.response)
    stopifnot( all(predictors %in% opts.predictors) )
    stopifnot( all(c(predictors,response) %in% data.names) )
    stopifnot( class(Predictor)=="SpatialPixelsDataFrame" )
    library(raster)
    #library(rgdal)
    library(mgcv)
    library(randomForest)
    library(rJava)
    library(dismo)
    Model<-switch(model,
    			  glm={
    			  	deps<-NULL
    			  	for(i in predictors){
    			  		deps<-switch(i,
    			  					 depth={c(deps,"s(depth)")},
    			  					 mgs={c(deps,"s(log(mgs))")},
    			  					 gsd={c(deps,"s(I(gsd^.5))")},
    			  					 mud={c(deps,"s(I(mud^.25))")},
    			  					 ssCav={c(deps,"s(log(ssCav))")},
    			  					 ssCmax={c(deps,"s(log(ssCmax))")},
    			  					 ssWav={c(deps,"s(log(ssWav))")},
    			  					 ssWmax={c(deps,"s(log(ssWmax))")},
    			  					 wtemp={c(deps,"s(I((1/wtemp)^.25))")},
    			  					 atemp={c(deps,"s(I((1/atemp)^.25))")},
    			  					 sal={c(deps,"s(I(log(1/sal)))")},
    			  					 pp={c(deps,"s(log(pp))")},
    			  					 trawlAgg={c(deps,"s(I(trawlAgg^.4))")},
    			  					 trawlAv={c(deps,"s(I(trawlAv^.4))")},
    			  					 dist={c(deps,"s(dist)")},
    			  					 ship={c(deps,"s(ship)")}
    			  		)            	
    			  	}
    			  	form<-as.formula(paste(response,"~",paste(deps, collapse= "+")))
    			  	fam<-switch(response,
    			  				mass={Tweedie(p=1.3,link=log)},
    			  				presence={poisson}
    			  	)
    			  	function(data)
    			  	{
    			  		glm(formula=form,family=fam,data=data)
    			  	}
    			  },
                  gam={
                      deps<-NULL
                      for(i in predictors){
                          deps<-switch(i,
                                       depth={c(deps,"s(depth)")},
                                       mgs={c(deps,"s(log(mgs))")},
                                       gsd={c(deps,"s(I(gsd^.5))")},
                                       mud={c(deps,"s(I(mud^.25))")},
                                       ssCav={c(deps,"s(log(ssCav))")},
                                       ssCmax={c(deps,"s(log(ssCmax))")},
                                       ssWav={c(deps,"s(log(ssWav))")},
                                       ssWmax={c(deps,"s(log(ssWmax))")},
                                       wtemp={c(deps,"s(I((1/wtemp)^.25))")},
                                       atemp={c(deps,"s(I((1/atemp)^.25))")},
                                       sal={c(deps,"s(I(log(1/sal)))")},
                                       pp={c(deps,"s(log(pp))")},
                                       trawlAgg={c(deps,"s(I(trawlAgg^.4))")},
                                       trawlAv={c(deps,"s(I(trawlAv^.4))")},
                                       dist={c(deps,"s(dist)")},
                                       ship={c(deps,"s(ship)")}
                          )            	
                      }
                      form<-as.formula(paste(response,"~",paste(deps, collapse= "+")))
                      fam<-switch(response,
                                  mass={Tweedie(p=1.3,link=log)},
                                  presence={poisson}
                      )
                      method<-"REML"
                      sel<-TRUE
                      function(data)
                      {
                          gam(formula=form,family=fam,data=data,method=method,select=sel)
                      }
                  },
                  rf={
                      obs<-switch(response,
                                  mass={response},
                                  presence={paste("as.factor(",response,")")}
                      )
                      form<-as.formula(paste(obs,"~",paste(predictors, collapse= "+")))
                      function(data)
                      {
                          randomForest(formula=form, data=data, importance=TRUE)
                      }
                  },
                  gbm.step={
                      gbmx<-as.numeric(match(predictors,data.names))
                      gbmy<-as.numeric(match(response,data.names))
                      
                      fam<-switch(response,
                                  mass={"gaussian"},
                                  presence={"bernoulli"}
                      )
                      function(data)
                      {
                          dat<-data.frame(data)
                          gbm.step(data=dat,gbm.x=gbmx,gbm.y=gbmy,family=fam,silent=TRUE, plot.main = FALSE)
                      }
                  },
    			  gbm={ 
    			  	fam<-switch(response,
    			  				mass={"gaussian"},
    			  				presence={"bernoulli"}
    			  	)
    			  	form<-as.formula(paste(response,"~",paste(predictors, collapse= "+")))
    			  	function(data)
    			  	{
    			  		suppressMessages(gbm(formula=form, data=data,distribution=fam,n.trees=gbm.trees,verbose=FALSE))
    			  	}
    			  },
                  max={
                      map<-as(Predictor[,match(predictors,names(Predictor))],"SpatialGridDataFrame")
                      function(data)
                      {
                          dat<-data[which(data[[response]]>0),]
                          maxent(map,dat)
                      }
                  }
    )
    return(Model)
}


#
#  predict method accordingly
#
#
makefn.predict<-function(response,model,environ,maxargs)
{
    library(raster)
    #library(rgdal)    
    library(mgcv)
    library(randomForest)
    library(dismo)
    library(gbm)
    Predict<-switch(model,
                    gam={
                        
                        function(fit,data)
                        {
                            predict(object=fit, newdata=data, type="response")
                        }
                    },
                    rf={
                        if(response=="presence")
                        {
                            function(fit,data)
                            {
                                predict(object=fit, newdata=data, type="prob")[,"1"]
                            }
                            
                        } else {
                            function(fit,data)
                            {
                                predict(object=fit, newdata=data, type="response")
                            }
                        }
                    },
    				gbm.step={
    					function(fit,data)
    					{
    						predict.gbm(object=fit, n.trees=fit$gbm.call$best.trees, newdata=data,const=add, type="response")
    					}
    				},
    				gbm={
    					function(fit,data)
    					{
    						best.iter<-suppressWarnings(gbm.perf(fit,method="OOB",plot.it=FALSE))
    						predict.gbm(object=fit, n.trees=best.iter, newdata=data, const=add,type="response")
    					}
    				},
    				max={
                        function(fit,data)
                        {
                            p<-predict(object=environ,model=fit, args=maxargs)
                            extract(p, data)
                        }
                    }
    )
    return(Predict)
}





makefn.free.model<-function(method,response,responsetype,predictors,Predictor,data.names,gbm.trees,check.names=TRUE, ...)
{
	#uses: null.device()
	opts.responsetype<-c("count","continuous" , "presence" )
	stopifnot( responsetype %in% opts.responsetype)
	opts.method<-c("gam" , "rf" , "gbm", "max", "gbm.step")
	stopifnot( all(c(predictors,response) %in% data.names) )
	stopifnot( class(Predictor)=="SpatialPixelsDataFrame" )
	library(raster)
	#library(rgdal)   
	library(mgcv)
	library(randomForest)
	library(rJava)
	library(dismo)
	Method<-switch(method,
				glm={
				   	deps<-NULL
				   	if(check.names){
				   	for(i in predictors){
				   		deps<-switch(i,
				   					 depth={c(deps,"depth")},
				   					 mgs={c(deps,"log(mgs)")},
				   					 gsd={c(deps,"I(gsd^.5)")},
				   					 mud={c(deps,"I(mud^.25)")},
				   					 ssCav={c(deps,"log(ssCav)")},
				   					 ssCmax={c(deps,"log(ssCmax)")},
				   					 ssWav={c(deps,"log(ssWav)")},
				   					 ssWmax={c(deps,"log(ssWmax)")},
				   					 wtemp={c(deps,"I((1/wtemp)^.25)")},
				   					 atemp={c(deps,"I((1/atemp)^.25)")},
				   					 sal={c(deps,"I(log(1/sal))")},
				   					 pp={c(deps,"log(pp)")},
				   					 trawlAgg={c(deps,"I(trawlAgg^.4)")},
				   					 trawlAv={c(deps,"I(trawlAv^.4)")},
				   					 dist={c(deps,"dist")},
				   					 ship={c(deps,"ship")},
                                     {c(deps,i)}
				   		)            	
				   	}	
				   	} else {
				   		deps<-predictors
				   	}
					form<-as.formula(paste(response,"~",paste(deps, collapse= "+")))
					fam<-	switch(responsetype,
								count={poisson},
								continuous={Tweedie(p=1.5,link=log)},
								presence={poisson}
							)
					function(data)
					{
						glm(formula=form,family=fam,data=data)
					}
				   	},
				gam={
					deps<-NULL
					if(check.names){
					for(i in predictors){
						deps<-	switch(i,
					 				depth={c(deps,"s(depth)")},
					 				mgs={c(deps,"s(log(mgs))")},
					 				gsd={c(deps,"s(I(gsd^.5))")},
					 				mud={c(deps,"s(I(mud^.25))")},
					 				ssCav={c(deps,"s(log(ssCav))")},
					 				ssCmax={c(deps,"s(log(ssCmax))")},
					 				ssWav={c(deps,"s(log(ssWav))")},
					 				ssWmax={c(deps,"s(log(ssWmax))")},
					 				wtemp={c(deps,"s(I((1/wtemp)^.25))")},
					 				atemp={c(deps,"s(I((1/atemp)^.25))")},
					 				sal={c(deps,"s(I(log(1/sal)))")},
					 				pp={c(deps,"s(log(pp))")},
					 				trawlAgg={c(deps,"s(I(trawlAgg^.4))")},
					 				trawlAv={c(deps,"s(I(trawlAv^.4))")},
					 				dist={c(deps,"s(dist)")},
					 				ship={c(deps,"s(ship)")},
									{c(deps,i)}
						)            	
					}
					} else {
						deps<-predictors
					}
					form<-as.formula(paste(response,"~",paste(deps, collapse= "+")))
					fam<-	switch(responsetype,
								count={poisson},
								continuous={Tweedie(p=1.3,link=log)},
								presence={poisson}
							)
					method<-"REML"
					sel<-TRUE
					function(data)
					{
						gam(formula=form,family=fam,data=data,method=method,select=sel)
					}
				},
				rf={
					obs<-	switch(responsetype,
								count={response},
								continuous={response},
								presence={paste("as.factor(",response,")")}
							)
					form<-as.formula(paste(obs,"~",paste(predictors, collapse= "+")))
					function(data)
					{
						randomForest(formula=form, data=data, importance=TRUE)
					}
				},
				gbm.step={
					gbmx<-as.numeric(match(predictors,data.names))
					gbmy<-as.numeric(match(response,data.names))
					fam<-	switch(responsetype,
								count={"poisson"},
								continuous={"gaussian"},
								presence={"bernoulli"}
							)
					function(data)
					{
						dat<-data.frame(data)
						gbm.step(data=dat,gbm.x=gbmx,gbm.y=gbmy,family=fam,silent=TRUE, plot.main = FALSE)
					}
				},
				gbm={ 
					fam<-	switch(responsetype,
								count={"poisson"},
								continuous={"gaussian"},
								presence={"bernoulli"}
							)
					form<-as.formula(paste(response,"~",paste(predictors, collapse= "+")))
					function(data)
					{
						suppressMessages(gbm(formula=form, data=data,distribution=fam,n.trees=gbm.trees,verbose=FALSE))
					}
				},
				max={
					map<-as(Predictor[,match(predictors,names(Predictor))],"SpatialGridDataFrame")
					function(data)
					{
						dat<-data[which(data[[response]]>0),]
						maxent(map,dat)
					}
				}
			)	
			return(Method)
}

makefn.free.predict<-function(responsetype,method,environ,maxargs)
{ 
  library(raster)
  #library(rgdal)   
  library(mgcv)
	library(randomForest)
	library(dismo)
	library(gbm)
	Predict<-switch(method,
					glm={
						
						function(fit,data)
						{
							predict(object=fit, newdata=data, type="response")
						}
					},
					gam={
						
						function(fit,data)
						{
							predict(object=fit, newdata=data, type="response")
						}
					},
					rf={
						if(responsetype=="presence")
						{
							function(fit,data)
							{
								predict(object=fit, newdata=data, type="prob")[,"1"]
							}
							
						} else {
							function(fit,data)
							{
								predict(object=fit, newdata=data, type="response")
							}
						}
					},
					gbm.step={
						function(fit,data)
						{
							predict.gbm(object=fit, n.trees=fit$gbm.call$best.trees, newdata=data,const=add, type="response")
						}
					},
					gbm={
						function(fit,data)
						{
							best.iter<-suppressWarnings(gbm.perf(fit,method="OOB",plot.it=FALSE))
							predict.gbm(object=fit, n.trees=best.iter, newdata=data, const=add,type="response")
						}
					},
					max={
						function(fit,data)
						{
							p<-predict(object=environ,model=fit, args=maxargs)
							extract(p, data)
						}
					}
	)
	return(Predict)
}
