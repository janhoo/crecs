metrics <- function(response,prediction,secondary=NULL,threshold.def="spec_sens",digits=6,pval=TRUE,...)
{
    # calculates various model metrics
	
	opts.threshold.def<-c("kappa","spec_sens" , "no_omission","prevalence","equal_sens_spec","sensitivty" )
	library(dismo)   # evaluate()
	library(minerva) # Maximal Information-Based Nonparametric Exploration (MINE) statistics
    if(is.null(secondary)){
    	secondary<-response
    }
    
    p <- as.numeric(prediction[presence.absence(response,logical=T)])
    a <- as.numeric(prediction[presence.absence(response,logical=T,inverse=T)])
    e<- evaluate(p=p,a=a)
    
    true.presence<-presence.absence(response,logical=T)
    threshold= as.numeric(threshold(e)[threshold.def])
    predicted.presence<-presence.absence(x=prediction,t=threshold,logical=TRUE)
    #predicted.presence<-prediction>threshold(e)[[threshold.def]]
    paretoopt<-which.max(e@TPR + e@TNR)
    # write results
	#browser()
    ###############
    result<-NULL
    result<-c(result, auc      = slot(e,'auc'))
    result<-c(result, kappa    = e@kappa[paretoopt])
    result<-c(result, overdp   = e@ODP[paretoopt])
    result<-c(result, sensi    = e@TPR[paretoopt])
    result<-c(result, spezi    = e@TNR[paretoopt])
    result<-c(result, tss      = e@TPR[paretoopt]+e@TNR[paretoopt]-1)
    result<-c(result, misclass = e@MCR[paretoopt])
       var.residuals <- (response-prediction)^2
       var.response  <- (response-mean(response))^2
    result<-c(result, mse      = mean( var.residuals ))
    result<-c(result, rmse      = sqrt(mean( var.residuals )))
    result<-c(result, rmse.trivi      = sqrt(mean( var.response )))
    result<-c(result, rmse.trivii      = sqrt(mean((response-ifelse(presence.absence(prediction,threshold,logical=TRUE),mean(response[presence.absence(response,logical=T)]),0))^2)))
    result<-c(result, rsqr      = 1-sum(var.residuals)/sum(var.response)  )
    result<-c(result, mape      = mean( abs(response-prediction) ))
	result<-c(result, mic      = mine(secondary,prediction)$MIC)
	result<-c(result, mic.predicted.presence      = mine(secondary[predicted.presence],prediction[predicted.presence])$MIC)
	result<-c(result, mic.true.presence      = mine(secondary[true.presence],prediction[true.presence])$MIC)
	result<-c(result, threshold=threshold)
	if(pval){
		spp<-secondary[predicted.presence]
		ppp<-prediction[predicted.presence]
		stp<-secondary[true.presence]
		ptp<-prediction[true.presence]
		
		
		C<-cor.test(secondary,prediction,method="pearson")
		result<-c(result, pearson.all                =  C$estimate)
		result<-c(result, pearson.all.p                =  C$p.value)
		C<-cor.test(spp,ppp,method="pearson")
		result<-c(result, pearson.predicted.presence                =  C$estimate)
		result<-c(result, pearson.predicted.presence.p                =  C$p.value)
		C<-cor.test(stp,ptp,method="pearson")
		result<-c(result, pearson.true.presence                =  C$estimate)
		result<-c(result, pearson.true.presence.p                =  C$p.value)		
		C<-cor.test(log(secondary+1),log(prediction+1),method="pearson")
		result<-c(result, pearson.log.all                =  C$estimate)
		result<-c(result, pearson.log.all.p                =  C$p.value)					
		C<-cor.test(log(stp+1),log(ptp+1),method="pearson")
		result<-c(result, pearson.log.true.presence                =  C$estimate)
		result<-c(result, pearson.log.true.presence.p                =  C$p.value)
		C<-cor.test(sqrt(secondary),sqrt(prediction),method="pearson")
		result<-c(result, pearson.sqrt.all                =  C$estimate)
		result<-c(result, pearson.sqrt.all.p                =  C$p.value)
		C<-cor.test(sqrt(stp),sqrt(ptp),method="pearson")
		result<-c(result, pearson.sqrt.true.presence                =  C$estimate)
		result<-c(result, pearson.sqrt.true.presence.p                =  C$p.value)		
		C<-cor.test(secondary,prediction,method="spearman")
		result<-c(result, spearman.all                =  C$estimate)
		result<-c(result, spearman.all.p                =  C$p.value)		
		C<-cor.test(spp,ppp,method="spearman")
		result<-c(result, spearman.predicted.presence                =  C$estimate)
		result<-c(result, spearman.predicted.presence.p                =  C$p.value)				
		C<-cor.test(stp,ptp,method="spearman")
		result<-c(result, spearman.true.presence                =  C$estimate)
		result<-c(result, spearman.true.presence.p                =  C$p.value)		
		C<-cor.test(secondary,prediction,method="kendall")
		result<-c(result, kendall.all                =  C$estimate)
		result<-c(result, kendall.all.p                =  C$p.value)		
		C<-cor.test(spp,ppp,method="kendall")
		result<-c(result, kendall.predicted.presence                =  C$estimate)
		result<-c(result, kendall.predicted.presence.p                =  C$p.value)				
		C<-cor.test(stp,ptp,method="kendall")
		result<-c(result, kendall.true.presence                =  C$estimate)
		result<-c(result, kendall.true.presence.p                =  C$p.value)			
	} else {
		result<-c(result, pearson.all                =  cor(secondary,prediction,method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.predicted.presence =  cor(secondary[predicted.presence],prediction[predicted.presence],method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.true.presence      =  cor(secondary[true.presence],prediction[true.presence],method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.log.all            =  cor(log(secondary+1),log(prediction+1),method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.log.true.presence  =  cor(log(secondary[true.presence]+1),log(prediction[true.presence]+1),method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.sqrt.all            =  cor(sqrt(secondary),sqrt(prediction),method="pearson",use="pairwise.complete.obs"))
		result<-c(result, pearson.sqrt.true.presence  =  cor(sqrt(secondary[true.presence]),sqrt(prediction[true.presence]),method="pearson",use="pairwise.complete.obs"))
		result<-c(result, spearman.all               =  cor(secondary,prediction,method="spearman",use="pairwise.complete.obs"))
		result<-c(result, spearman.predicted.presence= cor(secondary[predicted.presence],prediction[predicted.presence],method="spearman",use="pairwise.complete.obs"))
		result<-c(result, spearman.true.presence     =  cor(secondary[true.presence],prediction[true.presence],method="spearman",use="pairwise.complete.obs"))
		result<-c(result, kendall.all                =  cor(secondary,prediction,method="kendall",use="pairwise.complete.obs"))
		result<-c(result, kendall.predicted.presence =  cor(secondary[predicted.presence],prediction[predicted.presence],method="kendall",use="pairwise.complete.obs"))
		result<-c(result, kendall.true.presence      =  cor(secondary[true.presence],prediction[true.presence],method="kendall",use="pairwise.complete.obs"))		
	}	
    # if(any(sapply(result,is.na))){stop("NA in metrics")}
    return(result)
}

