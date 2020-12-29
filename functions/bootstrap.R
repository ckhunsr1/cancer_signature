#=================bootstrapping===========================================
resample.ie<-function(np,mx,exposure,aage,ppyear,bbyear,type.p="Auto",weights=NULL,order=NULL,type.c="Auto",order.c=NULL){
         
         #np represents the period years that need projections
         #mx is the mortality rates 
         #exposure is the exposure population corresponding to the mx
         #aage is the age category 
         #ppyear is the period category
         #bbyear is the cohort category
         #type.p indicates the method that is used for period coefficients extrapolation."uwp" means using unweighted liner projection
         #and "wp" means using linear projection with weights. If type.p=="wp", a vector indicates the weight function for the linear 
         #projection has to be specified here in the weights parameter.If type.p="Auto", ARIMA method is used for extending period coefficients 
         #and the algorithm would automatically choose the best ARIMA model. If type.p="ARIMA", ARIMA method is also implemented, but a user 
         #specific order, such as (1,0,1), has to assign to the parameter order to guide the ARIMA model. The default is type.p="Auto".     
         #type.c indicates the method that is used for cohort coefficients extrapolation."tc" means using all points for the cohort coefficients 
         #linear projection; "10c" means using the last 10 points; and "5c" means using the last 5 points. The ARIMA options for cohort coefficients
         #are the same as for the period coefficients: "Auto" and "ARIMA" with specified order.c. The default here is "Auto".     

         temp.mx<-ie.pred(np[1],mx,exposure,aage,ppyear,bbyear,type.p,weights,order,type.c,order.c,type.re="all")
         fitted.mx<-log(temp.mx[(1:(length(temp.mx)-length(aage)))])
         rs.mx<-log(mx)-fitted.mx
         mx.new<-exp(sample(rs.mx,replace=T)+fitted.mx)
         pre.temp<-ie.pred(np,mx.new,exposure,aage,ppyear,bbyear,type.p,weights,order,type.c,order.c,type.re="prediction")
         return(pre.temp)
} 



