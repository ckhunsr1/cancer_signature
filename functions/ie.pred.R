source("/gpfs/group/dxl46/default/private/poom/apc/functions/apcmat.fun.r")
source("/gpfs/group/dxl46/default/private/poom/apc/functions/constr.fun.r")
source("/gpfs/group/dxl46/default/private/poom/apc/functions/apcrdg.fun.r")
source("/gpfs/group/dxl46/default/private/poom/apc/functions/apcrdgse.fun.r")

library("forecast")

ie.pred<-function(np,mx,exposure,aage,ppyear,bbyear,type.p="Auto",weights=NULL,order=NULL,type.c="Auto",order.c=NULL,type.re="prediction"){

         #np represents the period years for projection
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
         #are the same as for the period coefficients: "Auto" and "ARIMA" with specified order.c. The default here is "Auto"     
         #type.re indicates the return value. "all" means return both fitted and predicted values."predition" indicates only return predicted values.
         #    The default is "prediction"

         #===================================obtain IE coefficients=======================================================================
         n.a<-length(aage)
         n.p<-length(ppyear)
         n.c<-length(bbyear)

         m6.mx<-matrix(mx,byrow=F,ncol=n.p)
         m6.exposure<-matrix(exposure,byrow=F,ncol=n.p)
         fit.6<-apcrdgse(r=m6.mx,nrisk=m6.exposure,fam="qlik",offst=T,lam=-1,stderr=F)
         
         #====================================period coefficient extrapolation=============================================================
         pp.df<-data.frame(ppyear=np)
         if (type.p=="uwp"){ 
              par.pp<-c(fit.6$pp,predict(lm(fit.6$pp~ppyear),pp.df))
         }
         if (type.p=="wp"){
              w.p<-weights
              par.pp<-c(fit.6$pp,predict(lm(fit.6$pp~ppyear,weights=w.p),pp.df))
         }
         if (type.p=="qp"){
             par.pp<-c(fit.6$pp,predict(lm(fit.6$pp~ppyear+I(ppyear^2)),pp.df))
         }
         if (type.p=="ep"){
             par.pp<-c(fit.6$pp,fit.6$pp[length(fit.6$pp)]*(0.98)^(1:length(np)))
         }  
         if (type.p=="ARIMA"){
             arima.pp<-arima(fit.6$pp,order)
             par.pp<-c(fit.6$pp,forecast(arima.pp,length(np))$mean)
         }
         if (type.p=="Auto"){
             arima.pp<-auto.arima(fit.6$pp)
             par.pp<-c(fit.6$pp,forecast(arima.pp,length(np))$mean)
         }       
         #====================================cohort coefficient extrapolation=============================================================
         if (type.c=="tc"){
             bb.df<-data.frame(bbyear=seq(from=bbyear[length(bbyear)]+5,by=5,length.out=length(np)))   
             par.bb<-c(fit.6$coh,predict(lm(fit.6$coh~bbyear),bb.df))
         }
         if (type.c=="10c"){
             bbbyear<-bbyear[-(1:(length(bbyear)-10))]
             bbb.df<-data.frame(bbbyear=seq(from=bbyear[length(bbyear)]+5,by=5,length.out=length(np)))   
             par.bb<-c(fit.6$coh,predict(lm(fit.6$coh[-(1:(length(bbyear)-10))]~bbbyear),bbb.df))         }
         if (type.c=="5c"){
             bbbyear<-bbyear[-(1:(length(bbyear)-5))]
             bbb.df<-data.frame(bbbyear=seq(from=bbyear[length(bbyear)]+5,by=5,length.out=length(np)))   
             par.bb<-c(fit.6$coh,predict(lm(fit.6$coh[-(1:(length(bbyear)-5))]~bbbyear),bbb.df))
         }
         if (type.c=="ARIMA"){
             arima.bb<-arima(fit.6$coh,order.c)
             par.bb<-c(fit.6$coh,forecast(arima.bb,length(np))$mean)
         }
         if (type.c=="Auto"){
             arima.bb<-auto.arima(fit.6$coh)
             par.bb<-c(fit.6$coh,forecast(arima.bb,length(np))$mean)
         }
         par.ie<-c(fit.6$aa,par.pp,par.bb)
     
         #====================================prepare the variable matrix xx==================================================================== 
         # xx represents the prediction vector or matrix. The first n.a column represents age indicator, the next n.p column is the period indicator, while
         #    the last n.a+n.p-1 is the cohort inidcator. Each XX cell represents the single categore of age, period or cohort

         xx.aa<-matrix(0,n.a*length(par.pp),n.a)
         for (i in 1:length(par.pp)){
             xx.aa[(((i-1)*n.a+1):(i*n.a)),]<-diag(x=1,n.a,n.a)
         }
         xx.aa<-rbind(aage,xx.aa)
         
         xx.pp<-matrix(0,n.a*length(par.pp),length(par.pp))
         for (i in 1:length(par.pp)){
                xx.pp[(((i-1)*n.a+1):(i*n.a)),i]<-rep(1,n.a)
         }
         xx.pp<-rbind(c(ppyear,np),xx.pp)

         xx.bb<-matrix(0,n.a*length(par.pp),length(par.bb))
         xx.bb<-rbind(c(bbyear,seq(from=bbyear[length(bbyear)]+5,by=5,length.out=length(np))),xx.bb)
         for (i in 2:(n.a*length(par.pp)+1)){
              b.temp<-xx.pp[1,xx.pp[i,]==1]-xx.aa[1,xx.aa[i,]==1]
              xx.bb[i,xx.bb[1,]==b.temp]=1
         }
         xx<-cbind(xx.aa[-1,],xx.pp[-1,],xx.bb[-1,])
         
        #=====================================function return (Mx)=================================================================================           
         if (type.re=="all"){
             return(exp(fit.6$inter+xx%*%par.ie))
         }
         if (type.re=="prediction"){
             return(exp(fit.6$inter+xx%*%par.ie)[-(1:(n.a*n.p))])
        }
}
             
        
        





