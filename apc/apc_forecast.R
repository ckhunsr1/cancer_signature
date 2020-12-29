args = commandArgs(trailingOnly=TRUE)

##################################################################################################################################
                                                ##Initiation##
##################################################################################################################################
source("/gpfs/group/dxl46/default/private/poom/apc/functions/apc_func.R")
library(segmented)
library(dplyr)

type = args[1] ##Either Incidence or Mortality##
pheno = args[2] ##Trait number according to the mapping file##

load(paste("/gpfs/group/dxl46/default/private/poom/apc/input/", type, "/", type, "_", pheno, ".RData", sep = "")) ##input to apc##
output = apc2(input)

n_pred = 5 ##Number of points to predict in the future##
n_bs = 1000 ##Number of bootstraps to be performed

##Define parameters according to input data##
age_range = output$Inputs$D$a
period_range = output$Inputs$D$p
cohort_range = output$Inputs$D$c
interval = age_range[2] - age_range[1]

##################################################################################################################################
						##Prediction according to Rosenberg et al.##
##################################################################################################################################
source("/gpfs/group/dxl46/default/private/poom/apc/functions/apc_pred.R")
##Find predictions from our original data##
m.rose = apc_pred(input, output)
colnames(m.rose)[4] = "rate"

fitted = log(as.vector(output$FittedRates$events/output$FittedRates$offset))
observed = as.vector(input$events/input$offset)
residual = log(observed) - fitted

for ( bs in 1:n_bs ) {
        print(bs)
	set.seed(bs)
        new_data = exp(sample(residual, replace = T) + fitted)
        new_data = matrix(new_data, byrow=F, ncol = length(period_range))
        input$events = input$offset * new_data
        output = apc2(input)

        m_temp = apc_pred(input, output)
        m.rose = cbind(m.rose, m_temp[,4])
        colnames(m.rose)[4+bs] = paste("rate_", bs, sep = "")
}


##################################################################################################################################
                                                ##Prediction according to Dr. Wenjiang J. Fu##
##################################################################################################################################
rate = as.vector(input$events/input$offset)
exp = as.vector(input$offset)
period_future = seq(from=interval + period_range[length(period_range)], to = (interval*n_pred) + period_range[length(period_range)], by = interval)

source("/gpfs/group/dxl46/default/private/poom/apc/functions/ie.pred.r")

##ARIMA projection
mx.arima <- ie.pred(period_future, rate, exp, age_range, period_range, cohort_range
                 ,type.p="Auto",type.c="Auto",type.re="all")
mx.arima.obs = mx.arima[1:(length(mx.arima) - length(age_range)*length(period_future)), ]
mx.arima.pred = mx.arima[(length(mx.arima) - length(age_range)*length(period_future) +1):length(mx.arima), ]

source("/gpfs/group/dxl46/default/private/poom/apc/functions/bootstrap.r")

##Generate bootstrap samples and perform arima prediction##
temp.arima<-sapply(1:n_bs, function(np,mx,exposure,age,pyear,byear,type.p,weights,order,type.c,order.c) 
		   resample.ie(period_future, rate, exp, age_range, period_range, cohort_range, type.p="Auto", weights=NULL, order=NULL, type.c="Auto", order.c=NULL))     
colnames(temp.arima) = paste("rate", 1:n_bs, sep = "_")

##Create table##
m.ie.obs = as.data.frame(cbind("Year" = rep(period_range, each = length(age_range)), "Age" = rep(age_range, length.out = length(mx.arima.obs)), rate = 100000*mx.arima.obs))
m.ie.obs$Cohort = m.ie.obs$Year - m.ie.obs$Age
m.ie.obs = m.ie.obs %>% select(Year, Cohort, Age, rate)
m.ie.obs = cbind(m.ie.obs, matrix( rep(m.ie.obs$rate, length.out = n_bs*length(m.ie.obs$rate)), ncol = n_bs, byrow=FALSE) )
colnames(m.ie.obs)[5:ncol(m.ie.obs)] = paste("rate", 1:n_bs, sep = "_")

m.ie.pred = as.data.frame(cbind("Year" = rep(period_future, each = length(age_range)), "Age" = rep(age_range, length.out = length(mx.arima.pred)), "rate" = 100000*mx.arima.pred))
m.ie.pred$Cohort = m.ie.pred$Year - m.ie.pred$Age
m.ie.pred = m.ie.pred %>% select(Year, Cohort, Age, rate)
m.ie.pred = cbind(m.ie.pred, 100000*temp.arima)

m.ie = rbind(m.ie.obs, m.ie.pred)

##################################################################################################################################
                                                ##Re-order both tables##
##################################################################################################################################
m.rose$id = paste(m.rose$Age, m.rose$Year, m.rose$Cohort, sep = "_")
m.ie$id = paste(m.ie$Age, m.ie$Year, m.ie$Cohort, sep = "_")
m.rose = m.rose[match(m.ie$id, m.rose$id), ]

m.rose = m.rose %>% select(-c("id"))
m.ie = m.ie %>% select(-c("id"))

##################################################################################################################################
                                                ##Find overall rates across age-group##
##################################################################################################################################
m.rose.all = data.frame()
pop = input$offset
future_period = setdiff(unique(m.rose$Year), period_range)

pop_pred <- function(period_range, pop_ss) {
	data = data.frame("period" = period_range, "pop" = pop_ss)
	lm = lm(pop ~ period, data = data)
	out = tryCatch({
                        temp = segmented(lm, seg.z = ~period, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE, seed = 123))
                        predict.segmented(temp, data.frame("period" = future_period))
                        }, error = function(err) {
			message(" Joinpoint fails and linear regression is used instead.")
                        temp = predict.lm(lm, data.frame("period" = future_period))
                        return(temp)
                        })
	return(out)
}

##Predict future population at risk using Joinpoint regresssion##
pop_future = data.frame()
for ( a in 1:length(age_range) ) {
	print(a)
	pop_ss = pop[a, ]
	temp = pop_pred(period_range, pop_ss)
	pop_future = rbind(pop_future, temp)
}
pop = cbind(pop, pop_future)
colnames(pop) = c(period_range, future_period)
rownames(pop) = age_range
pop = pop/100000

for ( y in unique(m.rose$Year) ) {
	m_temp = m.rose %>% dplyr::filter(Year == y)
	m_temp = m_temp[order(m_temp$Age), ]
	rownames(m_temp) = m_temp$Age
	m_temp = m_temp %>% select(-c(Year, Cohort, Age))	

	##Pooling all age group together##
	m_inc = m_temp * rep(pop[,as.character(y), drop = FALSE], each = ncol(m_temp))
	m_inc = colSums(m_inc)/sum(pop[,as.character(y), drop = FALSE])

	m_inc = as.data.frame(t(as.matrix(m_inc)))
	m_inc = cbind("Year" = y, "Cohort" = "All", "Age" = "All", m_inc)
	m.rose.all = rbind(m.rose.all, m_inc)
	rm(m_temp, m_inc)
}

m.ie.all = data.frame()
for ( y in unique(m.ie$Year) ) {
        m_temp = m.ie %>% dplyr::filter(Year == y)
        m_temp = m_temp[order(m_temp$Age), ]
        rownames(m_temp) = m_temp$Age
        m_temp = m_temp %>% select(-c(Year, Cohort, Age))

        ##Pooling all age group	together##
        m_inc = m_temp * rep(pop[,as.character(y), drop = FALSE], each = ncol(m_temp))
        m_inc =	colSums(m_inc)/sum(pop[,as.character(y), drop = FALSE])

        m_inc = as.data.frame(t(as.matrix(m_inc)))
        m_inc = cbind("Year" = y, "Cohort" = "All", "Age" = "All", m_inc)
        m.ie.all = rbind(m.ie.all, m_inc)
        rm(m_temp, m_inc)

}

save(m.rose, m.ie, m.rose.all, m.ie.all, file = paste("/gpfs/group/dxl46/default/private/poom/apc/output/", type, "/", type, "_", pheno, "_output.RData", sep = ""))
