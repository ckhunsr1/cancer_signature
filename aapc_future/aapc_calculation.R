args = commandArgs(trailingOnly=TRUE)

#####################################################################################################################################################################
aapc_calc <- function(seg) {
	if ( length(seg$coefficients) > 2 ) {
		out = 100*aapc(seg, exp.it = TRUE)[1]
	} else {
		out = 100*(exp(seg$coefficients[2]) - 1)
	}
	out
}
#####################################################################################################################################################################
library(dplyr)
library(segmented)

#type = "Incidence"
#pheno_list = c("2", "14", "15", "17", "32", "33", "37", "44", "50", "52", "54", "56", "58", "63", "64", "68", "69", "77", "80", "83", "87")

#type = "Mortality"
#pheno_list = c("3", "15", "16", "18", "22", "23", "27", "34", "40", "42", "44", "46", "48", "53", "54", "58", "59", "65", "68", "69", "71")

type = args[1] ##Either Incidence or Mortality##
pheno = args[2] ##Trait number according to the mapping file##

load(paste("/gpfs/group/dxl46/default/private/poom/apc/output/", type, "/", type, "_", pheno, "_output.RData", sep = ""))
m.rose.total = rbind(m.rose, m.rose.all)
m.rose.total[,1] = m.rose.total[,1] - 1977.5

m.ie.total = rbind(m.ie, m.ie.all)
m.ie.total[,1] = m.ie.total[,1] - 1977.5

rose_output = data.frame()
ie_output = data.frame()
for (a in as.character(unique(m.rose.total$Age))){
	print(a)
	m.rose.all = m.rose.total %>% filter(Age == a)
	m.ie.all = m.ie.total %>% filter(Age == a)

	##Calculate AAPC from actual estimate##
	lm.rose = lm(log(rate) ~ Year, data = m.rose.all)
	my.rose.seg = tryCatch({
			segmented(lm.rose, seg.z = ~Year, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE))
			}, error = function(err) {
			message(" Joinpoint fails and linear regression is used instead.")
			return(lm.rose)
			})
	aapc.rose.est = aapc_calc(my.rose.seg)

	lm.ie = lm(log(rate) ~ Year, data = m.ie.all)
	my.ie.seg = tryCatch({
			segmented(lm.ie, seg.z = ~Year, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE))
			}, error = function(err) {
                        message(" Joinpoint fails and linear regression is used instead.")
			return(lm.ie)
                        })
	aapc.ie.est =	aapc_calc(my.ie.seg)

	##Calculate AAPC from bootstrap results##
	aapc.rose.dist = c()
	aapc.ie.dist = c()
	for ( i in 5:ncol(m.rose.all) ) {
		m.rose.ss = m.rose.all[, c(1,i)]
		colnames(m.rose.ss)[2] = "rate"
		lm.rose.ss = lm(log(rate) ~ Year, data = m.rose.ss)
	        my.rose.seg.ss = tryCatch({
				segmented(lm.rose.ss, seg.z = ~Year, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE))
				}, error = function(err) {
				message(" Joinpoint fails and linear regression is used instead.")
				return(lm.rose.ss)
				})
		aapc.rose.est.ss = aapc_calc(my.rose.seg.ss)	
		aapc.rose.dist = c(aapc.rose.dist, aapc.rose.est.ss)
	
		m.ie.ss = m.ie.all[, c(1,i)]
	        colnames(m.ie.ss)[2] = "rate"
	        lm.ie.ss = lm(log(rate) ~ Year, data = m.ie.ss)
	        my.ie.seg.ss = tryCatch({
				segmented(lm.ie.ss, seg.z = ~Year, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE))
				}, error = function(err) {
				message(" Joinpoint fails and linear regression is used instead.")
                                return(lm.ie.ss)
                                })
		aapc.ie.est.ss = aapc_calc(my.ie.seg.ss)
		aapc.ie.dist = c(aapc.ie.dist, aapc.ie.est.ss)

		rm(m.rose.ss, lm.rose.ss, my.rose.seg.ss, m.ie.ss, lm.ie.ss, my.ie.seg.ss)
	}
	aapc.rose.dist = sort(aapc.rose.dist)
	aapc.ie.dist = sort(aapc.ie.dist)

	#rose_temp = data.frame("type" = "rose", "idx" = pheno, "Age" = a, "est" = aapc.rose.est, "lower" = aapc.rose.dist[0.025*(ncol(m.rose.all) - 4)], 
	#		"upper" = aapc.rose.dist[0.975*(ncol(m.rose.all) - 4)])

	rose_temp = data.frame("type" = "rose", "idx" = pheno, "Age" = a, "est" = aapc.rose.est, "lower" = aapc.rose.est - (1.96*sd(aapc.rose.dist)), 
			"upper" = aapc.rose.est + (1.96*sd(aapc.rose.dist)))
	rose_output = rbind(rose_output, rose_temp)

	#ie_temp = data.frame("type" = "ie", "idx" = pheno, "Age" = a, "est" = aapc.ie.est, "lower" = aapc.ie.dist[0.025*(ncol(m.ie.all) - 4)], 
	#		"upper" = aapc.ie.dist[0.975*(ncol(m.ie.all) - 4)])
	ie_temp = data.frame("type" = "ie", "idx" = pheno, "Age" = a, "est" = aapc.ie.est, "lower" = aapc.ie.est - (1.96*sd(aapc.ie.dist)), 
                        "upper" = aapc.ie.est + (1.96*sd(aapc.ie.dist)))
	ie_output = rbind(ie_output, ie_temp)
}

save(aapc.rose.est, aapc.rose.dist, aapc.ie.est, aapc.ie.dist, file = paste("/gpfs/group/dxl46/default/private/poom/apc/output/AAPC_", type, "/AAPC_", type, "_", pheno, "_output.RData", sep = ""))
write.table(rose_output, paste("/gpfs/group/dxl46/default/private/poom/apc/output/AAPC_", type, "/AAPC_", type, "_", pheno, "_rose_output.txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(ie_output, paste("/gpfs/group/dxl46/default/private/poom/apc/output/AAPC_", type, "/AAPC_", type, "_", pheno, "_ie_output.txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

