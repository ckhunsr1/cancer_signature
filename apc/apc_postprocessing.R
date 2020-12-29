library(dplyr)

#type = "Incidence"
#pheno_list = c("2", "14", "15", "17", "32", "33", "37", "44", "50", "52", "54", "56", "58", "63", "64", "68", "69", "77", "80", "83", "87")

type = "Mortality"
pheno_list = c("3", "15", "16", "18", "22", "23", "27", "34", "40", "42", "44", "46", "48", "53", "54", "58", "59", "65", "68", "69", "71")

for ( pheno in pheno_list ) {

	load(paste("/gpfs/group/dxl46/default/private/poom/apc/output/", type, "/", type, "_", pheno, "_output.RData", sep = ""))

	##Create desired output tables##
	m.rose.total = rbind(m.rose, m.rose.all)
	m.ie.total = rbind(m.ie, m.ie.all)
	m.rose.total$Age = as.character(m.rose.total$Age)
	m.ie.total$Age = as.character(m.ie.total$Age)

	##Re-order data##
	m.rose.total$id = paste(m.rose.total$Age, m.rose.total$Year, m.rose.total$Cohort, sep = "_")
	m.ie.total$id = paste(m.ie.total$Age, m.ie.total$Year, m.ie.total$Cohort, sep = "_")
	m.ie.total = m.ie.total[match(m.rose.total$id, m.ie.total$id), ]

	m.rose.total = m.rose.total %>% select(-c("id"))
	m.ie.total = m.ie.total %>% select(-c("id"))

	result = data.frame()
	for ( n in 1:nrow(m.rose.total) ) {
		if ( m.rose.total$Year[n] > 2015 ) {
			print(n)
			##Calculate confidence interval via percentile interval (1000 bootstraps)##
			m.rose.sort = sort(m.rose.total[n, 5:ncol(m.rose.total)])
			m.ie.sort = sort(m.ie.total[n, 5:ncol(m.ie.total)])
			
			result_temp = as.data.frame(cbind("Period" = m.rose.total[n,1], "Cohort" = m.rose.total[n,2], "Age" = m.rose.total[n,3], 
					    "Rate_rose" = m.rose.total[n,4], "Rate_upper_rose" = m.rose.sort[0.975*(ncol(m.rose.total) - 4)], 
					    "Rate_lower_rose" = m.rose.sort[0.025*(ncol(m.rose.total) - 4)],
					    "Rate_arima" = m.ie.total[n,4], "Rate_upper_arima" = m.ie.sort[0.975*(ncol(m.rose.total) - 4)], 
					    "Rate_lower_arima" = m.ie.sort[0.025*(ncol(m.rose.total) - 4)]))
			colnames(result_temp) = c("Period", "Cohort", "Age", "Rate_rose", "Rate_upper_rose", "Rate_lower_rose",
						  "Rate_arima", "Rate_upper_arima", "Rate_lower_arima")
			result = rbind(result, result_temp)
		} else {
			print(n)
			result_temp = as.data.frame(cbind("Period" = m.rose.total[n,1], "Cohort" = m.rose.total[n,2], "Age" = m.rose.total[n,3],
	                                    "Rate_rose" = m.rose.total[n,4], "Rate_upper_rose" = NA, "Rate_lower_rose" = NA,
	                                    "Rate_arima" = m.ie.total[n,4], "Rate_upper_arima" = NA, "Rate_lower_arima" = NA), stringsAsFactors = FALSE)
	                colnames(result_temp) =	c("Period", "Cohort", "Age", "Rate_rose", "Rate_upper_rose", "Rate_lower_rose",
	                                          "Rate_arima", "Rate_upper_arima", "Rate_lower_arima")
			result = rbind(result, result_temp)
		}
	}

	result$Period = as.numeric(as.character(result$Period))
	result$Age = factor(result$Age, levels = c("7.5", "12.5", "17.5", "22.5", "27.5", "32.5", "37.5", "42.5", "47.5", "52.5", "57.5", "62.5", "67.5", "72.5", "77.5", "82.5", "All"))
	result = result[order(result$Period, result$Age), ]

	write.table(result, paste("/gpfs/group/dxl46/default/private/poom/apc/output/", type, "/", type, "_", pheno, "_output.txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

}
