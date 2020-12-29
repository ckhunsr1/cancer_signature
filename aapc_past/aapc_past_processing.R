library(dplyr)
library(readxl)

##Incidence AAPC calculation for old time##
data <- read_xlsx(path = "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/Poom_analysis/SEER_database.xlsx", col_names = FALSE, sheet = 1)
data = as.data.frame(data)
colnames(data) = paste("V", 1:ncol(data), sep = "")

##Filter for cancer types##
name_data = (data[rowSums(is.na(data)) == ncol(data) - 1, ] %>% filter(V1 != "~"))$V1

##Filter for clean rows##
clean_data = data[rowSums(is.na(data)) == 0, ]
clean_data = clean_data %>% filter(V1 != "1975-2016")
gap = nrow(clean_data)/length(name_data)

name = c("2", "14", "15", "17", "32", "33", "37", "44", "50", "52", "54", "56", "58", "63", "64", "68", "69", "77", "80", "83", "87")
for ( idx in name ) {
	idx = as.numeric(as.character(idx))
	my_data = clean_data[(gap*(idx-1) + 1):(gap*idx), ] %>% select(-c("V1"))
	my_data = mutate_all(my_data, function(x) as.numeric(as.character(x)))

	#Calculate rate for all age group##
	result = data.frame()
	for ( i in 1:nrow(my_data) ) {
		result_temp = data.frame("Sex" = 0, "Year of diagnosis - 1975-2016" = i-1, "Crude Rate" = 100000*(sum(my_data[i, seq(1, ncol(my_data) - 1, by = 2)]))/(sum(my_data[i, seq(2, ncol(my_data), by = 2)])) )
		result = rbind(result, result_temp)
	}
	write.table(result, paste("/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/Poom_analysis/input/Incidence_aapc_old/Incidence_aapc_", idx, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

####################################################################################################################################################################################################################################################################
##Mortality AAPC calculation for old time##
data <- read_xlsx(path = "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/Poom_analysis/SEER_database.xlsx", col_names = FALSE, sheet = 3)
data = as.data.frame(data)
colnames(data) = paste("V", 1:ncol(data), sep = "")

##Filter for cancer types##
name_data = (data[rowSums(is.na(data)) == ncol(data) - 1, ] %>% filter(V1 != "^"))$V1

##Filter for clean rows##
clean_data = data[rowSums(is.na(data)) == 0, ]
clean_data = clean_data %>% filter(V1 != "1975-2016")
gap = nrow(clean_data)/length(name_data)

name = c("3", "15", "16", "18", "22", "23", "27", "34", "40", "42", "44", "46", "48", "53", "54", "58", "59", "65", "68", "69", "71")
for ( idx in name ) {
	idx = as.numeric(as.character(idx))
	my_data = clean_data[(gap*(idx-1) + 1):(gap*idx), ] %>% select(-c("V1"))
	my_data = data.frame(lapply(my_data, function(x) gsub("\\^", "0", x)))
	my_data = mutate_all(my_data, function(x) as.numeric(as.character(x)))

	#Calculate rate for all age group##
	result = data.frame()
	for ( i in 1:nrow(my_data) ) {
		result_temp = data.frame("Sex" = 0, "Year of diagnosis - 1975-2016" = i-1, "Crude Rate" = 100000*(sum(my_data[i, seq(1, ncol(my_data) - 1, by = 2)]))/(sum(my_data[i, seq(2, ncol(my_data), by = 2)])) )
		result = rbind(result, result_temp)
	}
	write.table(result, paste("/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/Poom_analysis/input/Mortality_aapc_old/Mortality_aapc_", idx, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

####################################################################################################################################################################################################################################################################
