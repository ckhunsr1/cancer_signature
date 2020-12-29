is.sequential <- function(x){
  all(diff(x) == rep(1,length(x)-1))
}    
library(readxl)

################################################################################################################################################################
##For incidence rate##
################################################################################################################################################################
data <- read_xlsx(path = "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/SEER_database.xlsx", col_names = FALSE, sheet = 1)
data = as.data.frame(data)
colnames(data) = paste("V", 1:ncol(data), sep = "")

##Filter for cancer types##
name_data = (data[rowSums(is.na(data)) == ncol(data) - 1, ] %>% filter(V1 != "~"))$V1

##Filter for clean rows##
clean_data = data[rowSums(is.na(data)) == 0, ]
clean_data = clean_data %>% filter(V1 != "1975-2016")

##Example input file to apc function (we will use it as a template)##
input = readRDS("/Users/chachritkhun/Downloads/APC_analysis_20201001_121449.495461_input.rds") ##input to apc##
age_range = input$ages
gap = nrow(clean_data)/length(name_data)

for ( idx in 1:length(name_data) ) {
	input$name = paste("US", name_data[idx], "cancer" ,sep = " ")
  	my_data = clean_data[(gap*(idx-1) + 1):(gap*idx), ] %>% select(-c("V1"))
	my_data = mutate_all(my_data, function(x) as.numeric(as.character(x)))
    
	##Prepare for the input file##
	interval = 5
	new_data = data.frame()
	##Reorganize data into format that apc can accept as input##
	for ( i in 1:floor(nrow(my_data)/interval) ) {
		data_ss = my_data[(interval*(i-1) + 1):(interval*(i)), ]
		temp = colSums(data_ss)
		new_data = rbind(new_data, temp)
	}

	format_data = data.frame(matrix(NA, nrow = ncol(new_data)/2, ncol = 2*nrow(new_data)))
	for ( j in 1:nrow(new_data) ) {
		for ( k in 1:(ncol(new_data)/2) ) {
			format_data[k, (2*j) - 1] = new_data[j, (2*k) - 1]
			format_data[k, (2*j)] = new_data[j, (2*k)]
		}
	}
    
	input$events = as.matrix(format_data[, seq(1, 15, by = 2)])
        row_include = which(rowSums(is.na(input$events)) == 0)

        if ( length(row_include) > 0 ) {
                if ( is.sequential(row_include) ) {

                        input$ages = seq(interval*row_include[1], interval*(1+row_include[length(row_include)]), by = interval)
                        input$events = input$events[row_include, ]
                        colnames(input$events) = NULL

                        input$offset = as.matrix(format_data[, seq(2, 16, by = 2)])
                        input$offset = input$offset[row_include, ]
                        colnames(input$offset) = NULL
                        save(input, file = paste("/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/input/Incidence/Incidence_", idx, ".RData", sep = ""))
                        print(paste(idx, "-Success with ", length(input$ages) - 1, " groups", sep = ""))
                } else {
                        print(paste(idx, "-", "Not sequential", sep = ""))
                }
        } else {
                print(paste(idx, "-", "No value", sep = ""))
        }
}

map = as.data.frame(cbind("Filename" = paste("Incidence", 1:length(name_data), sep = "_"), "Cancer" = name_data))
write.table(map, "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/input/Incidence_map.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

################################################################################################################################################################
##For mortality rate##
################################################################################################################################################################

data <- read_xlsx(path = "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/SEER_database.xlsx", col_names = FALSE, sheet = 3)
data = as.data.frame(data)
colnames(data) = paste("V", 1:ncol(data), sep = "")

##Filter for cancer types##
name_data = (data[rowSums(is.na(data)) == ncol(data) - 1, ] %>% filter(V1 != "^"))$V1

##Filter for clean rows##
clean_data = data[rowSums(is.na(data)) == 0, ]
clean_data = clean_data %>% filter(V1 != "1975-2016")

##Example input file to apc function (we will use it as a template)##
input = readRDS("/Users/chachritkhun/Downloads/APC_analysis_20201001_121449.495461_input.rds") ##input to apc##
age_range = input$ages
gap = nrow(clean_data)/length(name_data)

for ( idx in 1:length(name_data) ) {
	input$name = paste("US", name_data[idx], "cancer" ,sep = " ")
	my_data = clean_data[(gap*(idx-1) + 1):(gap*idx), ] %>% select(-c("V1"))
	my_data	= data.frame(lapply(my_data, function(x) gsub("\\^", "0", x)))
	my_data = mutate_all(my_data, function(x) as.numeric(as.character(x)))
  
	##Prepare for the input file##
	interval = 5
	new_data = data.frame()
  
	##Reorganize data into format that apc can accept as input##
	for ( i in 1:floor(nrow(my_data)/interval) ) {
		data_ss = my_data[(interval*(i-1) + 1):(interval*(i)), ]
		temp = colSums(data_ss)
		new_data = rbind(new_data, temp)
	}
  
	format_data = data.frame(matrix(NA, nrow = ncol(new_data)/2, ncol = 2*nrow(new_data)))
	for ( j in 1:nrow(new_data) ) {
		for ( k in 1:(ncol(new_data)/2) ) {
			format_data[k, (2*j) - 1] = new_data[j, (2*k) - 1]
			format_data[k, (2*j)] = new_data[j, (2*k)]
		}
	}
  
	input$events = as.matrix(format_data[, seq(1, 15, by = 2)])
	row_include = which(rowSums(is.na(input$events)) == 0)

	if ( length(row_include) > 0 ) {
	        if ( is.sequential(row_include) ) {
        
		        input$ages = seq(interval*row_include[1], interval*(1+row_include[length(row_include)]), by = interval)
	                input$events = input$events[row_include, ]
	                colnames(input$events) = NULL

	                input$offset = as.matrix(format_data[, seq(2, 16, by = 2)])
	                input$offset = input$offset[row_include, ]
	                colnames(input$offset) = NULL
	                save(input, file = paste("/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/input/Mortality/Mortality_", idx, ".RData", sep = ""))
	                print(paste(idx, "-Success with ", length(input$ages) - 1, " groups", sep = ""))
		} else {
			print(paste(idx, "-", "Not sequential", sep = ""))
		}
        } else {
                print(paste(idx, "-", "No value", sep = ""))
        }
}

map = as.data.frame(cbind("Filename" = paste("Mortality", 1:length(name_data), sep = "_"), "Cancer" = name_data))
write.table(map, "/Users/chachritkhun/Desktop/Liu_Labwork/cancer/cancer_signature/input/Mortality_map.txt",
            col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)


