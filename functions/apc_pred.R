apc_pred <- function(input, output) {
        ##Extract CRR from output##
        output$CohortRR = as.data.frame(output$CohortRR)
        output$CohortRR$logRR = log(output$CohortRR$`Rate Ratio`)

        ##Plot original CRR data##
        #ggplot(output$CohortRR, aes(x = Cohort, y = logRR)) + geom_line()

        ##Create linear model and perform joinpoint linear piecewise regression##
        lm = lm(logRR ~ Cohort, data = output$CohortRR)

	##Predict CRR into the future##
	year_start = cohort_range[length(cohort_range)] + interval
        year_end = cohort_range[length(cohort_range)] + (interval*n_pred)

	my.seg = tryCatch({
			segmented(lm, seg.z = ~Cohort, psi = NA, control = seg.control(n.boot = 50, it.max = 100, fix.npsi=FALSE, seed = bs))
			}, error = function(err) {
			message(" Joinpoint fails and linear regression is used instead.")
			return(lm)
                        })
	my.fitted <- fitted(my.seg)
	my.model <- data.frame(Cohort = cohort_range, logRR = my.fitted)

	if ( length(my.seg$coefficients) > 2 ){
	        crr_pred = data.frame("Cohort" = seq(year_start, year_end, by = interval),
        	                      "crr" = exp(predict.segmented(my.seg, data.frame("Cohort" = seq(year_start, year_end, by = interval)))))
	} else {
		crr_pred = data.frame("Cohort" = seq(year_start, year_end, by = interval),
                              "crr" = exp(predict.lm(my.seg, data.frame("Cohort" = seq(year_start, year_end, by = interval)))))
	}

        ##Create CRR dataframe that spans into the future##
        crr_old = output$CohortRR[,1:2]
        colnames(crr_old) = c("Cohort", "crr")
        crr_total = rbind(crr_old, crr_pred)

        ##Incidence rate prediction in future period##
        result = data.frame()
        period_start = period_range[length(period_range)] + interval
        period_end = period_range[length(period_range)] + (interval*n_pred)

        for ( year in seq(period_start, period_end, by = interval) ) {
                for ( i in 1:length(age_range) ) {
                        temp = data.frame("Year" = year, "Cohort" = year - age_range[i], "Age" = age_range[i],
                                        "expected_rate" = output$LongAge[i,2] * crr_total$crr[which(crr_total$Cohort == year - age_range[i])])
                        result = rbind(result, temp)
                }
        }

	##Fitted observed incidence rate according to input##
        result_old = data.frame()
        ##Process input data##
        for ( i in 1:length(age_range) ) {
                age = age_range[i]
                for (j in 1:ncol(input$events) ) {
                        year = period_range[j]
                        temp = data.frame("Year" = year,  "Cohort" = year - age, "Age" = age,
                                          "expected_rate" = 100000*output$FittedRates$events[i,j]/output$FittedRates$offset[i,j] )
                        result_old = rbind(result_old, temp)
                }
        }

	result_total = rbind(result_old, result)
        result_total$Age = as.factor(result_total$Age)
        return(result_total)
}
