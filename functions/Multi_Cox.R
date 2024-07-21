Multi_Cox <- function(time, event, vars, data, ...){
    require(stringr)
    
    fit <- coxph(as.formula(paste("Surv(", time, ", ", event, ") ~", str_c(vars, collapse = "+"), sep = "")), data = data)
    fit_summ = summary(fit)
    
    df <- data.frame(
        Variable = rownames(fit_summ$conf.int), 
        mean = round(fit_summ$coefficients[, 2, drop = T], 3),
        pvalue = as.character(round(fit_summ$coefficients[, 5, drop = T], 3)), 
        lower = round(fit_summ$conf.int[, 3, drop = T], 3), 
        upper = round(fit_summ$conf.int[, 4, drop = T], 3)
    )
    return(df)
}