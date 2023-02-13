plot.cov.betas <- function(m = m, output_wd = Output_wd, modelType = modelType){
  n.cov <- length(m$covNames) # Number of covariates without the intercept
  var.code <- vector()
  for (i in 1:n.cov){
    var.code[i] <- paste0("C", i)
  }
  
  var.name <- as.vector(m$covNames[1:n.cov])
  
  predictors <- as.data.frame(cbind(var.code, var.name))
  
  for (i in 1:nrow(predictors)){
    pdf(paste0(Output_wd, "/Betas_covariates_coef_plot_", 
               var.name[i], "_",modelType, '.pdf'))
    MCMCplot(mpost$Beta,
             params = predictors[i,1],
             ISB = FALSE,
             ref_ovl = TRUE,
             rank = FALSE,
             xlab = 'ESTIMATE',
             main = predictors[i,2],
             labels_sz = 0.5,
             med_sz = 1,
             thick_sz = 1,
             thin_sz = 1,
             ax_sz = 1,
             main_text_sz = 1)
    
    dev.off()
  }
}
