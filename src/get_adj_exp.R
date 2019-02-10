args = commandArgs(trailingOnly = TRUE)

fn_exp = args[1] # '../processed_data/transposed_expression/Adipose_experssion.txt'
fn_covariate = args[2] #'../processed_data/transposed_covariate/Adipose_covariate.txt'
fn_save = args[3] #'../processed_data/adj_expression/Adipose_adj-experssion.txt'


adjust_for_covariates <- function(expression_vec, cov_df) {
## adapted from https://github.com/hakyimlab/PredictDB_Pipeline_GTEx_v7/blob/master/model_training/scripts/gtex_v7_elasticnet.R#L46
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}


df_expr <- read.table(fn_exp, header = T, row.names = 1, sep = ',')
df_cov <- read.table(fn_covariate, header = TRUE, row.names = 1, sep = ',')



n_genes = ncol(df_expr)
n_sample = nrow(df_expr)


res = matrix(0, n_sample, n_genes)
for (i in 1:n_genes){
  expression_vec <- df_expr[ ,i]
  adj_expression <- adjust_for_covariates(expression_vec, df_cov)
  res[, i] <- adj_expression
}


write.table(res, file = fn_save, quote = FALSE, row.names=rownames(df_expr), sep = ',', col.names = colnames(df_expr))