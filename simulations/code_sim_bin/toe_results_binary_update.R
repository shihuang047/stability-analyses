#####################################################################################
### run all methods on simulated data with Toeplitz correlation ##################
#####################################################################################
args = commandArgs(trailingOnly=TRUE)
print(args)
dir = args[1]

library(FSA)
source('../../code_method/cv_method_binary_update.R')
source('../../code_method/getStability.R')
source('cv_sim_apply_binary_update.R')


dim.list = list()
size = c(50, 100, 500, 1000)
idx = 0
for (P in size){
    for (N in size){
        idx = idx + 1
        dim.list[[idx]] = c(P=P, N=N)
    }
}

## correlation strength
rou.list = seq(0.1, 0.9, 0.2)

files = NULL
for (rou in rou.list){
  for (dim in dim.list){
    p = dim[1]
    n = dim[2]
    files = cbind(files, paste0(dir, '/sim_toeplitz_corr', rou, paste('P', p, 'N', n, sep='_'), '.RData', sep=''))
  }
}


#------------------------
# compare methods
#------------------------
i <- as.numeric(Sys.getenv("PBS_ARRAYID"))

# ###################
# #### Lasso ########
# ###################
print('Lasso')
file_name = files[[i]]
results_toe_lasso = sim_evaluate_cv(sim_file=file_name, method='lasso', family='binomial')

when i = 61, sim_idx = 16, ROC cannot be calculated as y.test has all 0s
thus change seed number to avoid this extreme example
results_toe_lasso = sim_evaluate_cv(sim_file=file_name, method='lasso', family='binomial', seednum=5231)

save(file=paste0(dir, '/binary_update/toe_Lasso_binary_', i, '.RData'), results_toe_lasso, file_name)

# # #########################
# # #### Elastic Net ########
# # #########################
print('Elnet')
file_name = files[[i]]
results_toe_elnet = sim_evaluate_cv(sim_file=file_name, method='elnet', family='binomial')

when i = 61, sim_idx = 16, ROC cannot be calculated as y.test has all 0s
thus change seed number to avoid this extreme example
results_toe_elnet = sim_evaluate_cv(sim_file=file_name, method='elnet', family='binomial', seednum=5231)

save(file=paste0(dir, '/binary_update/toe_Elnet_binary_', i, '.RData'), results_toe_elnet, file_name)

# ############################
# #### Random Forests ########
# ############################
print('Random Forests')
file_name = files[[i]]
results_toe_rf = sim_evaluate_cv(sim_file=file_name, method='RF')

when i = 61, sim_idx = 16, ROC cannot be calculated as y.test has all 0s
thus change seed number to avoid this extreme example
results_toe_rf = sim_evaluate_cv(sim_file=file_name, method='RF', seednum=5231)

save(file=paste0(dir, '/binary_update/toe_RF_binary_', i, '.RData'), results_toe_rf, file_name)

# #########################################################
# #### Generalized compositional Lasso ########
# #########################################################
print('Generalized Compositional Lasso')
file_name = files[[i]]
results_toe_GenCompLasso = sim_evaluate_cv(sim_file=file_name, method='GenCompLasso')

save(file=paste0(dir, '/binary_update/toe_GenCompLasso_binary_', i, '.RData'),
     results_toe_GenCompLasso, file_name)

results_toe_GenCompLasso_dataSplit = sim_evaluate_cv(sim_file=file_name,
                                                     method='GenCompLasso',
                                                     data.split=TRUE)

# when i = 61, sim_idx = 16, ROC cannot be calculated as y.test has all 0s
# thus change seed number to avoid this extreme example
# results_toe_GenCompLasso_dataSplit = sim_evaluate_cv(sim_file=file_name,
#                                                      method='GenCompLasso',
#                                                      data.split=TRUE, seednum=5231)

