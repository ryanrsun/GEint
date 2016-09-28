# Run some tests on GEint
setwd('/users/ryansun/documents/research/paper2/software/geint/R')
source('GE_bias_normal_squaredmis.R')
source('GE_bias.R')
source('GE_enumerate_inputs.R')
source('GE_nleqslv.R')
source('GE_scoreeq_sim.R')
source('GE_test_moment_calcs.R')
source('GE_translate_inputs.R')


num_Z <- rpois(n=1, lambda=1)
if (num_Z==0) {num_Z <- 1}
num_W <- rpois(n=1, lambda=1)
if (num_W==0) {num_W <- 1}
beta_list <- list( runif(n=1, min=0.1, max=1),
				runif(n=1, min=0.1, max=1),
				runif(n=1, min=0.1, max=1),
				runif(n=1, min=0.1, max=1),
				runif(n=num_Z, min=0.1, max=1),
				runif(n=num_W, min=0.1, max=1)	)

rho_list <- list( runif(n=1, min=0.02, max=0.3),
					 runif(n=num_Z, min=0.02, max=0.3),
					  runif(n=num_Z, min=0.02, max=0.3),
					   runif(n=num_W, min=0.02, max=0.3),
					    runif(n=num_W, min=0.02, max=0.3),
					     runif(n=num_W*num_Z, min=0.02, max=0.3) )

prob_G = runif(n=1, min=0.05, max=0.95)

if (num_Z > 1) {
	temp <- num_Z*(num_Z-1) / 2
	cov_Z <- runif(n=temp, min=0.02, max=0.2)
} else {
	cov_Z <- NULL
}

if (num_W > 1) {
	temp <- num_W*(num_W-1) / 2
	cov_W <- runif(n=temp, min=0.02, max=0.2)
} else {
	cov_W <- NULL
}

# Test first
GE_test_moments_calcs(beta_list, rho_list, prob_G, cov_Z, cov_W)
sim_results <- GE_scoreeq_sim(beta_list=beta_list, prob_G=prob_G, rho_list=rho_list, cov_Z=cov_Z, cov_W=cov_W)
solve_results <- GE_bias_normal_squaredmis(beta_list, rho_list, prob_G, cov_Z, cov_W)

mu_list <- list( solve_results$mu_list[[1]], solve_results$mu_list[[2]],	
				solve_results$mu_list[[3]],rep(0, num_Z), rep(0, num_W) )

cov_list <- solve_results$cov_list
HOM_list <- solve_results$HOM_list
cov_mat_list <- list( MU_ZZ=solve_results$MU_ZZ, MU_WW=solve_results$MU_WW, MU_ZW=solve_results$MU_ZW,
					MU_WZ=solve_results$MU_WZ, MU_ZM=solve_results$MU_ZM, MU_WM=solve_results$MU_WM)

nleqslv_results <- GE_nleqslv(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)

ge_bias_results <- GE_bias(beta_list, cov_list, cov_mat_list, mu_list, HOM_list)


# Should all match
nleqslv_results$x
sim_results
solve_results$alpha_list

