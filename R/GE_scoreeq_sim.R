#' GE_scoreeq_sim.R
#'
#' Here we perform simulation to verify that we have solved for
#' the correct alpha values.
#' 
#' @param num_sims The number of repetitions to run, default is 5000
#' @param num_sub The number of subjects to simulate, default is 2000
#' @param beta_vec A vector of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' @param mu_vec A vector of the means for M, W, Z, in that order.
#' @param rho_vec A vector of the covariances in the order:
#' (rho_GE, rho_GZ, rho_GM, rho_GW, rho_EZ, rho_EM, rho_EW, rho_ZM, rho_ZW, rho_WM)
#'
#' @keywords simulation
#' @export
#' @examples 
#' library(mvtnorm)
#' temp_sig <- matrix(data=0.5, nrow=5, ncol=5)
#' diag(temp_sig) <- 1
#' sig_mat <- cor( rmvnorm(n=10, sigma=temp_sig) )
#' rho_vec <- sig_mat[lower.tri(sig_mat)]
#' mu_vec <- runif(3)
#' beta_vec <- runif(6)
#' GE_bias_results <- GE_bias_normal_squaredmis(beta_vec, mu_vec, rho_vec)
#' GE_nleqslv(GE_bias_results$beta_vec, GE_bias_results$mu_vec, GE_bias_results$cov_vec,
#' GE_bias_results$MU_ZZ, GE_bias_results$MU_ZM, GE_bias_results$MU_ZW, GE_bias_results$MU_WZ,
#' GE_bias_results$MU_WM, GE_bias_results$MU_WW, GE_bias_results$HOM_vec)
#' GE_scoreeq_sim(beta_vec=beta_vec, mu_vec=mu_vec, rho_vec=rho_vec)

GE_scoreeq_sim <- function(num_sims=5000, num_sub=2000, beta_vec, mu_vec, rho_vec)
{
	results <- matrix(data=NA, nrow=num_sims, ncol=6)
	sig_mat <- matrix(data=1, nrow=5, ncol=5)
	sig_mat[1,2] = sig_mat[2,1] = rho_vec[1]
	sig_mat[1,3] = sig_mat[3,1] = rho_vec[2]
	sig_mat[1,4] = sig_mat[4,1] = rho_vec[3]
	sig_mat[1,5] = sig_mat[5,1] = rho_vec[4]
	sig_mat[2,3] = sig_mat[3,2] = rho_vec[5]
	sig_mat[2,4] = sig_mat[4,2] = rho_vec[6]
	sig_mat[2,5] = sig_mat[5,2] = rho_vec[7]
	sig_mat[3,4] = sig_mat[4,3] = rho_vec[8]
	sig_mat[3,5] = sig_mat[5,3] = rho_vec[9]
	sig_mat[4,5] = sig_mat[5,4] = rho_vec[10]
	
	beta_0 <- beta_vec[1]
	beta_G <- beta_vec[2]
	beta_E <- beta_vec[3]
	beta_I <- beta_vec[4]
	BETA_Z <- beta_vec[5]
	BETA_M <- beta_vec[6]
	
	for(i in 1:num_sims)
	{
		sim_data = mvtnorm::rmvnorm(n=num_sub, mean=c(0,0,mu_vec[1], mu_vec[2], mu_vec[3]), sigma=sig_mat)
		G <- sim_data[,1]
		E <- sim_data[,2]
		Z <- sim_data[,3]
		M <- sim_data[,4]
		W <- sim_data[,5]
	
		Y = beta_0 + beta_G*G + beta_E*E^2 + beta_I*G*E^2 + BETA_Z*Z + BETA_M*M + rnorm(n=num_sub)
	
		temp_mod <- lm(Y~G+E+G*E+Z+W)
		results[i,1:3] <- summary(temp_mod)$coefficients[1:3,1]
		results[i,4] <- summary(temp_mod)$coefficients[6,1]
		results[i,5:6] <- summary(temp_mod)$coefficients[4:5,1]
	}
	colnames(results) <- c('alpha_0', 'alpha_G', 'alpha_E', 'alpha_I', 'ALPHA_Z', 'ALPHA_W')
	sim_alpha <- apply(results, 2, mean)
	
	return(sim_alpha)
}

