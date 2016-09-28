#' GE_nleqslv
#'
#' Deprecated version as of Sep 27, 2016.  Designed only to work with scalar Z and W/M.
#'
#' Here we use the package nleqslv to get a numerical solution which
#' we can use to check our direct solution.
#' 
#' @param beta_vec A vector of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' @param mu_vec A vector of the means for f, h, M, W, Z, in that order.
#' @param cov_vec A vector of the covariances in the order
#' (mu_GG, mu_GE, mu_Gf, mu_Gh, MU_GZ, MU_GM, MU_GW, mu_EE,
#' mu_Ef, MU_EZ, MU_EM, MU_EW, MU_fZ, MU_fW)
#' @param MU_ZZ Matrix mean
#' @param MU_ZM Matrix mean
#' @param MU_ZW Matrix mean
#' @param MU_WZ Matrix mean
#' @param MU_WM Matrix mean
#' @param MU_WW Matrix mean
#' @param HOM_vec A vector of the higher order moments in the order
#' (mu_GGE, mu_GGh, mu_GEE, mu_GEF, mu_GEh, MU_GEZ, MU_GEM, MU_GEW,
#' MU_GhW, MU_GhZ, mu_GGEE, mu_GGEf, mu_GGEh)
#'
#' @keywords nonlinear equation
#' @export
#' @examples 
#' GE_bias_results <- GE_bias_normal_squaredmis(runif(6), runif(3), runif(10))
#' GE_nleqslv(GE_bias_results$beta_vec, GE_bias_results$mu_vec, GE_bias_results$cov_vec,
#' GE_bias_results$MU_ZZ, GE_bias_results$MU_ZM, GE_bias_results$MU_ZW, GE_bias_results$MU_WZ,
#' GE_bias_results$MU_WM, GE_bias_results$MU_WW, GE_bias_results$HOM_vec)

scalar_GE_nleqslv <- function(beta_vec, mu_vec, cov_vec, MU_ZZ, MU_ZM, MU_ZW, MU_WM, MU_WW, HOM_vec)
{

	# Here we extract the relevant parameters from the inputs
	beta_0 <- beta_vec[1]
	beta_G <- beta_vec[2]
	beta_E <- beta_vec[3]
	beta_I <- beta_vec[4]
	BETA_Z <- beta_vec[5]
	BETA_M <- beta_vec[6]

	mu_f <- mu_vec[1]
	mu_h <- mu_vec[2]
	MU_M <- mu_vec[3]
	MU_Z <- 0
	MU_W <- 0
	
	mu_GG <- cov_vec[1]
	mu_GE <- cov_vec[2]
	mu_Gf <- cov_vec[3]
	mu_Gh <- cov_vec[4]
	MU_GZ <- cov_vec[5]
	MU_GM <- cov_vec[6]
	MU_GW <- cov_vec[7]
	mu_EE <- cov_vec[8]
	mu_Ef <- cov_vec[9]
	MU_EZ <- cov_vec[10]
	MU_EM <- cov_vec[11]
	MU_EW <- cov_vec[12]
	MU_fZ <- cov_vec[13]
	MU_fW <- cov_vec[14]
	
	mu_GGE <- HOM_vec[1]
	mu_GGh <- HOM_vec[2]
	mu_GEE <- HOM_vec[3]
	mu_GEf <- HOM_vec[4]
	mu_GEh <- HOM_vec[5]
	MU_GEZ <- HOM_vec[6]
	MU_GEM <- HOM_vec[7]
	MU_GEW <- HOM_vec[8]
	MU_GhW <- HOM_vec[9]
	MU_GhZ <- HOM_vec[10]
	mu_GGEE <- HOM_vec[11]
	mu_GGEf <- HOM_vec[12]
	mu_GGEh <- HOM_vec[13]
	
	#################################################################			
	# Define the set of score equations we will be solving
	score_eqs <- function(x)
	{
		alpha_0 <- x[1]
		alpha_G <- x[2]
		alpha_E <- x[3]
		alpha_I <- x[4]
		ALPHA_Z <- x[5]
		ALPHA_W <- x[6]
		y <- numeric(6)
		y[1] = alpha_0 + alpha_I*mu_GE + t(MU_Z) %*% ALPHA_Z + t(MU_W) %*% ALPHA_W - beta_0 - 
			beta_E*mu_f - beta_I*mu_Gh - t(MU_Z) %*% BETA_Z - t(MU_M) %*% BETA_M
	
		y[2] = alpha_G*mu_GG + alpha_E*mu_GE + alpha_I*mu_GGE + t(MU_GZ) %*% ALPHA_Z + t(MU_GW) %*% ALPHA_W - 
			beta_G*mu_GG - beta_E*mu_Gf - beta_I*mu_GGh - t(MU_GZ) %*% BETA_Z - t(MU_GM) %*% BETA_M
	
		y[3] = alpha_G*mu_GE + alpha_E*mu_EE + alpha_I*mu_GEE + t(MU_EZ) %*% ALPHA_Z + t(MU_EW) %*% ALPHA_W - 
			beta_G*mu_GE - beta_E*mu_Ef - beta_I*mu_GEh - t(MU_EZ) %*% BETA_Z - t(MU_EM) %*% BETA_M
	
		y[4] = alpha_0*mu_GE + alpha_G*mu_GGE + alpha_E*mu_GEE + alpha_I*mu_GGEE + t(MU_GEZ) %*% ALPHA_Z + 
			t(MU_GEW) %*% ALPHA_W - beta_0*mu_GE - beta_G*mu_GGE - beta_E*mu_GEf - beta_I*mu_GGEh -
			t(MU_GEZ) %*% BETA_Z - t(MU_GEM) %*% BETA_M
	
		y[5] = alpha_0*MU_Z + alpha_G*MU_GZ + alpha_E*MU_EZ + alpha_I*MU_GEZ + MU_ZZ %*% ALPHA_Z + 
			MU_ZW %*% ALPHA_W - beta_0*MU_Z - beta_G*MU_GZ - beta_E*MU_fZ - beta_I*MU_GhZ - 
			MU_ZZ %*% BETA_Z - MU_ZM %*% BETA_M
	
		y[6] = alpha_0*MU_W + alpha_G*MU_GW + alpha_E*MU_EW + alpha_I*MU_GEW + MU_ZW %*% ALPHA_Z + 
			MU_WW %*% ALPHA_W - beta_0*MU_W - beta_G*MU_GW - beta_E*MU_fW - beta_I*MU_GhW -
			MU_ZW %*% BETA_Z - MU_WM %*% BETA_M
	
		y
	}
	

	solved_scoreeqs = nleqslv::nleqslv(x=c(0,0,0,0,0,0), fn=score_eqs)
	return(solved_scoreeqs)
}


