#' GE_enumerate_inputs.R
#'
#' Call this function to display the necessary inputs for GE_bias_set. For terms in 
#' cov_mat_list like E[ZZ], these should be matrices where the (i,j) element is E[Z_i*Z_j].
#'
#' @return Nothing
#'
#' @export
#' @examples 
#' GE_enumerate_inputs()

GE_enumerate_inputs <- function()
{
	beta_help <- paste("beta_list should contain the following 6 inputs (in order):",
					"beta_0",
					"beta_G",
					"beta_E",
					"beta_I",
					"beta_Z",
					"beta_M", sep='\n')
	
	mu_help <- paste("mu_list should contain the following 5 inputs (in order):",
					"mu_f",
					"mu_h",
					"mu_Z",
					"mu_M",
					"mu_W", sep='\n')
	
	cov_help <- paste("cov_list should contain the following 14 inputs (in order):",
					"mu_GE",
					"mu_Gf",
					"mu_Gh",
					"mu_EE",
					"mu_Ef",
					"mu_EZ",
					"mu_EM",
					"mu_EW",
					"MU_fZ",
					"MU_fW", sep='\n')
					
	mat_help <- paste("cov_mat_list should be a list of the matrices:",
	        "mu_GG",
	        "mu_GZ",
	        "mu_GM",
	        "mu_GW",
					"mu_ZZ",
					"mu_ZW",
					"mu_WM",
					"mu_WW",
					"mu_WM", sep='\n')
					
	HOM_help <- paste("HOM_list should be a list of the following higher order moments:",
					"mu_GEG",
					"mu_GhG",
					"mu_GEE",
					"mu_GEf",
					"mu_GEh",
					"MU_GEZ",
					"MU_GEM",
					"MU_GEW",
					"MU_GhW",
					"MU_GhZ",
					"mu_GEEG",
					"mu_GEfG",
					"mu_GEhG", sep='\n')
	
	cat(beta_help, '\n\n')
	cat(mu_help, '\n\n')
	cat(cov_help, '\n\n')
	cat(mat_help, '\n\n')
	cat(HOM_help, '\n\n')
}

