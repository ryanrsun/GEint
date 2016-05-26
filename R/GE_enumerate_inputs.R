#' GE_enumerate_inputs.R
#'
#' Call this function only to display the necessary inputs for GE_bias.
#'
#' @keywords inputs
#' @export
#' @examples 
#' GE_enumerate_inputs()

GE_enumerate_inputs <- function()
{
	beta_help <- paste("beta_vec should contain the following 6 inputs (in order):",
					"beta_0",
					"beta_G",
					"beta_E",
					"beta_I",
					"BETA_Z",
					"BETA_M", sep='\n')

	mu_help <- paste("mu_vec should contain the following 5 inputs (in order):",
					"mu_f",
					"mu_h",
					"MU_Z",
					"MU_M",
					"MU_W", sep='\n')
					
	cov_help <- paste("cov_vec should contain the following 14 inputs (in order):",
					"mu_GG",
					"mu_GE",
					"mu_Gf",
					"mu_Gh",
					"MU_GZ",
					"MU_GM",
					"MU_GW",
					"mu_EE",
					"mu_Ef",
					"MU_EZ",
					"MU_EM",
					"MU_EW",
					"MU_fZ",
					"MU_fW", sep='\n')
					
	mat_help <- paste("You will also need to know the expectations of the following random matrices:",
					"MU_ZZ = Z%*%t(Z)",
					"MU_ZM",
					"MU_ZW",
					"MU_WZ",
					"MU_WM",
					"MU_WW", sep='\n')
					
	HOM_help <- paste("You will need to know the following higher order moments:",
					"mu_GGE",
					"mu_GGh",
					"mu_GEE",
					"mu_GEf",
					"mu_GEh",
					"MU_GEZ",
					"MU_GEM",
					"MU_GEW",
					"MU_GhW",
					"MU_GhZ",
					"mu_GGEE",
					"mu_GGEf",
					"mu_GGEh", sep='\n')
					
	cat(beta_help, '\n\n')
	cat(mu_help, '\n\n')
	cat(cov_help, '\n\n')
	cat(mat_help, '\n\n')
	cat(HOM_help, '\n\n')
}

