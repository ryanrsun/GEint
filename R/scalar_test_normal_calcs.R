#' test_normal_calcs
#'
#' Deprecated as of Sep 27, 2016.  Designed to work only with scalar Z and W/M.
#'
#' Only a test function to ensure the quantities calculated in GE_bias_normal_squaredmis
#' are correct.
#'
#' @param rho_vec A vector of the 10 pairwise covariances between the
#' covariates.  These should be in the order (1) cov_GE (2) cov_GZ (3) cov_EZ
#' (4) cov_GW (5) cov_EW (6) cov_ZW.
#' @param prob_G Probability that each allele is equal to 1.  Since each SNP has
#' two alleles, the expectation of G is 2*prob_G.
#' @param test_threshold How much margin for error on tests?


scalar_GE_normal_test <- function(rho_vec, prob_G, test_threshold=0.01)
{
	surv <- function(x) {1-pnorm(x)}
	
	# Translate inputs directly
  	rho_GE <- rho_vec[1]
  	rho_GZ <- rho_vec[2]
  	rho_EZ <- rho_vec[4]
  	rho_GW <- rho_vec[3]
  	rho_EW <- rho_vec[5]
  	rho_ZW <- rho_vec[6]
  	
  	# Calculated values for generation of data (given to rmvnorm), 
  	# different from inputs for terms involving G
  	w <- qnorm(1-prob_G)				
  	r_GE <- rho_GE / (2*dnorm(w))			
    	r_GZ <- rho_GZ / (2*dnorm(w))
  	r_GW <- rho_GW / (2*dnorm(w))
  
  	# Now generate the data (assume we already know we can)
  	temp_sig <- matrix(data=c(1, 0, r_GE, r_GZ, r_GW,
								0, 1, r_GE, r_GZ, r_GW,
								r_GE, r_GE, 1, rho_EZ, rho_EW,
								r_GZ, r_GZ, rho_EZ, 1, rho_ZW,
								r_GW, r_GW, rho_EW, rho_ZW, 1), nrow=5)
  	test_data <- tryCatch(rmvnorm(n=1000000, sigma=temp_sig), 
  				warning=function(w) w, error=function(e) e)
  	if (class(test_data)[1] != 'matrix') {stop('You specified an impossible covariance matrix!')}
  	G1 <- as.numeric(test_data[,1] > w)
  	G2 <- as.numeric(test_data[,2] > w)
  	G <- G1 + G2 - 2*prob_G
	E <- test_data[,3]
	Z <- test_data[,4]
	W <- test_data[,5]
	  	
  	# Start testing with calculated G covariances
  	temp_test <- rho_GE - cov(G,E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with rho_GE')}
  	
  	temp_test <- rho_GZ - cov(G,Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with rho_GZ')}
  	
  	temp_test <- rho_GW - cov(G,W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with rho_GW')}
  	
  	# Calculate more covariances
  	mu_GE <- rho_GE
 	mu_Gf <- 2*r_GE^2*w*dnorm(w) + 2*surv(w) - 2*prob_G
 	mu_Gh <- mu_Gf
 	mu_GG <- 2*prob_G*(1-prob_G)
  	MU_GZ <- rho_GZ  	# Vector
 	MU_GW <- rho_GW		# Vector
 	MU_GM <- 	2*r_GW^2*w*dnorm(w) + 2*surv(w) - 2*prob_G	# Vector, see gen_cor_bin_normal for explanation
 	MU_EM <- 	0				# Vector, in particular because third moment of W is 0
 	MU_EZ <- rho_EZ			# Vector
 	MU_EW <- rho_EW			# Vector
 	mu_EE <- 1
 	mu_Ef <- 0
 	MU_fZ <- 	0	# Vector
 	MU_fW <- 	0		# Vector
 	 
 	# Test the above
  	temp_test <- mu_GE - mean(G*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GE')}
  	
  	temp_test <- mu_Gf - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gf')}
  	
  	temp_test <- mu_Gh - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_Gh')}
  	
  	temp_test <- mu_GG - mean(G*G)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GG')}
  	
  	temp_test <- MU_GZ - mean(G*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GZ')}
  	
  	temp_test <- MU_GW - mean(G*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GW')}
  	
  	temp_test <- MU_GM - mean(G*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GM')}
  	
  	temp_test <- MU_EM - mean(E*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_EM')}
  	
  	temp_test <- MU_EZ - mean(E*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_EZ')}
  	
  	temp_test <- MU_EW - mean(E*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_EW')}
  	
  	temp_test <- MU_fZ - mean(E*E*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_fZ')}
  	
  	temp_test <- MU_fW - mean(E*E*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_fW')}
  	
   	########################
 	# Matrix covariances
 	MU_ZW <- rho_ZW		# Matrix	 
 	MU_ZM <- 0 		# Matrix
 	MU_WM <- 0 		# Matrix
 	MU_ZZ <- 1 		# Matrix
  	MU_WW <- 1 		# Matrix
  	
  	# Check these
  	temp_test <- MU_ZW - mean(Z*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZW')}
  	
  	temp_test <- MU_ZM - mean(Z*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_ZM')}
  	
  	temp_test <- MU_WM - mean(W*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_WM')}
  	
  	########################
 	# Higher order moments
  	mu_G1_E <- r_GE*dnorm(w)
  	mu_G1_EE <- r_GE^2*w*dnorm(w) + surv(w)
  	mu_G1_EEE <- r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w)
  
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_E <- function(x,w,r_GE) {
  		x*dnorm(x)*pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_E <- quadinf(f=f_G1_G2_E, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_EE <- function(x,w,r_GE) {
  		x^2*dnorm(x)*pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_EE <- quadinf(f=f_G1_G2_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	temp_sig <- matrix(data=c(1-r_GE^2, -r_GE^2, -r_GE^2, 1-r_GE^2), nrow=2)
  	f_G1_G2_EEE <- function(x,w,r_GE) {
  		x^3*dnorm(x)*pmvnorm(lower=c(w,w), upper=c(Inf,Inf), mean=c(r_GE*x, r_GE*x), sigma=temp_sig)
  	}
  	mu_G1_G2_EEE <- quadinf(f=f_G1_G2_EEE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE)$Q[1]
  
  	f_G1_E_Z <- function(x, w, r_EZ, r_GE, r_GZ) {
  		( r_EZ * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  			(r_GZ-r_GE*r_GZ) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  	}
  	mu_G1_E_Z <- quadinf(f= f_G1_E_Z, xa=-Inf, xb=Inf, w=w, r_EZ=rho_EZ, r_GE=r_GE, r_GZ=r_GZ)$Q
  
  	f_G1_E_W <- function(x, w, r_EW, r_GE, r_GW) {
  		( r_EW * x * surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / sqrt(1-r_GE^2) ) * 
  			(r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x* dnorm(x)
  	}
  	mu_G1_E_W <- quadinf(f= f_G1_E_W, xa=-Inf, xb=Inf, w=w, r_EW=rho_EW, r_GE=r_GE, r_GW=r_GW)$Q

  	f_G1_E_WW <- function(x, w, r_GE, r_GW, r_EW) {
  		( r_EW * x* surv( (w-x*r_GW) / sqrt(1-r_GW^2) ) + dnorm( (w-r_GW*x) / 
  			sqrt(1-r_GW^2) ) * (r_GE-r_GW*r_EW) / sqrt(1-r_GW^2) ) * x^2 * dnorm(x)
  	}
  	mu_G1_E_WW <- quadinf(f=f_G1_E_WW, xa=-Inf, xb=Inf, w=w , r_GE=r_GE, r_GW=r_GW, r_EW=rho_EW)$Q
  
  	f_G1_W_EE <- function(x, w, r_GE, r_GW, r_EW) {
  		( r_EW * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  			sqrt(1-r_GE^2) ) * (r_GW-r_GE*r_EW) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  	}
  	mu_G1_W_EE <- quadinf(f=f_G1_W_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GW=r_GW, r_EW=rho_EW)$Q
  
  	f_G1_Z_EE <- function(x, w, r_GE, r_GZ, r_EZ) {
  		( r_EZ * x* surv( (w-x*r_GE) / sqrt(1-r_GE^2) ) + dnorm( (w-r_GE*x) / 
  			sqrt(1-r_GE^2) ) * (r_GZ-r_GE*r_EZ) / sqrt(1-r_GE^2) ) * x^2 * dnorm(x)
  	}	
  	mu_G1_Z_EE <- quadinf(f=f_G1_Z_EE, xa=-Inf, xb=Inf, w=w, r_GE=r_GE, r_GZ=r_GZ, r_EZ=rho_EZ)$Q  
  	
  	# Check these	
  	temp_test <- mu_G1_E - mean(G1*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E')}
  	
  	temp_test <- mu_G1_EE - mean(G1*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_EE')}
  	
  	temp_test <- mu_G1_EEE - mean(G1*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_EEE')}
  	
  	temp_test <- mu_G1_G2_E - mean(G1*G2*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_E')}
  	
  	temp_test <- mu_G1_G2_EE - mean(G1*G2*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_EE')}
  	
  	temp_test <- mu_G1_G2_EEE - mean(G1*G2*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_G2_EEE')}
  	
  	temp_test <- mu_G1_E_Z - mean(G1*E*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_Z')}
  	
  	temp_test <- mu_G1_E_W - mean(G1*E*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_W')}
  	
  	temp_test <- mu_G1_E_WW - mean(G1*E*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_E_WW')}
  	
  	temp_test <- mu_G1_W_EE - mean(G1*W*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_W_EE')}
  	
  	temp_test <- mu_G1_Z_EE - mean(G1*Z*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_G1_Z_EE')}
  	
  	
  	# Last higher order terms
  	mu_GGE <- 2*mu_G1_E + 2*mu_G1_G2_E - 8*prob_G*mu_G1_E
  	mu_GGh <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  	mu_GEE <- mu_Gf
  	mu_GEf <- 2*(r_GE^3*w^2*dnorm(w) - r_GE^3*dnorm(w) + 3*r_GE*dnorm(w))
  	mu_GEh <- mu_GEf

  	MU_GEZ <- 2*mu_G1_E_Z - 2*prob_G*rho_EZ			# Vector
  	MU_GEW <- 2*mu_G1_E_W	- 2*prob_G*rho_EW		# Vector
  	MU_GEM <-	2*mu_G1_E_WW				# Vector
  	MU_GhW <- 2*mu_G1_W_EE
  	MU_GhZ <- 2*mu_G1_Z_EE
  
  	mu_GGEE <- 2*mu_G1_EE + 2*mu_G1_G2_EE + 4*prob_G^2*1 - 8*prob_G*mu_G1_EE
  	mu_GGEf <- 2*mu_G1_EEE + 2*mu_G1_G2_EEE + 4*prob_G^2*0 - 8*prob_G*mu_G1_EEE
  	mu_GGEh <- mu_GGEf

  	# Test
  	temp_test <- mu_GGE - mean(G*G*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGE')}
  	
  	temp_test <- mu_GGh - mean(G*G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGj')}
  	
  	temp_test <- mu_GEE - mean(G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEE')}
  	
  	temp_test <- mu_GEf - mean(G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEf')}
  	
  	temp_test <- mu_GEh - mean(G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GEh')}
  	
  	temp_test <- MU_GEZ - mean(G*E*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEZ')}
  	
  	temp_test <- MU_GEW - mean(G*E*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEW')}
  	
  	temp_test <- MU_GEM - mean(G*E*W*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with MU_GEM')}
  	
  	temp_test <- MU_GhW - mean(G*E*E*W)
  	if (abs(temp_test) > test_threshold) {warning('Problem with MU_GhW')}
  	
  	temp_test <- MU_GhZ - mean(G*E*E*Z)
  	if (abs(temp_test) > test_threshold) {warning('Problem with MU_GhZ')}
  	
  	temp_test <- mu_GGEE - mean(G*G*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGEE')}
  	
  	temp_test <- mu_GGEf - mean(G*G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGEf')}
  	
  	temp_test <- mu_GGEh - mean(G*G*E*E*E)
  	if (abs(temp_test) > test_threshold) {warning('Problem with mu_GGEj')}
  	
  	cat('Done')
}




