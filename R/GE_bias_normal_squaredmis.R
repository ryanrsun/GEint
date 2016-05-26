#' GE_bias_normal_squaredmis
#'
#' A function to calculate the bias in testing for GxE interaction.
#' Here we make the following assumptions:
#' (1) All covariates are marginally normally distributed with 
#' unit variance.  In addition, G and E have mean 0.
#' (2) The misspecification is of the form f(E)=h(E)=E^2
#' (3) All covariate terms exist and are scalars
#' 
#' @param beta_vec A vector of the effect sizes in the true model.
#' Use the order beta_0, beta_G, beta_E, beta_I, beta_Z, beta_M.
#' @param means_vec A vector of the means for M, W, Z, in that order.
#' @param rho_vec A vector of the 10 pairwise correlations between the
#' covariates.  Suppose that the order of the rows/columns is G,E,Z,M,W.
#' Then the correlation in row i, column j should be the (i-1)*(i-2)/2 + j
#' element in the vector.
#' 
#' @keywords bias normal squared_misspecification
#' @export
#' @examples 
#' GE_bias_normal_squaredmis( runif(6), runif(3), runif(10) )

GE_bias_normal_squaredmis <- function(beta_vec, means_vec, rho_vec)
{
  # First, make share we got good inputs
  if (!is.numeric(beta_vec) | !is.numeric(means_vec) | !is.numeric(rho_vec))
  {
    stop('Nonnumeric inputs!')
  }
  
  if (length(beta_vec) != 6 | length(means_vec) != 3 | length(rho_vec) != 10)
  {
    stop('Input vectors not the right size!')
  }
  
  # Translate inputs to quantities for calculation
  rho_GE <- rho_vec[1]
  rho_GZ <- rho_vec[2]
  rho_GM <- rho_vec[3]
  rho_GW <- rho_vec[4]
  rho_EZ <- rho_vec[5]
  rho_EM <- rho_vec[6]
  rho_EW <- rho_vec[7]
  rho_ZM <- rho_vec[8]
  rho_ZW <- rho_vec[9]
  rho_WM <- rho_vec[10]
  
  beta_0 <- beta_vec[1]
  beta_G <- beta_vec[2]
  beta_E <- beta_vec[3]
  beta_I <- beta_vec[4]
  BETA_Z <- beta_vec[5]  	
  BETA_M <- beta_vec[6]
  
  MU_M <- means_vec[1]  		
  MU_W <- means_vec[2]			
  MU_Z <- means_vec[3]			
  mu_f <- 1
  mu_h <- 1
  
  # Done translating inputs
  # Now calculate other necessary terms that have 
  # been determined by out assumptions + inputs
  
  ########################
  # Covariances
  mu_GE <- rho_GE
  mu_Gf <- 0
  mu_Gh <- 0
  mu_GG <- 1
  MU_GZ <- rho_GZ  	# Vector
  MU_GW <- rho_GW		# Vector
  MU_GM <- rho_GM		# Vector
  MU_EM <- rho_EM		# Vector
  MU_EZ <- rho_EZ			# Vecotr
  MU_EW <- rho_EZ			# Vector
  mu_EE <- 1
  mu_Ef <- 0
  MU_fZ <- MU_Z			# Vector
  MU_fW <- MU_W			# Vector
  
  ########################
  # Matrix covariances
  MU_ZW <- rho_ZW + MU_Z*MU_W		# Matrix	 
  MU_WZ <- rho_ZW + MU_Z*MU_W		# Matrix
  MU_ZM <- rho_ZM + MU_Z*MU_M		# Matrix
  MU_WM <- rho_WM + MU_W*MU_M		# Matrix
  MU_ZZ <- 1 + MU_Z^2			# Matrix
  MU_WW <- 1 + MU_W^2		# Matrix
  
  
  ########################
  # Higher order moments
  mu_GGE <- 0
  mu_GGh <- 1 - rho_GE^2 + 3*rho_GE^2
  mu_GEE <- 0
  mu_GEf <- rho_GE^3*3 + 3*rho_GE*(1-rho_GE^2)
  mu_GEh <- rho_GE^3*3 + 3*rho_GE*(1-rho_GE^2)
  MU_GEZ <- ( (rho_GE-rho_EZ*rho_GZ)*MU_Z + MU_Z*rho_GZ*(rho_EZ-rho_GE*rho_GZ) ) / (1-rho_GZ^2)			# Vector
  MU_GEW <- ( (rho_GE-rho_EW*rho_GW)*MU_W + MU_W*rho_GW*(rho_EW-rho_GE*rho_GW) ) / (1-rho_GW^2)			# Vector
  MU_GEM <- ( (rho_GE-rho_EM*rho_GM)*MU_M + MU_M*rho_GM*(rho_EM-rho_GE*rho_GM) ) / (1-rho_GM^2)			# Vector
  
  MU_GhW <- (1/(1-rho_GW^2))^2 * ( 3*rho_GW*(rho_GE-rho_EW*rho_GW)^2 + 
                                     (3*MU_W^2*rho_GW+3*rho_GW^3+3*rho_GW*(1-rho_GW^2)) * (rho_EW - rho_GE*rho_GW)^2 -
                                     4*MU_W^2*rho_GW*(rho_EW-rho_GE*rho_GW)^2 + MU_W^2*rho_GW*(rho_EW-rho_GE*rho_GW)^2 +
                                     2*(1-rho_GW^2+MU_W^2+3*rho_GW^2)*(rho_GE-rho_EW*rho_GW)*(rho_EW-rho_GE*rho_GW) - 
                                     2*MU_W^2*(rho_GE-rho_EW*rho_GW)*(rho_EW-rho_GE*rho_GW) ) + rho_GW - 
    (rho_GW/(1-rho_GW^2)) * (rho_GE^2+rho_EW^2-2*rho_GE*rho_EW*rho_GW)
  
  MU_GhZ <- (1/(1-rho_GZ^2))^2 * ( 3*rho_GZ*(rho_GE-rho_EZ*rho_GZ)^2 + 
                                     (3*MU_Z^2*rho_GZ+3*rho_GZ^3+3*rho_GZ*(1-rho_GZ^2)) * (rho_EZ - rho_GE*rho_GZ)^2 -
                                     4*MU_Z^2*rho_GZ*(rho_EZ-rho_GE*rho_GZ)^2 + MU_Z^2*rho_GZ*(rho_EZ-rho_GE*rho_GZ)^2 +
                                     2*(1-rho_GZ^2+MU_Z^2+3*rho_GZ^2)*(rho_GE-rho_EZ*rho_GZ)*(rho_EZ-rho_GE*rho_GZ) - 
                                     2*MU_Z^2*(rho_GE-rho_EZ*rho_GZ)*(rho_EZ-rho_GE*rho_GZ) ) + rho_GZ - 
    (rho_GZ/(1-rho_GZ^2)) * (rho_GE^2+rho_EZ^2-2*rho_GE*rho_EZ*rho_GZ)
  
  mu_GGEE <- 1 - rho_GE^2 + 3*rho_GE^2
  mu_GGEf <- 0
  mu_GGEh <- 0
  
  
  ########################
  # Some shortcut quantities
  A <- (mu_GE * MU_GZ / mu_GG - MU_EZ) / (mu_EE - mu_GE^2/mu_GG)
  B <- (mu_GE * MU_GW / mu_GG - MU_EW) / (mu_EE - mu_GE^2/mu_GG)
  
  O <- MU_Z%*%t(MU_Z) + MU_GZ%*%t(MU_GZ)/mu_GG - MU_ZZ - A %*% t(MU_EZ - MU_GZ*mu_GE/mu_GG)
  
  C <- (B %*% t(MU_EZ - MU_GZ*mu_GE/mu_GG) - MU_W%*%t(MU_Z) - MU_GW%*%t(MU_GZ)/mu_GG + MU_WZ) %*% solve(O)
  
  Q <- MU_W%*%t(MU_W) + MU_GW%*%t(MU_GW)/mu_GG - MU_WW + B %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + 
    C %*% ( MU_Z%*%t(MU_W) + t(MU_GZ)%*%t(MU_GW)/mu_GG - MU_ZW + A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) )
  
  D <- (mu_GE * mu_GGE / mu_GG - mu_GEE) / (mu_EE - mu_GE^2 / mu_GG)
  E <- t(MU_GEZ - MU_Z*mu_GE - MU_GZ*mu_GGE/mu_GG + D*(MU_EZ - MU_GZ*mu_GE/mu_GG)) %*% solve(O)
  EFF <- ( t(MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + D*(MU_GW * mu_GE / mu_GG - MU_EW)) + 
             E %*% (A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + MU_Z%*%t(MU_W) + MU_GZ%*%t(MU_GW)/mu_GG - MU_ZW) ) %*% solve(Q)
  
  
  # Solve for \alpha_I
  alpha_I_num <- beta_E * (-mu_f*mu_GE - mu_Gf*mu_GGE/mu_GG + mu_GEf + D * (mu_Ef - mu_Gf*mu_GE/mu_GG)) +
    beta_E * E %*% (-mu_f*MU_Z - MU_GZ*mu_Gf/mu_GG + MU_fZ + A * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (-mu_Gf*mu_GE - mu_GGh*mu_GGE/mu_GG + mu_GGEh + D * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    beta_I * E %*% (-mu_Gh*MU_Z - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    t(MU_GEM - MU_M*mu_GE - MU_GM*mu_GGE/mu_GG + D * (MU_EM - MU_GM*mu_GE/mu_GG)) %*% BETA_M + 
    E %*% (A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M - 
    beta_E * EFF %*% (-mu_f*MU_W - MU_GW*mu_Gf/mu_GG + MU_fW + B %*% (mu_GEf - mu_GGh*mu_GE/mu_GG)) - 
    beta_E * EFF %*% C %*% (-mu_f*MU_Z - mu_Gf*MU_GZ/mu_GG + MU_fZ + A %*% (mu_Ef - mu_Gf*mu_GE/mu_GG)) - 
    beta_I * EFF %*% (-mu_Gh*MU_W - MU_GW*mu_GGh/mu_GG + MU_GhW + B %*% (mu_GEh - mu_GGh*mu_GE/mu_GG)) -
    beta_I * EFF %*% C %*% (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A %*% (mu_GEh - mu_GGh*mu_GE/mu_GG)) - 
    EFF %*% ( -MU_W%*%t(MU_M) - MU_GW%*%t(MU_GM)/mu_GG + MU_WM + B %*% t(MU_EM - MU_GM*mu_GE/mu_GG) ) %*% BETA_M - 
    EFF %*% C %*% ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M
  
  alpha_I_denom <- EFF %*% ( MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + B * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) +
    EFF %*% C %*% ( MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) - 
    ( mu_GE^2 + mu_GGE^2/mu_GG - mu_GGEE + D * (mu_GGE*mu_GE/mu_GG - mu_GEE) ) - 
    E %*% ( MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG + A * (mu_GGE*mu_GE/mu_GG - mu_GEE) )
  
  alpha_I <- alpha_I_num / alpha_I_denom
  
  R <- beta_E * (-MU_W*mu_f - MU_GW*mu_Gf/mu_GG + MU_fW + B * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_E * C %*% (-mu_f*MU_Z - MU_GZ*mu_Gf/mu_GG + MU_fZ + A %*% (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (-MU_W*mu_Gh - MU_GW*mu_GGh/mu_GG + MU_GhW + B * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    beta_I * C %*% (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    ( B %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_W%*%t(MU_M) - MU_GW%*%t(MU_GM)/mu_GG + MU_WM) %*% BETA_M + 
    C %*% ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M + 
    alpha_I * (MU_W*mu_GE + MU_GW*mu_GGE/mu_GG - MU_GEW + B * (mu_GGE*mu_GE/mu_GG - mu_GEE)) + 
    alpha_I * C * (MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A*(mu_GGE*mu_GE/mu_GG - mu_GEE))
  
  ALPHA_W <- - solve(Q) %*% R
  
  P <- beta_E * (-MU_Z*mu_f - MU_GZ*mu_Gf/mu_GG + MU_fZ + A * (mu_Ef - mu_Gf*mu_GE/mu_GG)) + 
    beta_I * (MU_GhZ - MU_Z*mu_Gh - MU_GZ*mu_GGh/mu_GG + A * (mu_GEh - mu_GGh*mu_GE/mu_GG)) + 
    alpha_I * (MU_Z*mu_GE + MU_GZ*mu_GGE/mu_GG - MU_GEZ + A * (mu_GGE*mu_GE/mu_GG - mu_GEE)) + 
    ( A %*% t(MU_EM - MU_GM*mu_GE/mu_GG) - MU_Z%*%t(MU_M) - MU_GZ%*%t(MU_GM)/mu_GG + MU_ZM) %*% BETA_M + 
    ( A %*% t(MU_GW*mu_GE/mu_GG - MU_EW) + MU_Z%*%t(MU_W) + MU_GZ%*%t(MU_GW)/mu_GG - MU_ZW) %*% ALPHA_W
  
  Bz_Az <- solve(O) %*% P
  ALPHA_Z <- BETA_Z - solve(O) %*% P
  
  alpha_E <- ( beta_E * (mu_Ef - mu_Gf*mu_GE/mu_GG) + beta_I * (mu_GEh - mu_GGh*mu_GE/mu_GG) +  
                 alpha_I * (mu_GGE*mu_GE/mu_GG - mu_GEE) + t(MU_EZ - MU_GZ*mu_GE/mu_GG) %*% Bz_Az + 
                 t(MU_EM - MU_GM*mu_GE/mu_GG) %*% BETA_M + t(MU_GW*mu_GE/mu_GG - MU_EW) %*% ALPHA_W ) /
    (mu_EE - mu_GE^2/mu_GG)
  
  Bg_Ag <- ( alpha_E*mu_GE - beta_E*mu_Gf + alpha_I*mu_GGE - beta_I*mu_GGh - t(MU_GZ) %*% Bz_Az + 
               t(MU_GW) %*% ALPHA_W - t(MU_GM) %*% BETA_M ) / mu_GG
  alpha_G <- beta_G - Bg_Ag
  
  
  alpha_0 <- beta_0 + beta_E*mu_f + beta_I*mu_Gh - alpha_I*mu_GE + t(MU_Z) %*% Bz_Az +
    t(MU_M) %*% BETA_M - t(MU_W) %*% ALPHA_W
  
  # Return 
  return(list(alpha_0=alpha_0, 
              alpha_G=alpha_G, 
              alpha_E=alpha_E,
              alpha_I=alpha_I,
              ALPHA_Z=ALPHA_Z,
              ALPHA_W=ALPHA_W,
              beta_vec = c(beta_0, beta_G, beta_E, beta_I, BETA_Z, BETA_M),
              mu_vec = c(mu_f, mu_h, MU_Z, MU_M, MU_W),
              cov_vec = c(mu_GG, mu_GE, mu_Gf, mu_Gh, MU_GZ, MU_GM, MU_GW, mu_EE,
              mu_Ef, MU_EZ, MU_EM, MU_EW, MU_fZ, MU_fW),
              MU_ZZ=MU_ZZ, MU_ZM=MU_ZM, MU_ZW=MU_ZW, MU_WZ=MU_WZ, MU_WM=MU_WM, MU_WW=MU_WW,
              HOM_vec = c(mu_GGE, mu_GGh, mu_GEE, mu_GEf, mu_GEh, MU_GEZ, MU_GEM, MU_GEW,
              MU_GhW, MU_GhZ, mu_GGEE, mu_GGEf, mu_GGEh)))
}