#' @importFrom packagename functionname #sert à déclarer les fonctions provenant d'autres packages utilisées dans le code

# Compute confidence ellipses for ASCA-E score plots


confidence_ellipse_LimmPCA <- function(
    Thetahat, #: matrix of estimated coefficients
    Pf, #: matrix of PCA loadings (PCA on the effect matrix, not backtransformed)
    X, #: model matrix
    varParameters, #: list of var-cov matrix for each PC
    effect_cols,
    effect_name = "effect",
    effect_levels = NULL,
    alpha = 0.05,
    rf = 2,
    npoints = 100
  ) {
  
  # ============================================================
  # [0] VERIFICATION DES ENTREES
  # ============================================================
  
  Thetahat <- as.matrix(Thetahat)
  Pf <- as.matrix(Pf)
  X <- as.matrix(X)
  
  m <- ncol(Thetahat) #: number of transformed responses
  p <- ncol(X) #: number of parameters
  
  if (is.null(colnames(X))) {
    stop("X doit posseder des noms de colonnes.")
  }
  
  if (!all(effect_cols %in% colnames(X))) {
    stop("Certaines colonnes de effect_cols sont absentes de X.")
  }
  
  if (length(varParameters) != m) {
    stop("Le nombre de matrices varParameters doit correspondre au nombre de reponses transformées.")
  }
  
  if (rf != 2) {
    stop("Cette fonction génère des ellipses en deux dimensions : rf doit être égal a 2.")
  }
  
  # ============================================================
  # [1] CONSTRUCTION DE LA MATRICE DE CONTRASTES L
  # ============================================================
  
  # Extraire les combinaisons uniques de niveaux 
  X_effect <- X[, effect_cols]
  
  X_effect_unique <- unique(as.data.frame(X_effect))
  
  # Number of levels
  K <- nrow(X_effect_unique) 
  
  # Initialiser L (nb levels x nb parameters) : matrice de 0
  L <- matrix(
    0,
    nrow = K,
    ncol = p,
    dimnames = list(NULL, colnames(X))
  )
  
  # Index des colonnes dans X qui codent les levels 
  col_idx <- match(effect_cols, colnames(X))
  
  # Remplace les colonnes des paramètres concernés par cet effet par leur codage 
  L[, col_idx] <- as.matrix(X_effect_unique)
  
  # Final contrast matrix of dimension p x K
  L <- t(L) 
  
  # ============================================================
  # [2] PREPARATION DES LOADINGS PCA
  # ============================================================
  
  if (nrow(Pf) != m) {
    stop("Dimensions incompatibles entre Pf et Thetahat.")
    }
  
  Pf <- Pf[, seq_len(rf), drop = FALSE] #conserver uniquement les rf premières colonnes
  
  # ============================================================
  # [3] PROJECTION DES CENTRES
  # ============================================================
  
  centers_proj <- t(L) %*% Thetahat %*% Pf
  
  # ============================================================
  # [4] COMPUTATION OF SIGMA_Z MATRICES FOR EACH LEVEL
  # ============================================================
  
  compute_Sigma_z <- function(varParameters, L, Pf) {
    
    # Nombre de modèles et de contrastes
    m <- length(varParameters)
    K <- ncol(L)
    
    # Initialisation
    Sigma_z <- vector("list", K)
    
    # Calcul pour chaque contraste
    for (k in seq_len(K)) {
      
      l <- L[, k] # vecteur de contraste du contraste k 
      
      # Variance pour chaque modèle
      v <- numeric(m) 
      
      for (j in seq_len(m)) {
        
        VarB <- as.matrix(varParameters[[j]]) # variance de hat_beta_j
        
        v[j] <- drop(t(l) %*% VarB %*% l)
      }
      
      s11 <- sum(Pf[, 1]^2 * v)
      s12 <- sum(Pf[, 1] * Pf[, 2] * v)
      s22 <- sum(Pf[, 2]^2 * v)
      
      Sigma_z[[k]] <- matrix(
        c(s11, s12, s12, s22),
        nrow = 2
      )
    }
    
    names(Sigma_z) <- paste0("Contrast_", seq_len(K))
    
    return(Sigma_z)
  }
  
  sigma_z <- compute_Sigma_z(varParameters, L, Pf)
  
  # ============================================================
  # [5] SEUIL CHI-CARRE
  # ============================================================
  
  critical_value <- qchisq(1 - alpha, df = rf)
  
  # ============================================================
  # [6] GENERATION DES ELLIPSES
  # ============================================================
  
  generate_ellipse_points <- function(
    center,
    Sigma_proj,
    scale,
    npoints = 100
  ) {
    
    # center : vecteur contenant les coordonnées du centre de l'ellipse
    # Sigma_proj : matrice de covariance \(2\times2\) associée à l'estimation projetée
    # scale : facteur qui détermine la taille de l'ellipse, quantile chi carré 
    # npoints : nombre de points utilisés pour tracer le contour de l'ellipse
    
    angles <- seq(0, 2 * pi, length.out = npoints)
    
    unit_circle <- cbind(
      cos(angles),
      sin(angles)
    )
    
    eig <- eigen(Sigma_proj)
    
    if (min(eig$values) < -1e-8) {
      stop("Une matrice Sigma_z n'est pas semi-definie positive.")
    }
    
    Lambda = diag(eig$values)  # mettre les valeurs propresdans une matrice diagonale
    
    ellipse <- sqrt(scale) *
      unit_circle %*% sqrt(Lambda) %*% t(eig$vectors) +
      matrix(center, nrow = npoints, ncol = 2, byrow = TRUE)
    
    data.frame(
      PC1 = ellipse[, 1], #vérifier ordre PC ! 
      PC2 = ellipse[, 2]
    )
  }
  
  ellipses_list <- lapply(1:nrow(centers_proj), function(i) {
    
    df <- generate_ellipse_points(
      center = centers_proj[i, ],
      Sigma_proj = sigma_z[[i]],
      scale = critical_value,
      npoints = npoints
    )
    
    df$group <- paste0(effect_name, "_", i)
    
    df
  })
  
  ellipses_df <- do.call(rbind, ellipses_list) #verifier le résultat ! 
  
  
  
  # ============================================================
  # [7] ATTRIBUTION DES NIVEAUX DE L'EFFET
  # ============================================================
  
  if (is.null(effect_levels)) {
    effect_levels <- paste0("Level_", seq_len(K))
  }
  
  if (length(effect_levels) != K) {
    stop("Le nombre de niveaux doit correspondre au nombre de contrastes.")
  }
  
  ellipses_df[[effect_name]] <- factor( #verifier le resultat ! 
    rep(effect_levels, each = npoints),
    levels = effect_levels
  )
  
  # ============================================================
  # [8] DATAFRAME DES CENTRES
  # ============================================================
  
  centers_df <- data.frame(
    PC1 = centers_proj[, 1], # a verifier que PC1 ! 
    PC2 = centers_proj[, 2],
    group = paste0(effect_name, "_", seq_len(K))
  )
  
  centers_df[[effect_name]] <- factor(
    effect_levels,
    levels = effect_levels
  )
  
  # ============================================================
  # [9] RESULTATS
  # ============================================================
  
  return(list(
    L = L,
    P = Pf,
    centers_proj = centers_proj,
    sigma_z = sigma_z,
    critical_value = critical_value,
    ellipses_df = ellipses_df,
    centers_df = centers_df,
    effect_name = effect_name,
    effect_levels = effect_levels
  ))
}
