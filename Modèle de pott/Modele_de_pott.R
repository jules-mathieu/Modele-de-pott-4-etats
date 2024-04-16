###############
####### -- MODELE DE POTTS A 4 ETAT -- Objetcif : Retrouver la température de config.obs
##############

## init.configPotts
## -- Utilisation : init.configPotts initialise une configuration avec des valeurs aléatoires pour z.
## -- Arguments : nx et ny sont le nombre de colonnes et de lignes de la matrice
## -- Valeur : Objet de classe config

init.configPotts <- function(nx, ny){
  
  if (nx <= 0){stop("nx doit etre positif")}
  if (ny <= 0){stop("nx doit etre positif")}
  x <- c(1:nx)
  y <- c(1:ny)
  z <- matrix(data = sample(c(2,3,4,5),nx*ny,replace = TRUE), ncol = ny, nrow = nx, byrow = FALSE)
  res <- list(x,y,z)
  class(res) <- "configPotts"
  
  return(res)
}


## plot.configPotts
## -- Usage : plot.configPotts permet d’afficher une configuration de points.
## -- Arguments : un objet conf de la classe configPotts.
## -- Value : Rien.

plot.configPotts <- function(conf){
  if (class(conf)!= "configPotts"){
    stop("the argument conf is not a configPotts")
  }
  rouge <- rgb(0.55, 0, 0)  
  jaune <- rgb(1, 0.80, 0)   
  bleu_fonce <- rgb(0, 0.1, 0.58)  
  bleu_ciel <- rgb(0, 0.9, 0.99)
  
  image(conf[[3]], col = c(rouge,jaune,bleu_ciel,bleu_fonce))
  
}


## print.configPotts
## -- Usage : print.configPotts permet d’afficher dans la console une configuration de points. L’affichage dans
#             la console ne se fait que si la configuration a au moins 30 points en lignes et 40 en colonnes. 
#             Dans le cas inverse, elle affiche uniquement les dimensions de la configuration.
## -- Arguments : un objet conf de la classe configPotts.
## -- Value : Rien.

print.configPotts <- function(conf){
  if (class(conf)!= "configPotts"){
    stop("the argument conf is not a configPotts")
  }
  
  if (length(conf[[1]]) < 40 | length(conf[[2]]) < 30 ){
    print(dim(conf[[3]]))
  }
  else {
    print(conf[3]) 
  }
}


## get.pairs
## -- Usage : get.pairs permet de calculer le nombre de combinaisons de couleurs suivantes : (Bleu fonce,Bleu fonce), (Bleu fonce,Jaune)
#             (Bleu fonce,Bleu ciel), (Bleu fonce,Rouge), (Jaune,Jaune), (Jaune,Bleu ciel), (Jaune,Rouge), (Bleu ciel,Bleu ciel), (Bleu ciel, Rouge), (Rouge, Rouge)
## -- Arguments : un objet confPotts de la classe config.
## -- Value : La fonction retourne une liste comprenant le nombre d'occurence de chaque couple dans la config

get.pairs <- function(conf){
  if (class(conf)!= "configPotts"){stop("the argument conf is not a configPotts")}
  if (length(conf[[1]]) == 1 | length(conf[[2]]) == 1){stop("matrice trop petite")}
  
  x_temp <- length(conf[[1]])
  y_temp <- length(conf[[2]])
  
  neighbour_matrixH <- conf[[3]][,c(1:y_temp - 1)] * conf[[3]][,c(2:y_temp)]
  neighbour_matrixV <- conf[[3]][c(1:x_temp - 1),] * conf[[3]][c(2:x_temp),]  
  
  return (list( sum(neighbour_matrixH == 4) + sum(neighbour_matrixV == 4), 
                sum(neighbour_matrixH == 6) + sum(neighbour_matrixV == 6),
                sum(neighbour_matrixH == 8) + sum(neighbour_matrixV == 8),
                sum(neighbour_matrixH == 9) + sum(neighbour_matrixV == 9),
                sum(neighbour_matrixH == 10) + sum(neighbour_matrixV == 10),
                sum(neighbour_matrixH == 12) + sum(neighbour_matrixV == 12),
                sum(neighbour_matrixH == 15) + sum(neighbour_matrixV == 15),
                sum(neighbour_matrixH == 16) + sum(neighbour_matrixV == 16),
                sum(neighbour_matrixH == 20) + sum(neighbour_matrixV == 20),
                sum(neighbour_matrixH == 25) + sum(neighbour_matrixV == 25)))
}


## get.Hamiltonien
## -- Usage : get.Hamiltonien permet de calculer l’Hamiltonien d’une configuration.
## -- Arguments : un objet confPotts de la classe configuration, un ŕeel J et un reel h
## -- Value : La fonction retourne la valeur de l’Hamiltonien.

get.Hamiltonien <- function(conf, t){
  res1 <- 1/T*(get.pairs(conf)[[2]] + get.pairs(conf)[[3]] + get.pairs(conf)[[5]] + get.pairs(conf)[[6]] 
               + get.pairs(conf)[[7]] + get.pairs(conf)[[9]] - (get.pairs(conf)[[1]] + get.pairs(conf)[[4]]
                                                                + get.pairs(conf)[[8]] + get.pairs(conf)[[10]]))
  return(res1)
}


## iteration.MH
## -- Usage : iteration.MH permet d’effectuer une iteration de l’algorithme de Metropolis-Hastings.
## -- Arguments : un objet confPotts de la classe config, un reell β, un reel J et un reel h.
## -- Value : La fonction renvoir la nouvelle configuration modifiee ou non dans un objet de la classe configPotts.

iteration.MH <- function(conf, beta, t){
  
  conf2 <- conf
  P <- 0
  
  indice.random <- sample((length(conf2[[1]])*(length(conf2[[2]]))), size = 1)
  if (conf2[[3]][indice.random] == 2){ conf2[[3]][indice.random] <- sample(c(3,4,5),size = 1)}
  if (conf2[[3]][indice.random] == 3){ conf2[[3]][indice.random] <- sample(c(2,4,5),size = 1)}
  if (conf2[[3]][indice.random] == 4){ conf2[[3]][indice.random] <- sample(c(2,3,5),size = 1)}
  if (conf2[[3]][indice.random] == 5){ conf2[[3]][indice.random] <- sample(c(2,3,4),size = 1)}
  
  H1 <- get.Hamiltonien(conf,t)
  H2 <- get.Hamiltonien(conf2,t)
  
  if (H2 - H1 >= 0){P <- 1}
  else {P <- exp(-beta*(H2-H1))}
  
  if (P >= runif(n=1, min=0, max =1)){
    return(conf2)
  }
  return(conf)
}


## simul.Potts
## -- Usage : simul.Potts effecrtue n.sim iterations de Metropolis-Hastings pour simuler une configuration du modele de Potts.
## -- Arguments : — nx et ny:la taille de la grille
#                 — T et beta pour definir le modele a simuler
#                 — n.sim : un entier correspond au nombre d’iterations a effectuer
## -- Value : La fonction renvoit une liste avec la configuration finale (sous forme d’un objet de la classe configPotts)
#             et un vecteur stockant l'évolution de l'énergie

simul.Potts <- function(nx,ny,beta,t,n.sim){
  
  new_conf <- init.configPotts(nx,ny)
  Hamiltonien <- rep(NA, times = n.sim)
  
  for (k in 1:(n.sim)){
    new_conf <- iteration.MH(new_conf, beta, t)
    Hamiltonien[k] <- get.Hamiltonien(new_conf,t)
  }
  
  plot.configPotts(new_conf)
  return(list(new_conf, Hamiltonien))
}


## recherche.temp
## -- Usage : recherche.temp 
##
##


recherche.temp <- function(image_origine, beta, n.sim, n.precis, tmin, tmax, pas ){
  
  nx <- length(image_origine[[1]])
  ny <- length(image_origine[[2]])
  
  list_compare <- list()
  list_temp <- list()
  
  couple_ref <- get.pairs(image_origine)
  
  for (temp in seq(from = tmin, to = tmax, by = pas)) {
    matrice_compare <- matrix(ncol =10 , nrow = n.precis)
    
    for (k in 1:n.precis){
      
      conf_temporaire <- simul.Potts(nx,ny,beta,t,n.sim)[[1]]
      nouvelle_ligne <- unlist(get.pairs(conf_temporaire))
      matrice_compare[k,]<- nouvelle_ligne
    }
    couple_mean <- apply(matrice_compare, MARGIN = 2, mean)
    list_compare <- append(list_compare, list(couples = couple_mean))
  }
  vec_ref <- unlist(get.pairs(config.obs))
  
  for (i in seq(from = 1, to = (tmax - tmin)/pas +1)) {
    norm_temp <- sqrt(sum((unlist(list_compare[i]) - vec_ref)^2))
    list_temp <- append(list_temp, norm_temp)
    vec_temp <- unlist(list_temp)
    
  }
  indice <- which(vec_temp == min(vec_temp))
  return((indice-1)*pas + tmin)
}

recherche.temp(config.obs, beta =1, n.sim = 100, n.precis = 10, tmin = 0.1, tmax = 0.5, pas = 0.05)

