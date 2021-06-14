### KLAR KOD, inläsning av fMRI ### 

## Sätt upp working directory ## 
setwd("search path")

## Installerar paket ## 
#install.packages("oro.nifti")
#install.packages("R.matlab")
#install.packages("matlab")
library("oro.nifti")
library("R.matlab")
library("matlab")

## LÄSER IN ALLA FILER MED .nii.gz i dir()
temp <- list.files(pattern="*.nii.gz") 
### laddar in alla filer i enviorment
a <- list2env(lapply(setNames(temp, make.names(gsub("*.nii.gz$", "", temp))),readNIfTI))

### Roterar data, annars läses in bakvänt detta ger samma resultat som i matlab
MDGA <- function(Rdata){
  Rdata <- Rdata@.Data
  for(i in 1:dim(Rdata)[3]){
    for(j in 1:dim(Rdata)[4]){
      Rdata[,,i,j]<- flipud(Rdata[,,i,j])
    }
  }
  return(Rdata)
  
}

## erhåller alla subjekt + GICA som ligger i dir()
fMRI_data <- lapply(a,MDGA)

## MASKNING AV DATA ##
# Viktigt att GICA och Subjectdata är samma voxeldim. 

# Välj ut subjekt av intresse, tex. Subject 5, 10, 15
sub <- fMRI_data$X6_filtered_func_data_Subject5
# Välj ut hur många komponenter GICA ska ha 1:20 tex. som används i uppsatsen
GICA20 <- fMRI_data$X6_melodic_IC[,,,1:20]

# Börjar maskningen där alla värden under medelvärdet för hela individens aktivitet tas bort.
# Beräknar medelvärdet över alla tidsförlopp där voxlarna inte är 0, OBS för utvalda subjektet
mask <- sub < mean(sub[sub != 0])
# skapar mask, baserat på första tidsförloppet för att undvika dimensionsproblem 
mask <- mask[,,,1]

## Behåller sub för att använda som mall i visualising 
sub_mall <- sub

# maskar Subjektet samt GICA 
sub[mask] <- 0
GICA20[mask] <- 0

# Dimensionerna för individ-data och GICA
dim_T <- dim(sub)
dim_G <- dim(GICA20)

## omstrukturerar array i dim. (X*Y*Z*T) till en matris NXT, N=voxlar, T=tid

Y_N_T <- matlab::reshape(sub, c(dim_T[1]*dim_T[2]*dim_T[3],dim_T[4]))

## omstrukturerar array i dim. (X*Y*Z*M) till en matris NXM, N=voxlar, M=Komponenter 
X_N_M <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))

#Sätter TRUE på alla voxlar som inte är 0
row.sub <- apply(Y_N_T,1, function(row)all(row!=0))

# Väljer ut dessa (OBS, Y används som mask) 
Y_N_T <- Y_N_T[row.sub,]
X_N_M <- X_N_M[row.sub,]

## Standardisering ##
### Enligt Mejia 2019 ###

Standardisering <- function(ORGDATA){
  ## Centera data för varje voxel
  C1  <-   ORGDATA - apply(ORGDATA,1, mean)
  
  ## Centera data för varje tidspunkt
  S1 <- matrix(apply(C1, 2,mean),ncol=dim(C1)[2],nrow=dim(C1)[1],byrow=T)
  C2 <- C1 - S1
  
  ## Skapar matris med standardavvikelse från orginaldata innan centeringarna har gjorts
  MAT <- matrix(ncol =dim(C1)[2],nrow=dim(C1)[2],0) 
  diag(MAT) <- apply(ORGDATA,2, sd)
  # Standardiserar den centerade matrisen
  C3  <- C2 %*% solve(MAT)
  # Retunerar den färdiga centeringen
  return(C3)  
}
## Använder standardiseringsfunktionen för både GICA och Y 
Y <- Standardisering(Y_N_T)
X <- Standardisering(X_N_M)

### Frekventistisk anpassning 
# Paket för Psuedoinvers
# install.packages("pracma")
library("pracma")

BTC <-  pinv(X) %*% Y
BSM <- Y %*% pinv(BTC)

###################### Bayesianska Dual Regression - Modell 3 ######################
#install.packages("LaplacesDemon")
library("LaplacesDemon")
### BTC, hyperparameterar samt skapar matriser  ### 
B0_TC <- matrix(0, nrow = 20, ncol = 119)
V0_TC <- 1 * diag(119)
LAMBDA0_TC <- 1 * diag(20)
BN_TC <- solve(t(X)%*%X+ LAMBDA0_TC) %*% (t(X) %*%Y + LAMBDA0_TC%*%B0_TC)
LAMBDAN_TC <- t(X)%*%X + LAMBDA0_TC
Vn_TC <- V0_TC + t(Y - X%*%BN_TC) %*% (Y- X%*%BN_TC) + t(BN_TC-B0_TC)%*% LAMBDA0_TC %*%(BN_TC-B0_TC)

n <- dim(Y)[1]
nu0_TC <- 0
nu_n_TC  <- nu0_TC + n

### Definerar hyperparametrar för BSM ###
# Transponera Y
YSM <- t(Y)
B0_SM <- matrix(0, nrow = 20, ncol = n)
V0_SM <- 1 * diag(n)
LAMBDA0_SM <- 1 * diag(20)

nu_n_SM  <-  n + 2

setwd("search path")


## OBS TAR LÅNG TID ATT KÖRA. 
# Bör köras stegvis, varje dragning sparas i .RData fil.
for( i in 1:100){
  
  # Posterior för BTC
  set.seed(i)
  wish_BTC <- rinvwishart(nu = nu_n_TC  , S = (Vn_TC))
  # kovariansmatrisen för OLS
  set.seed(i)
  Bayes_BTC <- rmatrixnorm(M = BN_TC, U = as.symmetric.matrix(solve(LAMBDAN_TC)), V = (wish_BTC))
  
   # Hyperparametrar för BSM + transponera BTC 
  XSM <-  t(Bayes_BTC)
  BN_SM <- pinv(XSM) %*% YSM
  LAMBDAN_SM <- t(XSM)%*%XSM + LAMBDA0_SM
  Vn_SM <- V0_SM + t(YSM - XSM%*%BN_SM)%*% (YSM - XSM%*%BN_SM) + t(BN_SM-B0_SM) %*% LAMBDA0_SM %*% (BN_SM-B0_SM)
#  # Posterior för BSM
  set.seed(i)
  wish_BSM <- rinvwishart(nu = nu_n_SM  , S = (Vn_SM))
  # kovariansmatrisen för OLS sätts som fixa.
  
  set.seed(i)
  Bayes_BSM <- rmatrixnorm(M = BN_SM, U = as.symmetric.matrix(solve(LAMBDAN_SM)), V = (wish_BSM))

  save(Bayes_BTC,Bayes_BSM,wish_BTC, file = paste("Sub5_Sample",i,".RData",sep=""), compress = F)
  save(wish_BSM, file = paste("Sub5_Wish",i,".RData",sep=""), compress = F)
  print(i)
}



###################### Bayesianska Dual Regression - Modell 2 ######################

## BERÄKNING AV KOVARIANSMATRISERNA MED FREKVENTISTISK STRUKTUR ## 
## BTC 

Sb <- t(Y - X %*% BTC) %*% (Y - X %*% BTC)
s_kvad<-  Sb/nrow(Y)
Ko_BTC <- diag(ncol=119,nrow=119)

diag(Ko_BTC) <- diag(s_kvad)


## BSM 

Sb_SM <- t(t(Y) -  t(BTC) %*% t(BSM)) %*% (t(Y) - t(BTC) %*% t(BSM))
s_kvad_SM <-  Sb_SM/nrow(t(Y))

Ko_BSM <- diag(ncol=8349,nrow=8349)

diag(Ko_BSM) <- diag(s_kvad_SM)


### BTC, hyperparameterar samt skapar matriser  ### 
B0_TC <- matrix(0, nrow = 20, ncol = 119)
LAMBDA0_TC <- 1 * diag(20)
V0_TC <- 1 * diag(119)
BN_TC <- solve(t(X)%*%X+ LAMBDA0_TC) %*% (t(X) %*%Y + LAMBDA0_TC%*%B0_TC)
LAMBDAN_TC <- t(X)%*%X + LAMBDA0_TC
Vn_TC <- V0_TC + t(Y - X%*%BN_TC) %*% (Y- X%*%BN_TC) + t(BN_TC-B0_TC)%*% LAMBDA0_TC %*%(BN_TC-B0_TC)

n <- dim(Y)[1]

### Definerar hyperparametrar för BSM ###
# Transponera Y
YSM <- t(Y)
B0_SM <- matrix(0, nrow = 20, ncol = n)
LAMBDA0_SM <- 1 * diag(20)
V0_SM <- 1 * diag(8349)
setwd("search path")


for( i in 1:100){
  

  # kovariansmatrisen för OLS
  set.seed(i)
  Bayes_BTC <- rmatrixnorm(M = BN_TC, U = as.symmetric.matrix(solve(LAMBDAN_TC)), V = (Ko_BTC))
  
  # Hyperparametrar för BSM + transponera BTC 
  XSM <-  t(Bayes_BTC)
  BN_SM <- pinv(XSM) %*% YSM
  LAMBDAN_SM <- t(XSM)%*%XSM + LAMBDA0_SM
  Vn_SM <- V0_SM + t(YSM - XSM%*%BN_SM)%*% (YSM - XSM%*%BN_SM) + t(BN_SM-B0_SM) %*% LAMBDA0_SM %*% (BN_SM-B0_SM)
  
  # kovariansmatrisen för OLS sätts som fixa.
  set.seed(i)
  Bayes_BSM <- rmatrixnorm(M = BN_SM, U = as.symmetric.matrix(solve(LAMBDAN_SM)), V = (Ko_BSM))
  save(Bayes_BTC,Bayes_BSM, file = paste("FixKorSub5_Sample",i,".RData",sep=""), compress = F)
  print(i)
}




