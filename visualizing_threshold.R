s <- 1

sample <- c("Sub5_Sample","Sample","Sub15_Sample")
fix_sample <- c("FixKorSub5_Sample",
                "FixKorSub10_Sample",
                "FixKorSub15_Sample")

komp <- c("Sub5_Komp_","Sub10_Komp_","Sub15_Komp_")

sampleFixMap <- c("Sub_5_Kovarians","Sub_10_Kovarians","Sub_15_Kovarians")
sampleMap <- c("Subjekt5","Subjekt10","Subjekt15")
voxlar <- c(8349,7997,8075)


BSM_array_kor <- array(0, dim=c( 20,voxlar[s],100))
BTC_array_kor <- array(0, dim=c( 20,119 ,100))

BSM_array_wish <- array(0, dim=c( 20,voxlar[s] ,100))
BTC_array_wish <- array(0, dim=c( 20,119 ,100))


if( s == 1 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject5
}
if( s == 2 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject10
}
if( s == 3 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject15
}



GICA20 <- fMRI_data$X6_melodic_IC[,,,1:20]
mask <- sub < mean(sub[sub != 0])
mask <- mask[,,,1]
sub_mall <- sub
sub[mask] <- 0
GICA20[mask] <- 0

dim_T <- dim(sub)
dim_G <- dim(GICA20)


### 
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

###


c <- 1 
for(i in c(51:150)){
  setwd(paste("/Users/jacobwnilsson/Desktop/KOD3/SAMPLE/",sampleFixMap[s],sep=""))
  load(paste(fix_sample[s],i, ".RData", sep="") )
  print(i)
  BSM_array_kor[,,c] <- Bayes_BSM
  rm(Bayes_BSM)
  setwd(paste("/Users/jacobwnilsson/Desktop/KOD3/SAMPLE/",sampleMap[s],sep=""))
  load(paste(sample[s],i, ".RData", sep="") )
  print(i)
  BSM_array_wish[,,c] <- Bayes_BSM
  rm(Bayes_BSM)
  
  c <- c+1
}






## FREK ## 
X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))


X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- BSM
Frek <- reshape(X_N_M1, c(30,36,30,20))

## Bayesiansk Korrelation ##

setwd("/Users/jacobwnilsson/Desktop/KOD3/BILDER")

for(t in c(50,67)){

Greater_0_korr <- t(apply(BSM_array_kor,c(1,2), function(x) sum(x > 0) ))

Aktiv_Mall_korr<- Greater_0_korr[,] > t
Posterior_Mean_BSM_korr <- t(apply(BSM_array_kor,c(1,2), mean))
Posterior_Mean_BSM_korr[!Aktiv_Mall_korr] <- 0
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM_korr
Mean_greater_0_BSM_korr <- reshape(X_N_M1, c(30,36,30,20))
## Bayesiansk Wishart ## 
Greater_0 <- t(apply(BSM_array_wish,c(1,2), function(x) sum(x > 0) ))

Aktiv_Mall<- Greater_0[,] > t
Posterior_Mean_BSM <- t(apply(BSM_array_wish,c(1,2), mean))
Posterior_Mean_BSM[!Aktiv_Mall] <- 0
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM
Mean_greater_0_BSM <- reshape(X_N_M1, c(30,36,30,20))






MALL <- matlab::reshape(sub_mall, c(30*36*30,119))
M <- reshape(MALL, c(30,36,30,119))
M <- M[,,,1]





### SKRIV UT FÖR INDIVID 1  ### 
pdf(file=paste(c("Sub5_BSM_post",t,".pdf"),collapse = ""))
ortho2(M,Mean_greater_0_BSM[,,,6],mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,25,length.out = 65),col.y = rev(jet.colors(64)))
dev.off()


### SKRIV UT FÖR INDIVID 2  ### 
#pdf(file=paste(c("Sub10_BSM_post",t,".pdf"),collapse = ""))
#ortho2(M,Mean_greater_0_BSM[,,,1],mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,25,length.out = 65),col.y = rev(jet.colors(64)))
#dev.off()


### SKRIV UT FÖR INDIVID 2  ### 
#pdf(file=paste(c("Sub15_BSM_post",t,".pdf"),collapse = ""))
#ortho2(M,Mean_greater_0_BSM[,,,4],mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,25,length.out = 65),col.y = rev(jet.colors(64)))
#dev.off()

}





brainz1 <- image_read("Sub15_4_BSM_post50.pdf")
brainz2 <- image_read("Sub15_4_BSM_post67.pdf")

brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))


img <- c(brainz1,brainz2)

imz <- image_append(img, stack = TRUE)
image_info(imz)
colorbar <- function(lut, min, max=-min, nticks=11, title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4, seq(min, max, len=nticks),las=1,col="white",col.ticks="white",col.axis="white",cex.axis=1.3)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}

pdf(file="colorbarjet.pdf")
par(mar=c(1,26,0,6),bg ="black")
colorbar(lut=rev(jet.colors(64)), min=0,max=25, nticks = 5)
dev.off()


colorbar <- image_read("colorbarjet.pdf")
colorbar <- image_chop(colorbar,geometry = geometry_area(300,5,5,0))
colorbar <- image_scale(colorbar,geometry = geometry_area(0,380,0,0))

brain_plot <-image_append(c(imz,colorbar), stack = F)


brain_plot <- image_annotate(brain_plot, "Tröskelnivå: 50%", size = 20, color = "white",
                             degrees = 0, location = "+10+1")

brain_plot <- image_annotate(brain_plot, "Tröskelnivå: 67%", size = 20, color = "white",
                             degrees = 0, location = "+10+185")



pdf(file=paste("sub3_4_5067.pdf"), width = 5.5, height = 3.5)
par(mar=c(0,3,0,1),bg ="black")
plot(brain_plot)
dev.off()


