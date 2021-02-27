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




MALL <- matlab::reshape(sub_mall, c(30*36*30,119))
M <- reshape(MALL, c(30,36,30,119))
M <- M[,,,1]


c <- 1 
for(i in c(51:150)){
  setwd(paste("search path",sampleFixMap[s],sep=""))
  load(paste(fix_sample[s],i, ".RData", sep="") )
  print(i)
  BSM_array_kor[,,c] <- Bayes_BSM
  rm(Bayes_BSM)
  setwd(paste("search path",sampleMap[s],sep=""))
  load(paste(sample[s],i, ".RData", sep="") )
  print(i)
  BSM_array_wish[,,c] <- Bayes_BSM
  rm(Bayes_BSM)
  
  c <- c+1
}



BSMWISH <- t(apply(BSM_array_wish, c(1,2), mean))
BSMKOR <- t(apply(BSM_array_kor, c(1,2), mean))


colorbar <- function(lut, min, max=-min, nticks=11, title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4, seq(min, max, len=nticks),las=1,col="white",col.ticks="white",col.axis="white",cex.axis=1.3)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}
tmatrix <- matrix(0, nrow=20, ncol=8075)
x= BSM_array_wish
y= BSM_array_kor

for(j in 1:20){
  for(i in 1:8075){
    test <- t.test( x[j,i,] , y[j,i,])
    tmatrix[j,i] <- test$statistic[[1]]  
  }
}


setwd("search path")


pdf("tkvotcolorr.pdf")
par(mar=c(1,26,0,6),bg ="black")
colorbar(lut=rev(jet.colors(64)), min=-5,max=5, nticks = 5)
dev.off()

pdf("diffcolorr.pdf")
par(mar=c(1,26,0,6),bg ="black")
colorbar(lut=rev(jet.colors(64)), min=-20,max=20, nticks = 5)
dev.off()



setwd("search path")

for(k in c(1,4,6)){
  
X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- t(tmatrix)


diff <- reshape(X_N_M1, c(30,36,30,20))



pdf(file=paste("Sub10_Tkvot","_",k,".pdf",sep = ""))


ortho2(M,diff[,,,k],mfrow=c(1,3),ybreaks = seq(-5,5,length.out = 65),NA.y=T,NA.x=T,crosshairs = F,
       col.y = rev(jet.colors(64)))
dev.off()



X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- BSMWISH - BSMKOR


diff <- reshape(X_N_M1, c(30,36,30,20))

pdf(file=paste("Sub1_Diff","_",k,".pdf",sep = ""))

ortho2(M,diff[,,,k],mfrow=c(1,3),ybreaks = seq(-20,20,length.out = 65),NA.y=T,NA.x=T,crosshairs = F,
       col.y = rev(jet.colors(64)))
dev.off()

}

## DIFF ## 



brainz1 <- image_read("Sub15_Tkvot_1.pdf")
brainz2 <- image_read("Sub15_Tkvot_4.pdf")
brainz3 <- image_read("Sub15_Tkvot_6.pdf")

brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
brainz3 <- image_crop(brainz3, geometry = geometry_area(550,190,0,155))

img <- c(brainz1,brainz2,brainz3)

imz <- image_append(img, stack = TRUE)

colorbar <- image_read("tkvotcolorr.pdf")
colorbar <- image_chop(colorbar,geometry = geometry_area(300,5,5,0))
colorbar <- image_scale(colorbar,geometry = geometry_area(0,570,0,0))

brain_plot <-image_append(c(imz,colorbar), stack = F)


brain_plot <- image_annotate(brain_plot, "Nätverkskarta 1", size = 20, color = "white",
                             degrees = 0, location = "+10+2")

brain_plot <- image_annotate(brain_plot, "Nätverkskarta 4", size = 20, color = "white",
                             degrees = 0, location = "+10+185")

brain_plot <- image_annotate(brain_plot, "Nätverkskarta 6", size = 20, color = "white",
                             degrees = 0, location = "+10+375")
pdf("Sub15_Tkvot_146.pdf", width = 5.5, height = 3.5)
par(mar=c(0,3,0,1),bg ="black")
plot(brain_plot)
dev.off()



####################################

setwd("search path")



pdf("diffcolorr.pdf")
par(mar=c(1,26,0,6),bg ="black")
colorbar(lut=rev(jet.colors(64)), min=-10,max=10, nticks = 5)
dev.off()



for(k in c(1,4,6)){

X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- BSMWISH - BSMKOR


diff <- reshape(X_N_M1, c(30,36,30,20))

pdf(file=paste("Sub15_Diff","_",k,".pdf",sep = ""))

ortho2(M,diff[,,,k],mfrow=c(1,3),ybreaks = seq(-10,10,length.out = 65),NA.y=T,NA.x=T,crosshairs = F,
       col.y = rev(jet.colors(64)))
dev.off()

}



brainz1 <- image_read("Sub15_Diff_1.pdf")
brainz2 <- image_read("Sub15_Diff_4.pdf")
brainz3 <- image_read("Sub15_Diff_6.pdf")

brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
brainz3 <- image_crop(brainz3, geometry = geometry_area(550,190,0,155))

img <- c(brainz1,brainz2,brainz3)

imz <- image_append(img, stack = TRUE)

colorbar <- image_read("diffcolorr.pdf")
colorbar <- image_chop(colorbar,geometry = geometry_area(300,5,5,0))
colorbar <- image_scale(colorbar,geometry = geometry_area(0,570,0,0))

brain_plot <-image_append(c(imz,colorbar), stack = F)


brain_plot <- image_annotate(brain_plot, "Nätverkskarta 1", size = 20, color = "white",
                             degrees = 0, location = "+10+2")

brain_plot <- image_annotate(brain_plot, "Nätverkskarta 4", size = 20, color = "white",
                             degrees = 0, location = "+10+185")

brain_plot <- image_annotate(brain_plot, "Nätverkskarta 6", size = 20, color = "white",
                             degrees = 0, location = "+10+375")

pdf("Sub15_Diff_146.pdf", width = 5.5, height = 3.5)
par(mar=c(0,3,0,1),bg ="black")
plot(brain_plot)
dev.off()



image_negate(brain_plot)







####################################
setwd("search path")

library("neurobase")
pdf(file="slice.pdf")

ortho2(M,mfrow=c(1,3),ybreaks = seq(-20,20,length.out = 65),NA.y=T,NA.x=T,crosshairs = F,
       col.y = rev(jet.colors(64)))
dev.off()


library("magick")
slice <- image_read("slice.pdf")
slice <- image_crop(slice, geometry = geometry_area(550,220,0,120))

slice <- image_annotate(slice, "x - Koronal", size = 20, color = "white",
                             degrees = 0, location = "+20+2")

slice <- image_annotate(slice, "y - Sagittal", size = 20, color = "white",
                        degrees = 0, location = "+190+2")

slice <- image_annotate(slice, "z - Transversell", size = 20, color = "white",
                        degrees = 0, location = "+340+2")



pdf("tvärsnittxyz.pdf", width = 5.5, height = 3.5)
par(mar=c(0,1,0,1),bg ="black")
plot(slice)
dev.off()



setwd("search path")

### TA BORT KONVERGENS för voxlar ##
konv <- image_read("Sub15_konv6.pdf")
konv <- image_chop(konv,geometry = geometry_area(0,150,0,975))

image_write(konv, path = "Sub15_konv6_ny.pdf", format = "pdf",quality=100)

### TA BORT KONVERGENS komp 1-6##

konv <- image_read("Sub15_Komp_Konvergens_bayes.pdf")
konv <- image_chop(konv,geometry = geometry_area(0,20,0,0))
image_write(konv, path = "Sub15_Komp_Konvergens_bayes_ny.pdf", format = "pdf",quality=100)

