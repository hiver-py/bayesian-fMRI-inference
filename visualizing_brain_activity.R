## Visualiseringskod utgår från att alla dragningar ligger seperata i dir() 
# Här läses varje dragning in stegvis för att kunna göra all visualisering, 

library("oro.nifti")
library("neurobase")
library("magick")
library("png")
library("ggplot2")
library("reshape")
# Subjekt 5: 8349
# Subjekt 10:7997
# Subjekt 15: 8075
 #väljer ut individen som ska undersökas, s = 1 = individ 5,s = 2 = individ 10, s = 3 = individ 15
s <- 1

for(s in s){
  
# matchar filnamn för att hantera inläsning, Sub*_sample är från modell 3
sample <- c("Sub5_Sample","Sample","Sub15_Sample")
# matchar filnamn för att hantera inläsning, FixKor*_sample är från modell 2
fix_sample <- c("FixKorSub5_Sample",
"FixKorSub10_Sample",
"FixKorSub15_Sample")
# Förbereder namn för utvalda komponenter
komp <- c("Sub5_Komp_","Sub10_Komp_","Sub15_Komp_")

sampleFixMap <- c("Sub_5_Kovarians","Sub_10_Kovarians","Sub_15_Kovarians")
sampleMap <- c("Subjekt5","Subjekt10","Subjekt15")
voxlar <- c(8349,7997,8075)

# skapar alla matriser för BTC och BSM för modell 2 och 3.
BSM_array_kor <- array(0, dim=c( 20,voxlar[s],100))
BTC_array_kor <- array(0, dim=c( 20,119 ,100))

BSM_array_wish <- array(0, dim=c( 20,voxlar[s] ,100))
BTC_array_wish <- array(0, dim=c( 20,119 ,100))

# väljer ut individ alltså s i början
if( s == 1 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject5
}
if( s == 2 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject10
}
if( s == 3 ){
  sub <- fMRI_data$X6_filtered_func_data_Subject15
}
# Väljer ut GICA samt hur många antal komponenter
# Gör en maskning likt innan
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

### skapar mall för visualisering
MALL <- matlab::reshape(sub_mall, c(30*36*30,119))
M <- reshape(MALL, c(30,36,30,119))
M <- M[,,,1]

# laddar in sample nr: 51:150
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


## FREK ## 
X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))


X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- BSM
Frek <- reshape(X_N_M1, c(30,36,30,20))

## Bayesiansk Korrelation ##

Greater_0_korr <- t(apply(BSM_array_kor,c(1,2), function(x) sum(x > 0) ))

Aktiv_Mall_korr<- Greater_0_korr[,] > 67
Posterior_Mean_BSM_korr <- t(apply(BSM_array_kor,c(1,2), mean))
Posterior_Mean_BSM_korr[!Aktiv_Mall_korr] <- 0
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM_korr
Mean_greater_0_BSM_korr <- reshape(X_N_M1, c(30,36,30,20))
## Bayesiansk Wishart ## 
Greater_0 <- t(apply(BSM_array_wish,c(1,2), function(x) sum(x > 0) ))

Aktiv_Mall<- Greater_0[,] > 50
Posterior_Mean_BSM <- t(apply(BSM_array_wish,c(1,2), mean))
Posterior_Mean_BSM[!Aktiv_Mall] <- 0
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM
Mean_greater_0_BSM <- reshape(X_N_M1, c(30,36,30,20))

# sätter wd() för bilder
setwd("search path")
## Range för alla plottar ## 
maxrange <- matrix(0,ncol=3,nrow=6)
for(k in 1:6){
  maxrange[k,1] <- max(range(BSM[BSM[,k] > 0  ,k]))
  maxrange[k,2] <- max(range(Mean_greater_0_BSM_korr[,,,k]))
  maxrange[k,3] <- max(range(Mean_greater_0_BSM[,,,k]))
}
ymin <- 0
ymax <- round(max(maxrange))

# beräknar varje komponent/nätverkskarta k för varje individ 1 till 3. (5,10,15)
for(k in 1:6){
  
for(i in 1:3){
  pdf(file=paste(c(komp[s],k,"_",i,".pdf"),collapse = ""))
  if(i == 1){
    
  ortho2(M,Frek[,,,k],zlim.y=range(BSM[BSM[,k] > 0  ,k]),mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(ymin,ymax,length.out = 65),col.y = rev(jet.colors(64)))
    
    
  }
  if(i == 2){
    
  ortho2(M,Mean_greater_0_BSM_korr[,,,k],mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(ymin,ymax,length.out = 65),col.y = rev(jet.colors(64)))
    
    
  }
  if(i == 3){
    
    ortho2(M,Mean_greater_0_BSM[,,,k],mfrow = c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(ymin,ymax,length.out = 65),col.y = rev(jet.colors(64)))
    
  }

  
  dev.off()
}


brainz1 <- image_read(paste(komp[s],k,"_1.pdf",sep=""))
brainz2 <- image_read(paste(komp[s],k,"_2.pdf",sep=""))
brainz3 <- image_read(paste(komp[s],k,"_3.pdf",sep=""))

brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
brainz3 <- image_crop(brainz3, geometry = geometry_area(550,190,0,155))

img <- c(brainz1,brainz2,brainz3)

imz <- image_append(img, stack = TRUE)

# Skapar en colorbar för alla nätverksbilder
colorbar <- function(lut, min, max=-min, nticks=11, title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4, seq(min, max, len=nticks),las=1,col="white",col.ticks="white",col.axis="white",cex.axis=1.3)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }	
}
#
pdf(file="colorbarjet.pdf")
par(mar=c(1,26,0,6),bg ="black")
colorbar(lut=rev(jet.colors(64)), min=0,max=25, nticks = 5)
dev.off()

colorbar <- image_read("colorbarjet.pdf")
colorbar <- image_chop(colorbar,geometry = geometry_area(300,5,5,0))
colorbar <- image_scale(colorbar,geometry = geometry_area(0,570,0,0))

brain_plot <-image_append(c(imz,colorbar), stack = F)


brain_plot <- image_annotate(brain_plot, "Frekventistisk - Modell 1", size = 20, color = "white",
               degrees = 0, location = "+10+2")

brain_plot <- image_annotate(brain_plot, "Bayesiansk - Modell 2", size = 20, color = "white",
                             degrees = 0, location = "+10+185")

brain_plot <- image_annotate(brain_plot, "Bayesiansk - Modell 3", size = 20, color = "white",
               degrees = 0, location = "+10+375")
pdf(file=paste(c(komp[s],k,".pdf"),collapse = ""), width = 5.5, height = 3.5)
par(mar=c(0,3,0,1),bg ="black")
plot(brain_plot)
dev.off()


}
}



#### För jämförelse av modell 2 och 3 #### 


library("magick")
setwd("search path")

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

setwd("search path")


Greater_0_korr <- t(apply(BSM_array_kor,c(1,2), function(x) sum(x > 0) ))

Greater_0 <- t(apply(BSM_array_wish,c(1,2), function(x) sum(x > 0) ))


#### Visualiserar sannolikheten för att bsm > 0 för alla voxlar per nätverkskarta

library("reshape2")
Greater_Bayes2<- melt(Greater_0_korr, na.rm = TRUE)
Greater_Bayes3 <- melt(Greater_0, na.rm = TRUE)
Greater_Bayes2$Var1 <- "Modell 2"
Greater_Bayes3$Var1 <- "Modell 3"

colnames(Greater_Bayes2)[1] <- "Modell"
colnames(Greater_Bayes3)[1] <- "Modell"

Greater <- rbind(Greater_Bayes2,Greater_Bayes3)
names(Greater)[2] <- "Nätverkskarta"


ggplot(subset(Greater, Nätverkskarta %in% 1:6), aes(x = value)) +
  stat_ecdf(aes(color = Modell,linetype = Modell), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+
  labs(y = "f(weight)")


h <- ggplot((subset(Greater, Nätverkskarta %in% 1:6)), aes(y=(..count..)/(length(Greater$Nätverkskarta)/40),x=value,color=Modell,fill=Modell)) + geom_histogram(position="identity", alpha=0.8,bins=100,binwidth = 1,size=0) + 
facet_wrap(vars(Nätverkskarta), ncol=2,labeller = label_both) + theme_minimal()  + scale_fill_manual(values=c("darkred", "steelblue")) +  scale_color_manual(values = c("darkred", "steelblue"))  


  h <- h + ggtitle("Sannolikheten för att Bsm > 0 \n över alla voxlar per nätverkskarta ") + labs(y="Sannolikhet", x="Antal posteriordragningar") + theme(plot.title = element_text(hjust = 0.5)) 

  
  pdf(file=paste(c(komp[s],"histogram_bayes",".pdf"),collapse = ""), width = 5.5, height = 3.5)
  h
  dev.off()  
  
  
### Beräknar ackumulerat posteriormedelvärde
  
  # model 2 
  mean_matrix <- matrix(nrow=100, ncol=6)
  for(k in 1:6){
    for(i in 1:100){
      mean_matrix[i,k] <-   mean(BSM_array_kor[k,,1:i])
      
    }
    
  }
  
  # model 3
  mean_matrix2 <- matrix(nrow=100, ncol=6)
  for(k in 1:6){
    for(i in 1:100){
      mean_matrix2[i,k] <-   mean(BSM_array_wish[k,,1:i])
      
    }
    
  }
  
  
Konv_Bayes2<- melt(mean_matrix, na.rm = TRUE)
Konv_Bayes3 <- melt(mean_matrix2, na.rm = TRUE)
Konv_Bayes2$Var1 <- "Modell 2"
Konv_Bayes3$Var1 <- "Modell 3"
colnames(Konv_Bayes2)[1] <- "Modell"
colnames(Konv_Bayes3)[1] <- "Modell"

Konv <- rbind(Konv_Bayes2,Konv_Bayes3)
Konv$dragning <- rep(1:100,12)
Konv

 
k <- ggplot(Konv, aes(y=value, x=dragning, color=Modell)) + geom_line(size=0.8) + facet_wrap(vars(Var2), ncol=2) + theme_minimal() +  scale_color_manual(values = c("darkred", "steelblue"))
k <- k + ggtitle("")+ labs(y="Ack. posteriormedelvärde", x="Antal posteriordragningar") + theme(plot.title = element_text(hjust = 0.5)) 


pdf(file=paste(c(komp[s],"Konvergens_bayes",".pdf"),collapse = ""), width = 7.5, height = 3.5)
k
dev.off()  


#### Visualiserar inv. wishart matrisers korrelationer i ett histogram.  

setwd("search path")
load("Sub15MeanWish.RData")
wish_mean  <- mov10


get_upper_tri <- function(correlation_matrix){
  correlation_matrix[lower.tri(correlation_matrix)]<- NA
  return(correlation_matrix)
}

upper_tri <- get_upper_tri(cov2cor(wish_mean))
library("reshape")
melted_cormat <- melt(upper_tri, na.rm = TRUE)
library("ggplot2")
h <- ggplot(melted_cormat, aes(x=value)) + geom_histogram(aes(y = (..count..)/sum(..count..)), fill="steelblue",position="identity", alpha=0.5) + theme_minimal()
h <- h + ggtitle("Histogram")+ labs(y="Relativ frekvens", x="Korrelation") + theme(plot.title = element_text(hjust = 0.5)) 

pdf(file= "Sub15KorrHist.pdf", width = 7.5, height = 3.5)
h
dev.off()

load("sub15_ack_sd.RData")
wish_sd <- sd11
dim(wish_mean)
wish_tkvot <- wish_mean/sqrt(wish_sd/99)


upper_tri <- get_upper_tri(wish_tkvot)
library("reshape")
melted_cormat <- melt(upper_tri, na.rm = TRUE)
library("ggplot2")
k <- ggplot(melted_cormat, aes(x=value)) + geom_histogram(aes(y = (..count..)/sum(..count..)),fill="steelblue",position="identity", alpha=0.5) + theme_minimal()
k <- k + ggtitle("Histogram")+ labs(y="Relativ frekvens", x="Bayesiansk T-Kvot") + theme(plot.title = element_text(hjust = 0.5)) 

pdf(file= "Sub15TkvotHist.pdf", width = 7.5, height = 3.5)
k
dev.off()

### 

