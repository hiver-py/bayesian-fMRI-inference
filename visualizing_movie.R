# Skapar film med ackumuleratmedelvärde

# Subjekt 5: 8349
# Subjekt 10:7997
# Subjekt 15: 8075
s <- 1
k <- 1


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



X_N_M1 <- matlab::reshape(GICA20, c(dim_G[1]*dim_G[2]*dim_G[3],dim_G[4]))


Posterior_Mean_BSM <- t(apply(BSM_array_wish,c(1,2), mean))
pos <- Posterior_Mean_BSM>0
Posterior_Mean_BSM[!pos] <- 0 

X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM
aktiv_mall <- reshape(X_N_M1, c(30,36,30,20))

aktiv_pos <- aktiv_mall[,,15,k] > 0


ackmedel <- matrix(0,ncol=sum(sum(aktiv_pos)), nrow=100)
library("neurobase")

setwd("search path")
i <- 1
for(i in 1:100){
  BSMWISH <- t(apply(BSM_array_wish[,,1:i], c(1,2), mean))
  X_N_M1[X_N_M1 != 0] <- 0
  X_N_M1[row.sub,] <- BSMWISH
  
  
  BSMwish <- reshape(X_N_M1, c(30,36,30,20))
  ackmedel[i,]<- BSMwish[,,15,k][aktiv_pos]

  png(file=paste(c("Brain",i,".png"),collapse = ""), width = 450, height = 250)
try(  ortho2(M,BSMwish[,,,k],mfrow = c(1,1),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,25,length.out = 65),col.y = rev(jet.colors(64))), silent = TRUE )
  dev.off()

  df <- melt(ackmedel)
  
  line <- ggplot(data.frame(subset(df,Var1 %in% 1:i )), aes(x = Var1,group=Var2, y = value)) + 
    geom_line(size=0.2, alpha=0.4) + theme_minimal() + ggtitle("Konvergens för ackumulerade posteriormedelvärdet")+ labs(y="Ack. posteriormedelvärde", x="Antal posteriordragningar") + theme(plot.title = element_text(hjust = 0.5))  + scale_y_continuous(breaks= pretty_breaks())
  png(file=paste(c("Mean",i,".png"),collapse = ""), width = 450, height = 250)
  plot(line)
  dev.off()
  

}


for( p in 1:100){
  
  MEANS <- readPNG(paste(c("Mean",p,".png"),collapse = ""))
  BRAINS <- readPNG(paste(c("Brain",p,".png"),collapse = ""))
  IMG <- list(BRAINS,MEANS)
  png(file=paste(c("MOV",p,".png"),collapse = ""))
  
  layout(matrix(1:2, ncol=1, byrow=TRUE))
  for(i in 1:2) {
    par(mar=c(0,0,0,0),bg="black")
    plot(NA,xlim=0:1,ylim=0:1,xaxt="n",yaxt="n",bty="n",ann=FALSE)
    rasterImage(IMG[[i]],0,0,1,1)
  }
  dev.off()
  

  
}

png_files <- sprintf("MOV%d.png",1:100)
av::av_encode_video(png_files , 'sub15_NK6.mp4', framerate = 3.5)



## Skapar stillbild för sista dragningen att ha i uppsats ##


k <- 6


X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- Posterior_Mean_BSM
aktiv_mall <- reshape(X_N_M1, c(30,36,30,20))

aktiv_pos <- aktiv_mall[,,15,k] > 0
ackmedel <- matrix(0,ncol=sum(sum(aktiv_pos)), nrow=100)

for(i in 1:100){
  BSMWISH <- t(apply(BSM_array_wish[,,1:i], c(1,2), mean))
  X_N_M1[X_N_M1 != 0] <- 0
  X_N_M1[row.sub,] <- BSMWISH
  BSMwish <- reshape(X_N_M1, c(30,36,30,20))
  ackmedel[i,]<- BSMwish[,,15,k][aktiv_pos]
}

i <- 100


BSMWISH <- t(apply(BSM_array_wish[,,1:i], c(1,2), mean))
X_N_M1[X_N_M1 != 0] <- 0
X_N_M1[row.sub,] <- BSMWISH


BSMwish <- reshape(X_N_M1, c(30,36,30,20))
#ackmedel[i,]<- BSMwish[,,15,1][aktiv_pos]
setwd("search path")

pdf(file=paste(c("Brain",i,".pdf"),collapse = ""), width = 10, height = 3.5)
try(  ortho2(M,BSMwish[,,,k],mfrow = c(1,1),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,25,length.out = 65),col.y = rev(jet.colors(64))), silent = TRUE )
dev.off()

df <- melt(ackmedel)

line <- ggplot(data.frame(subset(df,Var1 %in% 1:i )), aes(x = Var1,group=Var2, y = value)) + 
  geom_line(size=0.2, alpha=0.4) + theme_minimal() + labs(y="Ack. posteriormedelvärde", x="Antal posteriordragningar") + theme(plot.title = element_text(hjust = 0.5))  + scale_y_continuous(breaks= pretty_breaks())
pdf(file=paste(c("Mean",i,".pdf"),collapse = ""), width = 10, height = 3.5)
par(mar=c(0,0,0,0),bg ="black")
plot(line)
dev.off()


B <- image_read_pdf("Brain100.pdf",pages=3)

M <- image_read_pdf("Mean100.pdf")


konv <- c(B,M)
konv <- image_append(konv, stack = T)

#konv <- image_resize(konv, geometry = geometry_area(750,350,0,0))

image_write(konv, path = "Sub15_konv6.pdf", format = "pdf",quality=100)

image_read_pdf("Brain100.pdf")




