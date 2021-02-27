# Visualisering av alla GICA 

GICA20 <- fMRI_data$X6_melodic_IC[,,,1:20]
mask <- sub < mean(sub[sub != 0])
mask <- mask[,,,1]
sub_mall <- sub
sub[mask] <- 0
GICA20[mask] <- 0

MALL <- matlab::reshape(sub_mall, c(30*36*30,119))
M <- reshape(MALL, c(30,36,30,119))
M <- M[,,,1]

### VISUALISERING AV GICA ###
dev.off()
hist(GICA20)
setwd("search path")

for(k in 1:20){
  pdf(file=paste(c("GICA",k,".pdf"),collapse = ""))
  
ortho2(M,GICA20[,,,k],mfrow=c(1,3),NA.y=T,NA.x=T,crosshairs = F,ybreaks=seq(0,max(GICA20),length.out = 65)
       ,col.y = rev(jet.colors(64)))

dev.off()
}

###
brainz1 <- image_read(paste("GICA",1,".pdf",sep=""))
brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
img <- c(brainz1)
for(i in 2:4){
  brainz2 <- image_read(paste("GICA",i,".pdf",sep=""))
  brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
  
  img <- c(img,brainz2)
  img <- image_append(img, stack = F)
  
}
### 

###
brainz1 <- image_read(paste("GICA",5,".pdf",sep=""))
brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
img2 <- c(brainz1)
for(i in 6:8){
  brainz2 <- image_read(paste("GICA",i,".pdf",sep=""))
  brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
  
  img2 <- c(img2,brainz2)
  img2 <- image_append(img2, stack = F)
  
}
### 

###
brainz1 <- image_read(paste("GICA",9,".pdf",sep=""))
brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
img3 <- c(brainz1)
for(i in 10:12){
  brainz2 <- image_read(paste("GICA",i,".pdf",sep=""))
  brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
  
  img3 <- c(img3,brainz2)
  img3 <- image_append(img3, stack = F)
  
}
### 

###
brainz1 <- image_read(paste("GICA",13,".pdf",sep=""))
brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
img4 <- c(brainz1)
for(i in 14:16){
  brainz2 <- image_read(paste("GICA",i,".pdf",sep=""))
  brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
  
  img4 <- c(img4,brainz2)
  img4 <- image_append(img4, stack = F)
  
}
### 

brainz1 <- image_read(paste("GICA",17,".pdf",sep=""))
brainz1 <- image_crop(brainz1, geometry = geometry_area(550,190,0,155))
img5 <- c(brainz1)
for(i in 18:20){
  brainz2 <- image_read(paste("GICA",i,".pdf",sep=""))
  brainz2 <- image_crop(brainz2, geometry = geometry_area(550,190,0,155))
  
  img5 <- c(img5,brainz2)
  img5 <- image_append(img5, stack = F)
  
}
### 



gica <- c(img,img2,img3,img4,img5)
gica <- image_append(gica, stack = T)

gica <- image_annotate(gica, "Komponent 1", size = 20, color = "white",
                             degrees = 0, location = "+10+1")

gica <- image_annotate(gica, "Komponent 2", size = 20, color = "white",
                       degrees = 0, location = "+510+1")

gica <- image_annotate(gica, "Komponent 3", size = 20, color = "white",
                       degrees = 0, location = "+1010+1")

gica <- image_annotate(gica, "Komponent 4", size = 20, color = "white",
                       degrees = 0, location = "+1510+1")

##
image_info(gica)
190 * 4
gica <- image_annotate(gica, "Komponent 5", size = 20, color = "white",
                       degrees = 0, location = "+10+190")

gica <- image_annotate(gica, "Komponent 6", size = 20, color = "white",
                       degrees = 0, location = "+510+190")

gica <- image_annotate(gica, "Komponent 7", size = 20, color = "white",
                       degrees = 0, location = "+1010+190")

gica <- image_annotate(gica, "Komponent 8", size = 20, color = "white",
                       degrees = 0, location = "+1510+190")

##

##

gica <- image_annotate(gica, "Komponent 9", size = 20, color = "white",
                       degrees = 0, location = "+10+380")

gica <- image_annotate(gica, "Komponent 10", size = 20, color = "white",
                       degrees = 0, location = "+510+380")

gica <- image_annotate(gica, "Komponent 11", size = 20, color = "white",
                       degrees = 0, location = "+1010+380")

gica <- image_annotate(gica, "Komponent 12", size = 20, color = "white",
                       degrees = 0, location = "+1510+380")

##

gica <- image_annotate(gica, "Komponent 13", size = 20, color = "white",
                       degrees = 0, location = "+10+570")

gica <- image_annotate(gica, "Komponent 14", size = 20, color = "white",
                       degrees = 0, location = "+510+570")

gica <- image_annotate(gica, "Komponent 15", size = 20, color = "white",
                       degrees = 0, location = "+1010+570")

gica <- image_annotate(gica, "Komponent 16", size = 20, color = "white",
                       degrees = 0, location = "+1510+570")

##

gica <- image_annotate(gica, "Komponent 17", size = 20, color = "white",
                       degrees = 0, location = "+10+760")

gica <- image_annotate(gica, "Komponent 18", size = 20, color = "white",
                       degrees = 0, location = "+510+760")

gica <- image_annotate(gica, "Komponent 19", size = 20, color = "white",
                       degrees = 0, location = "+1010+760")

gica <- image_annotate(gica, "Komponent 20", size = 20, color = "white",
                       degrees = 0, location = "+1510+760")
##



pdf(file=paste(c("all_gica",".pdf"),collapse = ""), width = 5.5, height = 2.5)
par(mar=c(0,0,0,0),bg ="black")
plot(gica)
dev.off()




