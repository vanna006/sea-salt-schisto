install.packages(vegan)
library(vegan)

#set to your wd 
setwd("F:/PhD/Undergrads/Anna_Ao_Yu")
salt <- read.csv("saltNMDS.csv")

best <- metaMDS(salt[,2:22], distance = "euclidean", k=2, trymax=100)		## Ordinate in 2 dim, 20 tries for best solution
plot(best, type = 't')
text(best$point[,1], best$points[,2], labels = salt$Substance)

dist(cbind(best$point[,1], best$points[,2]))
