
#findpeaks and findvalleys are derived from findPeaks and findValleys, respectively,
#from quantmod R package

findpeaks<- function(x){
  which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 
          0) + 1
}

findvalleys<-function(x){
  which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) > 
          0) + 1
}


peaks<- function(x,y){
  list<- c()
for (i in 1:7){
  g<- (sign((x[i+1]+y[i+1]) -  (x[i]-y[i])))
  #print(g)
  list[[length(list)+1]]<- g
}



#which(diff(as.numeric(list))<0) +2 -1 #valleys
which(diff(as.numeric(list))<0) +2 -1 #peaks
}

valleys<- function(x,y){
  list<- c()
  for (i in 1:7){
    g<- (sign((x[i+1]+y[i+1]) -  (x[i]-y[i])))
    #print(g)
    list[[length(list)+1]]<- g
  }
  
 
  
  #which(diff(as.numeric(list))<0) +2 -1 #valleys
  which(diff(as.numeric(list))>0) +2 -1 #peaks
}
