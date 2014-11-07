############ Statistics of the times measured #################

setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/")
gpu <- dir()[file.info(dir())$isdir]

for(i in gpu) {
  files <- list.files(paste("./", i, sep=""), pattern=".txt")
  
  for(j in files) {    
    temp <- as.matrix(read.table(paste("./", i, "/", j, sep=""), sep="\t", header=F, fill = TRUE))
    
    for(k in 1:length(temp[,1])) {
      print(paste("GeForce ", i, ", App ", j , ", Size No: ", k, sep=""))
      GPUTime <- temp[k,1:10]
      print(summary(GPUTime))
      print(try(t.test(GPUTime, alternative = "two.sided", conf.level = 0.95)))
    }
  }
}
