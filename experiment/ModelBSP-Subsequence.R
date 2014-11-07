############ Maximum Subsequence GLobal Memory #################
### BSP MULTILEVEL MEMORY

# The palette of color
cbbPalette <- c("red", "blue", "darkgray", "orange","black","lightblue", "lightblue","violet")

Nn <- 17:29;
N <- 2^Nn;

gridsize <- 32;
blocksize <- 128

numberthreads <- gridsize * blocksize;

N_perBlock <- N/gridsize;
N_perThread <- N_perBlock/blocksize;

# Miliseconds
globalLoadTrans <- numberthreads*N_perThread/32;
globalStoreTrans <- numberthreads*5;

# Miliseconds
sharedLoadTrans <- numberthreads*N_perThread;
sharedStoreTrans <- numberthreads*N_perThread;

latencySharedMemory <- 5; #Cycles per processor
latencyGlobalMemory <- 100 * latencySharedMemory; #Cycles per processor

#Cycles operations per Thread operations

timeComputationKernel <- 1000 * numberthreads * N_perThread;
latencyCommunication <- ((globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory)


SubSeqMax_gt630 <- 1:12;
SubSeqMax_gtx660 <- 1:12;
SubSeqMax_gtx680 <- 1:12;
SubSeqMax_Titan <- 1:12;
SubSeqMax_Tesla <- 1:12;

##### GeForce GT-630  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/gt-630/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gt630[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 96;
clockFrequency_GPU <- 1620; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_timeKernel_SubSeqMax_630 <- (timeComputationKernel + latencyCommunication)/(flopsTheoreticalpeakGPU*10^6)*10^3;

SpeedupSubSeqMax_630 <- timeKernel_timeKernel_SubSeqMax_630[1:12]/SubSeqMax_gt630;



##### GeForce GTX-660  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/gtx-660/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gtx660[i] <- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 960;
clockFrequency_GPU <- 1058; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

latencyCommunication <- (globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory;

timeKernel_timeKernel_SubSeqMax_660 <- (timeComputationKernel + latencyCommunication)/(flopsTheoreticalpeakGPU*10^6)*10^3;

SpeedupSubSeqMax_660 <- timeKernel_timeKernel_SubSeqMax_660[1:12]/SubSeqMax_gtx660;


##### GeForce GTX-680  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/gtx-680/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gtx680[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 1536;
clockFrequency_GPU <- 1006; # Mhz  

flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

latencyCommunication <- (globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory;

timeKernel_timeKernel_SubSeqMax_680 <- (timeComputationKernel + latencyCommunication)/(flopsTheoreticalpeakGPU*10^6)*10^3;

SpeedupSubSeqMax_680 <- timeKernel_timeKernel_SubSeqMax_680[1:12]/SubSeqMax_gtx680;


##### GeForce GTX-Titan  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/gtx-Titan/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_Titan[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 2688;
clockFrequency_GPU <- 876; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

latencyCommunication <- (globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory;

timeKernel_timeKernel_SubSeqMax_Titan <- (timeComputationKernel + latencyCommunication)/(flopsTheoreticalpeakGPU*10^6)*10^3;

SpeedupSubSeqMax_Titan <- timeKernel_timeKernel_SubSeqMax_Titan[1:12]/SubSeqMax_Titan;



##### Tesla K-20  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/Tesla-k20/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_Tesla[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 2496;
clockFrequency_GPU <- 706; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

latencyCommunication <- ((globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory);

timeKernel_timeKernel_SubSeqMax_Tesla <- (timeComputationKernel + latencyCommunication)/(flopsTheoreticalpeakGPU*10^6)*10^3;

SpeedupSubSeqMax_Tesla <- timeKernel_timeKernel_SubSeqMax_Tesla[1:12]/SubSeqMax_Tesla;


############ Difference between predicted and measured SubSequenceMaximum ##############

setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/")
png(filename="./SubSeMax-GPU.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:12], SpeedupSubSeqMax_630, type="l",  log="x", lty = 1, lwd=c(5,5), xaxt="n",
     ylim=c(0.25, 3.25), xlim=c(131072, 268435456),
     col=cbbPalette[1], ylab = "Time Predicted Vs. Measured", cex.axis=1.5, cex.lab=2,cex.main=2,
     xlab = "Size of the sequence ",  main = paste("Speedup Maximum subarray problem \n ", "in GPUs with Kepler Architectures", sep=""));
points(N[1:12], SpeedupSubSeqMax_630, col = cbbPalette[1], type = "p", pch=20,cex = 2)

lines(N[1:12], SpeedupSubSeqMax_660[1:12], col = cbbPalette[2],lwd=c(5,5));
points(N[1:12], SpeedupSubSeqMax_660, col = cbbPalette[2], pch=21,cex = 2);

lines(N[1:12], SpeedupSubSeqMax_680[1:12], col = cbbPalette[3],lwd=c(5,5));
points(N[1:12], SpeedupSubSeqMax_680, col = cbbPalette[3], pch=22,cex = 2);

lines(N[1:12], SpeedupSubSeqMax_Titan[1:12], col = cbbPalette[4],lwd=c(5,5));
points(N[1:12], SpeedupSubSeqMax_Titan, col = cbbPalette[4], pch=23,cex = 2);

lines(N[1:12], SpeedupSubSeqMax_Tesla[1:12], col = cbbPalette[5],lwd=c(5,5));
points(N[1:12], SpeedupSubSeqMax_Tesla, col = cbbPalette[5], pch=24,cex = 2);

#axis(1, 2^pow, cex.axis = 1.5,lwd=c(5,5));

axis(1, at = c(N[1:12]), labels = paste('2^',log2(c(N[1:12])),sep="") , cex.axis=1.5)
#axis(2, at=c(SpeedupMatSumGmUnSP_GT630));

grid()

legend('topright', 
       lty=1, col=c(cbbPalette[1:5]), pch=c(20:24), lwd=c(5,5), cex=2, 
       legend=c("GT-630", "GTX-660", "GT-680","GTX-Titan", "Tesla-K20"), title="Speedup")

dev.off()



setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/")
png(filename="./SubSeMax-GPU-Times.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:12], SubSeqMax_gt630, type="l",  log="xy", lty = 1, lwd=c(5,5), xaxt="n",
     ylim=c(0.1, max(SubSeqMax_gt630)), xlim=c(131072, 268435456),
     col=cbbPalette[1], ylab = "Time Predicted Vs. Measured ( % )", cex.axis=1.5, cex.lab=2,cex.main=2,
     xlab = "Size of the sequence ",  main = paste("Speedup Maximum subarray problem \n ", "in GPUs with Kepler Architectures", sep=""));
points(N[1:12], SubSeqMax_gt630, col = cbbPalette[1], type = "p", pch=20,cex = 2)

lines(N[1:12], SubSeqMax_gtx660[1:12], col = cbbPalette[2],lwd=c(5,5));
points(N[1:12], SubSeqMax_gtx660, col = cbbPalette[2], pch=21,cex = 2);

lines(N[1:12], SubSeqMax_gtx680[1:12], col = cbbPalette[3],lwd=c(5,5));
points(N[1:12], SubSeqMax_gtx680, col = cbbPalette[3], pch=22,cex = 2);

lines(N[1:12], SubSeqMax_Titan[1:12], col = cbbPalette[4],lwd=c(5,5));
points(N[1:12], SubSeqMax_Titan, col = cbbPalette[4], pch=23,cex = 2);

lines(N[1:12], SubSeqMax_Tesla[1:12], col = cbbPalette[5],lwd=c(5,5));
points(N[1:12], SubSeqMax_Tesla, col = cbbPalette[5], pch=24,cex = 2);

#axis(1, 2^pow, cex.axis = 1.5,lwd=c(5,5));

axis(1, at = c(N[1:12]), labels = paste('2^',log2(c(N[1:12])),sep="") , cex.axis=1.5)
#axis(2, at=c(SpeedupMatSumGmUnSP_GT630));

grid()

legend('bottomright', 
       lty=1, col=c(cbbPalette[1:5]), pch=c(20:24), lwd=c(5,5), cex=2, 
       legend=c("GT-630", "GTX-660", "GT-680","GTX-Titan", "Tesla-K20"), title="Speedup")

dev.off()

