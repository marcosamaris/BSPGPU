############ Maximum Subsequence GLobal Memory #################
### BSP MULTILEVEL MEMORY

# The palette of color
cbbPalette <- gray(1:6/ 8)#c("red", "blue", "darkgray", "orange","black","brown", "lightblue","violet")

Nn <- 17:29;
N <- 2^Nn;

gridsize <- 32;
blocksize <- 128

numberthreads <- gridsize * blocksize;

N_perBlock <- N/gridsize;
N_perThread <- N_perBlock/blocksize;

latencySharedMemory <- 5; #Cycles per processor
latencyGlobalMemory <- latencySharedMemory* 100; #Cycles per processor

latencyL1 <- latencySharedMemory; #Cycles per processor
latencyL2 <- latencyGlobalMemory*0.5; #Cycles per processor

SubSeqMax_gt630 <- 1:12;
SubSeqMax_gtx660 <- 1:12;
SubSeqMax_gtx680 <- 1:12;
SubSeqMax_Titan <- 1:12;
SubSeqMax_Tesla <- 1:12;
SubSeqMax_Tesla_k40 <- 1:12;

##### GeForce GT-630  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gt-630/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gt630[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 96;
clockFrequency_GPU <- 1620; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0.1*reads
#L2Effect <- c(0.1, 0.15, 0.1, 0.05, 0.05, 0.05, 0.1, 0.05, 0.025, 0 ,0,0,0)*reads
L2Effect <- 0.1*reads
W <- 1.1

CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_630 <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_630 <- timeKernel_timeKernel_SubSeqMax_630[1:12]/SubSeqMax_gt630;

##### GeForce GTX-660  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gtx-660/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gtx660[i] <- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 960;
clockFrequency_GPU <- 1058; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0
L2Effect <- 0.1*reads
W <- .65


CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_660 <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_660 <- timeKernel_timeKernel_SubSeqMax_660[1:12]/SubSeqMax_gtx660;


##### GeForce GTX-680  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gtx-680/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_gtx680[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 1536;
clockFrequency_GPU <- 1006; # Mhz  

flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0
L2Effect <- 0.1*reads
#  W <- .76

CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_680 <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_680 <- timeKernel_timeKernel_SubSeqMax_680[1:12]/SubSeqMax_gtx680;


##### GeForce GTX-Titan  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gtx-Titan/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_Titan[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 2688;
clockFrequency_GPU <- 876; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0
L2Effect <- 0.1*reads
# W <- .55

CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_Titan <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_Titan <- timeKernel_timeKernel_SubSeqMax_Titan[1:12]/SubSeqMax_Titan;



##### Tesla K-20  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/tesla-k20/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_Tesla[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 2496;
clockFrequency_GPU <- 706; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0
L2Effect <- 0.1*reads
#W <- 0.63

CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_Tesla <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_Tesla <- timeKernel_timeKernel_SubSeqMax_Tesla[1:12]/SubSeqMax_Tesla;

##### Tesla K-40  #####
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/tesla-k40/")
temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax_Tesla_k40[i]<- mean(temp[i,1:10])
}

PCoresNumber_GPU <- 2880;
clockFrequency_GPU <- 745; # Mhz  

#bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeakGPU <- clockFrequency_GPU * PCoresNumber_GPU ; # Mflops/second

#Cycles operations per Thread operations
timeComputationKernel <-  100*numberthreads * N_perThread;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
reads <-  numberthreads*N_perThread;

L1Effect <- 0
L2Effect <- 0.1*reads
#W <- 0.52

CommGM <- ((numberthreads*N_perThread - L1Effect - L2Effect + numberthreads*5)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

CommSM <- (numberthreads*N_perThread + numberthreads*5)*latencySharedMemory
timeKernel_timeKernel_SubSeqMax_Tesla_k40 <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeakGPU*10^6))*10^3;

SpeedupSubSeqMax_Tesla_k40 <- timeKernel_timeKernel_SubSeqMax_Tesla_k40[1:12]/SubSeqMax_Tesla_k40;

############ Difference between predicted and measured SubSequenceMaximum ##############
dataN <- 4:12
setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/")
png(filename="./SubSeqMax-GPU.png", width=800, height=600)
par(mar=c(1, 4, 2, 1) + 0.1)
layout(rbind(1,2), heights=c(15,1))  # put legend on bottom 1/8th of the chart
plot(N[dataN], SpeedupSubSeqMax_630[dataN], type="l",  log="x", lty = 1, lwd=c(7.5,7.5), xaxt="n",
     ylim=c(0.7, 1.2), xlim=c(1048576, 268435456),
     col=cbbPalette[1], ylab = " ", cex.axis=3.5, cex.lab=3.5,cex.main=3.5,
     xlab = " ",  main = paste(" ", sep=""));
points(N[dataN], SpeedupSubSeqMax_630[dataN], col = cbbPalette[1], type = "p", pch=20,cex = 3.5)

lines(N[dataN], SpeedupSubSeqMax_660[dataN], col = cbbPalette[2], lty = 2,lwd=c(7.5,7.5));
points(N[dataN], SpeedupSubSeqMax_660[dataN], col = cbbPalette[2], pch=21,cex = 3.5);

lines(N[dataN], SpeedupSubSeqMax_680[dataN], col = cbbPalette[3], lty = 3,lwd=c(7.5,7.5));
points(N[dataN], SpeedupSubSeqMax_680[dataN], col = cbbPalette[3], pch=22,cex = 3.5);

lines(N[dataN], SpeedupSubSeqMax_Titan[dataN], col = cbbPalette[4], lty = 4,lwd=c(7.5,7.5));
points(N[dataN], SpeedupSubSeqMax_Titan[dataN], col = cbbPalette[4], pch=23,cex = 3.5);

lines(N[dataN], SpeedupSubSeqMax_Tesla[dataN], col = cbbPalette[5], lty = 5,lwd=c(7.5,7.5));
points(N[dataN], SpeedupSubSeqMax_Tesla[dataN], col = cbbPalette[5], pch=24,cex = 3.5);

lines(N[dataN], SpeedupSubSeqMax_Tesla_k40[dataN], col = cbbPalette[6], lty = 6,lwd=c(5,5));
points(N[dataN], SpeedupSubSeqMax_Tesla_k40[dataN], col = cbbPalette[6], pch=25,cex = 3.5);
#axis(1, 2^pow, cex.axis = 1.5,lwd=c(5,5));

#   axis(1, at = c(N[dataN]), labels = paste('2^',log2(c(N[dataN])),sep="") , cex.axis=1)
#axis(2, at=c(SpeedupMatSumGmUnSP_GT630));

grid()
par(mar=c(2, 2, 0, 0))
plot.new()
legend('center', 
       col=c(cbbPalette[1:6]), lty=c(1:6),pch=c(20:25), lwd=c(5,5), cex=1, ncol=6,bty ="n",
       legend=c("GT-630", "GTX-660", "GTX-680","GTX-Titan", "Tesla-K20","Tesla-K0"))

dev.off()



setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/")
png(filename="./SubSeqMax-GPU-Times.png", width=800, height=600)
par(mar=c(1, 4, 2, 1) + 0.1)
layout(rbind(1,2), heights=c(15,1))  # put legend on bottom 1/8th of the chart
plot(N[dataN], SubSeqMax_gt630[dataN], type="l",  log="xy", lty = 1, lwd=c(5,5), xaxt="n",
     ylim=c(0.1, max(SubSeqMax_gt630[dataN])), xlim=c(1048576, 268435456),
     col=cbbPalette[1], ylab = " ", cex.axis=2.5, cex.lab=3,cex.main=3.5,
     xlab = " ",  main = paste(" ", sep=""));
points(N[dataN], SubSeqMax_gt630[dataN], col = cbbPalette[1], type = "p", pch=20,cex = 3.5)

lines(N[dataN], SubSeqMax_gtx660[dataN], col = cbbPalette[2], lty = 2,lwd=c(7.5,7.5));
points(N[dataN], SubSeqMax_gtx660[dataN], col = cbbPalette[2], pch=21,cex = 3.5);

lines(N[dataN], SubSeqMax_gtx680[dataN], col = cbbPalette[3], lty = 3,lwd=c(5,5));
points(N[dataN], SubSeqMax_gtx680[dataN], col = cbbPalette[3], pch=22,cex = 3.5);

lines(N[dataN], SubSeqMax_Titan[dataN], col = cbbPalette[4], lty =4,lwd=c(7.5,7.5));
points(N[dataN], SubSeqMax_Titan[dataN], col = cbbPalette[4], pch=23,cex = 3.5);

lines(N[dataN], SubSeqMax_Tesla[dataN], col = cbbPalette[5], lty = 5,lwd=c(7.5,7.5));
points(N[dataN], SubSeqMax_Tesla[dataN], col = cbbPalette[5], pch=24,cex = 3.5);

lines(N[dataN], SubSeqMax_Tesla_k40[dataN], col = cbbPalette[6], lty = 6,lwd=c(7.5,7.5));
points(N[dataN], SubSeqMax_Tesla_k40[dataN], col = cbbPalette[6], pch=25,cex = 3.5);

#axis(1, 2^pow, cex.axis = 1.5,lwd=c(5,5));

axis(1, at = c(N[dataN]), labels = paste('2^',log2(c(N[dataN])),sep="") , cex.axis=1)
#axis(2, at=c(SpeedupMatSumGmUnSP_GT630));

grid()
#   par(mar=c(2, 2, 0, 0))
#   plot.new()
#   legend('center', 
#          lty=1, col=c(cbbPalette[1:6]), pch=c(20:25), lwd=c(5,5), cex=2.2, ncol=6,bty ="n",
#          legend=c("GT-630", "GTX-660", "GT-680","GTX-Titan", "Tesla-K20", "Tesla-K40"))

dev.off()
