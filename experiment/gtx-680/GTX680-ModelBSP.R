setwd("/home/marcos/Dropbox/Doctorate/Results/Version3/gtx-680/")

GPU <- "Geforce GTX-680"

# The palette of color
cbbPalette <- c("red", "blue", "darkgray", "orange","black","lightblue", "lightblue","violet")

matMulGmUnSP <- 0
matMulGmCoSP <- 0
matMulSmUnSP <- 0
matMulSmCoSP <- 0
SubSeqMax <- 0

temp <- as.matrix(read.table("./matMul-Gm-Un-SP.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:7){
  matMulGmUnSP[i]<- mean(temp[i,1:10])
}

temp <- as.matrix(read.table("./matMul-Gm-SP.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:7){
  matMulGmCoSP[i]<- mean(temp[i,1:10])
}


temp <- as.matrix(read.table("./matMul-Sm-Un-SP.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:7){
  matMulSmUnSP[i]<- mean(temp[i,1:10])
}

temp <- as.matrix(read.table("./matMul-Sm-SP.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:7){
  matMulSmCoSP[i]<- mean(temp[i,1:10])
}


temp <- as.matrix(read.table("./SubSeqMax.txt", sep="\t", header=F,fill = TRUE))
for (i in 1:12){
  SubSeqMax[i]<- mean(temp[i,1:10])
}


Nn <- 7:14;
N <- 2^Nn;

syzeBytes <- N*2*4; # Bytes
numberMultiplication <- N;
numberSumas <- N;

tileWidth <- 16;
threadsPerBlock <- tileWidth*tileWidth;

gridsizes <- as.integer((N +  tileWidth -1)/tileWidth);
blocknumber <- gridsizes*gridsizes
numberthreads <- threadsPerBlock * blocknumber;


PCoresNumber_GT630 <- 1536;
clockFrequency_GT630 <- 1058; # Mhz	
cycleMultiplicationKepler <- 16;
cycleSumKepler <- 4;
cycleFMA <- 2

memoryBusWidth_GT630 <- 256/8; # Bytes
memoryClockRate_GT630 <- 3004; # Mhz
latencySharedMemory <- 5; #Cycles per processor
latencyGlobalMemory <- latencySharedMemory* 100; #Cycles per processor

bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeaK_GT630 <- clockFrequency_GT630 * PCoresNumber_GT630; # Mflops/second


GDDR5SizeMB_GT630 <- 2048; # MB
L2CacheSize_GT630 <- 512; # Kb
sharedMemorySize_GT630 <- 48;  # Kb
L1CacheSizeKB_GT630 <- 16; # Kb
registersPerBlock_GT630 <- 64; # Kb




### BSP MULTILEVEL MEMORY

#Cycles operations per Threadoperation
tempOperationcycles <- ((numberMultiplication * cycleFMA) + (0 * cycleSumKepler)) * numberthreads;
# Miliseconds
timeComputationKernel <-  (tempOperationcycles/ (flopsTheoreticalpeaK_GT630*10^6))*10^3; 
# Miliseconds
#communicationLatencyKernel <- (syzeBytes + 16)*numberthreads/(bandWidth_GT630*10^6)*10^3; 
# Miliseconds

globalLoadTrans <- numberthreads*N/2;
globalStoreTrans <- numberthreads;

processorLatencyGlobal <- (globalStoreTrans*latencyGlobalMemory  + globalLoadTrans*latencyGlobalMemory)/(flopsTheoreticalpeaK_GT630*10^6)*10^3;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Un <- timeComputationKernel + processorLatencyGlobal;

SpeedupMatMulGmUnSP_GT630 <- timeKernel_GM_Un[1:7]/matMulGmUnSP;

############ MatMul GLobal Memory UnCoalesced 

png(filename="./images/MatMul-Gm-Un-SP.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:7], matMulGmUnSP[1:7],type="l",log="xy", lty = 2, lwd=c(5,5), yaxt="n",xaxt="n", 
     ylim=c(min(matMulGmUnSP,timeKernel_GM_Un), c(max(matMulGmUnSP,timeKernel_GM_Un))), xlim=c(128, 8192), 
     col = "blue", ylab = "Time (milisecond)", xlab = "Size of the Matrices",  cex.lab=2,cex.main=2,
     main = paste("Matrix Multiplication - Global Memory Uncoalesced \n", GPU, sep =""), axes = FALSE);

lines(N, timeKernel_GM_Un, col = "red",lwd=c(5,5));
points(N, timeKernel_GM_Un, col = "red", type = "p",pch=20,cex = 2)

lines(N[1:7], matMulGmUnSP[1:7], col = "blue",lwd=c(5,5));
points(N[1:7], matMulGmUnSP, col = "blue", pch=21,cex = 2);

axis(1, at = c(N), labels = c(N), cex.axis = 1.5,lwd=c(5,5))

y <- floor(log10(range(matMulGmUnSP,timeKernel_GM_Un)))
pow <- seq(y[1], y[2]+1);
axis(2, 10^pow, cex.axis = 1.5,lwd=c(5,5));

legend('bottomright', c("Predicted time","Measured time"), 
       lty=1, col=c("red","blue"), pch=c(20,21), lwd=c(5,5), cex=2)

grid()
dev.off()



############ MatMul GLobal Memory Coalesced #################

globalLoadTrans <- numberthreads*N/16;
globalStoreTrans <- numberthreads/8;

processorLatencyGlobal <- (globalStoreTrans*latencyGlobalMemory  + globalLoadTrans*latencyGlobalMemory)/(flopsTheoreticalpeaK_GT630*10^6)*10^3;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Co <- timeComputationKernel + processorLatencyGlobal;
SpeedupMatMulGmCoSP_GT630 <- timeKernel_GM_Co[1:7]/matMulGmCoSP;

png(filename="./images/MatMul-Gm-Co-SP.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:7], matMulGmCoSP[1:7],type="l",log="xy", lty = 1, lwd=c(5,5), yaxt="n",xaxt="n", 
     ylim=c(min(matMulGmCoSP,timeKernel_GM_Co), c(max(matMulGmCoSP,timeKernel_GM_Co))), xlim=c(128, 8192), 
     col = "blue", ylab = "Time (milisecond)", xlab = "Size of the Matrices",  cex.lab=2,cex.main=2,
     main = paste("Matrix Multiplication - Global Memory Coalesced \n", GPU, sep =""), axes = FALSE);

lines(N, timeKernel_GM_Co, col = "red",lwd=c(5,5));
points(N, timeKernel_GM_Co, col = "red", type = "p",pch=20,cex = 2)

lines(N[1:7], matMulGmCoSP[1:7], col = "blue",lwd=c(5,5));
points(N[1:7], matMulGmCoSP, col = "blue", pch=21,cex = 2);

axis(1, at = c(N), labels = c(N),cex.axis = 1.5,lwd=c(5,5))

y <- floor(log10(range(matMulGmCoSP,timeKernel_GM_Co)))
pow <- seq(y[1], y[2]+1);
axis(2, 10^pow, cex.axis = 1.5,lwd=c(5,5));

legend('bottomright', c("Predicted time","Measured time"), 
       lty=1, col=c("red","blue"), pch=c(20,21), lwd=c(5,5), cex=2)
grid()
dev.off()


############ MatMul Shared Memory Un-Coalesced #################

globalLoadTrans <- numberthreads*N/16;
globalStoreTrans <- numberthreads/8;


sharedLoadTrans <- numberthreads*N;
sharedStoreTrans <- numberthreads;

processorLatencyGlobal <- ((globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory)/(flopsTheoreticalpeaK_GT630*10^6)*10^3;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Un <- timeComputationKernel + processorLatencyGlobal;
SpeedupMatMulSmUnSP_GT630 <- timeKernel_SM_Un[1:7]/matMulSmUnSP;

png(filename="./images/MatMul-Sm-Un-SP.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:7], matMulSmUnSP[1:7],type="l",log="xy", lty = 1, lwd=c(5,5), yaxt="n",xaxt="n", 
     ylim=c(min(matMulSmUnSP,timeKernel_GM_Un), c(max(matMulSmUnSP,timeKernel_GM_Un))), xlim=c(128, 8192), 
     col = "blue", ylab = "Time (milisecond)", xlab = "Size of the Matrices",  cex.lab=2,cex.main=2,
     main = paste("Matrix Multiplication - Shared Memory Uncoalesced \n", GPU, sep =""), axes = FALSE);

lines(N, timeKernel_SM_Un, col = "red",lwd=c(5,5));
points(N, timeKernel_SM_Un, col = "red", type = "p",pch=20,cex = 2)

lines(N[1:7], matMulSmUnSP[1:7], col = "blue",lwd=c(5,5));
points(N[1:7], matMulSmUnSP, col = "blue", pch=21,cex = 2);

axis(1, at = c(N), labels = c(N), cex.axis = 1.5,lwd=c(5,5))

y <- floor(log10(range(matMulSmUnSP,timeKernel_SM_Un)))
pow <- seq(y[1], y[2]+1);
axis(2, 10^pow, cex.axis = 1.5,lwd=c(5,5));

legend('bottomright', c("Predicted time","Measured time"), 
       lty=1, col=c("red","blue"), pch=c(20,21), lwd=c(5,5), cex=2)
grid()
dev.off()


############ MatMul Shared Memory Coalesced #################
globalLoadTrans <- numberthreads*N/32;
globalStoreTrans <- numberthreads/16;


sharedLoadTrans <- numberthreads*N/16;
sharedStoreTrans <- numberthreads/8;

processorLatencyGlobal <- ((globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory)/(flopsTheoreticalpeaK_GT630*10^6)*10^3;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Co <- timeComputationKernel + processorLatencyGlobal;
SpeedupMatMulSmCoSP_GT630 <- timeKernel_SM_Co[1:7]/matMulSmCoSP;


png(filename="./images/MatMul-Sm-Co-SP.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:7], matMulSmCoSP[1:7],type="l",log="xy", lty = 1, lwd=c(5,5), yaxt="n",xaxt="n", 
     ylim=c(min(matMulSmCoSP,timeKernel_SM_Co), c(max(matMulSmCoSP,timeKernel_SM_Co))), xlim=c(128, 8192), 
     col = "blue", ylab = "Time (milisecond)", xlab = "Size of the Matrices ",  cex.lab=2,cex.main=2,
     main = paste("Matrix Multiplication - Shared Memory Coalesced \n", GPU, sep =""), axes = FALSE);

lines(N, timeKernel_SM_Co, col = "red",lwd=c(5,5));
points(N, timeKernel_SM_Co, col = "red", type = "p",pch=20,cex = 2)

lines(N[1:7], matMulSmCoSP[1:7], col = "blue",lwd=c(5,5));
points(N[1:7], matMulSmCoSP, col = "blue", pch=21,cex = 2);

axis(1, at = c(N), labels = c(N), cex.axis = 1.5,lwd=c(5,5))

y <- floor(log10(range(matMulSmCoSP,timeKernel_SM_Un)))
pow <- seq(y[1], y[2]+1);
axis(2, 10^pow, cex.axis = 1.5,lwd=c(5,5));

legend('bottomright', c("Predicted time","Measured time"), 
       lty=1, col=c("red","blue"), pch=c(20,21), lwd=c(5,5), cex=2)
grid()
dev.off()


############ Difference between predicted and measured Time Matrix Multiplication ##############

png(filename="./images/speedup-MatMul.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:7], SpeedupMatMulGmUnSP_GT630, type="l",  log="x",axes = TRUE, xaxt="n", lty = 1, lwd=c(5,5), 
     ylim=c(0, 2), xlim=c(128, 8192),
     col=cbbPalette[1], ylab = "Time Predicted Vs. Measured", cex.axis=1.5, cex.lab=2, cex.main=2,
     xlab = "Size of the Matrices",  main = paste("Speedup Matrix Multiplication \n", GPU, sep=""));

points(N[1:7], SpeedupMatMulGmUnSP_GT630, col = cbbPalette[1], type = "p",pch=20,cex = 2)

lines(N[1:7], SpeedupMatMulGmCoSP_GT630, col = cbbPalette[2],lwd=c(5,5));
points(N[1:7], SpeedupMatMulGmCoSP_GT630, col = cbbPalette[2], type = "p",pch=21,cex = 2)

lines(N[1:7], SpeedupMatMulSmUnSP_GT630, col = cbbPalette[3],lwd=c(5,5));
points(N[1:7], SpeedupMatMulSmUnSP_GT630, col = cbbPalette[3], type = "p",pch=22,cex = 2)

lines(N[1:7], SpeedupMatMulSmCoSP_GT630, col = cbbPalette[4],lwd=c(5,5));
points(N[1:7], SpeedupMatMulSmCoSP_GT630, col = cbbPalette[4], type = "p",pch=23,cex = 2)

axis(1, at = c(N[1:7]), labels = c(N[1:7]), cex.axis = 1.5)
#axis(2, at=c(SpeedupMatMulGmUnSP_GT630), cex.axis = 1.5);

grid()
legend('topright', 
       lty=1, col=c(cbbPalette[1:4]), pch=c(20:23), lwd=c(5,5), cex=2, 
       legend=c("Global Uncoalesced","Global Coalesced","Shared Uncoalesced","Shared Coalesced"), title="Speedup Model Vs. Measured")

dev.off()

############ Maximum Subsequence GLobal Memory #################
### BSP MULTILEVEL MEMORY

Nn <- 17:29;
N <- 2^Nn;

gridsize <- 32;
blocksize <- 128

numberthreads <- gridsize * blocksize;

N_perBlock <- N/gridsize;
N_perThread <- N_perBlock/blocksize;

#Cycles operations per Threadoperation
tempOperationcycles <- 50 * numberthreads*N_perThread;
# Miliseconds
timeComputationKernel <-  (tempOperationcycles/ (flopsTheoreticalpeaK_GT630*10^6))*10^3; 
# Miliseconds
#communicationLatencyKernel <- (syzeBytes + 16)*numberthreads/(bandWidth_GT630*10^6)*10^3; 
# Miliseconds


globalLoadTrans <- numberthreads*N_perThread/32;
globalStoreTrans <- numberthreads*2*N_perThread;


sharedLoadTrans <- numberthreads*N_perThread*2;
sharedStoreTrans <- numberthreads*N_perThread*2;

processorLatencyGlobal <- ((globalStoreTrans + globalLoadTrans)*latencyGlobalMemory + (sharedLoadTrans + sharedStoreTrans)*latencySharedMemory)/(flopsTheoreticalpeaK_GT630*10^6)*10^3;

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_timeKernel_SubSeqMax <- timeComputationKernel + processorLatencyGlobal;

SpeedupSubSeqMax <- timeKernel_timeKernel_SubSeqMax[1:12]/SubSeqMax;

png(filename="./images/SubSeqMax.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:12], SubSeqMax[1:12],type="l",log="xy", lty = 1, lwd=c(5,5), yaxt="n",xaxt="n", 
     ylim=c(min(SubSeqMax,timeKernel_timeKernel_SubSeqMax), c(max(SubSeqMax,timeKernel_timeKernel_SubSeqMax))), xlim=c(131072, 536870912), 
     col = "blue", ylab = "Time (milisecond)", xlab = "Size of the Sequence ",  cex.lab=2,cex.main=2,
     main = paste("SubSeqMax - Shared Memory Coalesced \n", GPU, sep=""), axes = FALSE);

lines(N, timeKernel_timeKernel_SubSeqMax, col = "red",lwd=c(5,5));
points(N, timeKernel_timeKernel_SubSeqMax, col = "red", type = "p",pch=20,cex = 2)

lines(N[1:12], SubSeqMax[1:12], col = "blue",lwd=c(5,5));
points(N[1:12], SubSeqMax, col = "blue", pch=21,cex = 2);

axis(1, at = c(N), labels = c(N), cex.axis = 1.5,lwd=c(5,5))

y <- floor(log10(range(SubSeqMax,timeKernel_timeKernel_SubSeqMax)))
pow <- seq(y[1], y[2]+1);
axis(2, 10^pow, cex.axis = 1.5,lwd=c(5,5));

legend('bottomright', c("Predicted time","Measured time"), 
       lty=1, col=c("red","blue"), pch=c(20,21), lwd=c(5,5), cex=2)
grid()
dev.off()



############ Difference between predicted and measured SubSequenceMaximum ##############

png(filename="./images/speedup-SubSeMax.png", width=1200, height=800)
par(mar=c(6, 6, 4, 8) + 0.1)
plot(N[1:12], SpeedupSubSeqMax, type="l",  log="x",axes = TRUE, lty = 1, lwd=c(5,5), xaxt="n",
     ylim=c(0, 3), xlim=c(131072, 536870912),
     col=cbbPalette[1], ylab = "Time Predicted Vs. Measured ( % )", cex.axis=1.5, cex.lab=2,cex.main=2,
     xlab = "Size of the Matrices ",  main = paste("Speedup SubSeqMax \n ", GPU, sep=""));
points(N[1:12], SpeedupSubSeqMax, col = cbbPalette[1], type = "p",pch=20,cex = 2)

axis(1, at = c(N[1:12]), labels = c(N[1:12]), cex.axis=1.5)
#axis(2, at=c(SpeedupMatSumGmUnSP_GT630));

legend('topright', 
       lty=1, col=c(cbbPalette[1]), pch=c(20), lwd=c(5,5), cex=2, 
       legend=c("SubSeqMax"), title="Speedup")
grid()
dev.off()


