setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gtx-660")

GPU <- "Geforce GTX-660"

# The palette of color
cbbPalette <- gray(1:4/6)#c("red", "blue", "darkgray", "orange","black","lightblue", "lightblue","violet")


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


PCoresNumber_GT630 <- 960;
clockFrequency_GT630 <- 1058; # Mhz	
cycleMultiplicationKepler <- 16;
cycleSumKepler <- 4;
cycleFMA <- 2

memoryBusWidth_GT630 <- 192/8; # Bytes
memoryClockRate_GT630 <- 3004; # Mhz

latencySharedMemory <- 5; #Cycles per processor
latencyGlobalMemory <- latencySharedMemory* 100; #Cycles per processor

latencyL1 <- latencySharedMemory; #Cycles per processor
latencyL2 <- latencyGlobalMemory*0.5; #Cycles per processor

bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeaK_GT630 <- clockFrequency_GT630 * PCoresNumber_GT630 ; # Mflops/second


GDDR5SizeMB_GT630 <- 2048; # MB
L2CacheSize_GT630 <- 384; # Kb
sharedMemorySize_GT630 <- 48;  # Kb
L1CacheSizeKB_GT630 <- 16; # Kb
registersPerBlock_GT630 <- 64; # Kb


### BSP MULTILEVEL MEMORY

#Cycles operations per Threadoperation
tempOperationcycles <- ((numberMultiplication * cycleFMA) + (0 * cycleSumKepler)) * numberthreads;
# Miliseconds
timeComputationKernel <-  tempOperationcycles; 
# Miliseconds
#communicationLatencyKernel <- (syzeBytes + 16)*numberthreads/(bandWidth_GT630*10^6)*10^3; 
# Miliseconds

reads <- numberthreads*N*2

L1Effect <- 0
L2Effect <- 0.02*reads
W <- 4.8

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Un <- ( W^-1*(timeComputationKernel + CommGM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

SpeedupMatMulGmUnSP_GT630 <- timeKernel_GM_Un[1:7]/matMulGmUnSP;


############ MatMul GLobal Memory Coalesced #################
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 21

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Co <- ( W^-1*(timeComputationKernel + CommGM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

SpeedupMatMulGmCoSP_GT630 <- timeKernel_GM_Co[1:7]/matMulGmCoSP;


############ MatMul Shared Memory Un-Coalesced #################
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 20

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);
CommSM <- (numberthreads*N*2 + numberthreads)*latencySharedMemory

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Un <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

SpeedupMatMulSmUnSP_GT630 <- timeKernel_SM_Un[1:7]/matMulSmUnSP;


############ MatMul Shared Memory Coalesced #################
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 70

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);
CommSM <- ((numberthreads*N*2 + numberthreads)*latencySharedMemory)

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Co <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

SpeedupMatMulSmCoSP_GT630 <- timeKernel_SM_Co[1:7]/matMulSmCoSP;


############ Difference between predicted and measured Time Matrix Multiplication ##############
dataN <- 2:7
png(filename="./images/speedup-MatMul.png", width=800, height=600)
par(mar=c(4, 6, 2, 1) + 0.1)

plot(N[dataN], SpeedupMatMulGmUnSP_GT630[dataN], type="l",  log="x",axes = TRUE, xaxt="n", lty = 4, lwd=c(7.5,7.5), 
     ylim=c(0.9, 1.1), xlim=c(256, 8192),
     col=cbbPalette[1], ylab = " ", main=" ", cex.axis=3.5, cex.lab=3.5, cex.main=3.5,
     xlab = " ");

points(N[dataN], SpeedupMatMulGmUnSP_GT630[dataN], col = cbbPalette[1], pch=20,cex = 3.5)

lines(N[dataN], SpeedupMatMulGmCoSP_GT630[dataN], col = cbbPalette[2],lwd=c(7.5,7.5), lty = 3);
points(N[dataN], SpeedupMatMulGmCoSP_GT630[dataN], col = cbbPalette[2], pch=21,cex = 3.5)

lines(N[dataN], SpeedupMatMulSmUnSP_GT630[dataN], col = cbbPalette[3],lwd=c(7.5,7.5), lty = 2);
points(N[dataN], SpeedupMatMulSmUnSP_GT630[dataN], col = cbbPalette[3], pch=22,cex = 3.5)

lines(N[dataN], SpeedupMatMulSmCoSP_GT630[dataN], col = cbbPalette[4],lwd=c(7.5,7.5), lty = 1);
points(N[dataN], SpeedupMatMulSmCoSP_GT630[dataN], col = cbbPalette[4], pch=23,cex = 3.5)

axis(1, at = c(N[dataN]), labels = c(N[dataN]), cex.axis = 1.5)
#axis(2, at=c(SpeedupMatMulGmUnSP_GT630), cex.axis = 1.5);
grid()
# legend('bottom', 
#        lty = c(4:1), col=c(cbbPalette[1:4]), pch=c(20:23), lwd=c(7.5,7.5), cex=2.75, ncol=2,bty ="n",
#        legend=c("Optimization (#1)","Optimization (#2)","Optimization (#3)","Optimization (#4)"))
dev.off()
