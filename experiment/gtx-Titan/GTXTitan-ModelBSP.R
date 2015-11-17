setwd("/home/marcos/Dropbox/Doctorate/Results/2015/BSPGPU/experiment/gtx-Titan/")

GPU <- "Geforce GTX-Titan"

# The palette of color
cbbPalette <- gray(1:4/ 6)#c("red", "blue", "darkgray", "orange","black","lightblue", "lightblue","violet")


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


PCoresNumber_GT630 <- 2688;
clockFrequency_GT630 <- 876; # Mhz  
cycleMultiplicationKepler <- 16;
cycleSumKepler <- 4;
cycleFMA <- 2

memoryBusWidth_GT630 <- 384/8; # Bytes
memoryClockRate_GT630 <- 3004; # Mhz

latencySharedMemory <- 5; #Cycles per processor
latencyGlobalMemory <- latencySharedMemory* 100; #Cycles per processor

latencyL1 <- latencySharedMemory; #Cycles per processor
latencyL2 <- latencyGlobalMemory*0.5; #Cycles per processor

bandWidth_GT630 <-  memoryClockRate_GT630 * memoryBusWidth_GT630 *2; # MB/s
flopsTheoreticalpeaK_GT630 <- clockFrequency_GT630 * PCoresNumber_GT630; # Mflops/second


GDDR5SizeMB_GT630 <- 2048; # MB
L2CacheSize_GT630 <- 1536; # Kb
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
W <- 4.5

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Un <- ( W^-1*(timeComputationKernel + CommGM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

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
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 22

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_GM_Co <- ( W^-1*(timeComputationKernel + CommGM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

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
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 18

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);
CommSM <- (numberthreads*N*2 + numberthreads)*latencySharedMemory

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Un <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

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
L1Effect <- 0
L2Effect <- 0.02*reads
W <- 52

CommGM <- ((numberthreads*N*2 - L1Effect - L2Effect + numberthreads)*latencyGlobalMemory + L1Effect*latencyL1 + L2Effect*latencyL2);
CommSM <- ((numberthreads*N*2 + numberthreads)*latencySharedMemory)

# MODEL Multi_BSP_GPUModelPredictionTime Miliseconds
timeKernel_SM_Co <- ( W^-1*(timeComputationKernel + CommGM + CommSM)/(flopsTheoreticalpeaK_GT630*10^6))*10^3;

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
