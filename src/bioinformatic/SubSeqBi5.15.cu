#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <time.h> 

/*

Versão do Daniel com memória compartilhada e com 2 gpus
Melhor versão

*/

#define N_COL 2048 
#define N_LIN 2048

#define BLOCK_SIZE_SOMA_PREF 2
#define N_THREAD_PER_BLOCK_SOMA_PREF 1024

#define BLOCK_SIZE_CALC 16 
#define N_THREAD_PER_BLOCK_CALC 1024

#define N_GPU 1

__global__ void somaPrex(int* matriz, int nColuna, int nLin){ 

   int idThread = blockIdx.x * blockDim.x + threadIdx.x;

   int qtdColSoma = (nLin / (blockDim.x * gridDim.x));

   int comecDeTrabThr = ((qtdColSoma * nColuna) * idThread); 

   int i,j;
   for(j=0; j<qtdColSoma;j++){
	   for(i=1; i<nColuna; i++){
		   matriz[(comecDeTrabThr + i) + (j*nColuna)] += matriz[(j*nColuna) + (comecDeTrabThr + i) - 1];
	   }
   }
}

__global__ void calc(int* matriz, int nCol, int nLin, int qtdDeParesGH,  int *Subseqs, int *auxMatriz, int GPU){

	long int idThread = (blockIdx.x * blockDim.x + threadIdx.x);
	int tidAux        = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ int vetComp[49152];

	int i,j,k,l,auxComp,comecoVetThreadComp;
	
	//      M         t_m       S      suf
	int ini_M, fim_M, t_M, ini_S, fim_S, suf;
	ini_M = fim_M = ini_S = fim_S = -1;
 
	comecoVetThreadComp = (threadIdx.x * 12);
		
	auxComp = t_M = suf = 0;	

	for(l=0; l< ((qtdDeParesGH/(blockDim.x * gridDim.x))+1); l++){
	
		 float delta = ((4*(nCol*nCol)) - (4*nCol) + 1) - (8 * (idThread+1));
		 float auxG = ((((2*nCol) - 3) - sqrt(delta))/2);

		 if(auxG - (int)auxG){
		    auxG += 1;
         }

		 long int g = auxG;

		 float auxH= (idThread+1) - (((((2*nCol)-1)-g) * g)/2)+g;
    	 long int h = auxH;

		 if(g>=0 && g<=nCol && h>=0 && h<=nCol && idThread < qtdDeParesGH){
			 if(g==0){
			   for(j=0; j < nLin; j+=12){
				   for(k=0;k<12;k++){
					   vetComp[(threadIdx.x * 12)+k] = matriz[(h + (nCol * (k+j)))+(nCol*GPU)]; 
					}

					for(i = (comecoVetThreadComp -1); i < (comecoVetThreadComp + 12)-1; i++){
						if(i == fim_M){
							fim_S++;
							suf += vetComp[i+1];

							if(suf < 0){
								suf = 0;
								fim_S = -1;
							}

							ini_S = fim_S == 0 ? 0 : ini_S; // Inicio S

							  if(vetComp[i+1] > 0){
								   fim_M++;
								   t_M += vetComp[i+1];
								   ini_M = fim_M == 0 ? 0 : ini_M; // Inicio M
							   }
						}
						else{
							 if(suf + vetComp[i+1] > t_M){
								 fim_S++;
								 if(ini_M == -1){
								  fim_S = ini_S = i +1;
								 }

								 suf += vetComp[i+1];
								 ini_M = ini_S;
								 fim_M = fim_S;
								 t_M = suf;
							 }
							 else{
								    if(suf + vetComp[i+1] > 0){
								        fim_S++;
								        if(suf == 0){
								            ini_S = fim_S = i+1;
								        }

								        suf += vetComp[i+1];
								    }
								    else{
								        ini_S = fim_S = i + 2;
								        suf = 0;
								    }
							 }
						}
					}
			   }

				if(t_M > auxComp){																									
		  		   Subseqs[tidAux] = t_M;
				   auxComp = t_M;
			    }

				idThread += (blockDim.x * gridDim.x);
			 }
			 else{
				 for(j=0; j < nLin; j+=12){
				   for(k=0;k<12;k++){
					   vetComp[(threadIdx.x * 12)+k]  = matriz[(h + (nLin * (k+j)))+(nCol*GPU)] - matriz[((g-1) + (nLin * (k+j)))+(nCol*GPU)];
					}

					for(i = (comecoVetThreadComp -1); i < (comecoVetThreadComp + 12)-1; i++){
						if(i == fim_M){
							fim_S++;
							suf += vetComp[i+1];

							if(suf < 0){
								suf = 0;
								fim_S = -1;
							}

							ini_S = fim_S == 0 ? 0 : ini_S; // Inicio S

							  if(vetComp[i+1] > 0){
								   fim_M++;
								   t_M += vetComp[i+1];
								   ini_M = fim_M == 0 ? 0 : ini_M; // Inicio M
							   }
						}
						else{
							 if(suf + vetComp[i+1] > t_M){
								 fim_S++;
								 if(ini_M == -1){
								  fim_S = ini_S = i +1;
								 }

								 suf += vetComp[i+1];
								 ini_M = ini_S;
								 fim_M = fim_S;
								 t_M = suf;
							 }
							 else{
								    if(suf + vetComp[i+1] > 0){
								        fim_S++;
								        if(suf == 0){
								            ini_S = fim_S = i+1;
								        }

								        suf += vetComp[i+1];
								    }
								    else{
								        ini_S = fim_S = i + 2;
								        suf = 0;
								    }
							 }
						}
					}
			     }

				  if(t_M > auxComp){																									
		  		     Subseqs[tidAux] = t_M;
				     auxComp = t_M;
			      }
				  
		 		  idThread += (blockDim.x * gridDim.x);	
			 }
		 }
	}
}

int main(){

   float elapsedTime;    // Tempo
   cudaEvent_t start, stop; // Tempo

   int i; int qtdDeParesGH = (((N_COL/2)*((N_COL-1)/2)) / 2);
	
   printf("\n Pares G e H Calculados %d \n",qtdDeParesGH);

   //Alocando a matriz no host
   int *matriz_h = (int *)malloc(sizeof(int *) * (N_COL*N_LIN));
   int *subSeq_h = (int *)malloc(sizeof(int *) * ((BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC)*N_GPU));
   int *matriz_d; int *subSeq_d; int *auxMatriz_d;

   //Preenchendo a matriz no host
   for(i=0; i<(N_COL*N_LIN); i++){
       matriz_h[i] = -1;
   }

   for(i=0; i<(BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC); i++){
       subSeq_h[i] = 0;
   }

   matriz_h[228] = 580;   matriz_h[229] = 280;
		    

   for(i=0; i<N_GPU; i++){
		cudaSetDevice(i);

   		//Reservando espaco na GPU
   		cudaMalloc((void**)&matriz_d, (N_COL*N_LIN)  * sizeof(int)); 
   		cudaMalloc((void**)&auxMatriz_d, ((BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC) * N_LIN)  * sizeof(int)); 
   		cudaMalloc((void**)&subSeq_d, (BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC) * sizeof(int));

   		cudaMemcpy(matriz_d, matriz_h, (N_COL*N_LIN) * sizeof(int), cudaMemcpyHostToDevice);

	}

   cudaEventCreate(&start); // Contagem do tempo
   cudaEventCreate(&stop);
   cudaEventRecord(start, 0);

   for(i=0; i<N_GPU; i++){
		cudaSetDevice(i);
	    somaPrex<<<BLOCK_SIZE_SOMA_PREF, N_THREAD_PER_BLOCK_SOMA_PREF>>>(matriz_d,N_COL,N_LIN); 
   }

   cudaThreadSynchronize();

   for(i=0; i<N_GPU; i++){
      calc<<<BLOCK_SIZE_CALC, N_THREAD_PER_BLOCK_CALC>>>(matriz_d,N_COL/2,N_LIN, qtdDeParesGH,subSeq_d,auxMatriz_d,i);
   }

   for(i=0; i<N_GPU; i++){
   	   cudaMemcpy(subSeq_h + ((BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC)*i), subSeq_d, (BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC) * sizeof(int), cudaMemcpyDeviceToHost);
   }

   cudaEventRecord(stop, 0);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&elapsedTime, start, stop);

   //printf("\n Tempo do kernel (ms) = \%f\n\n", elapsedTime);


   //Encontrando a maior subSeq
   int maiorSubSeq = subSeq_h[0];
   for(i=0; i < ((BLOCK_SIZE_CALC*N_THREAD_PER_BLOCK_CALC)*N_GPU); i++){
	   maiorSubSeq = subSeq_h[i] > maiorSubSeq ? subSeq_h[i]:maiorSubSeq;
	   //printf("%d ", subSeq_h[i]);
   }

  printf("\n Maior SubSequencia encontrada \n\n %d \n\n",maiorSubSeq);

return 0;
}
