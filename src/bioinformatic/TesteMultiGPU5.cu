#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 50000//450000000//157370000//524288//157370000//1310720//262144//131072//262144//83886080

//Quantidade de threads por blocos
#define BLOCK_SIZE 1//1024//32//95536
#define nThreadsPerBlock 128//128//420 ou 416

#define nGPU 1
#define nVetorFinal 
#define nVetFinalGPUs ((BLOCK_SIZE * 128) * 5) 

/*
	6.14 versao com coalesced e thread's auxiliando na copia de memoria

	Instancia de 157370000

	com 420 threads

*/

__device__ void memoria(int *vetDados,int *vetComp, int ElemPorBlocos, int qtdProces){

	int aux = (qtdProces * 4096);

	int comecoBloco = blockIdx.x * ElemPorBlocos; // onde cada bloco irá comeca

	int idCompartilhada = threadIdx.x;
	int idGlobal = threadIdx.x + aux + comecoBloco;

	int i;
	for(i = 0; i < 4096; i += blockDim.x){
		vetComp[idCompartilhada] = vetDados[idGlobal];
		idCompartilhada += blockDim.x;
		idGlobal += blockDim.x;
	}
}

__global__ void subSeqMax(int *vet, int *vetFinal, int ElemPorThread, int n){	
		
	__shared__ int p[4096];
	
	//      M         t_m       S      suf
	int ini_M, fim_M, t_M, ini_S, fim_S, suf; //Variaveis do algoritmo
	t_M = suf = 0;

	int comecoThread = (threadIdx.x * 32);

	int j;
	for(j = 0; j < (n / 4096); j++){ // Quantas vezes terei que processa até chegar no n/blocos sendo que o vet compartilhado é de 4076

		memoria(vet,p,n,j);

		__syncthreads();

			if(threadIdx.x < 128){

				ini_M = fim_M = ini_S = fim_S = comecoThread -1;

				int i;
				for(i = comecoThread -1; i < comecoThread + 32; i++){
					if(i == fim_M){
				    	fim_S++;
				    	suf += p[i+1];

				    	if(suf < 0){
				        	suf = 0;
				        	fim_S = -1;
				   		 }
				    
						ini_S = fim_S == 0 ? 0 : ini_S; // Inicio S

				     	 if(p[i+1] > 0){
				           fim_M++;
				           t_M += p[i+1];
				           ini_M = fim_M == 0 ? 0 : ini_M; // Inicio M
				       	 }
					}
					else{
						 if(suf + p[i+1] > t_M){
						     fim_S++;
						     if(ini_M == -1){
						      fim_S = ini_S = i +1;
						     }

						     suf += p[i+1];
						     ini_M = ini_S;
						     fim_M = fim_S;
						     t_M = suf;
						 }
						 else{
						        if(suf + p[i+1] > 0){
						            fim_S++;
						            if(suf == 0){
						                ini_S = fim_S = i+1;
						            }

						            suf += p[i+1];
						        }
						        else{
						            ini_S = fim_S = i + 2;
						            suf = 0;
						        }
				     	}//else
					}//else
				}// 1* for
		}// If 128	
	}// 2* for

	if(threadIdx.x < 128){
		int idThread = blockIdx.x * blockDim.x + threadIdx.x;

		vetFinal[(idThread * 5)] =  vetFinal[(idThread * 5)+1] = vetFinal[(idThread * 5)+2] = vetFinal[(idThread * 5)+3] =
		vetFinal[(idThread * 5)+4] = -1;

		//Colocando o M
		vetFinal[(idThread * 5)+2] = t_M;

		//Calculando o Prefixo
		int pref_Max, soma_Pref;
		soma_Pref  = 0;
		pref_Max = 0;

		int i;
		if(ini_M > comecoThread -1){
		    for(i = 0; i < ini_M; i++){
		        soma_Pref += p[i];

		        if(soma_Pref > pref_Max){
		            pref_Max = soma_Pref;
		        }
		    }

		    if(pref_Max == 0){
				vetFinal[(idThread * 5)] = 0;
				vetFinal[(idThread * 5)+1] = soma_Pref;
		    }
		    else{
				vetFinal[(idThread * 5)] = pref_Max; //Prefixo
				vetFinal[(idThread * 5)+1] = soma_Pref - pref_Max; //Numeros negativos
		    }
		}

		//Calculo do sufixo
		int suf_Max, soma_Suf;
		soma_Suf = suf_Max = 0;

		if(fim_M < comecoThread + 32){
		    for(i = (comecoThread + 32)-1; i > fim_M; i--){
		        soma_Suf += p[i];

		        if(soma_Suf > suf_Max){
		            suf_Max = soma_Suf;
		        }
		    }

		    if(suf_Max == 0){
		        vetFinal[(idThread * 5)+3] = 0;	//Sufixo vazio
				vetFinal[(idThread * 5)+4] = suf_Max;//Os Numeros negativos

		    }
		    else{
		        vetFinal[(idThread * 5)+3] = suf_Max;	//Sufixo vazio
				vetFinal[(idThread * 5)+4] = soma_Suf - suf_Max;//Os Numeros negativos
		    }
		}
	}//if 128
}

void subSeqMaxFinal(int *vet, int n){

    //      M         t_m       S      suf
    int ini_M, fim_M, t_M, ini_S, fim_S, suf;
    ini_M = fim_M = ini_S = fim_S = -1;

    t_M = suf = 0;

	int start;
	int tmili;
	start = clock();


	int i;
    for(i = -1; i < n-1; i++){
        if(i == fim_M){
            fim_S++;
            suf += vet[i+1];

            if(suf < 0){
                suf = 0;
                fim_S = -1;
            }

            ini_S = fim_S == 0 ? 0 : ini_S; // Inicio S

              if(vet[i+1] > 0){
                   fim_M++;
                   t_M += vet[i+1];
                   ini_M = fim_M == 0 ? 0 : ini_M; // Inicio M
               }
        }
        else{
             if(suf + vet[i+1] > t_M){
                 fim_S++;
                 if(ini_M == -1){
                  fim_S = ini_S = i +1;
                 }

                 suf += vet[i+1];
                 ini_M = ini_S;
                 fim_M = fim_S;
                 t_M = suf;

             }
             else{
                    if(suf + vet[i+1] > 0){
                        fim_S++;
                        if(suf == 0){
                            ini_S = fim_S = i+1;
                        }

                        suf += vet[i+1];

                    }
                    else{
                        ini_S = fim_S = i + 2;
                        suf = 0;
                    }
             }
        }
    }

	tmili = (int)((clock()-start)*1000/CLOCKS_PER_SEC);

	printf(" \n\n A sub Sequencia deu %d  \n\n", t_M);

	printf("Tempo total do sequencial %d  \n\n",tmili);

}

int main(){

	float elapsedTime;    // Tempo
	cudaEvent_t start, stop; // Tempo

	//Vetores que serao usado no device 
	int *vet_d; // Vetor de dados device
	int *vetFinalGPUs_d; // Vetor final que as GPUs iram devolver

	//Vetores que serao usado no host
	int *vet_h = (int *) malloc(sizeof(int) * N); // Vetor Dados
	int *vetFinal_h = (int *) malloc (sizeof(int) * (nVetFinalGPUs * nGPU));// Vetor Final, as GPU iram devolver o seu vetor final dentro dle

	int i;
	for(i = 0; i < N; i++){ // Preenchimento dos dados
	     vet_h[i] = -1;
	}

	for(i = 0; i < (nVetFinalGPUs * nGPU); i++){ // Preenchimento dos dados
	     vetFinal_h[i] = -1;
	}

	vet_h[70] = 300;

//-----------------------------------------------Reservando espaço de memória e copiando os dados para o device
																													
	for(i = 0; i < nGPU; i++){
		cudaSetDevice(i);        //(N / nGPU) porque cada GPU vai cuida de intervalo da sequencia original

		cudaMalloc((void**)&vet_d, (N / nGPU) * sizeof(int)); //Vetor de dados em cada GPU
		cudaMalloc((void**)&vetFinalGPUs_d, nVetFinalGPUs * sizeof(int)); // Vetor final que cada GPU irá receber

    	cudaMemcpy(vet_d, (vet_h) + ((N / nGPU) * i),  (N / nGPU) * sizeof(int), cudaMemcpyHostToDevice);
						//(vet_h) + ((N / nGPU) * i) falo a onde esta cada intervalo das GPU
	}

//------------------------------------------------------------------------------------------------------Kernels

	int ElemPorGPU = N / nGPU; // Cada GPU fica responsavel por um intervalo de N/nGPU do vetor original
	int ElemPorBlocos = (ElemPorGPU / BLOCK_SIZE); //Cada bloco fica responsavel por um intervalo da sequencia origina
	int ElemPorThread = (ElemPorBlocos / nThreadsPerBlock); //Cada thread fica responsavel por um intervalo da sequencia original
	
	cudaEventCreate(&start); // Contagem do tempo
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for(i = 0; i < nGPU; i++){
		cudaSetDevice(i);
		subSeqMax<<<BLOCK_SIZE, nThreadsPerBlock>>>(vet_d, vetFinalGPUs_d, ElemPorThread, ElemPorBlocos);
	}


	cudaThreadSynchronize();

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);



//	printf("Primeiro kernel (ms) = \%f\n\n", elapsedTime);

	
//-----------------------------------------------------------------------------------------------------------------
	
	for(i = 0; i < nGPU; i++){
		cudaSetDevice(i);
		cudaMemcpy(vetFinal_h + (i * nVetFinalGPUs) , vetFinalGPUs_d, nVetFinalGPUs * sizeof(int), cudaMemcpyDeviceToHost); //Resposta Final
	}

//---------------------------------------------------------------------------------------Formando o vetor Final



	for(i = 0; i < 4096; i++){
		if(vetFinal_h[i] != 0 && vetFinal_h[i] != -1 )
			printf("%d ", vetFinal_h[i]);	
	}
		
	printf("\n\n");

	cudaFree(vetFinalGPUs_d);
	cudaFree(vet_d);

    	subSeqMaxFinal(vetFinal_h, (nVetFinalGPUs * nGPU));

	return 0;
}

