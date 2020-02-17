/**#################################################################################
#                                                                                  #
#        MODELAGEM ESTATÍSTICA DA LIBERAÇÃO DE FÁRMACOS NANO-ENCAPSULADOS  	   #
#                                Simulação - 3D                                    #
#                                                                                  #
#                                                                                  #
#    Orientador: Marco A. Alves                                                    #
#    Márcio Sampaio Gomes Filho                                                    #
#                                                                                  #
#    Criado: 01/06/2012 - UnB                                                      #
#                                                                                  #
#      Dscrição do Algoritmo:                                                      #
#           i)   Forma-se uma configuração inicial (rede cúbica) L^3               #
#           ii)  Sortea-se uma partícula aleatoriamente                            #
#           iii) Move-se a particula aleatoriamente                                #
#           iv)  Verifica-se se a partícula saiu da rede                           #
#           v)   Dado um passo de MC, grava-se o tempo-MC e o número de particulas #
#                remanescentes                                                     #
#                                                                                  #
#                                                                                  #
####################################################################################**/


/*  - Início: 01/06/2012
 *
 *  01/06/2012
 *    - criação da main
 *    - criação da função parametros
 *    - criação da matriz 3D
 *
 * 02/06/2012
 *   - criação da função salto
 *   - Implementação do passo de monte carlo
 *   - Nesta versão as particulas ficam saltando até sairem da rede
 *
 * 10/06/2012
 *  - correção do tempo (tempo de monte carlo 1/dt)
 *  - implementação de abrir e nomear os arquivos
 *
 * 15/06/2012
 *  - correção do tempo
 *
 * 20/04/2013
 *  - Implementção do passo de impressao
 *  - percentual da liberação a ser estuda -f
 *
 *  27/07/2014
 *   - Revisão codigo
 *   - Implementação da função de visualização
 *   - implementar funçao para particulas no tempo
 *
 * 06/05/2016
 *   - Revisão código
 *   - Implementação de uma função para a visualição no GNUPLOT
 *     - Função para a porosidade Ok
 *
 * 31/01/2017 META
 *  - função porosidade - ok
 *  - revisão do código
 *  - implementação do cálculo do tempo
 *  - funçoes debug e melhor função DEBUG
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#define COD_SAIDA   -2
#define COD_BORDA   -1

/* Cria arquivos para vizualização estrutural - Gnuplot -> compilação gcc    -D MOVIE */
//#define MOVIE

typedef struct
{
    int     tam_rede,               // tamanho da rede L x L x L
            numero_poros,           // numero de poros
            numero_simulacoes,      // numero de simulacoes
            numero_particulas,      // numero de particulas
            tempo_impressao;        // passos para a impressão
    float   fracao_liberacao;       // fração da liberação a ser estudada -f 0.6 = 60% da liberação
    double  PMC;                    // passo de monte carlo

} sistema;

sistema sis;

enum TIPO_SALTO
{
    p_direita = 1,
    p_esquerda,
    p_cima,
    p_baixo,
    p_tras,
    p_frente
};
/************ Funções ***************************************************************************************************************************************/
void parametros(int argc, char *argv[]);														/// recebe parâmetros via terminal
enum TIPO_SALTO salto_part(void);																/// sorteio do movimento, ex: p/direita
void grava_configuracao_arq_video (int ***y, FILE *arq1, FILE *arq2, FILE *arq3);				/// estrutura: particlas, membrana e poros
void particulas_video (int ***rede, FILE *arq3);												/// grava as posições das particlas em cada passo
int intervalo_medida(double t1, double t2);														/// verificar o intervalo de tempo
void dados_arquivo_saida(double *particulas_tempo_real, double *particulas_tempo_real2);        /// media, desvio quadratico medio
void escreve_membrana (int **membrana);															/// aloca os pontos da superfície da membrana
void forma_membrana(int ***rede, int **membrana);												/// forma configuração inicial - membrana, poros e partículas
void des_alloc(int ***rede, int  **membrana,double *vetor1, double *vetor2 );					/// libera os vetores
void imprime_rede(int x, int ***y);																/// imprime a rede-D
double second();																				/// função para calcular o tempo de execução
/*************************************************************************************************************************************************************/

int main (int argc, char *argv[])
{

    int i, j, k,                                         // contador
        ***rede,                                         // rede tridimensional
        **membrana,					                     // membrana
        x,y,z,                                           // poisições das partículas
        salto,                                           // salto
        salto_X, salto_Y, salto_Z,                       // saltos X, Y e Z
        cont_PART;                                       // conta o numero de particulas que sai da rede
    double tempo_mc, tempo_mc_antigo, *particulas_tempo_real, *particulas_tempo_real2,  tempo_inicial, tempo_final;

    /* início do tempo de execução */
    tempo_inicial = second();

    /* semente aleatória */;
    srand(time(NULL)*argc);
    parametros(argc, argv);

#ifdef MOVIE
    FILE *arq1=NULL, *arq2=NULL, *arq3=NULL;
    arq1 = fopen("membrana.xyz", "w+");
    arq2 = fopen("poros.xyz", "w+");
    arq3 = fopen("particulas.xyz", "w+");
#endif // MOVIE

///*############## alocando memória ######################

    /* WARNING:
    *   - ponteiro de ponteiro de ponteiro não é uma boa alternativa, pois a leitura na memoria é aleatória, - aqui podemos otimizar
         - além disso, um ponteiro aponta para um ponteiro de ponteiro (matriz) - tomar cuidado aqui
    *  Teste:  utilizaremos no máximo tam_rede=100, no entanto, p/ teste a criação do ponteiro funciona bem p/ tam_rede=1000 */
    rede =  (int***) calloc(sis.tam_rede + 1, sizeof(int**));
    for(i = 0; i < sis.tam_rede + 1; i++){                                             // começando em zero, por segurança!
        rede[i] = (int**) calloc(sis.tam_rede + 1, sizeof(int*));
        for(j = 0; j < sis.tam_rede + 1; j++)
            rede[i][j] = (int*) calloc(sis.tam_rede + 1, sizeof(int));
    }

    /* vetor guarda o nº de particulas em cada instante de tempo: a ser utilizado no cálculo da média */
    particulas_tempo_real  = (double *) calloc(sis.PMC + 1, sizeof(double));

    /* vetor guarda a o nº de particulas ao quadrado: a ser utilizado no cálculo da média do quadrado */
    particulas_tempo_real2 = (double *) calloc(sis.PMC + 1, sizeof(double));

    /* matriz representando a superfície da capsula - membrana */
    membrana = ( int** )   calloc (  6*(sis.tam_rede - 2)*(sis.tam_rede - 2) +  1, sizeof(  int * )); /// linhas na matriz - i -> Aij = Axy
        for(i = 0; i < 6*(sis.tam_rede - 2)*(sis.tam_rede - 2)  + 1; i++) /// colunas na matriz - j
            membrana[i] = ( int* )  calloc (3 + 1, sizeof( int ));       // 3 colunas

   /* Iniciamos escrevendo em um arquivo "membrana_inicial.dat"  os pares ordenados que representam a superfície (aqui fica independente da geometria da membrana).
    * Aqui poderia ser um outro programa que passa os pontos - independente da geometria . */
   escreve_membrana(membrana);


///*############## início da corrida ######################

    for(i = 1; i <= sis.numero_simulacoes; i++){

	 /*zera a rede - aloca a membrana - distribui os poros aleatoriamente - preenche de partículas*/
	 forma_membrana(rede, membrana);

#ifdef MOVIE
        grava_configuracao_arq_video(rede, arq1, arq2, arq3);
#endif

        tempo_mc = 0.0;
        cont_PART = 0;

        for(j = 1; j <= sis.PMC; j++) {

            for(k = 1; k <= sis.numero_particulas - cont_PART; k++){

                do{
                    x = rand() % sis.tam_rede + 1;
                    y = rand() % sis.tam_rede + 1;
                    z = rand() % sis.tam_rede + 1;
                } while(rede[x][y][z] <= 0);


                tempo_mc_antigo = tempo_mc;
                tempo_mc = tempo_mc + ( (double)   1.0   /  (sis.numero_particulas - cont_PART)  );

                if ( intervalo_medida(tempo_mc_antigo, tempo_mc) ) {
                    particulas_tempo_real[ (int) floor(tempo_mc) ] += (  (double)   sis.numero_particulas - cont_PART );
                    particulas_tempo_real2[(int) floor(tempo_mc) ] += (  (double)   sis.numero_particulas - cont_PART )* (  (double)  sis.numero_particulas - cont_PART ) ;
                }

                salto = salto_part();
                salto_X = salto_Y = salto_Z = 0;

                switch(salto){
                    case p_direita:
                        salto_X = 1;
                        break;
                    case p_esquerda:
                        salto_X = -1;
                        break;
                    case p_cima:
                        salto_Y = 1;
                        break;
                    case p_baixo:
                        salto_Y = -1;
						break;
                    case p_frente:
                        salto_Z = 1;
                        break;
					case p_tras:
                        salto_Z = -1;
                        break;
				}

                if(rede[ x + salto_X ][ y + salto_Y ][ z + salto_Z ] == COD_SAIDA){
                    rede[x][y][z] = 0;
                    cont_PART++;
#ifdef MOVIE
                fprintf(arq3,"C %d\t%d\t%d\n"
							 "C %d\t%d\t%d\n"
                             "C %d\t%d\t%d\n"
                             "C %d\t%d\t%d\n", x, y, z, x+salto_X, y+salto_Y , z+salto_Z, x+(2*salto_X), y+(2*salto_Y) , z+(2*salto_Z), x+(3*salto_X), y+(3*salto_Y) , z+(3*salto_Z));
#endif // MOVIE
                } else if(!rede[ x + salto_X ][ y + salto_Y ][ z + salto_Z ]) {
                    rede[ x + salto_X ][ y + salto_Y ][ z + salto_Z ] = rede[x][y][z];
                    rede[x][y][z] = 0;
                }
            } // fim de cada passo de monte carlo

#ifdef MOVIE
            particulas_video(rede, arq3);
#endif // MOVIE

		  if( ((int) floor(sis.numero_particulas * sis.fracao_liberacao) ) == cont_PART) break;
        } // fim do laço PMC
    }  // fim das simulacoes


    /* tempo final de execução */
    tempo_final = second();
    printf("# tempo de simulação: %f (segundos)  = %f (horas) ---> Parâmetros de entrada:  L: %d, N: %d, -p: %d, -t: %.3lf, -s: %d, -i: %d, -f: %f\n ", tempo_final - tempo_inicial, (tempo_final - tempo_inicial)*0.000277778, (sis.tam_rede-2), sis.numero_particulas, sis.numero_poros, sis.PMC, sis.numero_simulacoes, sis.tempo_impressao, sis.fracao_liberacao );

/// Gravando a media da liberação
    dados_arquivo_saida(particulas_tempo_real, particulas_tempo_real2);

#ifdef MOVIE
    fclose(arq1);
    fclose(arq2);
    fclose(arq3);
#endif // MOVIE

	/* liberando os vetores */
	des_alloc(rede, membrana, particulas_tempo_real, particulas_tempo_real2);

    return EXIT_SUCCESS;
}


///#################### FUNÇOES ######################
//---------------------------------------------------------------------------------------------
///* *** funçao parametros *** */
void parametros(int argc, char *argv[])
{
    register int i;

    sis.tempo_impressao = 1;
    sis.fracao_liberacao = 0.9999;
    sis.numero_poros = 1;

    for(i = 1; i < argc; i++){
        switch(argv[i][1]){

            case 'l': // tamanho da rede
                sis.tam_rede = atoi(argv[i+1]);
                break;
            case 'p':
                sis.numero_poros = atoi (argv[i+1]);   // número de buracos
                break;
            case 't': // passo de monte carlo
                sis.PMC = atoi(argv[i+1]);
                break;
            case 's': // numero de simulacoes
                sis.numero_simulacoes = atoi(argv[i+1]);
                break;
            case 'f':
                sis.fracao_liberacao = atof(argv[i+1]);
                break;
            case 'i': //tempo de impressao
                sis.tempo_impressao = atoi(argv[i+1]);
                break;
            case 'h': // help
                printf("\n  Modelagem estatística da liberação de fármacos\n"
                       "    Simulação em três-dimensões\n"
                       "    Orientador: Marco Aurélio A. Barbosa\n"
                       "    Márcio Sampaio Gomes filho - UnB\n"
                       "    Fevereiro 2017\n\n"
                       "    Opcoes de Enumeracao:\n"
                       "    -l         Tamanho da rede.\n"
                       "    -p         Número de poros. p <= 6*l^2 \n"
                       "    -n         Numero do particulas. Por padrão((N-2)²). \n"
                       "    -t         Passo de Monte Carlo. \n"
                       "    -s         Numero de simulações. \n"
                       "    -i         Tempo de impressão. \n"
                       "    -f         Fração da liberação a ser estudada. Ex:(-f 0.6 = 60 porcento). Por padrão -f = 0.9999\n"
                       "    -h         Imprime esta ajuda.\n\n");
                exit(1);
                break;

        } // fim do switch
    } // fim  do for

    /* condição de erro  */
    if(sis.numero_poros > 6*(sis.tam_rede - 2)*(sis.tam_rede - 2)) {
      printf("\nNúmero de Buracos invalido!\n");
      exit(1);
    }

    sis.numero_particulas = (sis.tam_rede - 2)*(sis.tam_rede - 2)*(sis.tam_rede - 2);
}
//---------------------------------------------------------------------------------------------

/* *** Função para alocar os pares ordenados da superfície da capsula em uma matriz *** */
//-------------------------------------------------------------------------
void escreve_membrana (int **membrana)
{
	int i,j,
        n1=0,
        n2=0,
		n3=0;

	/* arquivo com os pontos da superfície */
    FILE *arq=NULL;
    arq = fopen("membrana_inicial.dat", "w+");

    /* escrevendo um arquivo com os pontos que representam a superfície - rede cúbica */
    for(i = 2; i < sis.tam_rede; i++){
        for(j = 2; j < sis.tam_rede; j++){
			/* escrevendo os pontos da superfície no arquivo */
            fprintf(arq, "1\t%d\t%d\n"
			 "%d\t%d\t%d\n"
			 "%d\t1\t%d\n"
			 "%d\t%d\t%d\n"
			 "%d\t%d\t1\n"
			 "%d\t%d\t%d\n" , i,j,  sis.tam_rede,i,j,  i,j,   i, sis.tam_rede,  j,  i,j,  i,j, sis.tam_rede );
        }
    }

    /* fechando o arquivo */
    fclose(arq);

   /* Lendo o arquivo "membrana_inicial.dat" para aloca-lo como uma matriz na memória */
	do {
		arq = fopen("membrana_inicial.dat", "r");
        if(arq == NULL){
          printf("Erro ao abrir o arquivo!\n");
        }
   } while(arq == NULL);

	i=1;
	while (fscanf(arq, "%d %d %d\n", &n1, &n2, &n3) != EOF){
		membrana[i][1] = n1;
		membrana[i][2] = n2;
		membrana[i][3] = n3;
		i++;
    }

    /* fechando o arquivo e removendo o arquivo p/ não dar erro em simulações que serão rodadas ao mesmo tempo no cluster */
    fclose(arq);
    remove("membrana_inicial.dat");
}
//---------------------------------------------------------------------------------------------
/* *** função usada para formar a membrana do sistema ***  */
inline void forma_membrana(int ***rede, int **membrana)
{
	register int x, y, z, i, j, k;

	int	n_pontos    	    = 0,
		Npart		  		= 0,
		conta_numero_poros  = 0,
		coordenada_x        = 0;

    /* numero de pontos na membrana */
    n_pontos = 6*(sis.tam_rede - 2)*(sis.tam_rede - 2);

    /* zerando a rede */
    for(i = 1; i <= sis.tam_rede; i++){
       for(j = 1; j <= sis.tam_rede; j++){
            for(k = 1; k <= sis.tam_rede; k++){
                rede[i][j][k]=0;
            }
       }
    }

   /* formando a membrana */
	for(i = 1; i <= sis.tam_rede; i++){
		for(j = 1; j <= sis.tam_rede; j++){
			rede[1][i][j]       		 = COD_BORDA;   /// plano YZ, x = 1
            rede[sis.tam_rede][i][j]	 = COD_BORDA;   /// plano YZ, x = tam_rede
            rede[i][1][j]      		 	 = COD_BORDA;   /// plano XZ, Y = 1
            rede[i][sis.tam_rede][j] 	 = COD_BORDA;   /// plano XZ, Y = tam_rede
            rede[i][j][1]        		 = COD_BORDA;   /// plano XY, Y = 1
            rede[i][j][sis.tam_rede] 	 = COD_BORDA;   /// plano XY, Y = tam_rede
		}
   }

/* Descrição do Algoritmo
 *   1. sortea-se uma posição na membrana - sorteamos a coordenada x da matriz membrana
 *   2. verifica-se se tem um poro - caso não tenha - coloca um poro
 *   3. repete o processo até que tenha distribuido a quantidade de poros na superfície
 *
 * WARNING: Teste: L=52: 1 poro- tempo 0.041s e p/ 15000 poros- tempo 0.046s
*/

    do{
        coordenada_x = rand () % n_pontos + 1;

	if (  rede[ membrana[coordenada_x][1] ][ membrana[coordenada_x][2] ][ membrana[coordenada_x][3] ] !=  COD_SAIDA){ /*WARNING: na definição de membrana já excluimos os vértices */
	  rede[ membrana[coordenada_x][1] ][ membrana[coordenada_x][2] ][ membrana[coordenada_x][3] ] = COD_SAIDA;
	  conta_numero_poros++;
	}
    }while(conta_numero_poros < sis.numero_poros);

    /* preenchendo a matriz com as partículas */
    Npart=1;
    for (x = 2; x < sis.tam_rede; x++){
        for (y = 2; y < sis.tam_rede; y++){
            for (z = 2; z < sis.tam_rede; z++){
                rede[x][y][z] = Npart;
                Npart++;
            }
        }
    }
}

//-----------------------------------------------------------------------------------------
///* *** funçao para mostrar dados de impressão *** */
void dados_arquivo_saida(double *particulas_tempo_real, double *particulas_tempo_real2)
{
	int 	i,
			contador = 0;

	double  Nmedio 	= 0.0,
			Nmedio2	= 0.0,
			N2medio	= 0.0;

	printf("0	%d\n", sis.numero_particulas);

	for (i=1; i * sis.tempo_impressao < sis.PMC;i++){
		contador = i * sis.tempo_impressao;

		/* media */
		Nmedio = particulas_tempo_real[contador] /   (   sis.numero_simulacoes );
		/* quadrado da media */
		Nmedio2 = Nmedio*Nmedio;
		/* media do quadrado */
		N2medio = particulas_tempo_real2[contador] /   ( sis.numero_simulacoes );

		if ( Nmedio >= sis.numero_particulas * (1 - sis.fracao_liberacao) ){
            printf("%d	%lf %lf\n", contador, Nmedio, (double) sqrt((N2medio - Nmedio2)) );
		} else{
			break;
		}
	}
}
//-----------------------------------------------------------------------------------------
///* *** função para gravar a configuração em um arquivo.xyz *** */
void grava_configuracao_arq_video (int ***rede, FILE *arq1, FILE *arq2, FILE *arq3)
{
    int x, y, z;

    for (x = 1; x <= sis.tam_rede; x++){
        for (y = 1; y <= sis.tam_rede; y++){
            for (z = 1; z <= sis.tam_rede; z++){
                if(rede[x][y][z] == COD_BORDA){
                    fprintf(arq1, "%d\t%d\t%d\n", x, y, z);  ///gravando a membrana
                }else if(rede[x][y][z] == COD_SAIDA){
                    fprintf(arq2,"%d\t%d\t%d\n", x ,y, z);   /// gravando os poros
                } else if (rede[x][y][z] > 0){
                    fprintf(arq3,"%d\t%d\t%d\n", x ,y, z);    /// gravando as partículas
                }
            }
        }
    }
    fprintf(arq1,"\n");
    fprintf(arq2,"\n");
    fprintf(arq3,"\n");
}
//-----------------------------------------------------------------------------------------
///* *** funçao para salvar a posição das partículas *** */
inline void particulas_video (int ***rede, FILE *arq3)
{
  register int  i,j,k;

    for(i = 2; i < sis.tam_rede; i++){
        for(j = 2; j < sis.tam_rede; j++){
            for(k = 2; k < sis.tam_rede; k++){
                if(rede[i][j][k] > 0){
                    fprintf(arq3,"%d\t%d\t%d\n", i ,j, k);
                }
            }
        }
    }
    fprintf(arq3,"\n");
}
//-----------------------------------------------------------------------------------------
///* *** funçao para imprimir a rede *** */
void imprime_rede(int tam_rede, int ***rede)
{
    int i, j, k;

    printf("\n");
    for(i = 1; i <= tam_rede; i++){
       for(j = 1; j <= tam_rede; j++){
            for(k = 1; k <= tam_rede; k++){

                printf("%d\t", rede[i][j][k]);
            }
            printf("\n");
       }
        printf("\n");
    }

}
//-----------------------------------------------------------------------------------------
///* *** funçao para salto *** */
inline enum TIPO_SALTO salto_part(void)
{
    int salto;


    salto = rand () % 6;
    if(!salto) salto = 6;

    return (enum TIPO_SALTO) salto;
}

//-----------------------------------------------------------------------------------------
///* *** funçao desalocar os vetores *** */
void des_alloc(int ***rede, int **membrana, double *vetor1, double *vetor2 )
{
	int i, j;

	for(i = 0; i < sis.tam_rede + 1; i++){
		for(j = 0; j < sis.tam_rede + 1; j++){
			free(rede[i][j]);
        }
        free(rede[i]);
    }
    free(rede);

    for(i = 0; i < 6*(sis.tam_rede - 2)*(sis.tam_rede - 2) +  1; i++) free(membrana[i]);
    free(membrana);

	free(vetor1);
	free(vetor2);
}
//-----------------------------------------------------------------------------------------
///* *** funçao para verificar o intervalo de tempo *** */
int intervalo_medida(double t1, double t2)
{
    return (int) (floor(t2) - floor(t1));
}
//-----------------------------------------------------------------------------------------
/* Cálculo do tempo de execução do código */
double second() /* Returns elepsed seconds past from the last call to timer rest */
{

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}
////-----------------------------------------------------------------------------------------
/////* *** funçao para preencher a rede*** */
//inline int preenche_rede(int tam_rede,  int ***rede, FILE *arq3)
//{
//    register int x, y, z;
//    int Npart=1;
//
//    for (x = 2; x < tam_rede; x++){
//        for (y = 2; y < tam_rede; y++){
//            for (z = 2; z < tam_rede; z++){
//                rede[x][y][z] = Npart;
//               fprintf(arq3, "%d   %d  %d\n", x, y, z);
//                Npart++;
//            }
//        }
//    }
//
////    rede[3][3][3] = 10;
//
//    return (Npart-1);
//}


