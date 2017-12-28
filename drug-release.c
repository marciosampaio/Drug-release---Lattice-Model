/**#################################################################################
#                                                                                  #
#        MODELAGEM ESTATÍSTICA DA LIBERAÇÃO DE FÁRMACOS NANO-ENCAPSULADOS  	   #
#                                Simulação - 2D                                    #
#                                                                                  #
#                                                                                  #
#    Orientador: Marco A. Alves                                                    #
#    Márcio Sampaio Gomes Filho                                                    #
#                                                                                  #
#    Criado: 02/02/2012 - UnB                                                      #
#                                                                                  #
#      Descrição do Algoritmo:                                                     #
#           i)   Forma-se uma configuração inicial                                 #
#           ii)  Sortea-se uma partícula aleatoriamente                            #
#           iii) Move-se a particula aleatoriamente                                #
#           iv)  Verifica-se se a partícula saiu da rede                           #
#           v)   Dado um passo de MC, grava-se o tempo-MC e o número de particulas #
#                remanescentes                                                     #
#                                                                                  #
#                                                                                  #
####################################################################################**/

/*
*   - Início: 02/02/2012
*
*   23/02/12
*    - incluida enumeracao
*    - solicitacao doo nome dos arquivos
*
*   26/02/12
*    - incluida a função parametros com (argc e argv)
*    - incluida as structs
*
*  29/02/12 - nesta vesão a gravão é feita pelo terminal > tempo/cd1000.dat
*     - gravação pelo terminal
*
*  12/02/12
*    - melhoria de algumas funções do código(objetivo melhorar performace)
*    - função para preencher a matriz
*    - Usando ponteiro para a posição da particula
*
*  01/12/2012
*    - Correção final do TEMPO-DE-MONTE-CARLO
*
*  16/12/2012
*    - implementação da fração de poros -b
*    - sortea-se uma posição aleatoria na borda e coloca-se um buraco
*
*  24/10/2016
*    - revisão código
*    - desvio médio quadrático Delta_N
*    - tempo de simulação em CPU
*
*  26/10/2016
*   - inclusão das funções para vídeo
*   - shell-script para rodar e - se necessário gerar o vídeo de uma simulação.
*
*  28/10/2016
*   - Definindo a matriz retangular Ly por Lx , ou seja, utilizamos a notação de matrizes (ly linhas por lx colunas).
*   -  Debug e função da porosidade WARNING: a função de porosidade pode ser melhorada para casos mais gerais.
*   -  modificar a entrada do tamano da rede - gerar um vetor com todas as posições da membrana (pensar)
*
*  05/01/2017
*   - Nova Função distribuir os poros na superfície:
*      Função ou um arquivo  com os pares ordenados da membrana - aloca na memoria - então sortea-se uma posição para colocar um poro
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

#define COD_SAIDA -2
#define COD_BORDA -1

/* Cria arquivos para vizualização estrutural - Gnuplot -> compilação gcc    -D MOVIE */
//#define MOVIE

/* Abre funções de erro */
//#define DEBUG

typedef struct
{
    int lx,                       // tamanho da rede (representa o número de linhas em uma matriz A(lx x ly))
        ly,                       // tamanho da rede (representa o número de colunas em uma matriz A(lx x ly))
        numero_simulacoes,        // numero de simulações
        numero_particulas,        // n particulas
        tempo_impressao,          // tempo de impressao
        numero_poros;             // numero de poros na membrana
    float   fracao_liberacao;     // fração da liberação a ser
    double PMC;                   // passo de monte carlo

} sistema;

    sistema sis;

enum TIPO_SALTO
{
   direita=1,
   esquerda,
   acima,
   abaixo
};

int intervalo_medida(double t2, double t1);
void imprime_rede(int **y);
enum TIPO_SALTO salto_part (void);
void forma_membrana(int **y);
void parametros (int argc, char *argv[]);
void dados_arquivo_saida(double *particulas_tempo_real, double *particulas_tempo_real2);
void grava_configuracao_arq_video (int **rede, FILE *arq1, FILE *arq2, FILE *arq3);
void particulas_video (int **rede, FILE *arq3);
void membrana_escreve (int **membrana);
void forma_membrana_2(int **rede,  int **membrana);
double second();


int main (int argc, char *argv[])
{
    int i, j, k,                    // contador
        conta_particula,	         // conta as particulas que saem da rede
        x, y,                       // linha e colula
        **rede,                     // ponteiro de ponteiro para a rede 2D
        salto_x, salto_y,           // saltos
        **membrana;                 // matriz com os pares ordenados referentes a superfície da membrana (x y)

    double   tempo_mc, tempo_mc_antigo, *particulas_tempo_real,  *particulas_tempo_real2,  tempo_inicial, tempo_final;
    enum TIPO_SALTO salto;


#ifdef MOVIE
    FILE *arq1=NULL, *arq2=NULL, *arq3=NULL;
    arq1 = fopen("membrana.dat", "w+");
    arq2 = fopen("poros.dat"   , "w+");
    arq3 = fopen("particulas.dat", "w+");
#endif // MOVIE

    /* Començando a contar o tempo de simulação */
    tempo_inicial = second();

    /* paramentros  */
    parametros(argc, argv);
    srand(time(NULL)*argc);

    /* Condição de erro para o número de poros */
    if(sis.numero_poros > (2*sis.lx + 2*sis.ly)  ) {
        printf("\nNúmero de POROS invalido!\n");
        exit(1);
    }

///*############## alocando memória ######################

    /* WARNING:
    *   - ponteiro de ponteiro de ponteiro não é uma boa alternativa, pois a leitura na memoria é aleatória, - aqui podemos otimizar
    *   - aqui também mudamos para o usuário informar o apenas o tamanho interno da caixa de simulação - não informando a membrana
    *   - testado para tamanho extremamente grande de Lx:10000 e Ly:10000. */
    rede  = ( int** )   calloc ( sis.lx + 3, sizeof(  int * )); /// linhas na matriz - i -> Aij = Axy
    if (rede == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        exit(1);
    }
    for(i = 0; i < sis.lx + 3; i++) /// colunas na matriz - j
       rede[i] = (int*)  calloc (sis.ly + 3, sizeof( int ));


    /* vetor guarda o nº de particulas em cada instante de tempo: a ser utilizado no cálculo da média */
    particulas_tempo_real = (double*) calloc(sis.PMC + 1, sizeof(double));

    /* vetor guarda a o nº de particulas ao quadrado: a ser utilizado no cálculo da média do quadrado */
    particulas_tempo_real2 = (double *) calloc( sis.PMC + 1, sizeof( double ));

    membrana = ( int** )   calloc ( (2*sis.lx + 2*sis.ly) +  1, sizeof(  int * )); /// linhas na matriz - i -> Aij = Axy
    if (membrana == NULL) {
        printf ("** Erro: Memoria Insuficiente **");
        exit(1);
    }
    for(i = 0; i < (2*sis.lx + 2*sis.ly) + 1; i++) /// colunas na matriz - j
       membrana[i] = (int*)  calloc (2 + 1, sizeof( int ));


/* Iniciamos escrevendo em um arquivo "membrana_inicial.dat"  os pares ordenados que representam a superfície (aqui fica independente da geometria da membrana).
 * Aqui poderia ser um outro programa que passa os pontos.
 *  Neste caso, podemos ter tanto redes quadras quanto redes retangulares. */
    membrana_escreve(membrana);


/// ###### início da corrida ########

    for (i = 1; i <= sis.numero_simulacoes; i++){

        /* configuração inicial*/
        forma_membrana_2(rede, membrana); /* Função para qualquer geometria */
        //forma_membrana(rede, membrana); /* Função antiga - sortea-se uma aresta e posteriormente uma posição na borda */
        sis.numero_particulas = ( sis.lx ) * (sis.ly) ;
        conta_particula = tempo_mc = 0;

#ifdef MOVIE
   grava_configuracao_arq_video(rede, arq1, arq2, arq3);
#endif // MOVIE

#ifdef DEBUG
    printf("1) configuração inicial -> Simulação:%d ", i);
    imprime_rede(rede);
#endif

        for(j = 1; j <= sis.PMC; j++) {                             // inicio de cada passo de monte carlo

            for(k = 1; k <= sis.numero_particulas - conta_particula; k++){

                 /* sorteio da posição das partículas */

                /*WARNING: Aqui pela definição da rede precisamos nos preocupar como a rede é percorrida, aqui definimos a rede com a definição de matriz Aij = Axy*/
                 do {
                    x =  rand () % (sis.lx+1) + 1;
                    y  = rand () % (sis.ly+1) + 1;
                 } while(rede[x][y] <= 0);


               ///###### correção do tempo ######
               tempo_mc_antigo = tempo_mc;
               tempo_mc = tempo_mc + (   (double)  1  /  (sis.numero_particulas - conta_particula)   );

               if(intervalo_medida(tempo_mc_antigo, tempo_mc)) {
                   particulas_tempo_real[  (int) floor(tempo_mc) ] += (  (double)  (sis.numero_particulas - conta_particula) );
                   particulas_tempo_real2[ (int) floor(tempo_mc) ] += (  ( (double)  sis.numero_particulas - conta_particula) * (  (double) sis.numero_particulas - conta_particula) );
               }

               ///#### fim da correção ########

               /* sorteio do salto */
               salto = salto_part();
               salto_x = salto_y = 0;
               switch(salto) {
                   case direita:
                       salto_x = 1;
                       break;
                   case esquerda:
                       salto_x = -1;
                       break;
                   case acima:
                       salto_y = 1;
                       break;
                   case abaixo:
                       salto_y = -1;
                       break;
               }

               /* troca a posicão*/
               if( rede[ x + salto_x ][ y + salto_y ] == COD_SAIDA ) {   // sisdicao para p/ a saida da particula ** e sistar

#ifdef DEBUG
  printf("Sortea-se uma particula aleatoriamente: (ALERTA: X é o numero de linhas e Y: colunas)\n"
	"partcula: %d -> posicao (%d, %d)\n"
	"salto_x: %d  salto_y: %d, rede[ x + salto_x ][ y + salto_y ]: %d \n\n", rede[x][y], x, y, salto_x, salto_y, rede[ x + salto_x ][ y + salto_y ]);
#endif

#ifdef MOVIE
   fprintf(arq3,"%d    %d  \n"
                "%d    %d  \n"
                "%d    %d  \n"
                "%d    %d  \n", x , y, x+salto_x, y+salto_y ,  x+(2*salto_x), y+(2*salto_y), x+(3*salto_x), y+(3*salto_y));
#endif // MOVIE

                rede[x][y] = 0;
                conta_particula++;

                /* troca de posição aleatorimente - interação caroço duro ou volume excluído */
               } else if( !rede[ x + salto_x ][ y + salto_y ] ){       // codicao para buraco
                   rede[x + salto_x][ y + salto_y ] =  rede[x][y];
                   rede[x][y] = 0;
               }
           }  // fim de cada passo de MC

#ifdef MOVIE
    particulas_video(rede, arq3);
#endif // MOVIE

           if( ((int) floor(sis.numero_particulas  * sis.fracao_liberacao) ) == conta_particula) break;
        } // fim do passo de monte carlo
    } // fim de cada simulacao

    /*WARNING: tempo de simulação: não considera o tempo para calcular a média e o desvio  */
    tempo_final = second();
    printf("# tempo de simulação: %f (segundos)  = %f (horas) ---> Parâmetros de entrada:  -lx: %d, -ly: %d, -p: %d, -t: %.3lf, -s: %d, -i: %d, -f: %f\n ", tempo_final - tempo_inicial, (tempo_final - tempo_inicial)*0.000277778, sis.lx, sis.ly, sis.numero_poros, sis.PMC, sis.numero_simulacoes, sis.tempo_impressao, sis.fracao_liberacao );

    /* calculo das medias */
   dados_arquivo_saida(particulas_tempo_real, particulas_tempo_real2);



#ifdef MOVIE
    fclose(arq1);
    fclose(arq2);
    fclose(arq3);
#endif // MOVIE

    /*liberando os vetores e fechando os arquivos */
    for(i = 0; i < sis.lx + 3; i++) free(rede[i]);
    free(rede);
    for(i = 0; i < (2*sis.lx + 2*sis.ly) + 1; i++) free(membrana[i]);
    free(membrana);
    free(particulas_tempo_real);
    free(particulas_tempo_real2);


    return EXIT_SUCCESS;
}
//-----------------------------------------------------------------------------------------
///* *** funcao para tempo_real  *** */
int intervalo_medida(double t2, double t1)
{
    return (int) (floor(t2) - floor(t1));
}
//-----------------------------------------------------------------------------------------
///* *** funçao para mostrar dados de impressão *** */
void dados_arquivo_saida(double *particulas_tempo_real, double *particulas_tempo_real2)
{
    int contador,i;
    double Nmedio, Nmedio2,   N2medio;

    contador=0;
	Nmedio =  Nmedio2 =   N2medio = 0.0;

    printf("0 \t %d\n", sis.numero_particulas);
    for (i=1; i * sis.tempo_impressao < sis.PMC; i++){
        contador = i * sis.tempo_impressao;

        Nmedio = particulas_tempo_real[contador] /   (   (double) sis.numero_simulacoes );
        Nmedio2 = (double) Nmedio*Nmedio;

        N2medio = particulas_tempo_real2[contador] /  (   (double) sis.numero_simulacoes );

        if ( Nmedio >= sis.numero_particulas * (1 - sis.fracao_liberacao) ){
            printf("%d \t %lf \t %lf\n", contador, Nmedio, (double) sqrt((N2medio - Nmedio2)) );
            //	fprintf(arq4, "%d	%lf\n", contador, media);
       } else{
            break;
       }
    }
}

//-----------------------------------------------------------------------------------------
///* *** função para gravar a configuração em um arquivo.xyz *** */
void grava_configuracao_arq_video (int **rede, FILE *arq1, FILE *arq2, FILE *arq3)
{
    register int i, j;

    for (i = 1; i <= sis.lx + 2 ; i++){
        for (j = 1; j <= sis.ly + 2; j++){
            if(rede[i][j] == COD_BORDA){
                fprintf(arq1, "%d    %d \n", i, j);  ///gravando a membrana
            }else if(rede[i][j] == COD_SAIDA){
                fprintf(arq2,"%d    %d  \n", i, j);   /// gravando os poros
            } else if (rede[i][j] > 0){
                fprintf(arq3,"%d    %d  \n", i, j);    /// gravando as partículas
            }
        }
    }
    fprintf(arq1,"\n");
    fprintf(arq2,"\n");
    fprintf(arq3,"\n");
}
//-----------------------------------------------------------------------------------------
///* *** funçao para salvar a posição das partículas *** */
inline void particulas_video (int **rede, FILE *arq3)
{
  register int  i,j;

    for(i = 2; i <= sis.lx + 1; i++){
        for(j = 2; j <= sis.ly + 1; j++){
            if(rede[i][j] > 0){
                fprintf(arq3,"%d    %d  \n", i ,j);
            }
        }
    }
    fprintf(arq3,"\n");
}
//--------------------------------------------------------------------------
///* *** funçao para impressao da rede *** */
void imprime_rede(int **rede){
    int i, j;
    printf("\n");
    /*varrendo cada linha da matriz  */
    for(i = 1; i <= sis.lx + 2; i++) {
        for(j = 1; j <= sis.ly + 2; j++) {
           printf("%d\t", rede[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
 }
//--------------------------------------------------------------------------
///* *** funçao para o salto da particula *** */
inline enum TIPO_SALTO salto_part (void)
{
    int salto;

    salto = rand () % 4;

    // salto = (salto == 0 ) ? 4 : salto;

    if (!salto) salto = 4;

    return (enum TIPO_SALTO) salto;
}

/* *** Função para alocar os pares ordenados da superfície da capsula em uma matriz *** */
//-------------------------------------------------------------------------
void membrana_escreve (int **membrana)
{
    int i,
        n1=0,
        n2=0;

    /* arquivo com os pontos da superfície */
    FILE *arq=NULL;
    arq = fopen("membrana_inicial.dat", "w+");

    /* escrevendo um arquivo com os pontos que representam a superfície - rede retangular */
    /* escrevendo as laterais */
    for(i = 2; i < sis.lx + 2; i++){
     fprintf(arq, "%d\t1\n"
	     "%d\t\%d\n", i, i, sis.ly + 2 );
    }

    /* parte superior e inferior  */
    for(i = 2; i < sis.ly + 2; i++){
      fprintf(arq, "1\t%d\n"
	"%d\t\%d\n", i, sis.lx + 2, i );
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
    while (fscanf(arq, "%d %d\n", &n1, &n2) != EOF){
      membrana[i][1] = n1;
      membrana[i][2] = n2;
      i++;
    }

    /* fechando o arquivo e removendo o arquivo p/ não dar erro em simulações que serão rodadas ao mesmo tempo no cluster */
    fclose(arq);
    remove("membrana_inicial.dat");
}

/* *** função usada para formar a membrana do sistema ***  */
inline void forma_membrana_2(int **rede, int **membrana)
{
    register int i, j;

    int    n_linhas             = 0,
           conta_numero_poros   = 0,
           coordenada_x         = 0,
           k                    = 1;

    /* numero de linhas na matriz membrana */
    n_linhas = (2*sis.lx + 2*sis.ly);

    /* zerando a rede */
    for(i = 1; i <= sis.lx + 2; i++) {
        for(j = 1; j <= sis.ly + 2; j++) {
            rede[i][j] = 0;
        }
    }

    /* formando a  membrana */
    for(i = 1; i <= sis.ly + 2; i++){ /* parte superior e inferior  */
      rede[1][i] = COD_BORDA;
      rede[sis.lx + 2][i] = COD_BORDA;
    }

    for(i = 1; i <= sis.lx + 2; i++){ /* partes laterais */
      rede[i][1] = COD_BORDA;
      rede[i][sis.ly + 2] = COD_BORDA;
    }

/* Descrição do Algoritmo
 *   1. sortea-se uma posição na membrana - sorteamos a coordenada x da matriz membrana
 *   2. verifica-se se tem um poro - caso não tenha - coloca um poro
 *   3. repete o processo até que tenha distribuido a quantidade de poros na superfície
 *
 * WARNING: Teste p/ "forma_membrana_2" para Lx = 2000 e Ly = 2000  - totalmente preenchida com poros - e o tempo para a criação da estrutura inicial é bastante razoável  0.629502 (segundos).
*/

    do{
        coordenada_x = rand () % n_linhas + 1;    /// sorteando uma membrana (ou seja, aresta da rede)

        if (  rede[membrana[coordenada_x][1]][membrana[coordenada_x][2]] !=  COD_SAIDA){ /*WARNING: na definição de membrana já excluimos os vértices */
		rede[ membrana[coordenada_x][1] ][ membrana[coordenada_x][2] ] = COD_SAIDA;
	    	conta_numero_poros++;
	}
    }while(conta_numero_poros < sis.numero_poros);

    /* preenchendo a matriz com as partículas */
    for (i = 2; i <= sis.lx + 1 ; i++) {
        for(j = 2; j <= sis.ly + 1; j++){
            rede[i][j] = k;
            k++;
        }
    }

}
//--------------------------------------------------------------------------
///* *** funçao parametros *** */
void parametros (int argc, char *argv[])
{
    register int i;

    sis.tempo_impressao = 1;
    sis.numero_poros = 1;
    sis.fracao_liberacao = 0.9999;

    for(i = 1; i < argc; i++) {

        switch(argv[i][1]) {

            case 'x':
               sis.lx = atoi(argv[i+1]);           // tamanho da rede (direção x)
                break;
            case 'y':
               sis.ly = atoi(argv[i+1]);           // tamanho da rede (direção y)
                break;
            case 'p':
                sis.numero_poros = atoi (argv[i+1]);   // número de poros
                break;
            case 'n':
                sis.numero_particulas = atoi(argv[i+1]); // numero de particulas - inativo
                break;
            case 't':
                sis.PMC = atoi(argv[i+1]);               // passo de monte carlo
                break;
            case 's':
                sis.numero_simulacoes = atoi(argv[i+1]);
                break;
            case 'f':
                sis.fracao_liberacao = atof(argv[i+1]);
                break;
            case 'i':
                sis.tempo_impressao = atoi(argv[i+1]);
                break;
            case 'h':
                    printf("\n  Modelagem estatística da liberação de fármacos nano-encapsulados\n"
                           "    Simulação em duas-dimensões\n"
                           "    Orientador: Marco Aurélio A. Barbosa\n"
                           "    Márcio Sampaio Gomes filho - UnB\n"
                           "    Janeiro/2012\n\n"
                           "    Opcoes de Enumeracao:\n"
                           "    -x         Tamanho da rede  (direção x). Ou melhor, similar a matriz, temos Aij = Axy, x linhas e y colunas.\n"
                           "    -y         Tamanho da rede  (direção y).\n"
                           "    -p         Número de poros. (definição retangular) \n"
                           "    -n         Numero do particulas. Por padrão((N-2)²). \n"
                           "    -t         Passo de Monte Carlo. \n"
                           "    -s         Numero de simulações. \n"
                           "    -i         Tempo de impressão (por padrão i=1). \n"
                           "    -f         Fração da liberação a ser estudada. Ex:(-f 0.6 = 60 porcento). Por padrão -f = 0.9999\n"
                           "    -h         Imprime esta ajuda.\n\n");
                exit(1);
                break;
        } // fim do switch
    } // fim do for
}

/* Cálculo do tempo de execução do código */
double second() /* Returns elepsed seconds past from the last call to timer rest */
{

    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}
//--------------------------------------------------------------------------
///* *** funçao para formar as paredes da rede *** */
////// ESTA FUNCAO PODE SER OTIMIZADA ATRAVES DE UMA VARREDURA COBRINDO SOMENTE AS BORDAS!
//inline void parede_rede(int TAM_REDE, int **rede){
//    register int i,j;
//
//    for(i = 1; i <= sis.TAM_REDE; i++) {
//        for(j = 1; j <= sis.TAM_REDE; j++) {
//            if(i == 1 ||  j == 1 || TAM_REDE == i || TAM_REDE== j) { // sisdiçao para a borda da rede
//                rede[i][j] = COD_BORDA;
//            }
//           /* if(j == TAM_REDE) {      // buraco para efusão....  ok - testado -2
//                efusao =  TAM_REDE / 2;
//                rede[efusao][j] = -2;
//            }*/
//        }
//    }
//    rede[TAM_REDE/2][TAM_REDE]=COD_SAIDA;
//}
//--------------------------------------------------------------------------
///* *** distribuindo as particulas aleatoriamente na rede *** */
//inline void distribui_particulas (int TAM_REDE, int **rede, int numero_particulas)
//{
//    int x, y;
//    register int i;
//
//    for(i = 1; i <= sis.numero_particulas; i++) {
//        do {
//            x = rand () % sis.TAM_REDE + 1;
//            y = rand () % sis.TAM_REDE + 1;
//        } while( rede[x][y] != 0 );   // evita caroço duro
//
//        rede[x][y] = i;
//    }
//}
//
// /* *** função usada para formar a parede do sistema ***  */
// inline void forma_membrana(int **rede)
// {
//     register int i, j;
//
//
//     int num_aleatorio           = 0,
//     numero_poros_1              = 0,
//     numero_poros_2              = 0,
//     numero_poros_3              = 0,
//     numero_poros_4              = 0,
//     conta_numero_poros          = 0,
//     membrana                    = 0,
//     k                           = 1;
//
//     /* zerando a rede */
//     for(i = 1; i <= sis.lx + 2; i++) {
//         for(j = 1; j <= sis.ly + 2; j++) {
//             rede[i][j] = 0;
//         }
//     }
//
//     /* formando a  membrana */
//     for(i = 1; i <= sis.ly + 2; i++){ /* parte superior e inferior  */
//       rede[1][i] = COD_BORDA;
//       rede[sis.lx + 2][i] = COD_BORDA;
//     }
//     for(i = 1; i <= sis.lx + 2; i++){ /* partes laterais */
//       rede[i][1] = COD_BORDA;
//       rede[i][sis.ly + 2] = COD_BORDA;
//     }
//
// /* sorteando os poros aleatoriamente, excluindo os vertices!
// *    1. sortea-se uma aresta
// *    2. coloca-se um poro
// *
// * WARNING: é possível pensar em um único vetor (ou matriz) com as posições da membrana, aqui pode ser melhorado, no entanto, testei para Lx = 2000 e Ly = 2000  e o tempo para a criação da estrutura inicial é bastante razoável 5.80 (segundos).
// */
//     do{
//         membrana = rand () % 4 + 1;    /// sorteando uma membrana (ou seja, aresta da rede)
//
//         switch(membrana){
//             case 1:
//                 if(numero_poros_1 < sis.ly){ /// Retiramos as arestas
//                     do {
//                         num_aleatorio = rand() % sis.ly + 2;
//                      }while( rede[1][num_aleatorio ] != COD_BORDA);
//
//                      rede [1][num_aleatorio]  = COD_SAIDA;
//                      numero_poros_1++;
//                 }
//                 break;
//             case 2:
//                 if(numero_poros_2 < sis.ly){ /// Retiramos as arestas
//                     do {
//                         num_aleatorio = rand() % sis.ly + 2;
//                     }while( rede[sis.lx + 2][num_aleatorio] != COD_BORDA ) ;
//
//                     rede [sis.lx + 2][num_aleatorio]  = COD_SAIDA;
//                     numero_poros_2++;
//                 }
//                 break;
//             case 3:
//                  if(numero_poros_3 < sis.lx){ /// Retiramos as arestas
//                      do {
//                          num_aleatorio = rand() % sis.lx + 2;
//
//                      }while( rede[num_aleatorio][1] != COD_BORDA);
//
//                      rede[num_aleatorio][1]   = COD_SAIDA;
//                      numero_poros_3++;
//                  }
//                  break;
//             case 4:
//                 if(numero_poros_4 < sis.lx){ /// Retiramos as arestas
//                     do {
//                         num_aleatorio = rand() % sis.lx + 2;
//                      }while( rede [ num_aleatorio ][ sis.ly + 2 ] != COD_BORDA) ;
//
//                      rede [ num_aleatorio ][ sis.ly + 2 ]  = COD_SAIDA;
//                      numero_poros_4++;
//                  }
//                  break;
//         }
//
//         conta_numero_poros =  numero_poros_1 + numero_poros_2  + numero_poros_3 + numero_poros_4;
//     }while(conta_numero_poros < sis.numero_poros);
//
//     /* preenchendo a matriz com as partículas */
//     for (i = 2; i <= sis.lx + 1 ; i++) {
//         for(j = 2; j <= sis.ly + 1; j++){
//             rede[i][j] = k;
//             k++;
//         }
//     }
// }
//
