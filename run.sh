#!/bin/bash

#  DRUG RELEASE - 2D ( Márcio Sampaio, Marco Aurélio e Fernando Albuquerque) - 07/05/2016
#  	- Sricpt desenvolvido para rodar as simulacoes em 2D - poroso
# 	- Com opção de vídeo OBS: aqui o vídeo é para uma simulação -> Compilar com a opção de -D MOVIE

#- $MOVIE /// 0:OFF   1:ON -> abre arquivos para vizualização estrutural - Gnuplot 

# WARNING:vídeo exclusivo para uma simulação!!!
MOVIE=0

##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
#    			COMPILANDO O CÓDIGO
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
case "${MOVIE}" in   "1")
  gcc drug-release.c  -o nano -lm -O3  -D MOVIE  	# movie open 
  ;;
  
  # Roda o código semm a opção de gerar o vídeo da simulação
  "0")
  gcc drug-release.c  -o nano -lm -O3      			# movie off
  ;;
  esac

# Impressão do help
#./ nano2d -h 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# "\n  Modelagem estatística da liberação de fármacos nano-encapsulados\n"
# "    Opcoes de Enumeracao:\n"
#     "    Opcoes de Enumeracao:\n"
#     "    -x         Tamanho da rede  (direção x). Ou melhor, similar a matriz, temos Aij = Axy, x linhas e y colunas.\n"
#     "    -y         Tamanho da rede  (direção y).\n"
#     "    -p         Número de poros. (definição retangular) \n"
#     "    -n         Numero do particulas. Por padrão((N-2)²). \n"
#     "    -t         Passo de Monte Carlo. \n"
#     "    -s         Numero de simulações. \n"
#     "    -i         Tempo de impressão. \n"
#     "    -f         Fração da liberação a ser estudada. Ex:(-f 0.6 = 60 porcento). Por padrão -f = 0.9999\n"
#     "    -h         Imprime esta ajuda.\n\n");


##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
#    			RODANDO
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
lx=10;		     	 #número de linhas
ly=10;                   #número de colunas
N=$(( lx*ly ));          # número de partículas
poros=15;                # numero de poros (máximo pmax = 2*lx + 2*ly)
n_passos=1000;    	 # passos de MC (~ tempo) 
s=1000;			 # número de simulações 

echo ''
echo  Parâmetros de entrada: 
echo  Lx=$lx, Ly=$ly N=$N Poros=$poros s=$s, PMC=$n_passos 

   ./nano -x $lx  -y $ly -p $poros -t $n_passos -s $s -i 2 -f 0.9999 > ly$ly-lx$lx-p$poros-s$s.dat

 arq=ly$ly-lx$lx-p$poros-s$s.dat
#---------------------------------------------------------------
#			 VÍDEO: READ ME			       #
#---------------------------------------------------------------

#  1. Images of particles and membrane 		(created with gnuplot software)
#  2. Dynamical plot 				(gnuplot)
#  2.1 Dynamical plot of  release	 	(gnuplot)
#  3. Putting all figures together  		(montage)
#  4. Created movie 				(ffmpeg)
#---------------------------------------------------------------


case "${MOVIE}" in "1")


#Pengando o tempo final no arquivo
tfinal=$(awk 'END{print $1}' < $arq)

largura_imagemX=$lx+2;
largura_imagemY=$ly+2;

eixo_y=$N+10;
eixo_x=$tfinal+10;

## Criando pastas no diretorio
mkdir estrutura;
mkdir liberacao;
mkdir juntas;

#$($t_final+1)
#Montando as figuras da simulação
for ((i=1; i <= $tfinal+1; i++))
  do gnuplot -e "set terminal jpeg; 
	set xrange [*:$largura_imagemX];
	set yrange [*:$largura_imagemY];
	unset border;
	unset xtics;
	unset ytics; 
	set format x '';
	set format y '';
	set nokey;
	set view 70,96;
	set xyplane 0;
	plot 'particulas.dat' every :::$i::$i  w points pt 7 lc 3, 'membrana.dat'  w points pt 5 lc 1, 'poros.dat'   w points pt 5 lc -1" > estrutura/${i}.jpeg 
  done;
# Teste para a dinâmica dos poros 
#splot 'membrana.xyz' every :::$i::$i  w points pt 5, 'poros.xyz' every :::$i::$i  w points pt 5 " > membrana/membrana${i}.jpeg 


## Criando gráfico do número de particulas em função do tempo
for ((i=1; i <= $tfinal+1; i++))
  do gnuplot -e "set terminal jpeg; 
  set yrange [0:${eixo_y}];
  set xrange [0:${eixo_x}];
  set xlabel 'time[MCS]';
  set ylabel 'N(t)';
  set arrow from 20.10,3.8 to 20.10,6 nohead;
  plot '<head -n $i $arq' w lines lw 2 notitle" > liberacao/txSt$i.jpeg;
done;
	
#Formando uma figura só!
for ((i=1; i <= $tfinal+1; i++))
      do
	    montage estrutura/${i}.jpeg   liberacao/txSt$i.jpeg -geometry +2+2 juntas/juntos$i.jpeg;
      done;

#Criando o vídeo release-3D.mpeg
ffmpeg  -r 20 -i juntas/juntos%d.jpeg  ly$ly-lx$lx-p$poros-s$s.mpeg; # para simulações longas

#ffmpeg -r 15 -i juntas/juntos%d.jpeg - vcodec flv  L$L-P$poros.flv

#Fim do case
;;
# Roda o código normalmente 



  "0")
  echo ''
  echo "Video > off"
  ;;

  *)
  echo "Opcao invalida ... Digite 1 (sim) ou 0 (nao): video"
  ;;
  esac
 
#Apagando pastas!
rm -rf estrutura;
rm -rf liberacao;
rm -rf juntas;  # excluindo figuras p/não dar erro no vídeo!


