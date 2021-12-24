# Equazione_di_Poisson_2D
Corso di Metodi Computazionali della Fisica. Anno Accademico 2019/2020.

Prof. Paolo Umari

Soluzione dell'equazione di Poisson in 2D attraverso il metodo di Gauss.

## Istruzioni per la compilazione
compile:

c++ -Wall -I. -o main *.cc


Initial data: 10, 10, 100, 100, 4, 4, 2, 80

gnuplot:

 - set hidden3d
 - splot ’output.txt ’ w l
