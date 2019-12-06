compile:

c++ -Wall -I. -o main *.cc


Initial data: 10, 10, 100, 100, 4, 4, 2, 80

gnuplot:

set hidden3d
splot ’output.txt ’ w l