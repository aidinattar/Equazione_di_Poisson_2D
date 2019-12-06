#include <iostream>
#include "poisson.h"

int main( int argc, char* argv[] ){

    double Lx, Ly,
           px, py,
           Lc, rho;
    int Nx, Ny;
    std::string file;
    
    std::cout << "lati rettangolo:" 
              << std::endl;
    std::cin  >>  Lx  >> Ly;
    std::cout << "numero di punti su x e y:" 
              << std::endl;
    std::cin  >>  Nx  >> Ny;
    std::cout << "coordinate quadrato:" 
              << std::endl;
    std::cin  >>  px  >> py;
    std::cout << "lato quadrato:" 
              << std::endl;
    std::cin  >>  Lc;
    std::cout << "rho_0:" 
              << std::endl;
    std::cin  >>  rho;
    std::cout << "nome file output: ";
    std::cin  >> file; 

    poisson phi( Lx, Ly, Nx, Ny, px, py, Lc, rho );
    phi.calcH();
    phi.calcMatrix();
    phi.calcf();
    phi.calcTilde();
    phi.gaussmethod();
    phi.output( file );
    
    return 1;
}