#include <fstream>
#include "poisson.h"
#include <cmath>

poisson::poisson( double lx, double ly,
                  double nx, double ny,
                  double Px, double Py,
                  double lc, double Rho ):
    Lx( lx ), Ly( ly ),
    Nx( nx ), Ny( ny ),
    px( Px ), py( Py ),
    Lc( lc ), rho( Rho ){
        matrix = new double*[ Nx * Ny ];
        for( int i = 0; i < ( Nx * Ny ); ++i )
            matrix[ i ] = new double[ Nx * Ny ];
        
        matilde = new double*[ Nx * Ny ];
        for( int i = 0; i < ( Nx * Ny ); ++i )
            matilde[ i ] = new double[ Nx * Ny ];
        
        f      = new double[ Nx * Ny ];
        ftilde = new double[ Nx * Ny ];

        phi    = new double[ Nx * Ny ];
}

poisson::~poisson( void ){
    for( int i = 0; i < ( Nx * Ny ); ++i )
            delete matrix[ i ];
    
    delete matrix;
    delete      f;

    for( int i = 0; i < ( Nx * Ny ); ++i )
            delete matilde[ i ];

    delete matilde;
    delete ftilde;
    delete phi;
}

void poisson::calcH( void ){

    hx = Lx / Nx;
    hy = Ly / Ny;

    return;
}

void poisson::calcMatrix( void ){

    int l ,  m,
        i ,  j,
        i1, j1;

    for( l = 1; l <= ( Nx * Ny ); ++l ){
        j = round( l / Ny );
        i =    l - j * Nx  ;
        for( m = 1; m <= ( Nx * Ny ); ++m ){
            i1 = round( m / Ny );
            j1 =    m - j * Nx  ;

            if( i == 1 || i == l || j == 1 || j == l ){
                    if( l == m ) matrix[ l - 1 ][ m - 1 ] = 1;
                    else         matrix[ l - 1 ][ m - 1 ] = 0;
                }
            else{
                if( ( i1 == i + 1 || i1 == i - 1 ) && j1 == j ) 
                    matrix[ l - 1 ][ m - 1 ] = 1 / pow( hx, 2 );
                else
                    if( ( j1 == j + 1 || j1 == j - 1 ) && i1 == i )
                         matrix[ l - 1 ][ m - 1 ] = 1 / pow( hy, 2 );
                    else matrix[ l - 1 ][ m - 1 ] = 0;
            }
        }
    }
}

void poisson::calcf( void ){

    int xp  = px * Nx / Lx;
    int yp  = py * Ny / Ly;
    int Lcx = Lc * Nx / Lx;
    int Lcy = Lc * Ny / Ly; 

    int i, j; 
    for( int l = 1; l <= Nx * Ny; ++l ){
        j = round( l / Ny );
        i =    l - j * Nx  ;
        if( ( i >= xp && i <= xp + Lcx ) && ( j >= yp && j <= yp + Lcy ) )
             f[ l - 1 ] = rho;
        else f[ l - 1 ] =   0;
    }

    return;
}

void poisson::calcTilde( void ){
    
    for( int l = 0; l < Nx * Ny; ++l )
        for( int m = 0; m < Nx * Ny; ++m )
            if( l == m )    matilde[ l ][ m ] = 0;
            else            matilde[ l ][ m ] = matrix[ l ][ m ] /
                                                matrix[ l ][ l ];

    for( int l = 0; l < Nx * Ny; ++l )
        ftilde[ l ] = f[ l ] / matrix[ l ][ l ];

    return;
}

void poisson::gaussmethod( void ){

    double *phi0, *phi1;
    phi0 = new double[ Nx * Ny ];
    phi1 = new double[ Nx * Ny ];

    for( int i = 0; i < Nx * Ny; ++i )
        phi0[ i ] = ftilde[ i ];

    while( equal( phi0, phi1 ) == false ){
        phi0 = phi1;
        phi1 = multi( matilde, phi0 );
        for( int i = 0; i < Nx * Ny; ++i )
            phi1[ i ] += ftilde[ i ];    
    }

    phi = phi1;

    return;
}

double* poisson::multi( double** matrix, double* v ){

    double *product = new double[ Nx * Ny ];

    for( int i = 0; i < Nx * Ny; ++i )
        for( int k = 0; k < Nx * Ny; ++k )
            product[ i ] = matilde[ i ][ k ] * v[ k ];

    return product;
}


bool poisson::equal( double* a, double* b ){
    bool ans = true;

    for ( int i = 0; i < Nx * Ny; ++i )
        if( a[ i ] != b[ i ] )  ans = false;

    return ans;
}

void poisson::output( std::string nomefile ){

    int i, j;
    std::ofstream file( nomefile.c_str() );
    for( int l = 0; l < Nx * Ny; ++l ){
        j = round( l / Ny );
        i =    l - j * Nx  ;

        file << i << ' ' << j << ' ' << phi[ i ];
    }

    return;
}