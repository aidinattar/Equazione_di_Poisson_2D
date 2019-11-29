#ifndef POISSON_H
#define POISSON_H

class poisson{

    public:
        poisson( double, double,
                 double, double,
                 double, double,
                 double, double );
        ~poisson( void );

        void calcH      ( void );
        void calcMatrix ( void );
        void calcf      ( void );
        void calcTilde  ( void );
        void calcPhi    ( void );
        void gaussmethod( void );

        double* multi   ( double**, double* );
        bool equal      ( double*, double* );
        void output     ( std::string nomefile );

    private:
        double  Lx;
        double  Ly;
        int     Nx;
        int     Ny;
        double  px;
        double  py;
        double  Lc;
        double rho;

        double hx;
        double hy;

        double **matrix;
        double **matilde;
        double *f;
        double *ftilde;
        double *phi;
};

#endif