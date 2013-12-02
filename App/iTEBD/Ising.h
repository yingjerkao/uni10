#ifndef TFISING_H
#define TFISING_H

double Sx[ 2 * 2 ] = { 0.0, 1.0, 1.0, 0.0 };
double Sz[ 2 * 2 ] = { 1.0, 0.0, 0.0, -1.0 };

double * tfIsing( double transField, double longiField = 0.0 );


double * tfIsing( double transField, double longiField )
{
    const int size = 16;
    double * hamiltonian = new double [ size ];
    for ( int i = 0; i < size; ++i )
        hamiltonian[ i ] = 0.0;

    // H0 = -Sz * Sz
    hamiltonian[ 0 ] -= 1.0;
    hamiltonian[ 5 ] += 1.0;
    hamiltonian[ 10 ] += 1.0;
    hamiltonian[ 15 ] -= 1.0;

    if ( transField != 0.0 )
    {
        // H = H0 + transField * Sx
        hamiltonian[ 1 ] += transField / 2.0, hamiltonian[ 2 ] += transField / 2.0;
        hamiltonian[ 4 ] += transField / 2.0, hamiltonian[ 7 ] += transField / 2.0;
        hamiltonian[ 8 ] += transField / 2.0, hamiltonian[ 11 ] += transField / 2.0;
        hamiltonian[ 13 ] += transField / 2.0, hamiltonian[ 14 ] += transField / 2.0;
    }

    if ( longiField != 0.0 )
    {
        // H = H0 + longiField * Sz
        hamiltonian[ 0 ] += longiField, hamiltonian[ 15 ] -= longiField;
    }

    return hamiltonian;
}

#endif
