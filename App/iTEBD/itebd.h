#ifndef ITEBD_H
#define ITEBD_H

#include <mkl.h>
#include "SyTensor.h"

#include "Ising.h"

void prod( double * gamma, const double * data, int bondDim, int physicalDim, bool right = true );

void invData( double * data, int n );

void expm( double * outMatrix, int n, double alpha, const double * inMatrix, int order = 20 );

void truncate( double * dest, const double * data, int row, int col, int n, bool truncRow = true );

double energy( const std::vector<SyTensor_t> & gamma, const std::vector<SyTensor_t> & l, double * hermitonian, int physicalDim );

double expectionValue( const std::vector<SyTensor_t> & gamma, const std::vector<SyTensor_t> & l, double * oper, int physicalDim );

double entanglementEtropy( const SyTensor_t & l );


// ======================================================
void expm( double * outMatrix, int n, double alpha, const double * inMatrix, int order )
{
    int size = n * n;
    for ( int i = 0; i < size; ++i )
        outMatrix[ i ] = 0.0;

    for ( int i = 0; i < n; ++i )
        outMatrix[ i * n + i ] = 1.0;

    double * matrixPow = new double [ size ];
    double * tmpMatrix = new double [ size ];

    for ( int i = 0; i < size; ++i )
        matrixPow[ i ] = outMatrix[ i ];

    int inc = 1;
    double alpha2;
    double beta = 0.0;
    double * ptr;
    for ( int i = 1; i < order; ++i )
    {
        alpha2 = alpha / i;
        dgemm( "N", "N", &n, &n, &n, &alpha2, inMatrix, &n, matrixPow, &n, &beta, tmpMatrix, &n );

        ptr = tmpMatrix;
        tmpMatrix = matrixPow;
        matrixPow = ptr;

        alpha2 = 1.0;
        daxpy( &size, &alpha2, matrixPow, &inc, outMatrix, &inc );

    } 

    delete [] matrixPow;
    delete [] tmpMatrix;
}


void prod( double * gamma, const double * data, int bondDim, int physicalDim, bool right )
{
    int size = physicalDim * bondDim;
    int offset;
    if ( right )
    {
        for ( int i = 0; i < size; ++i )
        {
            offset = i * bondDim;
            for ( int j = 0; j < bondDim; ++j )
                gamma[ offset + j ] *= data[ j ];
        }
    }
    else
    {
        for ( int i = 0; i < bondDim; ++i )
        {
            offset = i * size;
            for ( int j = 0; j < size; ++j )
                gamma[ offset + j ] *= data[ i ];
        }
    }
}


void invData( double * data, int n )
{
    for ( int i = 0; i < n; ++i )
        data[ i ] = 1.0 / data[ i ];
}


void truncate( double * dest, const double * data, int row, int col, int n, bool truncRow )
{
    int size = truncRow ? n * col : row * n;
    if ( truncRow )
    {
        for ( int i = 0; i < size; ++i )
            dest[ i ] = data[ i ];
    }
    else
    {
        int offset1, offset2;
        for ( int i = 0; i < row; ++i )
        {
            offset1 = i * n;
            offset2 = i * col;
            for ( int j = 0; j < n; ++j )
                dest[ offset1 + j ] = data[ offset2 + j ];
        }
    }
}


double energy( const std::vector<SyTensor_t> & gamma, const std::vector<SyTensor_t> & l, double * hermitonian, int physicalDim )
{
    const int bondDim = l[ 0 ].getElemNum();
    int tensorBondLabel[ 4 ];
    std::vector<Qnum_t> qnums;
    std::vector<Bond_t> bonds;

    for ( int i = 0; i < physicalDim; ++i )
        qnums.push_back( Qnum_t() );

    bonds.push_back( Bond_t( BD_ROW, qnums ) ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
    bonds.push_back( Bond_t( BD_COL, qnums ) ); bonds.push_back( Bond_t( BD_COL, qnums ) );
    tensorBondLabel[ 0 ] = 111, tensorBondLabel[ 1 ] = 222, tensorBondLabel[ 2 ] = -1, tensorBondLabel[ 3 ] = 1;

    SyTensor_t H( bonds, tensorBondLabel, "hermitonian" );
    H.addRawElem( hermitonian );

    SyTensor_t theta, Htheta;

    int order[ 4 ] = { 11, 111, 222, 22 };
    int size = bondDim * physicalDim * physicalDim * bondDim;
    int inc = 1;
    double energy[ 2 ];
    for ( int A = 0; A < 2; ++A )
    {
        int B = ( A + 1 ) % 2;

        theta = gamma[ A ];
        prod( theta.elem, l[ A ].elem, bondDim, physicalDim, true );
        prod( theta.elem, l[ B ].elem, bondDim, physicalDim, false );

        Htheta = gamma[ B ];
        prod( Htheta.elem, l[ B ].elem, bondDim, physicalDim, true );

        tensorBondLabel[ 0 ] = 11, tensorBondLabel[ 1 ] = -1, tensorBondLabel[ 2 ] = 0;
        theta.addLabel( tensorBondLabel );

        tensorBondLabel[ 0 ] = 0, tensorBondLabel[ 1 ] = 1, tensorBondLabel[ 2 ] = 22;
        Htheta.addLabel( tensorBondLabel );

        theta *= Htheta;
        Htheta = theta * H;

        Htheta.reshape( order, 2 );

        energy[ A ] = ddot( &size, Htheta.elem, &inc, theta.elem, &inc );
    }

    return ( energy[ 0 ] + energy[ 1 ] ) / 2.0;
}


double expectionValue( const std::vector<SyTensor_t> & gamma, const std::vector<SyTensor_t> & l, double * oper, int physicalDim )
{
    const int bondDim = l[ 0 ].getElemNum();
    std::vector<Qnum_t> qnums;
    std::vector<Bond_t> bonds;
    int tensorBondLabel[ 3 ];

    for ( int i = 0; i < physicalDim; ++i )
        qnums.push_back( Qnum_t() );

    bonds.push_back( Bond_t( BD_ROW, qnums ) ); bonds.push_back( Bond_t( BD_COL, qnums ) );
    tensorBondLabel[ 0 ] = 111, tensorBondLabel[ 1 ] = 0;

    SyTensor_t Op( bonds, tensorBondLabel, "operator" );
    Op.addRawElem( oper );

    SyTensor_t theta, OpTheta;
    tensorBondLabel[ 0 ] = 1, tensorBondLabel[ 1 ] = 0, tensorBondLabel[ 2 ] = 2;
    int order[ 3 ] = { 1, 111, 2 };
    int size = bondDim * physicalDim * physicalDim;
    int inc = 1;
    double expVal[ 2 ];
    for ( int A = 0; A < 2; ++A )
    {
        int B = ( A + 1 ) % 2;

        theta = gamma[ A ];
        prod( theta.elem, l[ A ].elem, bondDim, physicalDim, true );
        prod( theta.elem, l[ B ].elem, bondDim, physicalDim, false );

        theta.addLabel( tensorBondLabel );
        OpTheta = theta * Op;

        OpTheta.reshape( order, 2 );
        expVal[ A ] = ddot( &size, OpTheta.elem, &inc, theta.elem, &inc );
    }

    return ( expVal[ 0 ] + expVal[ 1 ] ) / 2.0;
}


double entanglementEtropy( const SyTensor_t & l )
{
    int bondDim = l.getElemNum();
    const double * dataPtr = l.elem;

    double sEntangle = 0.0;
    for ( int i = 0; i < bondDim; ++i )
        sEntangle -= ( dataPtr[ i ] * dataPtr[ i ] ) * log( dataPtr[ i ] * dataPtr[ i ] );

    return sEntangle;
}

#endif
