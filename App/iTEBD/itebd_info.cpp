#include <iostream>

#include <string>

#include <fstream>

#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "itebd.h"


bool compare( double a, double b )
{
    return a > b;
}

namespace prog_opts = boost::program_options;

int main( int argc, char * argv[] )
{
    std::ifstream wavefunction;
    std::string waveFile;

    prog_opts::options_description description( "iTEBD_info Usage", 100 );

    description.add_options()
        ( "help,h", "Display this help message" )
        ( "entropy,e", "Show entanglement entropy" )
        ( "density,d", "Show eigenvalues of density matrix" )
        ( "corre,c", "Show correlation length" )
        ( "wavefunction", prog_opts::value<std::string>( &waveFile )->required(), "The wavefunction file" );

    prog_opts::positional_options_description positionalOptions;
    positionalOptions.add( "wavefunction", 1 ); 

    prog_opts::variables_map vm;

    try
    {
        prog_opts::store( prog_opts::command_line_parser( argc, argv ).options( description ).positional( positionalOptions ).run(), vm );

        if ( vm.count( "help" ) )
        {
            std::cerr << "usage: itebd_info [options] wavefunction\n";
            std::cerr << description << std::endl;
            return 1;
        }

        prog_opts::notify( vm );
    }
    catch ( prog_opts::required_option & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << "usage: itebd_info [options] wavefunction\n";
        std::cerr << description << std::endl;
        return 1;
    }
    catch( prog_opts::error & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << "usage: itebd_info [options] wavefunction\n";
        std::cerr << description << std::endl;
        return 1;
    }

    int bondDim;
    int physicalDim;
    double * dptr = NULL;
    std::vector<SyTensor_t> gamma;
    std::vector<SyTensor_t> l;

    std::vector<Qnum_t> qnums;
    std::vector<Bond_t> bonds;

    wavefunction.exceptions( std::ifstream::failbit | std::ifstream::badbit );

    try
    {
        wavefunction.open( waveFile.c_str(), std::ios_base::in | std::ios_base::binary );

        char identify[ 9 ] = { 0 };
        wavefunction.read( identify, 8 );
        if ( std::string( "MY_ITEBD" ).compare( identify ) != 0 )
        {
            std::cerr << "Invailed wavefunction file format!\n";
            return 2;
        }

        wavefunction.read( (char *)&bondDim, sizeof( int ) );
        wavefunction.read( (char *)&physicalDim, sizeof( int ) );

        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( physicalDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_COL, qnums ) );

        for ( int i = 0; i < 2; ++i )
            gamma.push_back( SyTensor_t( bonds ) );

        bonds.clear();
        bonds.push_back( Bond_t( BD_COL, qnums ) );

        for ( int i = 0; i < 2; ++i )
            l.push_back( SyTensor_t( bonds ) );


        for ( int i = 0; i < 2; ++i )
        {
            dptr = new double [ bondDim * physicalDim * bondDim ];
            wavefunction.read( (char *)dptr, bondDim * physicalDim * bondDim * sizeof( double ) );
            gamma[ i ].addRawElem( dptr );

            delete [] dptr;

            dptr = new double [ bondDim ];
            wavefunction.read( (char *)dptr, bondDim * sizeof( double ) );
            l[ i ].addRawElem( dptr );

            delete [] dptr;
        }

        wavefunction.close();
    }
    catch ( std::ifstream::failure & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        wavefunction.close();
        return 1;
    }

    qnums.clear();
    bonds.clear();

    std::cout << boost::format( "# physical dimension: %=5d\n"
            "# bond dimension: %=5d\n" )
        % physicalDim % bondDim;

    if ( vm.count( "entropy" ) )
    {
        boost::format showEntropyFmt( "%-25.16g" );
        std::cout << "# Entropy:\n";
        for ( int i = 0; i < 2; ++i )
            std::cout << showEntropyFmt % entanglementEtropy( l[ i ] ) << std::endl;
        std::cout << std::endl;
    }
    else if ( vm.count( "density" ) )
    {
        boost::format showEntropyFmt( "%-25.16g" );
        boost::format showDensityFmt( "  %-25.16g" );
        std::cout << "# Density matrix:\n";
        for ( int i = 0; i < 2; ++i )
        {
            std::cout << "# Site: " << i + 1 << std::endl;
            std::cout << "Entropy:  " << showEntropyFmt % entanglementEtropy( l[ i ] ) << std::endl;
            const double * dataPtr = l[ i ].elem;
            for ( int j = 0; j < bondDim; ++j )
                std::cout << j << showDensityFmt % ( dataPtr[ j ] * dataPtr[ j ] ) << std::endl;
            std::cout << std::endl;
        }
    }
    else if ( vm.count( "corre" ) )
    {
        int tensorBondLabel[ 4 ];
        int order[ 4 ] = { 11, 111, 222, 22 };

        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( physicalDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_COL, qnums ) );

        tensorBondLabel[ 0 ] = 11, tensorBondLabel[ 1 ] = -1, tensorBondLabel[ 2 ] = 0;
        SyTensor_t T1( bonds, tensorBondLabel, "T1" );

        qnums.resize( bondDim, Qnum_t() ); bonds[ 0 ] = Bond_t( BD_ROW, qnums );
        qnums.resize( physicalDim, Qnum_t() ); bonds[ 1 ] = Bond_t( BD_COL, qnums );
        qnums.resize( bondDim, Qnum_t() ); bonds[ 2 ] = Bond_t( BD_COL, qnums );

        tensorBondLabel[ 0 ] = 0, tensorBondLabel[ 1 ] = 1, tensorBondLabel[ 2 ] = 22;
        SyTensor_t T2( bonds, tensorBondLabel, "T2" );

        qnums.clear();
        bonds.clear();

        T1.addRawElem( gamma[ 0 ].elem );
        prod( T1.elem, l[ 0 ].elem, bondDim, physicalDim, true );

        T2.addRawElem( gamma[ 1 ].elem );
        prod( T2.elem, l[ 1 ].elem, bondDim, physicalDim, true );

        T1 *= T2;
        T2 = T1;

        tensorBondLabel[ 0 ] = 111, tensorBondLabel[ 1 ] = -1, tensorBondLabel[ 2 ] = 1, tensorBondLabel[ 3 ] = 222;
        T2.addLabel( tensorBondLabel ); 

        T1 *= T2;
        T1.reshape( order, 2 );

        int size = bondDim * bondDim;
        double * WR = new double [ size ];
        double * WI = new double [ size ];

        int SDIM;
        int lwork = -1;
        double * work = new double [ 1 ];
        int info;

        dgees( "N", "N", NULL, &size, T1.elem, &size, &SDIM, WR, WI, NULL, &size, work, &lwork, NULL, &info );

        if ( info != 0 )
            throw "DGEES info != 0";
        else
        {
            lwork = work[ 0 ];
            delete [] work;
        }

        work = new double [ lwork ];
        dgees( "N", "N", NULL, &size, T1.elem, &size, &SDIM, WR, WI, NULL, &size, work, &lwork, NULL, &info );

        std::vector<double> eigvals;
        for ( int i = 0; i < size; ++i )
        {
            if ( WI[ i ] == 0.0 )
                eigvals.push_back( abs( WR[ i ] ) );
        }

        delete [] WR;
        delete [] WI;
        delete [] work;

        std::sort( eigvals.begin(), eigvals.end(), compare );
/*
        for ( int i = 0; i < eigvals.size(); ++i )
            std::cout << eigvals[ i ] << std::endl;
*/
        std::cout << boost::format( "%=25.16g" ) % ( -1.0 / log( eigvals[ 1 ] ) ) << std::endl;
    }

    return 0;
}
