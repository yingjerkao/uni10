#include <iostream>

#include <string>

#include <fstream>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

#include "itebd.h"

namespace prog_opts = boost::program_options;

std::ostream * output = NULL;

int main( int argc, char * argv[] )
{
    int bondDim = 10;
    double deltaT = 0.05;
    int numIter = 100;

// ======================================================
    std::fstream wavefunction;
    std::string outputFile;
    std::string waveFile;

    prog_opts::options_description description( "iTEBD Usage", 100 );

    description.add_options()
        ( "help,h", "Display this help message" )
        ( "showIter", "Show information for each iteration" )
        ( "quiet", "In quiet mode" )
        ( "deltaT", prog_opts::value<double>( &deltaT )->default_value( 0.05 ), "The timestep per iteration" )
        ( "bondDim,d", prog_opts::value<int>( &bondDim )->default_value( 10 ), "The bond dimension of MPS" )
        ( "numIter,n", prog_opts::value<int>( &numIter )->default_value( 100 ), "Number of timesteps to perform" )
        ( "wavefunction,w", prog_opts::value<std::string>( &waveFile )->required(), "The file name to store wavefunction" )
        ( "output,o", prog_opts::value<std::string>( &outputFile )->implicit_value( "itebd.log" ), "The file name to store the output" );

// ======================================================
    const int physicalDim = 2;

    double transField = 0.0;
    double longiField = 0.0;
    double Jz = 0.0;

    description.add_options()
        ( "tf", prog_opts::value<double>( &transField )->default_value( 0.0 ), "Transversal field strength (for itf hamiltonian)" )
        ( "lf", prog_opts::value<double>( &longiField )->default_value( 0.0 ), "Longitudinal field strength (for itf hamiltonian)" );

// ======================================================
    prog_opts::variables_map vm;

    try
    {
        prog_opts::store( prog_opts::command_line_parser( argc, argv ).options( description ).run(), vm );

        if ( vm.count( "help" ) )
        {
            std::cerr << "usage: itebd [options] -w wavefunction\n";
            std::cerr << description << std::endl;
            return 1;
        }

        prog_opts::notify( vm );
    }
    catch ( prog_opts::required_option & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << "usage: itebd [options] -w wavefunction\n";
        std::cerr << description << std::endl;
        return 1;
    }
    catch( prog_opts::error & e )
    {
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr << "usage: itebd [options] -w wavefunction\n";
        std::cerr << description << std::endl;
        return 1;
    }

    if ( vm.count( "output" ) )
        output = new std::ofstream( outputFile.c_str(), std::ofstream::out );
    else
        output = &std::cout;

    boost::format showEntropyFmt( "%=25.16f" );

// ======================================================
    std::vector<Qnum_t> qnums;
    std::vector<Bond_t> bonds;
    int inc = 1;
    int size;
    double alpha;
    double * dptr = NULL;
    int tensorBondLabel[ 4 ];

// ======================================================
    std::vector<SyTensor_t> gamma;
    std::vector<SyTensor_t> l;

    wavefunction.open( waveFile.c_str(), std::ios_base::in | std::ios_base::binary );
    if ( wavefunction.is_open() )
    {
        char identify[ 9 ] = { 0 };
        wavefunction.read( identify, 8 );
        if ( std::string( "MY_ITEBD" ).compare( identify ) != 0 )
        {
            std::cerr << "Invailed wavefunction file format!\n";
            return 2;
        }

        int inPhysicalDim;
        wavefunction.read( (char *)&bondDim, sizeof( int ) );
        wavefunction.read( (char *)&inPhysicalDim, sizeof( int ) );
        if ( inPhysicalDim != physicalDim )
        {
            std::cerr << "Invailed wavefunction file format!\n";
            return 3;
        }

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
    else
    {
        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( physicalDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
        qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_COL, qnums ) );

        for ( int i = 0; i < 2; ++i )
        {
            gamma.push_back( SyTensor_t( bonds ) );
            gamma[ i ].randomize();
        }

        bonds.clear();
        bonds.push_back( Bond_t( BD_COL, qnums ) );

        for ( int i = 0; i < 2; ++i )
        {
            l.push_back( SyTensor_t( bonds ) );
            l[ i ].randomize();
        }
    }

    dptr = NULL;
    qnums.clear();
    bonds.clear();

// ======================================================
    qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
    qnums.resize( physicalDim, Qnum_t() ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
    qnums.resize( bondDim, Qnum_t() ); bonds.push_back( Bond_t( BD_COL, qnums ) );

    tensorBondLabel[ 0 ] = 11, tensorBondLabel[ 1 ] = -1, tensorBondLabel[ 2 ] = 0;
    SyTensor_t lbAla( bonds, tensorBondLabel, "lbAla" );

    qnums.resize( bondDim, Qnum_t() ); bonds[ 0 ] = Bond_t( BD_ROW, qnums );
    qnums.resize( physicalDim, Qnum_t() ); bonds[ 1 ] = Bond_t( BD_COL, qnums );
    qnums.resize( bondDim, Qnum_t() ); bonds[ 2 ] = Bond_t( BD_COL, qnums );

    tensorBondLabel[ 0 ] = 0, tensorBondLabel[ 1 ] = 1, tensorBondLabel[ 2 ] = 22;
    SyTensor_t Blb( bonds, tensorBondLabel, "Blb" );

    qnums.clear();
    bonds.clear();

// ======================================================
    double * hamiltonian = NULL;
    
    hamiltonian = tfIsing( transField, longiField );

    if ( !vm.count( "quiet" ) || vm.count( "output" ) )
    {
        *output << boost::format( "\n# timestep: %=10.6g\n# bond dimension: %=10d\n" ) % deltaT % bondDim;
        *output << boost::format( "\n# Ising model\n# Transversal field strength: %=10f\n# Longitudinal field strength: %=10f\n" )
            % transField % longiField;

        if ( vm.count( "showIter" ) )
            *output << boost::format( "\n#%=25s%=25s%=25s\n" ) % "energy" % "<Sz>" % "Entropy";
    }

    boost::format showIterFmt( " %=25.16f%=25.16f" );


// ======================================================
    for ( int i = 0; i < physicalDim; ++i )
        qnums.push_back( Qnum_t() );

    bonds.push_back( Bond_t( BD_ROW, qnums ) ); bonds.push_back( Bond_t( BD_ROW, qnums ) );
    bonds.push_back( Bond_t( BD_COL, qnums ) ); bonds.push_back( Bond_t( BD_COL, qnums ) );

    tensorBondLabel[ 0 ] = 111, tensorBondLabel[ 1 ] = 222, tensorBondLabel[ 2 ] = -1, tensorBondLabel[ 3 ] = 1;

    SyTensor_t U( bonds, tensorBondLabel, "U" );

    dptr = new double [ U.getElemNum() ];
    expm( dptr, physicalDim * physicalDim, -deltaT, hamiltonian );
    U.addRawElem( dptr );

    delete [] dptr;

    dptr = NULL;
    qnums.clear();
    bonds.clear();

// ======================================================
    size = bondDim * physicalDim;
    double * u = new double [ size * size ];
    double * s = new double [ size ];
    int lwork = -1;
    double * work = new double [ 1 ];
    int info;

    SyTensor_t theta;
    dgesvd( "O", "A", &size, &size, theta.elem, &size, s, NULL, &size, u, &size, work, &lwork, &info );

    if ( info != 0 )
        throw "DGESVD info != 0";
    else
    {
        lwork = work[ 0 ];
        delete [] work;
    }
    work = new double [ lwork ];

    int order[ 4 ] = { 11, 111, 222, 22 };
    int n = 0;
    while ( n < numIter )
    {
        for ( int A = 0; A < 2; ++A )
        {
            int B = ( A + 1 ) % 2;

            lbAla.addRawElem( gamma[ A ].elem );
            prod( lbAla.elem, l[ A ].elem, bondDim, physicalDim, true );
            prod( lbAla.elem, l[ B ].elem, bondDim, physicalDim, false );

            Blb.addRawElem( gamma[ B ].elem );
            prod( Blb.elem, l[ B ].elem, bondDim, physicalDim, true );

            theta = lbAla * Blb;
            theta *= U;

            theta.reshape( order, 2 );

            dgesvd( "O", "A", &size, &size, theta.elem, &size, s, NULL, &size, u, &size, work, &lwork, &info );

            alpha = 1.0 / dnrm2( &bondDim, s, &inc );
            dcopy( &bondDim, s, &inc, l[ A ].elem, &inc );
            dscal( &bondDim, &alpha, l[ A ].elem, &inc );

            truncate( gamma[ A ].elem, u, size, size, bondDim, false );
            truncate( gamma[ B ].elem, theta.elem, size, size, bondDim, true );

            dcopy( &bondDim, l[ B ].elem, &inc, s, &inc );
            invData( s, bondDim );
            prod( gamma[ A ].elem, s, bondDim, physicalDim, false );
            prod( gamma[ B ].elem, s, bondDim, physicalDim, true );
        }

        if ( vm.count( "showIter" ) && ( !vm.count( "quiet" ) || vm.count( "output" ) ) )
        {
            *output << showIterFmt % energy( gamma, l, hamiltonian, physicalDim ) % expectionValue( gamma, l, Sz, physicalDim );
            *output << std::endl;
        }

        ++n;
    }

// ======================================================
    if ( !vm.count( "quiet" ) || vm.count( "output" ) )
    {
        *output << boost::format( "\n\n# energy: %=25.16f\n# <Sz>: %=25.16f\n" )
            % energy( gamma, l, hamiltonian, physicalDim )
            % expectionValue( gamma, l, Sz, physicalDim );

        *output << "# Entropy\n";
        for ( int i = 0; i < 2; ++i )
            *output << showEntropyFmt % entanglementEtropy( l[ i ] ) << std::endl;
    }

// ======================================================
    wavefunction.open( waveFile.c_str(), std::ios_base::out | std::ios_base::binary );
    if ( wavefunction.is_open() )
    {
        wavefunction.write( "MY_ITEBD", 8 );
        wavefunction.write( (char *)&bondDim, sizeof( int ) );
        wavefunction.write( (char *)&physicalDim, sizeof( int ) );

        for ( int i = 0; i < 2; ++i )
        {
            wavefunction.write( (char *)( gamma[ i ].elem ), gamma[ i ].getElemNum() * sizeof( double ) );
            wavefunction.write( (char *)( l[ i ].elem ), l[ i ].getElemNum() * sizeof( double ) );
        }

        wavefunction.close();
    }

// ======================================================
    if ( vm.count( "output" ) )
        delete output;

    delete [] u;
    delete [] s;
    delete [] work;
    delete [] hamiltonian;

    return 0;
}
