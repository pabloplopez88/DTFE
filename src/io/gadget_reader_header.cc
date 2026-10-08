/*
 *  Copyright (c) 2013       Marius Cautun
 *
 *                           Institute for Computational Cosmology
 *                           Durham University, Durham, UK
 *
 *
 *  Modified in 2026 by Pablo Lopez (https://github.com/pabloplopez88/DTFE); see the git history for the changes.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
 
 
#include <cstring>
#include <vector>

#define SWAP_HEADER_ENDIANNESS(x1,x2,x3,x4) { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 ); x4.swapBytes();} }
#define SWAP_ENDIANNESS(x1,x2,x3)           { if( x1 ) {BYTESWAP( x2 ); BYTESWAP( x3 );} }
#define READ_DELIMETER \
    inputFile.seekg( offset, std::ios::cur ); \
    inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) ); \
    if ( swapEndian ) BYTESWAP( buffer1 )
#define DELIMETER_CONSISTANCY_CHECK(field) \
    inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) ); \
    if ( swapEndian ) BYTESWAP( buffer2 ); \
    if ( buffer1!=buffer2 ) \
        throwError( "The integers before and after the particle " field " data block in the GADGET file '" + fileName + "' did not match. The GADGET snapshot file is corrupt." )


// Header structure for reading Gadget snapshots
struct Gadget_header
{
    int      npart[6];
    double   mass[6];
    double   time;
    double   redshift;
    int      flag_sfr;
    int      flag_feedback;
    int      npartTotal[6];
    int      flag_cooling;
    int      num_files;
    double   BoxSize;
    double   Omega0;
    double   OmegaLambda;
    double   HubbleParam;
    char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */


    // return the file name for a Gadget snapshot saved in single or multiple files - note that the name must contain a '%i' or '%s' character
    std::string filename(std::string fileRoot, int const fileNumber, bool checkFileExists=true )
    {
        char buf[500];
        sprintf( buf, fileRoot.c_str(), fileNumber );
        std::string fileName( buf );
        if ( not bfs::exists(fileName) and checkFileExists )
            throwError( "The program could not open the input GADGET snapshot file/files: '" + fileName + "'. It cannot find the file/files." );
        return fileName;
    }

    // Function that prints the Gadget header.
    void print()
    {
        std::cout << "\nThe header of the Gadget file contains the following info:\n"
            << "npart[6]     =  " << npart[0] << "  " << npart[1] << "  " << npart[2] << "  " << npart[3] << "  " <<  npart[4] << "  " <<  npart[5] << "\n"
            << "mass[6]      =  " << mass[0] << "  " << mass[1] << "  " << mass[2] << "  " << mass[3] << "  " << mass[4] << "  " << mass[5] << "\n"
            << "time         =  " << time << "\n"
            << "redshift     =  " << redshift << "\n"
            << "flag_sfr     =  " << flag_sfr << "\n"
            << "flag_feedback=  " << flag_feedback << "\n"
            << "npartTotal[6]=  " << npartTotal[0] << "  " << npartTotal[1] << "  " << npartTotal[2] << "  " << npartTotal[3] << "  " << npartTotal[4] << "  " << npartTotal[5] << "  " << "\n"
            << "flag_cooling =  " << flag_cooling << "\n"
            << "num_files    =  " << num_files << "\n"
            << "BoxSize      =  " << BoxSize << "\n"
            << "Omega0       =  " << Omega0 << "\n"
            << "OmegaLambda  =  " << OmegaLambda << "\n"
            << "h            =  " << HubbleParam << "\n\n";
    }

    // Swap endianness
    void swapBytes()
    {
        ByteSwapArray( npart, 6 );
        ByteSwapArray( mass, 6 );
        BYTESWAP( time );
        BYTESWAP( redshift );
        BYTESWAP( flag_sfr );
        BYTESWAP( flag_feedback );
        ByteSwapArray( npartTotal, 6 );
        BYTESWAP( flag_cooling );
        BYTESWAP( num_files );
        BYTESWAP( BoxSize );
        BYTESWAP( Omega0 );
        BYTESWAP( OmegaLambda );
        BYTESWAP( HubbleParam );
    }

    /* Gadget-4 (unless compiled with the GADGET2_HEADER option) writes a different header to its binary snapshots (SnapFormat 1 and 2):
            long long npart[N];  long long npartTotal[N];  double mass[N];  double time, redshift, BoxSize;
            int num_files;  long long Ntrees, TotNtrees;
       with N the number of particle types (6 by default). Its size is 24*N+48 bytes (192 bytes for N=6).
       Returns true if 'blockSize' is the size of such a header and stores N in 'noTypes'. */
    static bool isGadget4HeaderSize(int const blockSize, int *noTypes)
    {
        if ( blockSize<=48 or (blockSize-48)%24!=0 or blockSize==256 ) return false;
        *noTypes = (blockSize-48) / 24;
        return (*noTypes)>=1 and (*noTypes)<=64;
    }

    // Returns true if 'blockSize' is the size of a header that can be read: Gadget-1/2 (256 bytes) or Gadget-4.
    static bool isKnownHeaderSize(int const blockSize)
    {
        int noTypes;
        return blockSize==256 or isGadget4HeaderSize( blockSize, &noTypes );
    }

    // Checks for the type of the Gadget file -> can detect Gadget file type 1 & 2, with a Gadget-1/2 or a Gadget-4 header. Returns true if it could identify the gadget file type.
    bool detectSnapshotType(int const bufferValue,
                            int *gadgetFileType,
                            bool *swapEndian)
    {
        int buffer1 = bufferValue;
        *swapEndian = false;

        if ( buffer1 == 8 )                         // gadget file format 2
            *gadgetFileType = 2;
        else if ( isKnownHeaderSize(buffer1) )      // gadget file format 1 (Gadget-1/2 or Gadget-4 header)
            *gadgetFileType = 1;
        else                                        // check for swapped endianness
        {
            BYTESWAP( buffer1 );
            *swapEndian = true;
            if ( buffer1 == 8 )                     // gadget file format 2
                *gadgetFileType = 2;
            else if ( isKnownHeaderSize(buffer1) )  // gadget file format 1
                *gadgetFileType = 1;
            else                                    // could not detect the file type
                return false;
        }
        return true;
    }

    /* Reads the header data block (including the integers before and after it) from a Gadget binary snapshot. The file must be positioned at the
       beginning of the header block (after the 16 bytes block label for Gadget file format 2). It first tries the Gadget-1/2 header (256 bytes),
       then the Gadget-4 header; if none of them matches it stops the program. The header is always stored in the Gadget-1/2 layout of this
       structure, with the native endianness. Returns 2 for a Gadget-1/2 header and 4 for a Gadget-4 header. */
    int readHeaderBlock(std::fstream &inputFile,
                        bool const swapEndian,
                        std::string const &fileName)
    {
        int buffer1, buffer2, noTypes, headerType;
        inputFile.read( reinterpret_cast<char *>(&buffer1), sizeof(buffer1) );
        if ( swapEndian ) BYTESWAP( buffer1 );

        if ( buffer1==256 )                                     // Gadget-1/2 header
        {
            inputFile.read( reinterpret_cast<char *>(this), 256 );
            if ( swapEndian ) this->swapBytes();
            headerType = 2;
        }
        else if ( isGadget4HeaderSize(buffer1, &noTypes) )      // Gadget-4 header
        {
            std::vector<char> raw( buffer1 );
            inputFile.read( &(raw[0]), buffer1 );
            this->setFromGadget4Header( &(raw[0]), noTypes, swapEndian, fileName );
            headerType = 4;
        }
        else
            throwError( "Could not read the header of the GADGET snapshot file '" + fileName + "'. Tried the Gadget-1/2 header (256 bytes) and the Gadget-4 header (24*N+48 bytes for N particle types), but the header block has a different size. The file is corrupt or is not a GADGET binary snapshot." );

        inputFile.read( reinterpret_cast<char *>(&buffer2), sizeof(buffer2) );
        if ( swapEndian ) BYTESWAP( buffer2 );
        if ( not inputFile or buffer1!=buffer2 )
            throwError( "The integers before and after the header of the GADGET snapshot file '" + fileName + "' do not match. The GADGET snapshot file is corrupt." );
        return headerType;
    }

    // Copies the values of a Gadget-4 header (stored in 'raw') to this structure.
    void setFromGadget4Header(char const *raw,
                              int const noTypes,
                              bool const swapEndian,
                              std::string const &fileName)
    {
        std::memset( this, 0, sizeof(*this) );
        size_t const N = noTypes;
        for (size_t i=0; i<N; ++i)
        {
            long long np    = readValue<long long>( raw + 8*i, swapEndian );
            long long npTot = readValue<long long>( raw + 8*(N+i), swapEndian );
            double    m     = readValue<double>( raw + 8*(2*N+i), swapEndian );
            if ( i>=6 )         // DTFE works with 6 particle types
            {
                if ( np!=0 ) throwError( "The GADGET snapshot file '" + fileName + "' has particles of type 6 or larger, but DTFE can only read the particle types 0 to 5." );
                continue;
            }
            if ( np>2147483647LL ) throwError( "The GADGET snapshot file '" + fileName + "' has more than 2^31 particles of one type in a single file. Please write the snapshot in several files." );
            npart[i] = int( np );
            npartTotal[i] = int( npTot & 0xffffffffLL );   // only used for information (as the low word in Gadget-2)
            mass[i] = m;
        }
        time      = readValue<double>( raw + 24*N,      swapEndian );
        redshift  = readValue<double>( raw + 24*N + 8,  swapEndian );
        BoxSize   = readValue<double>( raw + 24*N + 16, swapEndian );
        num_files = readValue<int>(    raw + 24*N + 24, swapEndian );
    }

    // Reads a value of type T from a byte array, swapping the endianness if needed.
    template <typename T>
    static T readValue(char const *p, bool const swapEndian)
    {
        T value;
        std::memcpy( &value, p, sizeof(T) );
        if ( swapEndian ) BYTESWAP( value );
        return value;
    }
};

