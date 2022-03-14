#ifndef FRINGE_ULONG_MASK_H
#define FRINGE_ULONG_MASK_H

#include <iostream>
#include <cassert>


//Check code assumptions
//For C98: Unsigned Long was 4 bytes
//For C11: Unsigned int is 4 bytes
//Compilation assumes C++11 for now 
inline void checkLongSetting()
{

//    std::cout << "Size of unsigned int  : " << sizeof(unsigned int) << " bytes \n";
//    std::cout << "Size of unsigned long : " << sizeof(unsigned long) << " bytes \n";
    
    static_assert( sizeof(unsigned int) == 4, "Code assumes that size of unsigned int is 32 bits");

}




/*struct Ulongmask
 *This structure is meant to capture neighborhood mask 
 * for a (2*Ny+1) x (2*Nx+1) sized window in a unsigned long array
 * Each pixel is assigned to a bit and the bit is set/unset
 * as needed to represent the local neighborhood.
 * The index runs from [-Ny,Ny] and [-Nx,Nx].
 * Ny and Nx represent half window sizes*/

struct Ulongmask
{
    int Ny;    //Half window size in height direction
    int Nx;    //Half window size in width direction

    //Function to set a particular bit for a neighborhood pixel
    void setbit(unsigned int* arr,  int ii, int jj, bool flag) const;

    //Function to get a particular bit for a neighborhood pixel
    bool getbit(unsigned int* arr, int ii, int jj) const;

    //Print local neighborhood matrix
    void print(unsigned int* arr) const;
};

//Constant unsigned long for one
const unsigned int U1=1;

/*Set a specific bit but with no bounds checking
 * ii in [-Ny,Ny]
 * jj in [-Nx,Nx]
 * arr is of size atleast floor((2*Ny+1)*(2*Nx+1)/32) */
inline void Ulongmask::setbit(unsigned int* arr, int ii, int jj, bool flag) const
{
    //Flattened index into 2D array
    int flat = (ii+Ny)*(2*Nx+1) + jj+Nx;

    //Byte number in arr to modify
    int num = flat / 32;

    //Bit number of specific byte to modify
    int bit = flat % 32;

    if(flag)  //Set the bit to 1
    {
        arr[num] |= (((unsigned int)1)<<bit);
    }
    else      //Set the bit to 0
    {
        arr[num] &= (~(((unsigned int)(1))<<bit));
    }
}


/*Get a specific bit but with no bounds checking
 * ii in [-Ny,Ny]
 * jj in [-Nx,Nx]
 * arr is of size atleast floor((2*Ny+1)*(2*Nx+1)/32) */
inline bool Ulongmask::getbit(unsigned int* arr, int ii, int jj) const
{
    //Flattened index into 2D array
    int flat = (ii+Ny)*(2*Nx+1) + jj+Nx;

    //Byte number in arr to modify
    int num = flat / 32;

    //Bit number of specific byte to modify
    int bit = flat % 32;

    //Return the bit value
    return  (U1 == ((arr[num] >> bit) & U1)); 
}

//Print mask for debugging
inline void Ulongmask::print(unsigned int* arr) const
{
    //For each row
    for(int ii=-Ny; ii <= Ny; ii++)
    {
        //For each col
        for(int jj=-Nx;jj<=Nx;jj++)
        {
            std::cout << (int) getbit(arr, ii, jj) <<" ";
        }
        std::cout << "\n";
    }
}


//Sum up all the set bits
//http://stackoverflow.com/questions/109023/how-to-count-the-number-of-set-bits-in-a-32-bit-integer
inline int numberOfSetBits(unsigned int i)
{
    i = i - ((i>>1) & 0x55555555);
    i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
    return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
}


#endif
