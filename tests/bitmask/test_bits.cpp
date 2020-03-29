#include <iostream>
#include "ulongmask.hpp"



void testOneLong()
{
    const Ulongmask msk = {2,2};

    unsigned int test = 0;
    for(int ii=-2; ii<=2; ii++)
    {
        for(int jj=-2; jj<=2; jj++)
        {
            test = 0;
            msk.setbit(&test, ii, jj, true);

            int flat = (ii+2)*(2*2+1) + jj+2;
            int num = flat / 32;
            unsigned int bit = flat % 32;

            std::cout << " flat = " <<  flat << "   "
                      <<  (((unsigned int)1)<<bit) << "  "
                      << test << "\n";
        }
    }
}



void testFiveLong()
{
    const Ulongmask msk = {7,7};

    unsigned int test[8] = {0,0,0,0,0,0,0,0};

    for(int ii=-7; ii<=7; ii++)
    {
        for(int jj=-7; jj<=7; jj++)
        {
            for(int kk=0; kk<8; kk++) test[kk] = 0;
            msk.setbit(test, ii, jj, true);

            int flat = (ii+7)*(2*7+1) + jj+7;
            int num = flat / 32;
            unsigned int bit = flat % 32;

            std::cout << " flat = " <<  flat << "   " << test[num] << "  " 
                      <<  (((unsigned int)1)<<bit) << "  "
                      << test[0] <<" " <<test[1] << " " 
                      << test[2] <<" " << test[3] << " " 
                      << test[4] <<" " << test[5] << " "
                      << test[6] <<" " << test[7] << " "
                      << msk.getbit(test,ii,jj) << "\n";
        }
    }
}





int main(int argc, char **argv)
{

    std::cout << "Testing single uint32 field \n";
    testOneLong();

    std::cout << "Testing with five uint32 fields \n";
    testFiveLong();

}

