#include <armadillo>
#include <iostream>
#include "EigenLapack.hpp"

int main(int argc, char **argv)
{
   /* arma::mat A = arma::randu<arma::mat>(50,50);
    arma::mat B = A.t()*A;  // generate a symmetric matrix

    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, B);


    // for matrices with complex elements

    arma::cx_mat C = arma::randu<arma::cx_mat>(50,50);
    arma::cx_mat D = C.t()*C;

    arma::vec eigval2;
    arma::cx_mat eigvec2;

    eig_sym(eigval2, eigvec2, D);

    std::cout << "Eigen values of random 50 x 50\n";
    eigval2.raw_print();*/
    
    
    arma::mat A1(3,3);
    A1(0,0) = 3;
    A1(0,1) = -2;
    A1(0,2) = 4;
    A1(1,0) = -2;
    A1(1,1) = 8;
    A1(1,2) = 2;
    A1(2,0) = 4;
    A1(2,1) = 2;
    A1(2,2) = 3;

    arma::vec val1;
    arma::mat vec1;

    eig_sym(val1, vec1, A1);

    std::cout << "Testing known 3 x 3 \n";
    A1.raw_print();
    std::cout << "Eigen values: \n";
    val1.raw_print();

    std::cout << "Eigen vectors: \n";
    vec1.raw_print();


    arma::cx_mat Acx(3,3);
    Acx(0,0) = 3;
    Acx(0,1) = -2;
    Acx(0,2) = 4;
    Acx(1,0) = -2;
    Acx(1,1) = 8;
    Acx(1,2) = 2;
    Acx(2,0) = 4;
    Acx(2,1) = 2;
    Acx(2,2) = 3;

    std::cout << "Input array: \n";
    Acx.raw_print();

    {
        arma::cx_mat Ac1 = Acx;
        EVWorker ev(3);
        std::cout << "Order = " << ev.size << "\n";
        ev.largestEigen( Ac1.memptr(), true);

        std::cout << "Largest eigen value: " << ev.eigval[0] << "\n";
        {
            arma::cx_vec vec(ev.eigvec, 3, false);
            vec.raw_print();
        }
    }

    {
        arma::cx_mat Ac1 = Acx;

        EVWorker ev(3);
        std::cout << "Order  " << ev.size << "\n";
        ev.smallestEigen( Ac1.memptr(), true);

        std::cout << "Smallest eigen value: " << ev.eigval[0] << "\n";
        {
            arma::cx_vec vec(ev.eigvec, 3, false);
            vec.raw_print();
        }
    }

    {
        arma::cx_mat Apos(3,3, arma::fill::eye);
        Apos *= 2;
        Apos += Acx;

        arma::cx_mat Ainv = Apos.i();

        std::cout << "Matrix inverse: \n";
        Ainv.raw_print();

        EVWorker ev(3);
        arma::cx_mat Ac1 = Apos + 0.0;

        int stat = ev.positiveDefiniteInverse(Ac1.memptr());

        std::cout << "Lapack inverse: " << stat << "\n";
        Ac1.raw_print();

        arma::cx_mat test = Apos * Ac1;
        test.raw_print();

    }


}
