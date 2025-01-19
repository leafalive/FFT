#include <iostream>
#include "fft.h"


int main() {
    //Eigen::MatrixXcd
    Eigen::VectorXd a = Eigen::VectorXd::LinSpaced(3,1,3);
    std::cout << a << std::endl;
    Eigen::MatrixXcd A(3,3);
    A << 2,1,2,  3,2,1,  1,2,3;
    std::cout << A << std::endl;
    std::cout << FFT::fft(A,-1,1);
//    std::cout << "Hello, World!" << std::endl;
//    for(int i = 0;i < 10;i++) {
//        std::cout << i << " ";
//    }
    return 0;
}
/*
 * input
(2,0) (1,0) (2,0)
(3,0) (2,0) (1,0)
(1,0) (2,0) (3,0)
 * output
(6,0)        (5,0)        (6,0)
(0,-1.73205)       (-1,0)  (0,1.73205)
(0,1.73205)       (-1,0) (0,-1.73205)
*/
