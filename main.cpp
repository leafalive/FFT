#include <iostream>
#include "fft.h"

using namespace std;
int main() {
    Eigen::VectorXd a = Eigen::VectorXd::LinSpaced(3,1,3);
    std::cout << a << std::endl;
    Eigen::MatrixXcd A(3,3);
    A << 2,1,2,  3,2,1,  1,2,3;
    Eigen::MatrixXcd u(3,1);
    u << 1 , 0 , 1;
    Eigen::MatrixXcd v(2,1);
    v << 2, 7;
    Eigen::MatrixXcd B = FFT::fft(u,3,1);
    std::cout << FFT::ifft(B,3,1) << std::endl;
    std::cout << u << std::endl;
    std::cout << v << std::endl;
    std::cout << FFT::conv(u,v) << std::endl;
    return 0;
}
/*
 * input
(2,0) (1,0) (2,0)
(3,0) (2,0) (1,0)
(1,0) (2,0) (3,0)
 1 0 1
 2 7
 * output
(6,0)        (5,0)        (6,0)
(0,-1.73205)       (-1,0)  (0,1.73205)
(0,1.73205)       (-1,0) (0,-1.73205)
 2 7 2 7
*/
