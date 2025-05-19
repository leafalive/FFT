#include <iostream>
#include "fft.h"
#include "arrayfire.h"
using namespace af;
//using namespace std;


int main() {
    af::setDevice(0);
    af::info();
//
//  int N = 1000;
//  af::array t = af::seq(0, N-1) / (double)N;
//  af::array x = af::sin(2 * af::Pi * 5 * t) + 0.5 * af::randn(N, f64);
//
//  af::array b = af::constant(0.2, 5, f64);


    std::vector<float> h_x = {0,0,0,0,0,1,1,1,1,1};
    af::array x(10, h_x.data());

    std::vector<float> h_b = {0.2f, 0.2f, 0.2f, 0.2f, 0.2f};
    std::vector<float> h_a = {1.0f};
    af::array b(5, h_b.data());
    af::array a(1, h_a.data());

    af::array y = filtfilt::filtfilt(b, a, x);
    af_print(y);

    return 0;
}
