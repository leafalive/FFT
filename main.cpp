#include <iostream>
#include "fft.h"
#include "arrayfire.h"
using namespace af;
//using namespace std;


int main() {
    af::setDevice(0);
    af::info();

    std::vector<float> h_x = {0,0,0,0,0,1,1,1,1,1};
    af::array x(10, h_x.data());

    std::vector<float> h_b = {0.2f, 0.2f, 0.2f, 0.2f, 0.2f};
    std::vector<float> h_a = {1.0f};
    af::array b(5, h_b.data());
    af::array a(1, h_a.data());


//    int n = 100;
//    std::vector<float> h_x(n);
//    for (int i=0; i<n; ++i)
//      h_x[i] = sin(2 * af::Pi * 0.1 * i);
//
//    af::array x(n, h_x.data());
//
//    std::vector<float> h_b = {0.1f, 0.2f, 0.4f, 0.2f, 0.1f};
//    std::vector<float> h_a = {1.0f, -0.5f, 0.2f};
//    af::array b(5, h_b.data());
//    af::array a(3, h_a.data());
    af::array y = filtfilt::filtfilt(b, a, x);
    af_print(y);

    return 0;
}
