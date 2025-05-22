#include <iostream>
#include "fft.h"
#include "arrayfire.h"
using namespace af;
//using namespace std;


int main() {
    af::setDevice(0);
    af::info();

    std::vector<double> h_x;
    for(int i = 0;i < 160;i++) {
      h_x.push_back(cos(af::Pi / 4 * i));
    }
    af::array x(160,h_x.data());
    std::vector<double> h_b = {1, 3, 3, 1};
    std::vector<double> h_a = {3, 0, 1, 0};
    for(int i =0 ;i < 4;i++) {
        h_b[i] /= double (6);
        h_a[i] /= double (3);
    }

    af::array b(4,h_b.data());
    af::array a(4,h_a.data());
    af::array y = filtfilt::filtfilt(b, a, x);
    af_print(y);

    return 0;
}
