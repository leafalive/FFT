//
// Created by Yanran LI on 2025/1/19.
//
#include "fft.h"
using namespace af;
namespace FFT {

    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A, int n,int dim ) {
        if(n == 0 || A.size() == 0) {
            return Eigen::MatrixXcd(0,0);
        }
        const long long Rows = A.rows();
        const long long Cols = A.cols();
        Eigen::MatrixXcd Temp = A;
        Eigen::FFT<double> F;

        if(n > 0) {
            if(dim == 1) {
                Temp.conservativeResize(n, Cols);
                if(n > Rows) { // 补0
                    for(long long k = Rows;k < n;k++) {
                        Temp.row(k).setZero();
                    }
                }
            }
            else if(dim == 2) {
                Temp.conservativeResize(Rows, n);
                if(n > Cols) { // 补0
                    for(long long k = Cols;k < n;k++) {
                        Temp.col(k).setZero();
                    }
                }
            }
            else {}
        }

        Eigen::MatrixXcd Result(Temp.rows(), Temp.cols());  // 结果
        if(dim == 1) { // 各列计算
            for(int k = 0;k < Temp.cols();k++) {
                Eigen::VectorXcd tmp(Temp.cols());
                F.fwd(tmp, Temp.col(k));
                Result.col(k) = tmp;
            }
        }
        else if(dim == 2) { // 各行计算
            for(int k = 0;k < Temp.rows();k++) {
                Eigen::VectorXcd tmp(Temp.rows());
                F.fwd(tmp, Temp.row(k));
                Result.row(k) = tmp;
            }
        }
        else {}
        return Result;
    }

    Eigen::MatrixXcd ifft(Eigen::MatrixXcd &A, int n,int dim ) {
      if(n == 0 || A.size() == 0) {
        return Eigen::MatrixXcd(0,0);
      }
      const long long Rows = A.rows();
      const long long Cols = A.cols();
      Eigen::MatrixXcd Temp = A;
      Eigen::FFT<double> F;

      if(n > 0) {
        if(dim == 1) {
          Temp.conservativeResize(n, Cols);
          if(n > Rows) { // 补0
            for(long long k = Rows;k < n;k++) {
              Temp.row(k).setZero();
            }
          }
        }
        else if(dim == 2) {
          Temp.conservativeResize(Rows, n);
          if(n > Cols) { // 补0
            for(long long k = Cols;k < n;k++) {
              Temp.col(k).setZero();
            }
          }
        }
        else {}
      }

      Eigen::MatrixXcd Result(Temp.rows(), Temp.cols());  // 结果
      if(dim == 1) { // 各列计算
        for(int k = 0;k < Temp.cols();k++) {
          Eigen::VectorXcd tmp(Temp.cols());
          F.inv(tmp, Temp.col(k));
          Result.col(k) = tmp;
        }
      }
      else if(dim == 2) { // 各行计算
        for(int k = 0;k < Temp.rows();k++) {
          Eigen::VectorXcd tmp(Temp.rows());
          F.inv(tmp, Temp.row(k));
          Result.row(k) = tmp;
        }
      }
      else {}
      return Result;
    }

    Eigen::MatrixXcd conv(Eigen::MatrixXcd &A,Eigen::MatrixXcd &B) {
        int result_length = A.size() + B.size() - 1;

        Eigen::MatrixXcd u_fft = fft(A,result_length,1);
        Eigen::MatrixXcd v_fft = fft(B,result_length,1);

        Eigen::MatrixXcd temp = u_fft.cwiseProduct(v_fft);

        Eigen::MatrixXcd result = ifft(temp,result_length, 1);

        return result;
    }
}

namespace filtfilt {

  af::array filtfilt(const af::array &b, const af::array &a, const af::array &x) {

    int len_a = a.dims(0);
    int len_b = b.dims(0);
    int filter_order = max(len_a, len_b) - 1;
    int n_pad = 3 * filter_order;


    int n = x.dims(0);


    if (n <= n_pad) {
      return x;
    }


    af::array x_pre = af::flip(x(af::seq(n_pad)), 0);

    af::array x_post = af::flip(x(af::seq(n - n_pad, n - 1)), 0);

    af::array x_padded = af::join(0, x_pre, x, x_post);

    // 前向滤波
    auto forward_filter = [](const af::array &x, const af::array &b,
                           const af::array &a) {
        int len_x = x.dims(0);
        int len_b = b.dims(0);
        int len_a = a.dims(0);


        std::vector<float> h_a(a.elements());
        a.host(h_a.data());
        std::vector<float> h_b(b.elements());
        b.host(h_b.data());


        std::vector<float> h_x(x.elements());
        x.host(h_x.data());

        std::vector<float> h_y(len_x, 0.0f);

        for (int i = 0; i < len_x; ++i) {
            float sum_b = 0.0f;
            for (int k = 0; k < len_b; ++k) {
                if (i - k >= 0) {
                    sum_b += h_b[k] * h_x[i - k];
                }
            }

            float sum_a = 0.0f;
            for (int k = 1; k < len_a; ++k) {
                if (i - k >= 0) {
                    sum_a += h_a[k] * h_y[i - k];
                }
            }

            h_y[i] = (sum_b - sum_a) / h_a[0];
        }

        return af::array(len_x, h_y.data());
    };

    af::array y_forward = forward_filter(x_padded, b, a);

    // 反转
    af::array y_reversed = af::flip(y_forward, 0);


    af::array y_reversed_filtered = forward_filter(y_reversed, b, a);

    // 反转
    af::array y_filtfilt = af::flip(y_reversed_filtered, 0);
    y_filtfilt = y_filtfilt(af::seq(n_pad, n_pad + n - 1));

    return y_filtfilt;
}
}
