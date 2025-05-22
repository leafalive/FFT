//
// Created by Yanran LI on 2025/1/19.
//
#include "fft.h"
using namespace af;
//namespace FFT {
//
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A, int n,int dim ) {
//        if(n == 0 || A.size() == 0) {
//            return Eigen::MatrixXcd(0,0);
//        }
//        const long long Rows = A.rows();
//        const long long Cols = A.cols();
//        Eigen::MatrixXcd Temp = A;
//        Eigen::FFT<double> F;
//
//        if(n > 0) {
//            if(dim == 1) {
//                Temp.conservativeResize(n, Cols);
//                if(n > Rows) { // 补0
//                    for(long long k = Rows;k < n;k++) {
//                        Temp.row(k).setZero();
//                    }
//                }
//            }
//            else if(dim == 2) {
//                Temp.conservativeResize(Rows, n);
//                if(n > Cols) { // 补0
//                    for(long long k = Cols;k < n;k++) {
//                        Temp.col(k).setZero();
//                    }
//                }
//            }
//            else {}
//        }
//
//        Eigen::MatrixXcd Result(Temp.rows(), Temp.cols());  // 结果
//        if(dim == 1) { // 各列计算
//            for(int k = 0;k < Temp.cols();k++) {
//                Eigen::VectorXcd tmp(Temp.cols());
//                F.fwd(tmp, Temp.col(k));
//                Result.col(k) = tmp;
//            }
//        }
//        else if(dim == 2) { // 各行计算
//            for(int k = 0;k < Temp.rows();k++) {
//                Eigen::VectorXcd tmp(Temp.rows());
//                F.fwd(tmp, Temp.row(k));
//                Result.row(k) = tmp;
//            }
//        }
//        else {}
//        return Result;
//    }
//
//    Eigen::MatrixXcd ifft(Eigen::MatrixXcd &A, int n,int dim ) {
//      if(n == 0 || A.size() == 0) {
//        return Eigen::MatrixXcd(0,0);
//      }
//      const long long Rows = A.rows();
//      const long long Cols = A.cols();
//      Eigen::MatrixXcd Temp = A;
//      Eigen::FFT<double> F;
//
//      if(n > 0) {
//        if(dim == 1) {
//          Temp.conservativeResize(n, Cols);
//          if(n > Rows) { // 补0
//            for(long long k = Rows;k < n;k++) {
//              Temp.row(k).setZero();
//            }
//          }
//        }
//        else if(dim == 2) {
//          Temp.conservativeResize(Rows, n);
//          if(n > Cols) { // 补0
//            for(long long k = Cols;k < n;k++) {
//              Temp.col(k).setZero();
//            }
//          }
//        }
//        else {}
//      }
//
//      Eigen::MatrixXcd Result(Temp.rows(), Temp.cols());  // 结果
//      if(dim == 1) { // 各列计算
//        for(int k = 0;k < Temp.cols();k++) {
//          Eigen::VectorXcd tmp(Temp.cols());
//          F.inv(tmp, Temp.col(k));
//          Result.col(k) = tmp;
//        }
//      }
//      else if(dim == 2) { // 各行计算
//        for(int k = 0;k < Temp.rows();k++) {
//          Eigen::VectorXcd tmp(Temp.rows());
//          F.inv(tmp, Temp.row(k));
//          Result.row(k) = tmp;
//        }
//      }
//      else {}
//      return Result;
//    }
//
//    Eigen::MatrixXcd conv(Eigen::MatrixXcd &A,Eigen::MatrixXcd &B) {
//        int result_length = A.size() + B.size() - 1;
//
//        Eigen::MatrixXcd u_fft = fft(A,result_length,1);
//        Eigen::MatrixXcd v_fft = fft(B,result_length,1);
//
//        Eigen::MatrixXcd temp = u_fft.cwiseProduct(v_fft);
//
//        Eigen::MatrixXcd result = ifft(temp,result_length, 1);
//
//        return result;
//    }
//}

namespace filtfilt {

af::array filtfilt(const af::array& b, const af::array& a, const af::array& x) {
  // 提取滤波器系数
  std::vector<double> b_vec(b.elements());
  std::vector<double> a_vec(a.elements());
  b.host(b_vec.data());
  a.host(a_vec.data());

  // 确保a[0]为1
  if (std::abs(a_vec[0] - 1.0f) > 1e-6) {
    float a0 = a_vec[0];
    for (auto& coef : b_vec) coef /= a0;
    for (auto& coef : a_vec) coef /= a0;
  }

  int n_b = b_vec.size();
  int n_a = a_vec.size();
  int maxNza = max(n_a, n_b);

  // 数据长度和边缘计算
  int nx = x.dims(0);
  int nfilt = max(n_a, n_b);
  int nfact = 3 * (nfilt - 1); // 边缘长度
  if (nfact > nx) nfact = nx - 1;

  // 提取并扩展数据
  std::vector<double> x_vec(nx);
  x.host(x_vec.data());

  std::vector<double> x_ext(nx + 2 * nfact);

  // 填充扩展数据
  // 左边界
  for (int i = 0; i < nfact; i++) {
    x_ext[i] = 2 * x_vec[0] - x_vec[nfact - i];
  }
  // 主体
  for (int i = 0; i < nx; i++) {
    x_ext[nfact + i] = x_vec[i];
  }
  // 右边界
  for (int i = 0; i < nfact; i++) {
    x_ext[nfact + nx + i] = 2 * x_vec[nx - 1] - x_vec[nx - 2 - i];
  }

  // 计算滤波器初始状态
  std::vector<double> zi(n_a - 1, 0.0f);
  if (n_a > 1) {
    // 解决方程 A*zi = B*x[0]
    // 这里简化为直接设置zi，实际上应该解方程
    for (int i = 0; i < n_a - 1; i++) {
      zi[i] = x_ext[0] * b_vec[i+1];
      for (int j = 0; j < i + 1; j++) {
        zi[i] -= a_vec[j+1] * zi[i-j];
      }
    }
  }

  // 正向滤波
  std::vector<double> y(nx + 2 * nfact);
  y[0] = b_vec[0] * x_ext[0];

  for (int i = 1; i < n_b; i++) {
    y[0] += b_vec[i] * x_ext[0];
  }

  for (int i = 1; i < n_a; i++) {
    y[0] -= a_vec[i] * zi[i-1];
  }

  for (int i = 1; i < nx + 2 * nfact; i++) {
    y[i] = b_vec[0] * x_ext[i];

    for (int j = 1; j < n_b; j++) {
      if (i - j >= 0)
        y[i] += b_vec[j] * x_ext[i - j];
    }

    for (int j = 1; j < n_a; j++) {
      if (i - j >= 0)
        y[i] -= a_vec[j] * y[i - j];
    }
  }

  // 反向滤波
  // 反转信号
  std::vector<double> y_rev(nx + 2 * nfact);
  for (int i = 0; i < nx + 2 * nfact; i++) {
    y_rev[i] = y[nx + 2 * nfact - 1 - i];
  }

  // 计算反向滤波的初始状态
  std::vector<double> zi_rev(n_a - 1, 0.0f);
  if (n_a > 1) {
    for (int i = 0; i < n_a - 1; i++) {
      zi_rev[i] = y_rev[0] * b_vec[i+1];
      for (int j = 0; j < i + 1; j++) {
        zi_rev[i] -= a_vec[j+1] * zi_rev[i-j];
      }
    }
  }

  // 执行反向滤波
  std::vector<double> z(nx + 2 * nfact);
  z[0] = b_vec[0] * y_rev[0];

  for (int i = 1; i < n_b; i++) {
    z[0] += b_vec[i] * y_rev[0];
  }

  for (int i = 1; i < n_a; i++) {
    z[0] -= a_vec[i] * zi_rev[i-1];
  }

  for (int i = 1; i < nx + 2 * nfact; i++) {
    z[i] = b_vec[0] * y_rev[i];

    for (int j = 1; j < n_b; j++) {
      if (i - j >= 0)
        z[i] += b_vec[j] * y_rev[i - j];
    }

    for (int j = 1; j < n_a; j++) {
      if (i - j >= 0)
        z[i] -= a_vec[j] * z[i - j];
    }
  }

  // 再次反转信号
  std::vector<double> result(nx);
  for (int i = 0; i < nx; i++) {
    result[i] = z[nx + 2 * nfact - 1 - (nfact + i)];
  }

  // 转换回ArrayFire数组
  return af::array(nx, result.data());
}





//    af::array filtfilt(const af::array &b, const af::array &a, const af::array &x) {
//
//        int len_a = a.dims(0);
//        int len_b = b.dims(0);
//        int n = x.dims(0);
//        int filter_order = max(len_a, len_b) - 1;
//        int n_pad = 3 * filter_order;
//
//
//        n_pad = min(n_pad, n - 1);
//        if (n_pad < 1) n_pad = 1;
//
//
//        af::filter();
//        af::array x_pre = af::flip(x(af::seq(n_pad)), 0);
//
//        af::array x_post = af::flip(x(af::seq(n - n_pad, n - 1)), 0);
//
//        af::array x_padded = af::join(0, x_pre, x, x_post);
//
//        // 前向滤波
//        auto forward_filter = [](const af::array &x, const af::array &b,
//                           const af::array &a) {
//            int len_x = x.dims(0);
//            int len_b = b.dims(0);
//            int len_a = a.dims(0);
//
//
//            std::vector<float> h_a(a.elements());
//            a.host(h_a.data());
//            std::vector<float> h_b(b.elements());
//            b.host(h_b.data());
//
//
//            std::vector<float> h_x(x.elements());
//            x.host(h_x.data());
//
//            std::vector<float> h_y(len_x, 0.0f);
//
//            for (int i = 0; i < len_x; ++i) {
//                float sum_b = 0.0f;
//                for (int k = 0; k < len_b; ++k) {
//                    if (i - k >= 0) {
//                        sum_b += h_b[k] * h_x[i - k];
//                    }
//                }
//
//                float sum_a = 0.0f;
//                for (int k = 1; k < len_a; ++k) {
//                    if (i - k >= 0) {
//                        sum_a += h_a[k] * h_y[i - k];
//                    }
//                }
//
//                h_y[i] = (sum_b - sum_a) / h_a[0];
//            }
//
//            return af::array(len_x, h_y.data());
//        };
//
//        af::array y_forward = forward_filter(x_padded, b, a);
//
//        // 反转
//        af::array y_reversed = af::flip(y_forward, 0);
//
//
//        af::array y_reversed_filtered = forward_filter(y_reversed, b, a);
//
//        // 反转
//        af::array y_filtfilt = af::flip(y_reversed_filtered, 0);
//        y_filtfilt = y_filtfilt(af::seq(n_pad, n_pad + n - 1));
//
//        return y_filtfilt;

}
