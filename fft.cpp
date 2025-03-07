//
// Created by Yanran LI on 2025/1/19.
//
#include "fft.h"
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
