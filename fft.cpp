//
// Created by Yanran LI on 2025/1/19.
//
#include "fft.h"
namespace FFT {
//    Eigen::VectorXd fft(Eigen::VectorXd & A) {
//
//    }
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A) {
//        if(A.size() == 0) {
//            return Eigen::MatrixXcd(0,0);
//        }
//        const long long Rows = A.rows();
//        const long long Cols = A.cols();
//        Eigen::MatrixXcd Temp = A;
//        Eigen::FFT<double> F;
//        Eigen::MatrixXcd Result(Temp.rows(), Temp.cols());  // 结果
//    }
//
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A, int n) {
//
//    }

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
}
