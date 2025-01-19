//
// Created by Yanran LI on 2025/1/19.
//
#include "fft.h"
namespace FFT {
//    Eigen::VectorXd fft(Eigen::VectorXd & A) {
//
//    }
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A) {
//        std::cout << A << std::endl;
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
        Eigen::MatrixXcd Result(Rows, Cols);  // 结果


        if(dim == 1) { // 各列计算
            for(int k = 0;k < Cols;k++) {
                Eigen::VectorXcd tmp(Cols);
                F.fwd(tmp, Temp.col(k));
                Result.col(k) = tmp;
            }
        }
        else if(dim == 2) { // 各行计算
            for(int k = 0;k < Rows;k++) {
                Eigen::VectorXcd tmp(Rows);
                F.fwd(tmp, Temp.row(k));
                Result.row(k) = tmp;
            }
        }
        else {}
        return Result;
    }
}
