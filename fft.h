//
// Created by Yanran LI on 2025/1/19.
//

#ifndef FFT_FFT_H
#define FFT_FFT_H
#include "Eigen/Dense"
#include "unsupported/Eigen/FFT"

namespace FFT {
//    Eigen::VectorXd fft(Eigen::VectorXd & A);
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A);
//    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A, int n);
    Eigen::MatrixXcd fft(Eigen::MatrixXcd &A, int n = 0,int dim = 1 );
}
#endif //FFT_FFT_H

