#include <iostream>
#include <vector>

#include <matrix_impl.hpp>
#include <macros.hpp>

int main()
{
    // Initialize a column vector with n x m elements initialized to zero
    int n = 3;
    int m = 4;
    Math::Matrix<double> mat(n, m, 0);
    check(mat.rows() == n, "");
    check(mat.cols() == m, "");

    // Initialize a row vector and its transpose
    std::vector<double> vec {2, 3, -1};
    Math::Matrix<double> mat2(vec, false);
    Math::Matrix<double> mat2T = mat2.getTranspose();
    check(mat2(0, 2) == mat2T(2, 0), "");

    // Multiply mat with matT and store result in res
    Math::Matrix<double> res;
    Math::Matrix<double>::multiplyMatrices(mat2, mat2T, &res);
    check(res(0,0) == 14, "");

    // Initialize an empty matrix
    Math::Matrix<double>mat3;
    // Resize a matrix to be n x m and contain all zeros
    mat3.resize(n, m, 0);
    for(int i = 0; i < n; ++i)
    {
        for(int j = 0; j < m; ++j)
        {
            check(mat3(0, 0) == 0, "");
        }
    }

    return 0;
}
