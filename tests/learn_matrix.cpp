#include <iostream>
#include <vector>

#include <matrix_impl.hpp>

int main()
{
    // Initialize a column vector with n x m elements initialized to zero
    int n = 3;
    int m = 4;
    Math::Matrix<double> mat(n, m, 0);
    std::vector<double> vec {2, 3, -1};
    // Initialize a row vector
    Math::Matrix<double> mat(vec, false);
    Math::Matrix<double> matT = mat.getTranspose();

    std::cout << mat(0, 2) << std::endl;
    std::cout << matT(2, 0) << std::endl;

    Math::Matrix<double> res;
    Math::Matrix<double>::multiplyMatrices(mat, matT, &res);
    res.writeIntoTextFile("result.txt");

    return 0;
}
