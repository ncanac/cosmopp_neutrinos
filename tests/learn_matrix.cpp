#include <iostream>
#include <vector>

#include <matrix_impl.hpp>

int main()
{
    std::vector<double> vec {2, 3, -1};
    Math::Matrix<double> mat(vec, false);
    Math::Matrix<double> matT = mat.getTranspose();

    std::cout << mat(0, 2) << std::endl;
    std::cout << matT(2, 0) << std::endl;

    Math::Matrix<double> res;
    Math::Matrix<double>::multiplyMatrices(mat, matT, &res);
    res.writeIntoTextFile("result.txt");

    return 0;
}
