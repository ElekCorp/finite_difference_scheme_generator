#include <iostream>
#include <gmpxx.h>

#include "MATRIX.h"

template <class T>
T pow(const T& x, int i)
{
    if (i < 0)
        return pow(T(1) / x, -i);
        //std::cout << "pow exponent negative\n";
    else if (i == 0)
        return T(1);
    else if (i == 1)
        return x;
    else if (i % 2 == 0)
        return pow(x, i / 2) * pow(x, i / 2);
    else//i%2==1
        return x * pow(x, i / 2) * pow(x, i / 2);
}
template<>
mpq_class pow(const mpq_class& x, int i)
{
    if (i < 0)
    {
        mpq_class inv;
        mpq_inv(inv.get_mpq_t(), x.get_mpq_t());
        return pow(inv, -i);
    }
    else if (i == 0)
        return mpq_class(1);
    else if (i == 1)
        return x;
    else if (i % 2 == 0)
        return pow(x, i / 2) * pow(x, i / 2);
    else//i%2==1
        return x * pow(x, i / 2) * pow(x, i / 2);
}

mpq_class matrix_coef(int i, int j, int k);
MATRIX<mpq_class> init_matrix(int p, int k);
MATRIX<mpq_class> init_vec(int m, int p);
MATRIX<mpq_class> get_coef_vec(int p, int m, int k);

int main(int argc, char** argv)
{
    
    int p = atoi(argv[1]);
    if (p <= 0) p = 1;
    /*int k = atoi(argv[2]);
    if (k <= 0) k = p;
    int m = atoi(argv[3]);
    if (m <= 0) m = 2;//*/
    

    //int p = 30;
    int k = p;
    int m = 2;


    MATRIX<mpq_class> cucc = get_coef_vec(p, m, k);

    std::cout << cucc;

    return 0;
}

mpq_class matrix_coef(int i, int j, int k)
{
    mpq_class gmpj_k = j-k;

    return pow(gmpj_k, i);

}

MATRIX<mpq_class> init_matrix(int p, int k)
{
    int num = 2 * p + 1;
    MATRIX<mpq_class> A(num, num, 0);
    for (int i = 0; i < num; ++i)
    {
        for (int j = 0; j < num; ++j)
        {
            A.ij(i, j, matrix_coef(i, j, k));
        }
    }

    return A;
}

MATRIX<mpq_class> init_vec(int m, int p)
{
    int num = 2 * p + 1;
    MATRIX<mpq_class> b(num, 1, 0);

    b.ij(m, 0, factorial(mpz_class(m)));

    return b;
}

MATRIX<mpq_class> get_coef_vec(int p, int m, int k)
{
    MATRIX<mpq_class> A = init_matrix(p, k);
    MATRIX<mpq_class> b = init_vec(m, p);

    /*std::cout << A << "\n\n";
    std::cout << b << "\n\n";//*/
    return (A.inverse_gj()) * b;
}
