///Elek Péter

#ifndef MATRIX_H
#define MATRIX_H
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <fstream>
#include <cstdarg>
#include <vector>
#include <omp.h>



template <class T>
class MATRIX
{
private:
    friend T;
    T* pt;
    int sor;///tetszoleges sorban levo elemek szama
    int oszlop;///tetszoleges pszlopban levo elemek szama

    template <class U> friend MATRIX<U> pow2n(MATRIX<U> const& A, unsigned int n);


public:
    T pol_ertek(T const& x) const;

    MATRIX(int oszlop, int sor, T const& x);
    MATRIX(int oszlop, int sor);
    ~MATRIX();
    MATRIX(const MATRIX& other);
    const MATRIX& operator=(MATRIX const& z)
    {
        if (this != &z)
            delete[] pt;
        else
            return z;
        oszlop = z.oszlop;
        sor = z.sor;
        pt = new T[oszlop * sor];


        for (int i = 0; i < oszlop * sor; ++i)
        {
            pt[i] = z.pt[i];
        }

        return z;
    }
    void ij(int i, int j, T const&);
    T ij_ki(int i, int j);

    void print(void) const;

    MATRIX operator-(MATRIX const& B) const;
    MATRIX operator+(MATRIX const& B) const;
    MATRIX operator*(MATRIX const& B) const;
    MATRIX operator*(T const& z) const;
    MATRIX operator/(T const& z) const;
    friend MATRIX operator*(T const& z, MATRIX const& A)
    {
        return A * z;
    }

    template <class U> friend MATRIX<U> pow(MATRIX<U> const& A, int n);
    template <class U> friend MATRIX<U> exp(MATRIX<U> const&);
    template <class U> friend MATRIX<U> sin(MATRIX<U> const&);///csak komplex T vel megy
    template <class U> friend MATRIX<U> cos(MATRIX<U> const&);///csak komplex T vel megy
    template <class U> friend MATRIX<U> sinh(MATRIX<U> const&);
    template <class U> friend MATRIX<U> cosh(MATRIX<U> const&);
    template <class U> friend MATRIX<U> syl(MATRIX<U> const& X, U(*f)(U));///T re megvalositott fv ertelmezheto diagonalizalhato matrixokra
    template <class U> friend MATRIX<U> syl(MATRIX<U> const& X, U(*f)(U const&));///Sylvester's formula komplex kell legyen a matrix illetve nem lehet 2 azonos sajaterteke
    template <class U> friend MATRIX<U> ftag(MATRIX<U> const& X, U(*f)(U));///elemenkenti uvelet elvegzese
    template <class U> friend MATRIX<U> ftag(MATRIX<U> const& X, U(*f)(U const&));




    void rand_re(void);
    void rand_kp(void);

    T trace(void) const;
    MATRIX charpol(void) const;
    T det(void) const;
    T det_gj(void) const;//Ez pontosabb és gyorsabb
    MATRIX inverse(void) const;
    MATRIX inverse_gj(void) const;//Ez pontosabb és gyorsabb
    MATRIX adj(void) const;
    MATRIX eigenvalue(void) const;
    MATRIX eigenvector(T const&) const;/// a megadott T sajat ertekhez legkozelebbi sajatvektort keresi meg
    MATRIX eigenvector_all(void) const;
    MATRIX eigenvector_check(void) const;
    MATRIX eigenvector_power_iteration(void) const;
    T vector_norm(void) const;

    T oszlopjoeleme(int& j, int n);
    T oszlopjoeleme_det(int& j, int n);//ez a det_gj() hez kell
    void sorcsere(int i, int j);
    void sorszorzasa(int j, T const& z);
    void sorkivonas(int i, int j, T const& z);




    template <class U> friend std::ostream& operator<<(std::ostream& os, MATRIX<U> const& z);
    template <class U> friend bool fajlbair1(const char* s, MATRIX<U> const& z);
    template <class U> friend bool fajlbair(const char* s, MATRIX<U> const& z);
    template <class U> friend bool fajlbeolvas(const char* s, MATRIX<U>& z);


    //Neuron hálóhoz
    template <class U> friend MATRIX<U> sigmoid(MATRIX<U> const& X);
    template <class U> friend MATRIX<U> next_layer(MATRIX<U> const& M, MATRIX<U> const& M_in, MATRIX<U> const& b);



};

template <class T>
std::ostream& operator<<(std::ostream& os, MATRIX<T> const& z)
{
    for (int i = 0; i < z.oszlop * z.sor; ++i)
    {
        os << z.pt[i] << " ";
        if (((i + 1) % z.sor) == 0)
        {
            os << "\n";
        }
    }
    //os<<";";

    return os;
}

template <class T>
bool fajlbair1(const char* s, MATRIX<T> const& z)
{
    std::ofstream file(s);

    file << z.sor << z.oszlop;


    for (int i = 0; i < z.oszlop * z.sor; ++i)
    {
        file << z.pt[i];

    }

    file.close();

    return true;
}

template <class T>
bool fajlbair(const char* s, MATRIX<T> const& z)
{
    std::ofstream file(s, std::ios::out | std::ios::binary | std::ios::trunc);

    file.write((char*)&z.sor, sizeof(int));
    file.write((char*)&z.oszlop, sizeof(int));



    file.write((char*)z.pt, size_t(z.oszlop) * size_t(z.sor) * sizeof(T));

    /*for (int i = 0; i < z.oszlop * z.sor; ++i)
    {
        file.write((char*)&z.pt[i], sizeof(T));

    }*/


    file.close();

    return true;
}

template <class T>
bool fajlbeolvas(const char* s, MATRIX<T>& z)
{
    delete[] z.pt;
    std::ifstream myFile(s, std::ios::in | std::ios::binary);

    //file.seekg(0, std::ios::beg);
    myFile.read((char*)&z.sor, sizeof(int));

    //std::cout << z.sor << "\n\n";

    myFile.read((char*)&z.oszlop, sizeof(int));

    z.pt = new T[size_t(z.sor) * size_t(z.oszlop)];

    myFile.read((char*)z.pt, size_t(z.sor) * size_t(z.oszlop) * sizeof(T));

    /*for (int i = 0; i < z.oszlop * z.sor; ++i)
    {
        myFile.read((char*)&z.pt[i], sizeof(T));

    }*/

    myFile.close();

    return true;
}


template <class T>
MATRIX<T>::MATRIX(int oszlop, int sor, T const& x)
{
    this->oszlop = oszlop;
    this->sor = sor;
    this->pt = new T[oszlop * sor];

    T z(0);
    //#pragma omp parallel for
    for (int i = 0; i < oszlop * sor; ++i)
    {
        this->pt[i] = z;
    }


    int a;

    if (oszlop < sor)
    {
        a = oszlop;
    }
    else
    {
        a = sor;
    }

    //#pragma omp parallel for
    for (int i = 0; i < a; ++i)
    {
        this->pt[i * sor + i] = x;
    }

}

template <class T>
MATRIX<T>::MATRIX(int oszlop, int sor)
{
    this->oszlop = oszlop;
    this->sor = sor;
    this->pt = new T[oszlop * sor];
    T x(0);
    //#pragma omp parallel for
    for (int i = 0; i < oszlop * sor; ++i)
    {
        this->pt[i] = x;
    }
}

template <class T>
MATRIX<T>::MATRIX(const MATRIX<T>& other)
{
    oszlop = other.oszlop;
    sor = other.sor;
    pt = new T[other.oszlop * other.sor];

    //#pragma omp parallel for
    for (int i = 0; i < other.oszlop * other.sor; ++i)
    {
        pt[i] = other.pt[i];
    }
}

template <class T>
MATRIX<T>::~MATRIX()
{
    delete[] pt;
}

template <class T>
void MATRIX<T>::print(void) const
{

    for (int i = 0; i < oszlop * sor; ++i)
    {
        std::cout << pt[i] << " ";
        if (((i + 1) % sor) == 0)
        {
            std::cout << "\n";
        }
    }

}

template <class T>
MATRIX<T> MATRIX<T>::operator+ (MATRIX<T> const& B) const
{
    MATRIX<T> A(oszlop, B.sor);

    if ((oszlop == B.oszlop) && (sor == B.sor))
    {
        //#pragma omp parallel for
        for (int i = 0; i < oszlop * sor; ++i)
        {
            A.pt[i] = pt[i] + B.pt[i];
        }
    }
    else
    {
        throw std::invalid_argument("HIBA nem osszeadhato matrixok\n");
    }



    return A;
}

template <class T>
MATRIX<T> MATRIX<T>::operator- (MATRIX<T> const& B) const
{
    MATRIX<T> A(oszlop, B.sor);

    if ((oszlop == B.oszlop) && (sor == B.sor))
    {
        //#pragma omp parallel for
        for (int i = 0; i < oszlop * sor; ++i)
        {
            A.pt[i] = pt[i] - B.pt[i];
        }
    }
    else
    {
        throw std::invalid_argument("HIBA nem kivonhato matrixok\n");
    }



    return A;
}

template <class T>
MATRIX<T> MATRIX<T>::operator* (MATRIX<T> const& B) const
{
    MATRIX<T> A(oszlop, B.sor);
    if (sor == B.oszlop)
    {
#pragma omp parallel for collapse(2)
        for (int i = 0; i < oszlop; ++i)
        {
            for (int j = 0; j < B.sor; ++j)
            {
                T z(0);
                for (int k = 0; k < sor; ++k)
                {
                    z = z + pt[i * sor + k] * B.pt[k * B.sor + j];
                }

                A.pt[i * A.sor + j] = z;
            }

        }


        return A;
    }
    else
    {
        throw std::invalid_argument("HIBA nem osszeszorozhato matrixok\n");
    }


}

template <class T>
MATRIX<T> MATRIX<T>::operator* (T const& z) const
{
    MATRIX<T> A(oszlop, sor);
    //#pragma omp parallel for
    for (int i = 0; i < sor * oszlop; ++i)
    {
        A.pt[i] = pt[i] * z;
    }



    return A;

}

template <class T>
MATRIX<T> MATRIX<T>::operator/ (T const& z) const
{
    MATRIX<T> A(oszlop, sor);
    //#pragma omp parallel for
    for (int i = 0; i < sor * oszlop; ++i)
    {
        A.pt[i] = pt[i] / z;
    }



    return A;

}

template <class T>
MATRIX<T> pow(MATRIX<T> const& A, int n)
{
    T z(1);
    MATRIX<T> I(A.oszlop, A.sor, z);

    if (n < 0)
        return pow(A.inverse_gj(), (-n));
    else if (n == 0)
        return  I;
    else if (n == 1)
        return  A;
    else if (n % 2 == 0)
        return pow((A * A), (n / 2));
    else /*if (n%2==1)*/
        return A * pow((A * A), ((n - 1) / 2));
}

template <class T>
MATRIX<T> pow2n(MATRIX<T> const& A, unsigned int n)
{
    T z(1);
    MATRIX<T> I(A.oszlop, A.sor, z);
    if (n == 1)
        return  A;
    else
        return pow((A * A), (n / 2));

}

template <class T>
T MATRIX<T>::trace(void) const
{
    T x(0);

    int a;

    if (oszlop < sor)
    {
        a = oszlop;
    }
    else
    {
        a = sor;
    }

    for (int i = 0; i < a; ++i)
    {
        x = x + (this->pt[i * sor + i]);
    }
    return x;
}
template <class T>
MATRIX<T> MATRIX<T>::charpol(void) const
{
    MATRIX Z(1, sor + 1);

    if (sor == oszlop)
    {
        T c(1);
        MATRIX M(sor, sor), I(sor, sor, c);

        Z.pt[sor] = c;
        T a(-1);
        T g(1);

        for (int k = 1; k <= sor; ++k)
        {
            g = a / T(k);
            M = (*this) * M + I * c;
            c = g * ((*this) * M).trace();
            Z.pt[sor - k] = c;

        }

    }

    else
    {
        throw std::invalid_argument("NEM n*n-es a matrix\n");
    }

    return Z;


}
template <class T>
T MATRIX<T>::det(void) const
{
    T x(-1);


    if (((this->oszlop) % 2) == 0)
    {
        x = (this->charpol().pt[0]);
        return x;
    }
    else
    {
        x = (this->charpol().pt[0]) * x;
        return x;
    }
}



template <class T>
MATRIX<T> MATRIX<T>::inverse(void) const
{

    if (oszlop == sor)
    {
        MATRIX pol(1, 1);
        pol = this->charpol();
        MATRIX a(sor, sor);
        T zero(0);

        if (pol.pt[0] == zero)
            throw std::invalid_argument("determinans 0\nnem invertalhato a matrix\n");
        ///#pragma omp parallel for reduction(+:a)///Ez így kajak mûködik????-NEM
        for (int i = 1; i <= (sor); ++i)
        {
            a = a + (pow((*this), (i - 1)) * pol.pt[i]);
        }

        if ((sor % 2) == 0)
        {
            T x(-1);
            a = a * (x / (this->det()));
            return a;
        }
        else
        {
            T x(1);
            a = a * (x / (this->det()));
            return a;
        }


    }
    else
    {
        throw std::invalid_argument("nem n*n es a matrix tehat az inverz nem egyertelmu\n");
    }


}


template <class T>
T MATRIX<T>::det_gj(void) const
{



    if (oszlop == sor)
    {
        int j;
        MATRIX A = *this;
        T x(1), e(1), ex(-1), zero(0), g;




        for (int k = 0; k < sor; ++k)
        {
            if (A.oszlopjoeleme_det(j, k) == zero)
            {
                return zero;
            }

            A.sorcsere(j, k);

            if (j != k)
            {
                x = x * ex;
            }


            g = e / (A.pt[k * sor + k]);
            x = x * g;
            A.sorszorzasa(k, g);

            //... Ezt nem lehet párhuzamosítani ...///#pragma omp parallel for
            for (int l = k; l < sor; ++l)
            {
                if (l != k)
                {
                    g = A.pt[l * sor + k];
                    A.sorkivonas(k, l, g);


                }
            }

        }





        return e / x;
    }
    else
    {
        throw std::invalid_argument("nem n*n es a matrix tehat determinanst nem ertelmezek\n");
    }
}


template <class T>
MATRIX<T> MATRIX<T>::inverse_gj(void) const
{


    if (oszlop == sor)
    {
        int j;
        MATRIX A = *this;
        T e(1), g;
        MATRIX I(sor, sor, e);



        for (int k = 0; k < sor; ++k)
        {
            A.oszlopjoeleme(j, k);
            A.sorcsere(j, k);
            I.sorcsere(j, k);

            g = e / (A.pt[k * sor + k]);
            A.sorszorzasa(k, g);
            I.sorszorzasa(k, g);
            //... Ezt nem lehet párhuzamosítani ...///#pragma omp parallel for
            for (int l = 0; l < sor; ++l)
            {
                if (l != k)
                {
                    g = A.pt[l * sor + k];
                    A.sorkivonas(k, l, g);
                    I.sorkivonas(k, l, g);

                }
            }

        }





        return I;
    }
    else
    {
        throw std::invalid_argument("nem n*n es a matrix tehat az inverz nem egyertelmu\n");
    }


}

template <class T>
T MATRIX<T>::oszlopjoeleme(int& j, int n)
{
    T zero(0);
    if ((0 <= n) && (n < sor))
    {

        for (int i = n; i < oszlop; ++i)
        {

            if (pt[i * sor + n] != zero)
            {
                j = i;
                return pt[i * sor + n];
            }

        }

        std::cout << *this << "\n";

        throw std::invalid_argument("nincs jo elem az egyik oszlopban\n");

    }
    else
    {
        throw std::invalid_argument("tul indexeles oszlopjoeleme hivasaban\n");
    }


}

template <class T>
T MATRIX<T>::oszlopjoeleme_det(int& j, int n)
{
    T zero(0);
    if ((0 <= n) && (n < sor))
    {

        for (int i = n; i < oszlop; ++i)
        {

            if (pt[i * sor + n] != zero)
            {
                j = i;
                return pt[i * sor + n];
            }

        }

        return zero;

    }
    else
    {
        throw std::invalid_argument("tul indexeles oszlopjoeleme hivasaban\n");
    }


}

template <class T>
void MATRIX<T>::sorszorzasa(int j, T const& z)
{

    if (j >= 0 && j < oszlop)
    {
        //#pragma omp parallel for
        for (int i = 0; i < sor; ++i)
        {
            pt[j * sor + i] = z * pt[j * sor + i];
        }

    }
    else
    {
        throw std::invalid_argument("tul indexeles sorszorzasaban\n");
    }


}
template <class T>
void MATRIX<T>::sorcsere(int i, int j)
{
    T t;

    if ((i >= 0) && (j >= 0) && (j < oszlop) && (i < oszlop))
    {


        if (i != j)
        {
            //#pragma omp parallel for
            for (int k = 0; k < sor; ++k)
            {
                t = pt[i * sor + k];
                pt[i * sor + k] = pt[j * sor + k];
                pt[j * sor + k] = t;
            }
        }
    }

    else
    {
        throw std::invalid_argument("tul indexeles sorcsereben\n");
    }
}

template <class T>
void MATRIX<T>::sorkivonas(int i, int j, T const& z)
{
    if ((i >= 0) && (j >= 0) && (j < oszlop) && (i < oszlop))
    {
        if (i != j)
        {
            //#pragma omp parallel for
            for (int k = 0; k < sor; ++k)
            {
                pt[j * sor + k] = pt[j * sor + k] - z * pt[i * sor + k];
            }
        }
    }

    else
    {
        throw std::invalid_argument("tul indexeles sorcsereben\n");
    }
}

template <class T>
T MATRIX<T>::pol_ertek(T const& x) const
{
    T a(0, 0);
    for (int i = 0; i < sor; ++i)
    {
        a = a + (pt[i] * pow(x, i));
    }
    return a;
}

template <class T>
MATRIX<T> MATRIX<T>::adj(void) const
{
    MATRIX<T> A(sor, oszlop);

    T z;
    //#pragma omp parallel for collapse(2)
    for (int j = 0; j < oszlop; ++j)
    {
        for (int i = 0; i < sor; ++i)
        {
            z = pt[i + sor * j];
            z.imag = -z.imag;
            A.pt[j + oszlop * i] = z;
        }
    }



    return A;
}

template <class T>
T MATRIX<T>::vector_norm(void) const ///oszlop vektor legyen
{
    if (sor != 1)
    {
        throw std::invalid_argument("nem oszlop vektorra lett hivva a vector_norm\n");
    }


    return sqrt(((this->adj()) * (*this)).pt[0]);
}

template <class T>
MATRIX<T> MATRIX<T>::eigenvalue(void) const
{
    if (sor != oszlop)
    {
        throw std::invalid_argument("eigenvalue hivas nem nxn es matrixra\n");
    }

    MATRIX<T> a = this->charpol();


    MATRIX<T> x(1, a.sor - 1), z(1, a.sor - 1);

    MATRIX c(1, a.sor - 1);
    for (int i = 0; i < a.sor - 1; ++i)
    {
        c.pt[i] = a.pt[i + 1] * (i + 1);
    }


    T f(0, 0), g(0, 0), one(1, 0);///itt mar kell hogy komplex legyen ezert 2 egészes a konstruktor hivasa ide dob errort ha nem komplex MATRIX-amire hivjuk
//#pragma omp parallel for
    for (int i = 0; i < a.sor - 1; ++i) ///kezdo ertekek iniciizalasa ez az alapmegoldas
    {
        x.pt[i].real = ((double)rand()) / RAND_MAX;
        x.pt[i].imag = ((double)rand()) / RAND_MAX;

    }
    for (int h = 0; h < 10000; ++h)
    {
        //#pragma omp parallel for
        for (int i = 0; i < a.sor - 1; ++i)
        {
            f = (a.pol_ertek(x.pt[i])) / (c.pol_ertek(x.pt[i]));
            g.real = 0;
            g.imag = 0;

            for (int k = 0; k < i; ++k)
            {
                g = g + one / (x.pt[i] - x.pt[k]);
            }
            for (int k = i + 1; k < a.sor - 1; ++k)
            {
                g = g + one / (x.pt[i] - x.pt[k]);
            }

            z.pt[i] = x.pt[i] - f / (one - f * g);

        }
        x = z;
        ///std::cout<<x<<"\n";
    }



    return x;
}

template <class T>
MATRIX<T> MATRIX<T>::eigenvector(T const& z) const///ha egzaktul vagy az abrazolasi pontossagon belul talalod el a sajaterteket akkor nem invertalhato az iteralando matrix ez a nincs jo elemet dobja mint hiba
{
    MATRIX<T> v(oszlop, 1);

    v.rand_re();
    v = v * (T(1) / v.vector_norm());

    if (sor != oszlop)
    {
        throw std::invalid_argument("eigenvector hivas nem nxn es matrixra\n");
    }
    MATRIX<T> I(oszlop, sor, T(1)), A(1, 1);
    A = *this;
    A = *this - I * z;
    A = A.inverse_gj();

    for (unsigned int i = 0; i < 100000; ++i)
    {
        v = A * v * (T(1) / (A * v).vector_norm());
    }

    return v;

}

template <class T>
MATRIX<T> MATRIX<T>::eigenvector_power_iteration(void) const
{
    MATRIX<T> v(oszlop, 1);

    v.rand_kp();
    v = v * (T(1) / v.vector_norm());

    if (sor != oszlop)
    {
        throw std::invalid_argument("eigenvector hivas nem nxn es matrixra\n");
    }

    for (unsigned int i = 0; i < 1000; ++i)
    {
        v = (*this) * v * (T(1) / ((*this) * v).vector_norm());
    }

    return v;

}

template <class T>
MATRIX<T> MATRIX<T>::eigenvector_all(void) const///oszlopaiba lesznek a sajat vektorai
{
    MATRIX<T> se = eigenvalue();
    MATRIX<T> v(oszlop, 1);

    MATRIX<T> sv(sor, sor);
    for (int i = 0; i < sor; ++i)
    {
        v = this->eigenvector(se.pt[i]);
        //#pragma omp parallel for
        for (int j = 0; j < sor; ++j)
        {
            sv.pt[j * sor + i] = v.pt[j];
        }
    }




    return sv;

}

template <class T>
MATRIX<T> MATRIX<T>::eigenvector_check(void) const///oszlopaiba lesznek a sajat vektorai
{
    MATRIX<T> se = eigenvalue();
    MATRIX<T> v(oszlop, 1);

    MATRIX<T> sv(sor, sor);
    for (int i = 0; i < sor; ++i)
    {
        v = this->eigenvector(se.pt[i] + T(0.01));
        v = (*this) * v - se.pt[i] * v;
        //#pragma omp parallel for
        for (int j = 0; j < sor; ++j)
        {
            sv.pt[j * sor + i] = v.pt[j];
        }
    }






    return sv;

}



template <class T>
void MATRIX<T>::rand_re(void)
{
    //#pragma omp parallel for
    for (int i = 0; i < oszlop * sor; ++i)
    {
        this->pt[i] = T(double(rand()) / RAND_MAX);
    }
}

template <class T>
void MATRIX<T>::rand_kp(void)
{
    //#pragma omp parallel for
    for (int i = 0; i < oszlop * sor; ++i)
    {
        this->pt[i] = T(double(rand()) / RAND_MAX, double(rand()) / RAND_MAX);
    }
}

template<class T>
void MATRIX<T>::ij(int i, int j, T const& z)
{
    pt[i * sor + j] = z;
}

template<class T>
T MATRIX<T>::ij_ki(int i, int j)
{
    return pt[i * sor + j];
}


template <class T>
MATRIX<T> exp(MATRIX<T> const& x)
{
    unsigned int n = 2147483648;
    MATRIX<T> b(x.sor, x.oszlop, T(1));
    T c = T(n) * T(2);

    return pow2n(b + x / c, n) * pow2n(b + x / c, n);
}

template <class T>
MATRIX<T> sin(MATRIX<T> const& x)
{
    return (exp(T(0, 1) * x) - exp(T(0, -1) * x)) / T(0, 2);
}

template <class T>
MATRIX<T> cos(MATRIX<T> const& x)
{
    return (exp(T(0, 1) * x) + exp(T(0, -1) * x)) / T(2, 0);
}

template <class T>
MATRIX<T> sinh(MATRIX<T> const& x)
{
    return (exp(x) - exp(T(-1) * x)) / T(2);
}

template <class T>
MATRIX<T> cosh(MATRIX<T> const& x)
{
    return (exp(x) + exp(T(-1) * x)) / T(2);
}

template <class T>
MATRIX<T> ftag(MATRIX<T> const& x, T(*f)(T))
{
    MATRIX<T> A(x.oszlop, x.sor);
    //#pragma omp parallel for
    for (int i = 0; i < x.sor * x.oszlop; ++i)
    {
        A.pt[i] = f(x.pt[i]);
    }

    return A;
}

template <class T>
MATRIX<T> ftag(MATRIX<T> const& x, T(*f)(T const&))
{
    MATRIX<T> A(x.oszlop, x.sor);
    //#pragma omp parallel for
    for (int i = 0; i < x.sor * x.oszlop; ++i)
    {
        A.pt[i] = f(x.pt[i]);
    }

    return A;
}



template <class T>
MATRIX<T> syl(MATRIX<T> const& x, T(*f)(T&))
{

    MATRIX<T> X(x.oszlop, x.sor), I(x.oszlop, x.sor, T(1));
    MATRIX<T> pol = x.eigenvalue();
    MATRIX<T> A(x.oszlop, x.sor);

    for (int i = 0; i < pol.sor; ++i)
    {
        A = I;

        for (int j = 0; j < pol.sor; ++j)
        {
            if (i != j)
            {
                A = A * (T(1) / (pol.pt[i] - pol.pt[j])) * (x - pol.pt[j] * I);
            }

        }
        X = X + f(pol.pt[i]) * A;
    }

    return X;
}

template <class T>
MATRIX<T> syl(MATRIX<T> const& x, T(*f)(T const&))
{

    MATRIX<T> X(x.oszlop, x.sor), I(x.oszlop, x.sor, T(1));
    MATRIX<T> pol = x.eigenvalue();
    MATRIX<T> A(x.oszlop, x.sor);


    for (int i = 0; i < pol.sor; ++i)
    {
        A = I;

        for (int j = 0; j < pol.sor; ++j)
        {
            if (i != j)
            {
                A = A * (T(1) / (pol.pt[i] - pol.pt[j])) * (x - pol.pt[j] * I);
            }

        }
        X = X + f(pol.pt[i]) * A;
    }

    return X;
}



template <class T>
MATRIX<T> next_layer(MATRIX<T> const& M, MATRIX<T> const& M_in, MATRIX<T> const& b)
{
    return sigmoid((M * M_in) - b);
}


template <class T>
MATRIX<T> sigmoid(MATRIX<T> const& x)
{
    MATRIX<T> A(x.oszlop, x.sor);
    //#pragma omp parallel for
    for (int i = 0; i < x.sor * x.oszlop; ++i)
    {
        A.pt[i] = 1.0 / (1 + exp(-x.pt[i]));
    }

    return A;
}


#endif // MATRIX_H






