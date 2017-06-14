#include <su3.h>
#include <complex.h>

void SU3BaseTests()
{
    /*
     * Basic tests of the complex class.
     */
    complex a = complex(1,1);
    complex b = complex(2,2);
    a *= b;
    cout << a << endl;
    SU3 A, B, C, D;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A.mat[(i*3+j)] = complex(i,j);
            B.mat[(i*3+j)] = complex(i*i,j*j);
            C.mat[(i*3+j)] = complex(0,0);
            D.mat[(i*3+j)] = complex(0,0);
        }
    }
    cout << endl;
    A.print();
    cout << endl;
    B.print();
    cout << A.get(1,1) << endl;
    C = A + B;
    D = A * B;
    D.print();
    std::cout << "exiting in SU3MatrixGenerator.cpp line 56" << std::endl;
}
