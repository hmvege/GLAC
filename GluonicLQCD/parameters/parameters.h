#ifndef PARAMETERS_H
#define PARAMETERS_H


class parameters
{
public:
    parameters();
    // Total lattice sizes
    static int NSpatial;
    static int NTemporal;
    static int latticeSize;
    // Sub lattice sizes
    static unsigned int N[4];
    static int subLatticeSize;
    // Beta value constant
    static double beta;
    // Lattice spacing
    static double a;

};

#endif // PARAMETERS_H
