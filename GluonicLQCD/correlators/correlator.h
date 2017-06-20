#ifndef CORRELATOR_H
#define CORRELATOR_H


class Correlator
{
private:
    int m_latticeSize;
public:
    Correlator(int latticeSize);
    ~Correlator();
    virtual double calculate(Links *lattice);
};

#endif // CORRELATOR_H
