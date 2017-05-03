#ifndef ACTION_H
#define ACTION_H

class Action
{
protected:
    double (*potential)(double x);
    double a;
    int N;
public:
    Action();
    Action(int NPathPoints, double new_a);
    virtual ~Action() {}
    virtual double getAction(double * x, int i);
    // Setters
    void setPotential(double (*pot)(double x)) { potential = pot; }
    void setLatticeSpacing(double new_a) { a = new_a; }
    void setNLatticePoints(int new_N) { N = new_N; }
    // Getters
    double getLatticeSpacing() { return a; }
    int getNLatticePoints() { return N; }
};

#endif // ACTION_H
