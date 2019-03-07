#include "taylorexp.h"

TaylorExp::TaylorExp(unsigned int N)
{
    m_N = N;
    if (m_N < 1) {
        printf("Error: invalid degree %d of Taylor expansion.", m_N);
        exit(1);
    }
    m_QMul.identity();
    m_QSum.identity();
}

SU3 TaylorExp::exp(SU3 Q)
{
    /*
     * Exponentiate using regular Taylor expansion.
     */
    m_QMul.identity();
    m_QSum.identity();
    m_taylorFactor = 1;

    for (unsigned int i = 1; i < m_N + 1; i++) {
        m_QMul *= Q;
        m_QSum += (m_QMul / m_taylorFactor);
        m_taylorFactor *= (double(i) + 1.0);

    }
    return m_QSum;
}

void TaylorExp::setTaylorDegree(unsigned int N) {
    m_N = N;
    m_QMul.identity();
    m_QSum.identity();
    m_taylorFactor = 1;
}
