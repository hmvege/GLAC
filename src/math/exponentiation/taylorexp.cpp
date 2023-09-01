#include "taylorexp.h"

#include <stdexcept>

TaylorExp::TaylorExp(unsigned int N)
{
    verifyTaylorDegree(N);
}

SU3 TaylorExp::exp(SU3 Q)
{
    SU3 QMul;
    SU3 QSum;
    
    QMul.identity();
    QSum.identity();

    double taylorFactor = 1;

    for (unsigned int i = 1; i < m_N + 1; i++) {
        QMul *= Q;
        QSum += (QMul / taylorFactor);
        taylorFactor *= (static_cast<double>(i) + 1.0);
    }
    return QSum;
}

/*!
 * \brief TaylorExp::setTaylorDegree
 * \param N taylor polynomial degree.
 *
 * // TODO: remove this function, since it is strictly redundant
 * 
 * Sets the taylor polynomial degree.
 */
void TaylorExp::setTaylorDegree(unsigned int N) {
    verifyTaylorDegree(N);
    m_N = N;
}

void TaylorExp::verifyTaylorDegree(const int n) {
    if (n < 1) {
        std::string error_message = "Error: invalid degree " + std::to_string(n) + " of Taylor expansion. Minimum degree allowed is 1.";
        throw std::invalid_argument(error_message);
    }
    if (n > 16) {
        std::string error_message = "Error: invalid degree " + std::to_string(n) + " of Taylor expansion. Maximum allowed degree is 16";
        throw std::invalid_argument(error_message);
    }
}
