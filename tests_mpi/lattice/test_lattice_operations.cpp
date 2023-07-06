// #include <actions/action.h>
#include <actions/actions.h>
#include <config/parameters.h>
#include <lattice/test_lattice_operations.h>
#include <parallelization/communicator.h>
#include <parallelization/index.h>

#include <cmath>

class BaseTest
{
public:
  BaseTest() {}
  ~BaseTest() = default;

private:
};

bool testLatticeOperations()
{
  bool passed = true;

  const size_t n_spatial = 8;
  const size_t n_temporal = n_spatial * 2;

  Parameters::setBeta(6.0);
  Parameters::setNSpatial(n_spatial);
  Parameters::setNTemporal(n_temporal);

  // The sub-lattice size is always initialized in the Communicator
  Parallel::Communicator::initializeSubLattice();

  const auto n_dim = Parameters::getN();

  WilsonGaugeAction S;

  Lattice<SU3>* L = new Lattice<SU3>[4];
  Lattice<SU3>* LNew = new Lattice<SU3>[4];
  for (int i = 0; i < 4; i++)
  {
    L[i].allocate(n_dim);
    LNew[i].allocate(n_dim);
  }

  for (int mu = 0; mu < 4; mu++)
  {
    L[mu].identity();
    LNew[mu].zeros();
  }

  // Tests the action derivative as used in flow.
  for (int mu = 0; mu < 4; mu++)
  {
    LNew[mu] = S.getActionDerivative(L, mu);
  }

  SU3 updateLink;
  updateLink.identity();
  updateLink *= 2;

  double dS = 0;

  for (unsigned int x = 0; x < L->m_dim[0]; x++)
  {
    for (unsigned int y = 0; y < L->m_dim[1]; y++)
    {
      for (unsigned int z = 0; z < L->m_dim[2]; z++)
      {
        for (unsigned int t = 0; t < L->m_dim[3]; t++)
        {
          for (int mu = 0; mu < 4; mu++)
          {
            S.computeStaple(L, x, y, z, t, mu);
            dS = S.getDeltaAction(L[mu][Parallel::Index::getIndex(x, y, z, t)],
                                  updateLink);

            if (std::fabs(dS + 36) > 1e-15)
            {
              cout << "Error in action(should equal 36): " << dS << "\n";
              passed = false;
              x = L->m_dim[0];
              y = L->m_dim[1];
              z = L->m_dim[2];
              t = L->m_dim[3];
              break;
            }
          }
        }
      }
    }
  }

  delete[] L;
  delete[] LNew;

  return passed;
}
