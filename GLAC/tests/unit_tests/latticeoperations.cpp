#include "latticeoperations.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"

// Needed for performing gauge invariance test
#include "io/fieldio.h"
#include "observables/plaquette.h"

LatticeOperations::LatticeOperations()
{

}

////////////////////////////////////
////////// LATTICE TESTS ///////////
////////////////////////////////////
bool LatticeOperations::testLatticeAddition() {
    /*
     * Tests the addition of lattice object of type SU3, complex and double.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice(m_dim);
    tempSU3Lattice.zeros();
    Lattice<complex>tempComplexLattice(m_dim);
    tempComplexLattice.zeros();
    Lattice<double>tempDoubleLattice(m_dim);
    tempDoubleLattice.zeros();
    tempSU3Lattice = latticeSU3_U1 + latticeSU3_U2;
    tempComplexLattice = latticeComplex_z1 + latticeComplex_z2;
    tempDoubleLattice = latticeDouble1 + latticeDouble2;
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],UAdd)) {
            passed = false;
            cout << "    FAILED: lattice SU3 addition is not correct." << endl;
            break;
        }
        if (!compareComplex(tempComplexLattice[iSite],zAdd)) {
            passed = false;
            cout << "    FAILED: lattice complex addition is not correct." << endl;
            break;
        }
        if (fabs(tempDoubleLattice[iSite] - (latticeDoubleValue1 + latticeDoubleValue2)) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: lattice double addition is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: lattice SU3, complex and double addition is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeSubtraction() {
    /*
     * Tests the subtraction of lattice object of type SU3, complex and double.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice(m_dim);
    tempSU3Lattice.zeros();
    Lattice<complex>tempComplexLattice(m_dim);
    tempComplexLattice.zeros();
    Lattice<double>tempDoubleLattice(m_dim);
    tempDoubleLattice.zeros();
    tempSU3Lattice = latticeSU3_U1 - latticeSU3_U2;
    tempComplexLattice = latticeComplex_z1 - latticeComplex_z2;
    tempDoubleLattice = latticeDouble1 - latticeDouble2;
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],USub)) {
            passed = false;
            cout << "    FAILED: lattice SU3 subtraction is not correct." << endl;
            break;
        }
        if (!compareComplex(tempComplexLattice[iSite],zSub)) {
            passed = false;
            cout << "    FAILED: lattice complex subtraction is not correct." << endl;
            break;
        }
        if (fabs(tempDoubleLattice[iSite] - (latticeDoubleValue1 - latticeDoubleValue2)) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: lattice double subtraction is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: lattice SU3, complex and double subtraction is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeMultiplication() {
    /*
     * Tests the multiplication of lattice object of type SU3, complex and double.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice(m_dim);
    tempSU3Lattice.zeros();
    Lattice<complex>tempComplexLattice(m_dim);
    tempComplexLattice.zeros();
    Lattice<double>tempDoubleLattice(m_dim);
    tempDoubleLattice.zeros();
    tempSU3Lattice = latticeSU3_U1*latticeSU3_U2;
    tempComplexLattice = latticeComplex_z1*latticeComplex_z2;
    tempDoubleLattice = latticeDouble1*latticeDouble2;
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],UMul)) {
            passed = false;
            cout << "    FAILED: lattice SU3 multiplication is not correct." << endl;
            break;
        }
        if (!compareComplex(tempComplexLattice[iSite],zMul)) {
            passed = false;
            cout << "    FAILED: lattice complex multiplication is not correct." << endl;
            break;
        }
        if (fabs(tempDoubleLattice[iSite] - (latticeDoubleValue1*latticeDoubleValue2)) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: lattice double multiplication is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: lattice SU3, complex and double multiplication is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeDivision() {
    /*
     * Tests the division of lattice object of type SU3, complex and double.
     */
    bool passed = true ;
    Lattice<complex>tempComplexLattice(m_dim);
    tempComplexLattice.zeros();
    Lattice<double>tempDoubleLattice(m_dim);
    tempDoubleLattice.zeros();
    tempComplexLattice = latticeComplex_z1 / z2;
    tempDoubleLattice = latticeDouble1 / latticeDoubleValue2;
    for (unsigned int iSite = 0; iSite < tempComplexLattice.m_latticeSize; iSite++) {
        if (!compareComplex(tempComplexLattice[iSite],zDiv)) {
            passed = false;
            cout << "    FAILED: lattice complex division is not correct." << endl;
            break;
        }
        if (fabs(tempDoubleLattice[iSite] - (latticeDoubleValue1/latticeDoubleValue2)) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: lattice double division is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: lattice complex and double division is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeRealTrace() {
    /*
     * Tests taking the real trace of a lattice object.
     */
    bool passed = true ;
    Lattice<double>tempDoubleLattice(m_dim);
    double tempDoubleTrace = U1.trace().z[0];
    tempDoubleLattice = realTrace(latticeSU3_U1);
    for (unsigned int iSite = 0; iSite < tempDoubleLattice.m_latticeSize; iSite++) {
        if (fabs(tempDoubleLattice[iSite] - tempDoubleTrace) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: realTrace() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: realTrace() of a SU3 lattice is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeImagTrace() {
    /*
     * Tests taking the imaginary trace of a lattice object.
     */
    bool passed = true ;
    Lattice<double>tempDoubleLattice(m_dim);
    double tempDoubleTrace = U1.trace().z[1];
    tempDoubleLattice = imagTrace(latticeSU3_U1);
    for (unsigned int iSite = 0; iSite < tempDoubleLattice.m_latticeSize; iSite++) {
        if (fabs(tempDoubleLattice[iSite] - tempDoubleTrace) > m_machineEpsilon) {
            passed = false;
            cout << "    FAILED: imagTrace() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: imagTrace() of a SU3 lattice is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeTrace() {
    /*
     * Tests taking the trace of a lattice object(complex return value).
     */
    bool passed = true ;
    Lattice<complex>tempComplexLattice(m_dim);
    complex tempComplexTrace = U1.trace();
    tempComplexLattice = trace(latticeSU3_U1);
    for (unsigned int iSite = 0; iSite < tempComplexLattice.m_latticeSize; iSite++) {
        if (!compareComplex(tempComplexLattice[iSite],tempComplexTrace)) {
            passed = false;
            cout << "    FAILED: trace() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: trace() of a SU3 lattice is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeSubtractReal() {
    /*
     * Tests subtracting from the real elements along the diagonal in a SU3 lattice object.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice  = latticeSU3_U1;;
    SU3 tempSU3;
    tempSU3 = U1;
    for (int iMat = 0; iMat < 18; iMat+=8) {
        tempSU3[iMat] -= latticeDoubleValue1;
    }
    tempSU3Lattice = subtractReal(tempSU3Lattice, latticeDouble1);
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],tempSU3)) {
            passed = false;
            cout << "    FAILED: subtractReal() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) {
        cout << "    SUCCESS: subtractReal() of a SU3 lattice is correct." << endl;
    }
    return passed;
}

bool LatticeOperations::testLatticeSubtractImag() {
    /*
     * Tests subtracting from the imaginary elements along the diagonal in a SU3 lattice object.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice = latticeSU3_U2;
    SU3 tempSU3 = U2;
    for (int iMat = 1; iMat < 18; iMat+=8) {
        tempSU3[iMat] -= latticeDoubleValue2;
    }
    tempSU3Lattice = subtractImag(tempSU3Lattice, latticeDouble2);
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],tempSU3)) {
            passed = false;
            cout << "    FAILED: subtractImag() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (passed && m_verbose) cout << "    SUCCESS: subtractImag() of a SU3 lattice is correct." << endl;
    return passed;
}

bool LatticeOperations::testLatticeSum() {
    /*
     * Tests summing the lattice into a single variable T.
     */
    bool passed = true ;
    SU3 tempSU3;
    tempSU3 = U1*double(latticeSU3_U1.m_latticeSize);
    if (!compareSU3(sum(latticeSU3_U1),tempSU3)) {
        passed = false;
        cout << "    FAILED: sum() of a SU3 lattice is not correct." << endl;
    }
    if (m_verbose && passed) cout << "    SUCCESS: sum() of a SU3 lattice is correct." << endl;
    return passed;
}

bool LatticeOperations::testLatticeSumRealTrace() {
    /*
     * Tests summing, tracing and taking the real part of the lattice.
     */
    bool passed = true ;
    double tempDouble = U1.trace().z[0]*double(latticeSU3_U1.m_latticeSize);
    double tempSU3SumRealTrace = sumRealTrace(latticeSU3_U1);
    if (fabs(tempSU3SumRealTrace - tempDouble) > m_machineEpsilon) {
        passed = false;
        cout << "    FAILED: sumRealTrace() of a SU3 lattice is not correct." << endl;
    }
    if (m_verbose && passed) cout << "    SUCCESS: sumRealTrace() of a SU3 lattice is correct." << endl;
    return passed;
}

bool LatticeOperations::testLatticeSumRealTraceMultiplication() {
    /*
     * Tests multiplying two lattices and then summing, tracing and taking the real part of the lattice.
     */
    bool passed = true ;
    double tempDouble = (U1*U2).trace().z[0]*double(latticeSU3_U1.m_latticeSize);
    double tempSU3SumRealTrace = sumRealTraceMultiplication(latticeSU3_U1,latticeSU3_U2);
    if (fabs(tempSU3SumRealTrace - tempDouble) > m_machineEpsilon) {
        passed = false;
        cout << "    FAILED: sumRealTraceMultiplication() of a SU3 lattice is not correct." << endl;
    }
    if (m_verbose && passed) cout << "    SUCCESS: sumRealTraceMultiplication() of a SU3 lattice is correct." << endl;
    return passed;
}

bool LatticeOperations::testLatticeInverse() {
    /*
     * Tests the inverse of the lattice.
     */
    bool passed = true ;
    Lattice<SU3>tempSU3Lattice(m_dim);
    tempSU3Lattice = inv(latticeSU3_U1);
    for (unsigned int iSite = 0; iSite < tempSU3Lattice.m_latticeSize; iSite++) {
        if (!compareSU3(tempSU3Lattice[iSite],UCT)) {
            passed = false;
            cout << "    FAILED: inv() of a SU3 lattice is not correct." << endl;
            break;
        }
    }
    if (m_verbose && passed) cout << "    SUCCESS: lattice inverse is correct." << endl;
    return passed;
}

bool LatticeOperations::testLatticeShift() {
    /*
     * Tests the lattice shift function for all possible directions.
     */
    bool passed = true;

    Parallel::Communicator::setBarrier();

    if (Parallel::ParallelParameters::active) {
        unsigned int position = 0;
        unsigned int N1,N2,N3;
        Lattice<SU3>L;
        L.allocate(m_dim);
        for (unsigned int dir = 0; dir < 8; dir++) {

            // Sets the lattice sites with the value of the processor they are running on.
            for (unsigned int iSite = 0; iSite < m_subLatticeSize; iSite++) {
                L[iSite] = double(m_processRank);
            }

            // Specifies if we are going backwards or farwards
            if (dir % 2 == 0) {
                // Backwards
                L = shift(L,BACKWARDS,dir / 2);
            } else {
                // Forwards
                L = shift(L,FORWARDS,dir / 2);
            }

            // Sets correct sublattice face to check.
            if (dir / 2 == 0) {
                N1 = L.m_dim[1];
                N2 = L.m_dim[2];
                N3 = L.m_dim[3];
            }
            else if (dir / 2 == 1) {
                N1 = L.m_dim[0];
                N2 = L.m_dim[2];
                N3 = L.m_dim[3];
            }
            else if (dir / 2 == 2) {
                N1 = L.m_dim[0];
                N2 = L.m_dim[1];
                N3 = L.m_dim[3];
            }
            else {
                N1 = L.m_dim[0];
                N2 = L.m_dim[1];
                N3 = L.m_dim[2];
            }

            // Loops over cube and checks that it has received the correct values
            for (unsigned int i = 0; i < N1; i++) {
                for (unsigned int j = 0; j < N2; j++) {
                    for (unsigned int k = 0; k < N3; k++) {

                        // Sets correct position
                        if (dir / 2 == 0) { // x direction
                            position = Parallel::Index::getIndex((L.m_dim[0] - 1) * (dir % 2),i,j,k);
                        } else if (dir / 2 == 1) { // y direction
                            position = Parallel::Index::getIndex(i,(L.m_dim[1] - 1) * (dir % 2),j,k);
                        } else if (dir / 2 == 2) { // z direction
                            position = Parallel::Index::getIndex(i,j,(L.m_dim[2] - 1) * (dir % 2),k);
                        } else if (dir / 2 == 3) { // t direction
                            position = Parallel::Index::getIndex(i,j,k,(L.m_dim[3] - 1) * (dir % 2));
                        } else {
                            Parallel::Communicator::MPIExit("Error in testLatticeShift");
                        }

                        // Compares matrices at position.
                        for (unsigned int iMat = 0; iMat < 18; iMat++) {

                            // Remember, we are shifting the lattice, such that we for say FORWARDS shift, we update the values at 0
                            // in a given face with the values from the preceding shift.
                            if (L.m_sites[position].mat[iMat] != double(Parallel::Neighbours::get((dir) % 2 + (dir / 2) * 2))) {
                                passed = false;

                                // Prints diagnostics by rank order.
                                for (int iRank = 0; iRank < Parallel::Communicator::getNumProc(); iRank++) {
                                    if (m_processRank==iRank) {
                                        printf("\nRank: %d Shift: %d Dir: %d Expected-neighbour from %d: %d Value received: %f\n",
                                               m_processRank, dir%2, dir, (dir + 1) % 2 + (dir / 2) * 2,
                                               Parallel::Neighbours::get((dir + 1) % 2 + (dir / 2) * 2),
                                               L.m_sites[position].mat[iMat]);
                                        Parallel::Neighbours::getNeighbours(m_processRank)->print();
                                        cout << "Processor position: ";
                                        for (int iDim = 0; iDim < 4; iDim++) {
                                            cout << Parallel::Neighbours::getProcessorDimensionPosition(iDim) << " ";
                                        }
                                        cout << endl;
                                    }
                                    Parallel::Communicator::setBarrier();
                                }
                                if (m_processRank == 0) {
                                    cout << "Error print done. Stopping test." << endl;
                                }

                                i = N1;
                                j = N2;
                                k = N3;
                                iMat = 18;
                                dir = 8;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    Parallel::Communicator::setBarrier();
    if (passed) {
        if (m_processRank == 0 && Parameters::getUnitTestingVerbose()) cout << "    SUCCESS: Lattice shift." << endl;
    } else {
        if (m_processRank == 0) cout << "    FAILED: Lattice shift." << endl;
    }

    return passed;
}

bool LatticeOperations::testFieldGaugeInvariance() {
    /*
     * Tests if the lattice remaince gauge invariant
     */
    bool passed = true;
    double obsBefore = 0, obsAfter = 0;

    // Temporary sets the number of observables we are to store to 2
    unsigned int tempNCf = Parameters::getNCf();
    Parameters::setNCf(2);

    // Temporary changes the input folder for the files, as the file is provided from the command line
    std::string tempInputFolderPath = Parameters::getInputFolder();
    Parameters::setInputFolder(std::string("/"));

    // Initializes the observable - the plaquette
    Plaquette P(false);
    P.setLatticeSize(m_subLatticeSize);

    // Initializes the lattice and allocates it with correct dimensions
    Lattice<SU3> L[4];
    for (int i = 0; i < 4; i++) L[i].allocate(m_dim);

    // Retrieves the gauge field to check the gauge invariance of
    std::string gaugeFieldFileName = Parameters::getGaugeFieldToCheck();
    IO::FieldIO::loadFieldConfiguration(gaugeFieldFileName,L);

    // Calculates the plaquette of the lattice
    P.calculate(L,0);
    obsBefore = P.getObservable(0);
    MPI_Allreduce(&obsBefore,&obsBefore,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    obsBefore /= double(m_numprocs);
    if (m_processRank == 0 && m_verbose) printf("\n    Plaquette before gauge transformation: %20.16f\n",obsBefore);

    // Sets up the gauge invariance lattice and populates it with random SU3 matrices.
    Lattice<SU3> Omega,OmegaNext;
    Omega.allocate(m_dim);
    for (unsigned int iSite = 0; iSite < Omega.m_latticeSize; iSite++) {
        Omega[iSite] = m_SU3Generator->generateRST();
    }

    // Loops over the different directions and performs the gauge transformation.
    for (int mu = 0; mu < 4; mu++) {
        // Left hand multiplies.
        OmegaNext = inv(shift(Omega,FORWARDS,mu));
        L[mu] = Omega*L[mu]*OmegaNext;
    }

    // Calculates the plaquette after the gauge transformation
    P.calculate(L,1);
    obsAfter = P.getObservable(1);
    MPI_Allreduce(&obsAfter,&obsAfter,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    obsAfter /= double(m_numprocs);
    if (m_processRank == 0 && m_verbose) printf("    Plaquette after gauge transformation:  %20.16f\n",obsAfter);
    if (fabs(obsAfter - obsBefore) > 1e-15) {
        passed = false;
        if (m_processRank == 0) cout << "    FAILED: Lattice " << gaugeFieldFileName << " is not invariant under a SU3 gauge transformation." << endl;
    } else {
        if (m_processRank == 0) cout << "    SUCCESS: Lattice " << gaugeFieldFileName << " is invariant under a SU3 gauge transformation." << endl;
    }
    Parameters::setNCf(tempNCf);
    Parameters::setInputFolder(tempInputFolderPath);
    return passed;
}

bool LatticeOperations::fullLatticeTests()
{
    /*
     * Tests that require an entire lattice(with sublattices) to run tests on.
     */
    bool passed = true;
    // Initializes
    if (m_processRank == 0 && m_verbose) {
        printf("Running Lattice parallel tests on sublattice of size %d^3 x %d.\n",m_N,m_NT);
    }
    // Runs tests
    if (Parameters::getCheckFieldGaugeInvariance()) {
        passed = testFieldGaugeInvariance();
    }
    if (passed && testLatticeShift()) {
        if (m_processRank == 0) cout << "PASSED: Parallel lattice tests." << endl;
    } else {
        passed = false;
        if (m_processRank == 0) cout << "FAILURE: Parallel lattice tests." << endl;
    }
    return passed;
}

bool LatticeOperations::runLatticeTests()
{
    /*
     * Function for testing all aspects of the lattice class:
     *  - addition
     *  - subtraction
     *  - multiplication
     *  - division(for complex and real)
     *  - real/imag trace
     *  - subtracting a real/imag
     *  - summing the lattice
     *  - summing and taking the real trace
     *  - multiplying two lattices and taking the real sum
     *  - finding the inverse
     *  - lattice shift
     *  - (if config is provided) gauge field invariance tests
     */
    if (m_verbose && m_processRank == 0) cout << "Running lattice tests." << endl;

    bool passed = false;
    if (m_processRank == 0)
    {
        passed = (testLatticeAddition() && testLatticeSubtraction() && testLatticeMultiplication()
                  && testLatticeDivision() && testLatticeRealTrace() && testLatticeImagTrace()
                  && testLatticeSubtractReal() && testLatticeSubtractImag() && testLatticeSum()
                  && testLatticeSumRealTrace() && testLatticeSumRealTraceMultiplication()
                  && testLatticeInverse() && testLatticeTrace());
    }

    MPI_Bcast(&passed, 1, MPI_BYTE, 0, MPI_COMM_WORLD);
    Parallel::Communicator::setBarrier();

    passed = passed && fullLatticeTests();

    if (passed) {
        if (m_processRank == 0) cout << "PASSED: Lattice operations and functions." << endl;
    } else {
        if (m_processRank == 0) cout << "FAILURE: Lattice operations and functions." << endl;
    }
    Parallel::Communicator::setBarrier();

    return passed;
}

