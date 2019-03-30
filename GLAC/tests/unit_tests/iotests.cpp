#include "iotests.h"
#include "parallelization/parallel.h"
#include "config/parameters.h"
#include "io/fieldio.h"

IOTests::IOTests()
{

}

// IO tests
bool IOTests::testIOLatticeWriteRead()
{
    if ((m_processRank == 0) && m_verbose) {
        printf("Running Lattice IO tests on lattice of size %d^3 x %d\n", Parameters::getNSpatial(), Parameters::getNTemporal());
    }

    using std::chrono::steady_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;


    bool passed = true;

    if (Parallel::ParallelParameters::active) {

        // Create lattice
        // Initializes the lattice and allocates it with correct dimensions
        Lattice<SU3> * LBefore = new Lattice<SU3>[4];
        Lattice<SU3> * LAfter = new Lattice<SU3>[4];
        for (int i = 0; i < 4; i++) LBefore[i].allocate(m_dim);
        for (int i = 0; i < 4; i++) LAfter[i].allocate(m_dim);

        for (int mu = 0; mu < 4; mu++)
        {
            for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
            {
                LBefore[mu][iSite] = m_SU3Generator->generateRandom(); // Fully random
            }
        }

        steady_clock::time_point t0_write;
        steady_clock::time_point t1_write;
        steady_clock::time_point t0_read;
        steady_clock::time_point t1_read;

        t0_write = steady_clock::now();

        // Save lattice
        IO::FieldIO::writeFieldToFile(LBefore, 0);

        t1_write= steady_clock::now();

        // Since we are by default writing to output, we temporarily set the input to be output.
        std::string tmpInputFolder = Parameters::getInputFolder();
        Parameters::setInputFolder("/output/" + Parameters::getBatchName() + "/field_configurations/");

        Parallel::Communicator::setBarrierActive();

        // Simple cp from writeFieldToFile()
        std::string cfg_name = Parameters::getBatchName() + "_b" + std::to_string(Parameters::getBeta())
                + "_N" + std::to_string(Parameters::getNSpatial())
                + "_NT" + std::to_string(Parameters::getNTemporal())
                + "_np" + std::to_string(Parallel::Communicator::getNumProc())
                + "_config" + std::string("00000") + ".bin";

        t0_read = steady_clock::now();

        // Load into new lattice
        IO::FieldIO::loadFieldConfiguration(cfg_name, LAfter);

        t1_read = steady_clock::now();

        // Resets the input folder
        Parameters::setInputFolder(tmpInputFolder);

        if (m_processRank == 0 && m_verbose) {
            printf("    Write time: %.16f seconds.\n    Read time:  %.16f\n",
                   duration_cast<duration<double>>(t1_write - t0_write).count(),
                   duration_cast<duration<double>>(t1_read - t0_read).count());
        }

        // Compare loaded lattice with created lattice
        for (int mu = 0; mu < 4; mu++)
        {
            for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
            {
                for (int i = 0; i < 18; i++) {
                    if (fabs(LBefore[mu][iSite][i] - LAfter[mu][iSite][i]) > 1e-15) {
                        cout << "Error in: " << LBefore[mu][iSite][i] << " " << LAfter[mu][iSite][i] << endl;
                        passed = false;
                        break;
                    }
                }
            }
        }

        delete [] LBefore;
        delete [] LAfter;
    }

    if (passed) {
        if (m_processRank == 0) cout << "PASSED: IO lattice write and load." << endl;
    } else {
        if (m_processRank == 0) cout << "FAILED: IO lattice write or load." << endl;
    }


    return passed;
}

bool IOTests::testIOWriteDoubles()
{
    if ((m_processRank == 0) && m_verbose) {
        printf("    Testing IO doubles writing for lattice of size %d^3 x %d\n", Parameters::getNSpatial(), Parameters::getNTemporal());
    }

    using std::chrono::steady_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;

    bool passed = true;

    if (Parallel::ParallelParameters::active) {

        // Create lattice
        // Initializes the lattice and allocates it with correct dimensions
        Lattice<double> LBefore;
        Lattice<double> LAfter;
        LBefore.allocate(m_dim);
        LAfter.allocate(m_dim);

        for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
        {
            LBefore[iSite] = m_uniform_distribution(m_generator); // Fully random
        }

        steady_clock::time_point t0_write;
        steady_clock::time_point t1_write;

        t0_write = steady_clock::now();

        std::string filenameIOTest = "ioDoublesIOTest";

        // Save lattice
        IO::FieldIO::writeDoublesFieldToFile(LBefore, 0, filenameIOTest);

        t1_write= steady_clock::now();

        // Since we are by default writing to output, we temporarily set the input to be output.
        std::string tmpInputFolder = Parameters::getInputFolder();
        Parameters::setInputFolder("/output/" + Parameters::getBatchName() + "/field_configurations/");

        Parallel::Communicator::setBarrierActive();

        // Simple cp from writeFieldToFile()
        std::string cfg_name = Parameters::getBatchName() + "_b" + std::to_string(Parameters::getBeta())
                + "_N" + std::to_string(Parameters::getNSpatial())
                + "_NT" + std::to_string(Parameters::getNTemporal())
                + "_np" + std::to_string(Parallel::Communicator::getNumProc())
                + "_config" + std::string("00000") + ".bin";

        std::string filenamePath = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName()
                + "/scalar_fields/" + filenameIOTest + "/" + cfg_name;

        // Load into new lattice
        MPI_File file;
        MPI_File_open(Parallel::ParallelParameters::ACTIVE_COMM, filenamePath.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

        MPI_Offset offset = 0;
        std::vector<unsigned int> N =Parameters::getN();
        long long nt = 0, nz = 0, ny = 0, nx = 0;

        for (unsigned long int t = 0; t < N[3]; t++) {
            nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * N[3] + t);
            for (unsigned long int z = 0; z < N[2]; z++) {
                nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * N[2] + z);
                for (unsigned long int y = 0; y < N[1]; y++) {
                    ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * N[1] + y);
                    for (unsigned long int x = 0; x < N[0]; x++) {
                        nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * N[0] + x);
                        offset = Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*sizeof(double);
                        MPI_File_read_at(file, offset, &LAfter[Parallel::Index::getIndex(x,y,z,t)], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
        MPI_File_close(&file);

        // Resets the input folder
        Parameters::setInputFolder(tmpInputFolder);

        if (m_processRank == 0 && m_verbose) {
            printf("    Write time: %.16f seconds.\n",
                   duration_cast<duration<double>>(t1_write - t0_write).count());
        }

        // Compare loaded lattice with created lattice
        for (unsigned long iSite = 0; iSite < m_subLatticeSize; iSite++)
        {
            for (int i = 0; i < 18; i++) {
                if (fabs(LBefore[iSite] - LAfter[iSite]) > 1e-15) {
                    cout << "Error in: " << LBefore[iSite] << " " << LAfter[iSite] << endl;
                    passed = false;
                    break;
                }
            }
        }
    }

    if (passed) {
        if (m_processRank == 0) cout << "PASSED: IO lattice doubles write and load." << endl;
    } else {
        if (m_processRank == 0) cout << "FAILED: IO lattice doubles write or load." << endl;
    }

    return passed;
}

bool IOTests::runIOTests()
{
    bool passed = (testIOWriteDoubles() && testIOLatticeWriteRead());

    if (m_processRank == 0) {
        if (passed) {
            cout << "PASSED: IO properties." << endl;
        } else {
            cout << "FAILURE: IO properties." << endl;
        }
    }
    Parallel::Communicator::setBarrier();

    return passed;
}
