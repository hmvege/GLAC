#include "fieldio.h"
#include "config/parameters.h"
#include "parallelization/index.h"
#include "parallelization/communicator.h"
#include <mpi.h>
#include <cmath>
#include <fstream>

using std::cout;
using std::endl;

const long long IO::FieldIO::m_SU3Doubles = 18;
const long long IO::FieldIO::m_SU3Size = m_SU3Doubles*sizeof(double);
const long long IO::FieldIO::m_linkDoubles = m_SU3Doubles*4;
const long long IO::FieldIO::m_linkSize = m_linkDoubles*sizeof(double);
std::vector<unsigned int> IO::FieldIO::m_N;

/*!
 * \namespace IO
 *
 * Holds the field input/output class FieldIO as well as methods for writing observables to file.
 */

/*!
 * \brief IO::FieldIO::FieldIO initializes the FieldIO object.
 *
 * Empty constructor.
 */
IO::FieldIO::FieldIO()
{
}

/*!
 * \brief IO::FieldIO::~FieldIO
 *
 * Empty destructor.
 */
IO::FieldIO::~FieldIO()
{
}

/*!
 * \brief IO::FieldIO::init updates the lattice dimensions.
 *
 * Used in unit testing.
 */
void IO::FieldIO::init()
{
    m_N = Parameters::getN();
}

/*!
 * \brief IO::FieldIO::writeFieldToFile writes a lattice of SU3 matrix objects to file.
 * \param lattice is a pointer of four lattice objects, one for each Lorentz index.
 * \param configNumber the configuration number to name the file with.
 *
 * Uses MPI_File_write_at to write the configuration to a file.
 *
 * If debug is true in the passed .json parameter file, it will perform a check for lattice corruption.
 */
void IO::FieldIO::writeFieldToFile(Lattice<SU3> *lattice, unsigned int configNumber)
{
    /*
     * C-method for writing out configuration to file.
     *
     * Arguments:
     *  configNumber   : (int) configuration number
     */
    MPI_File file;

    // Converting config number to a more machine friendly layout
    char cfg_number[6];
    sprintf(cfg_number, "%05d", configNumber + Parameters::getConfigStartNumber());

    std::string filename = Parameters::getBatchName() + "_b" + std::to_string(Parameters::getBeta())
            + "_N" + std::to_string(Parameters::getNSpatial())
            + "_NT" + std::to_string(Parameters::getNTemporal())
            + "_np" + std::to_string(Parallel::Communicator::getNumProc())
            + "_config" + std::string(cfg_number) + ".bin";

    std::string filenamePath = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + "/field_configurations/" + filename;

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(lattice, "Configuration is corrupt in IO::FieldIO write field before");
    }

    MPI_File_open(Parallel::ParallelParameters::ACTIVE_COMM, filenamePath.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    long long nt = 0, nz = 0, ny = 0, nx = 0;
    MPI_Offset offset = 0;

    for (long long mu = 0; mu < 4; mu++) {
        for (long long t = 0; t < m_N[3]; t++) {
            nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
            for (long long z = 0; z < m_N[2]; z++) {
                nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
                for (long long y = 0; y < m_N[1]; y++) {
                    ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                    for (long long x = 0; x < m_N[0]; x++) {
                        nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                        offset = Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size;
                        MPI_File_write_at(file, offset, &lattice[mu][Parallel::Index::getIndex(x,y,z,t)], m_SU3Doubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
    }
    MPI_File_close(&file);

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(lattice, "Configuration is corrupt in IO::FieldIO write field after");
    }

    if (Parallel::Communicator::getProcessRank() == 0 && !Parameters::getUnitTesting()) {
        printf("    %s written.", filename.c_str());
    }
}

/*!
 * \brief IO::FieldIO::writeDoublesFieldToFile writes a lattice of doubles to file.
 * \param lattice is a pointer of four lattice objects, one for each Lorentz index.
 * \param configNumber configuration number.
 * \param observable name of the observable.
 *
 * Useful if one wants to visualize the lattice.
 * \sa https://github.com/hmvege/LatViz for a configuration visualization method.
 *
 * Uses MPI_File_write_at to write the configuration to a file.
 *
 * If debug is true in the passed .json parameter file, it will perform a check for lattice corruption.
 */
void IO::FieldIO::writeDoublesFieldToFile(Lattice<double> lattice, unsigned int configNumber, std::string observable)
{
    /*
     * C-method for writing out single double lattice field to file. No lorentz indices
     *
     * Arguments:
     *  lattice: Lattice<double>, single field of doubles
     *  configNumber   : (int) configuration number
     */
    MPI_File file;

    // Converting config number to a more machine friendly layout
    char cfg_number[6];
    sprintf(cfg_number, "%05d", configNumber + Parameters::getConfigStartNumber());

    std::string filename = Parameters::getBatchName() + "_b" + std::to_string(Parameters::getBeta())
            + "_N" + std::to_string(Parameters::getNSpatial())
            + "_NT" + std::to_string(Parameters::getNTemporal())
            + "_np" + std::to_string(Parallel::Communicator::getNumProc())
            + "_config" + std::string(cfg_number) + ".bin";

    std::string filenamePath = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName()
            + "/scalar_fields/" + observable + "/" + filename;

    MPI_File_open(Parallel::ParallelParameters::ACTIVE_COMM, filenamePath.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    MPI_Offset offset = 0;
    long long nt = 0, nz = 0, ny = 0, nx = 0;

    for (unsigned long int t = 0; t < m_N[3]; t++) {
        nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned long int z = 0; z < m_N[2]; z++) {
            nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned long int y = 0; y < m_N[1]; y++) {
                ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned long int x = 0; x < m_N[0]; x++) {
                    nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                    offset = Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*sizeof(double);
                    MPI_File_write_at(file, offset, &lattice[Parallel::Index::getIndex(x,y,z,t)], 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                }
            }
        }
    }
    MPI_File_close(&file);

    if (Parallel::Communicator::getProcessRank() == 0) {
        printf("\n    %s written.", filenamePath.c_str());
    }
}

/*!
 * \brief IO::FieldIO::loadFieldConfiguration loads a regular configuration into memory.
 * \param filename
 * \param lattice is a pointer of four lattice objects, one for each Lorentz index, that will filled with the configuration read from file.
 *
 * If debug is true in the passed .json parameter file, it will perform a check for lattice corruption.
 */
void IO::FieldIO::loadFieldConfiguration(std::string filename, Lattice<SU3> *lattice)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     *
     * Arguments:
     * - filename
     * - lattice
     */
    MPI_File file;

    // Sets up file name
    std::string fname = Parameters::getFilePath() + Parameters::getInputFolder() + filename;

    // Checks if file we are trying to load exists or not
    if (!check_file_existence(fname.c_str())) {
        Parallel::Communicator::MPIExit("File " + fname + " does not exist");
    }

    MPI_File_open(Parallel::ParallelParameters::ACTIVE_COMM, fname.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    MPI_Offset offset = 0;
    long long nt = 0, nz = 0, ny = 0, nx = 0;

    for (long long mu = 0; mu < 4; mu++) {
        for (long long t = 0; t < m_N[3]; t++) {
            nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
            for (long long z = 0; z < m_N[2]; z++) {
                nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
                for (long long y = 0; y < m_N[1]; y++) {
                    ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                    for (long long x = 0; x < m_N[0]; x++) { // TODO: Maybe change this to read x at once
                        nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                        offset = Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size;
                        MPI_File_read_at(file, offset, &lattice[mu][Parallel::Index::getIndex(x,y,z,t)], m_SU3Doubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
    }

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(lattice, "Configuration is corrupt in IO::FieldIO load field after loaded");
    }

    MPI_File_close(&file);
    if (Parallel::Communicator::getProcessRank() == 0 && !Parameters::getUnitTesting()) {
        printf("\nConfiguration %s loaded", fname.c_str());
    }
}

/*!
 * \brief IO::FieldIO::loadChromaFieldConfiguration loads a configuration from Chroma into memory.
 * \param filename
 * \param lattice is a pointer of four lattice objects, one for each Lorentz index, that will filled with the configuration read from file.
 *
 * If debug is true in the passed .json parameter file, it will perform a check for lattice corruption.
 */
void IO::FieldIO::loadChromaFieldConfiguration(std::string filename, Lattice<SU3> *lattice)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    MPI_File file;

    // Sets up file name
    std::string fname = Parameters::getFilePath() + Parameters::getInputFolder() + filename;

    // Checks if file we are trying to load exists or not
    if (!check_file_existence(fname.c_str())) {
        Parallel::Communicator::MPIExit("File " + fname + " does not exist");
    }

    //    if (Parallel::Communicator::getProcessRank() == 0) printf("\nConfiguration to be loaded: %s", fname.c_str());

    MPI_File_open(Parallel::ParallelParameters::ACTIVE_COMM, fname.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    MPI_Offset offset = 0;
    long long nt = 0, nz = 0, ny = 0, nx = 0;


    double temp = 0;
    for (long long mu = 0; mu < 4; mu++) {
        for (long long t = 0; t < m_N[3]; t++) {
            nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
            for (long long z = 0; z < m_N[2]; z++) {
                nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
                for (long long y = 0; y < m_N[1]; y++) {
                    ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                    for (long long x = 0; x < m_N[0]; x++) {
                        nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                        for (long long i = 0; i < 18; i++) {
                            offset = Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size + i*sizeof(double);
                            MPI_File_read_at(file, offset, &temp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                            lattice[mu][Parallel::Index::getIndex(x,y,z,t)][i] = reverseDouble(temp);
                        }
                    }
                }
            }
        }
    }

    MPI_File_close(&file);

    if (Parameters::getDebug()) {
        Parallel::Communicator::checkLattice(lattice, "Configuration is corrupt in IO::FieldIO load chroma field after loaded.");
    }

    if (Parallel::Communicator::getProcessRank() == 0) printf("\nConfiguration %s loaded", fname.c_str());

}

inline bool IO::FieldIO::check_file_existence (const std::string fname) {
    std::ifstream infile(fname);
    return infile.good();
}

inline double IO::FieldIO::reverseDouble(const double inDouble)
{
    /*
    * Method for reversing the bytes in a double, since Chroma uses a reversed ordering for its doubles.
    */
    double retVal;
    char *doubleToConvert = ( char* ) & inDouble;
    char *returnDouble = ( char* ) & retVal;

    // swap the bytes into a temporary buffer
    returnDouble[0] = doubleToConvert[7];
    returnDouble[1] = doubleToConvert[6];
    returnDouble[2] = doubleToConvert[5];
    returnDouble[3] = doubleToConvert[4];
    returnDouble[4] = doubleToConvert[3];
    returnDouble[5] = doubleToConvert[2];
    returnDouble[6] = doubleToConvert[1];
    returnDouble[7] = doubleToConvert[0];
    //   delete [] doubleToConvert; // MEMORY LEAK HERE?
    return retVal;
}
