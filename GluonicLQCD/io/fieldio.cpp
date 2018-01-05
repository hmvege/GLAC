#include "fieldio.h"
#include "config/parameters.h"
#include "parallelization/index.h"
#include "parallelization/communicator.h"
#include <mpi.h>
#include <cmath>

using std::cout;
using std::endl;

const int IO::FieldIO::m_SU3Doubles = 18;
const int IO::FieldIO::m_SU3Size = m_SU3Doubles*sizeof(double);
const int IO::FieldIO::m_linkDoubles = m_SU3Doubles*4;
const int IO::FieldIO::m_linkSize = m_linkDoubles*sizeof(double);
std::vector<unsigned int> IO::FieldIO::m_N;

IO::FieldIO::FieldIO()
{
}

IO::FieldIO::~FieldIO()
{
}

void IO::FieldIO::init()
{
    m_N = Parameters::getN();
}

void IO::FieldIO::writeFieldToFile(Lattice<SU3> *lattice, int configNumber)
{
    /*
     * C-method for writing out configuration to file.
     * Arguments:
     *  configNumber   : (int) configuration number
     */
    MPI_File file;
    std::string filename = Parameters::getFilePath() + Parameters::getOutputFolder() + Parameters::getBatchName() + "/field_configurations/"
                                            + Parameters::getBatchName()
                                            + "_beta" + std::to_string(Parameters::getBeta())
                                            + "_spatial" + std::to_string(Parameters::getNSpatial())
                                            + "_temporal" + std::to_string(Parameters::getNTemporal())
                                            + "_threads" + std::to_string(Parallel::Communicator::getNumProc())
                                            + "_config" + std::to_string(configNumber + Parameters::getConfigStartNumber()) + ".bin";

    MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
            nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
                ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
                    nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                    for (int mu = 0; mu < 4; mu++) {
                        MPI_File_write_at(file, Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size, &lattice[mu][Parallel::Index::getIndex(x,y,z,t)], m_SU3Doubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
    }
    MPI_File_close(&file);
}

void IO::FieldIO::loadFieldConfiguration(std::string filename, Lattice<SU3> *lattice)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     * - lattice
     */
    MPI_File file;
    MPI_File_open(MPI_COMM_SELF, (Parameters::getFilePath() + Parameters::getInputFolder() + filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    for (unsigned int mu = 0; mu < 4; mu++) {
        for (unsigned int t = 0; t < m_N[3]; t++) {
            nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
            for (unsigned int z = 0; z < m_N[2]; z++) {
                nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
                for (unsigned int y = 0; y < m_N[1]; y++) {
                    ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                    for (unsigned int x = 0; x < m_N[0]; x++) {
                        nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                        MPI_File_read_at(file, Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size, &lattice[mu][Parallel::Index::getIndex(x,y,z,t)], m_SU3Doubles, MPI_DOUBLE, MPI_STATUS_IGNORE);
                    }
                }
            }
        }
    }
    MPI_File_close(&file);
    if (Parallel::Communicator::getProcessRank() == 0) printf("\nConfiguration %s loaded", (Parameters::getFilePath() + Parameters::getInputFolder() + filename).c_str());
}

void IO::FieldIO::loadChromaFieldConfiguration(std::string filename, Lattice<SU3> *lattice)
{
    /*
     * Method for loading a field configuration and running the plaquettes on them.
     * Arguments:
     * - filename
     */
    MPI_File file;
    MPI_File_open(MPI_COMM_SELF, (Parameters::getFilePath() + Parameters::getInputFolder() + filename).c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);
    MPI_Offset nt = 0, nz = 0, ny = 0, nx = 0;

    double temp = 0;
    for (unsigned int t = 0; t < m_N[3]; t++) {
        nt = (Parallel::Neighbours::getProcessorDimensionPosition(3) * m_N[3] + t);
        for (unsigned int z = 0; z < m_N[2]; z++) {
            nz = (Parallel::Neighbours::getProcessorDimensionPosition(2) * m_N[2] + z);
            for (unsigned int y = 0; y < m_N[1]; y++) {
                ny = (Parallel::Neighbours::getProcessorDimensionPosition(1) * m_N[1] + y);
                for (unsigned int x = 0; x < m_N[0]; x++) {
                    nx = (Parallel::Neighbours::getProcessorDimensionPosition(0) * m_N[0] + x);
                    for (unsigned int mu = 0; mu < 4; mu++) {
                        for (int i = 0; i < 18; i++) {
                            MPI_File_read_at(file, Parallel::Index::getGlobalIndex(nx,ny,nz,nt)*m_linkSize + mu*m_SU3Size + i*sizeof(double), &temp, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
                            lattice[mu][Parallel::Index::getIndex(x,y,z,t)][i] = reverseDouble(temp);
                            // Checking for corruption
                            if (std::isnan(lattice[mu][Parallel::Index::getIndex(x,y,z,t)][i]))
                            {
                                lattice[mu][Parallel::Index::getIndex(x,y,z,t)].print();
                                printf("\nConfiguration is corrupt.\n");
                                exit(1);
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_File_close(&file);
    if (Parallel::Communicator::getProcessRank() == 0) printf("\nConfiguration %s loaded", (Parameters::getInputFolder() + Parameters::getBatchName()).c_str());
}


inline double IO::FieldIO::reverseDouble(const double inDouble)
{
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
//   delete [] doubleToConvert; MEMORY LEAK HERE?
   return retVal;
}
