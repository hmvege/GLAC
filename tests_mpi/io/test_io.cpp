#include <config/parameters.h>
#include <io/fieldio.h>
#include <math/lattice.h>
#include <mpi.h>
#include <parallelization/communicator.h>
#include <parallelization/index.h>
#include <parallelization/neighbours.h>
#include <parallelization/parallelparameters.h>

#include <array>
#include <catch2/catch_all.hpp>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace
{
  const std::string TAGS = "[mpi][io]";
  constexpr unsigned int N_SPATIAL = 8;
  constexpr unsigned int N_TEMPORAL = 16;

  std::filesystem::path getTestRoot()
  {
    return std::filesystem::temp_directory_path() / "glac_mpi_tests";
  }

  std::filesystem::path getOutputRoot() { return getTestRoot() / "output"; }

  std::filesystem::path getBatchRoot(const std::string& batchName)
  {
    return getOutputRoot() / batchName;
  }

  std::string withTrailingSlash(const std::filesystem::path& path)
  {
    return path.string() + "/";
  }

  std::string getConfigurationFilename(const unsigned int configNumber)
  {
    char cfgNumber[6];
    std::snprintf(cfgNumber, sizeof(cfgNumber), "%05d",
                  configNumber + Parameters::getConfigStartNumber());

    return Parameters::getBatchName() + "_b" +
           std::to_string(Parameters::getBeta()) + "_N" +
           std::to_string(Parameters::getNSpatial()) + "_NT" +
           std::to_string(Parameters::getNTemporal()) + "_np" +
           std::to_string(Parallel::Communicator::getNumProc()) + "_config" +
           std::string(cfgNumber) + ".bin";
  }

  unsigned int getGlobalCoordinate(const unsigned int dim,
                                   const unsigned int localCoordinate)
  {
    return static_cast<unsigned int>(
             Parallel::Neighbours::getProcessorDimensionPosition(dim)) *
             Parameters::getN().at(dim) +
           localCoordinate;
  }

  std::size_t getGlobalScalarIndex(const unsigned int x, const unsigned int y,
                                   const unsigned int z, const unsigned int t)
  {
    const std::size_t n = Parameters::getNSpatial();
    return x + n * (y + n * (z + n * t));
  }

  double encodeScalarValue(const unsigned int x, const unsigned int y,
                           const unsigned int z, const unsigned int t)
  {
    return static_cast<double>(x + 10 * y + 100 * z + 1000 * t);
  }

  SU3 encodeGaugeMatrix(const unsigned int x, const unsigned int y,
                        const unsigned int z, const unsigned int t,
                        const unsigned int mu)
  {
    SU3 matrix;
    for (int i = 0; i < 18; ++i)
    {
      matrix[i] = static_cast<double>(i + 100 * mu + 1000 * x + 10000 * y +
                                      100000 * z + 1000000 * t);
    }
    return matrix;
  }

  void configureIOEnvironment(const std::string& batchName)
  {
    Parameters::setLatticeSize(1);
    Parameters::setSubLatticeSize(1);
    Parameters::setSubLatticePreset(false);
    Parameters::setNSpatial(N_SPATIAL);
    Parameters::setNTemporal(N_TEMPORAL);
    Parallel::Index::setNTot(N_SPATIAL, N_TEMPORAL);
    Parameters::setBeta(6.0);
    Parameters::setConfigStartNumber(0);
    Parameters::setUnitTesting(true);
    Parameters::setDebug(false);
    Parameters::setFilePath("");
    Parameters::setOutputFolder(withTrailingSlash(getOutputRoot()));
    Parameters::setInputFolder("");
    Parameters::setBatchName(batchName);

    Parallel::Communicator::initializeSubLattice();
    IO::FieldIO::init();
  }

  void prepareBatchDirectories(const std::string& batchName,
                               const std::string& observable = "")
  {
    Parallel::Communicator::setBarrierActive();

    if (Parallel::Communicator::getProcessRank() == 0)
    {
      const auto batchRoot = getBatchRoot(batchName);
      std::filesystem::remove_all(batchRoot);
      std::filesystem::create_directories(batchRoot / "field_configurations");

      if (!observable.empty())
      {
        std::filesystem::create_directories(batchRoot / "scalar_fields" /
                                            observable);
      }
    }

    Parallel::Communicator::setBarrierActive();
  }

  std::array<Lattice<SU3>, 4> allocateGaugeField()
  {
    std::array<Lattice<SU3>, 4> field;
    for (auto& muField : field)
    {
      muField.allocate(Parameters::getN());
    }
    return field;
  }

  void fillGaugeField(std::array<Lattice<SU3>, 4>& field)
  {
    const auto dims = Parameters::getN();

    for (unsigned int t = 0; t < dims[3]; ++t)
    {
      const unsigned int globalT = getGlobalCoordinate(3, t);
      for (unsigned int z = 0; z < dims[2]; ++z)
      {
        const unsigned int globalZ = getGlobalCoordinate(2, z);
        for (unsigned int y = 0; y < dims[1]; ++y)
        {
          const unsigned int globalY = getGlobalCoordinate(1, y);
          for (unsigned int x = 0; x < dims[0]; ++x)
          {
            const unsigned int globalX = getGlobalCoordinate(0, x);
            const auto site = Parallel::Index::getIndex(x, y, z, t);

            for (unsigned int mu = 0; mu < 4; ++mu)
            {
              field[mu][site] =
                encodeGaugeMatrix(globalX, globalY, globalZ, globalT, mu);
            }
          }
        }
      }
    }

    Parallel::Communicator::setBarrierActive();
  }

  Lattice<double> allocateScalarField()
  {
    Lattice<double> field;
    field.allocate(Parameters::getN());
    return field;
  }

  void fillScalarField(Lattice<double>& field)
  {
    const auto dims = Parameters::getN();

    for (unsigned int t = 0; t < dims[3]; ++t)
    {
      const unsigned int globalT = getGlobalCoordinate(3, t);
      for (unsigned int z = 0; z < dims[2]; ++z)
      {
        const unsigned int globalZ = getGlobalCoordinate(2, z);
        for (unsigned int y = 0; y < dims[1]; ++y)
        {
          const unsigned int globalY = getGlobalCoordinate(1, y);
          for (unsigned int x = 0; x < dims[0]; ++x)
          {
            const unsigned int globalX = getGlobalCoordinate(0, x);
            field[Parallel::Index::getIndex(x, y, z, t)] =
              encodeScalarValue(globalX, globalY, globalZ, globalT);
          }
        }
      }
    }

    Parallel::Communicator::setBarrierActive();
  }
}  // namespace

TEST_CASE("MPI gauge field write/read round-trip", TAGS + "[gauge]")
{
  if (!Parallel::ParallelParameters::active)
  {
    SUCCEED("Inactive MPI ranks do not participate in IO tests.");
    return;
  }

  const std::string batchName = "mpi_io_gauge_roundtrip";
  configureIOEnvironment(batchName);
  prepareBatchDirectories(batchName);

  auto before = allocateGaugeField();
  auto after = allocateGaugeField();
  fillGaugeField(before);

  for (auto& muField : after)
  {
    muField.zeros();
  }

  IO::FieldIO::writeFieldToFile(before.data(), 0);
  Parallel::Communicator::setBarrierActive();

  const auto outputFile = getBatchRoot(batchName) / "field_configurations" /
                          getConfigurationFilename(0);

  if (Parallel::Communicator::getProcessRank() == 0)
  {
    CAPTURE(outputFile.string());
    REQUIRE(std::filesystem::exists(outputFile));
    REQUIRE(std::filesystem::file_size(outputFile) ==
            static_cast<std::uintmax_t>(Parameters::getLatticeSize()) * 4 * 18 *
              sizeof(double));
  }

  Parameters::setInputFolder(
    withTrailingSlash(getBatchRoot(batchName) / "field_configurations"));
  IO::FieldIO::loadFieldConfiguration(getConfigurationFilename(0),
                                      after.data());

  const auto dims = Parameters::getN();
  for (unsigned int t = 0; t < dims[3]; ++t)
  {
    for (unsigned int z = 0; z < dims[2]; ++z)
    {
      for (unsigned int y = 0; y < dims[1]; ++y)
      {
        for (unsigned int x = 0; x < dims[0]; ++x)
        {
          const auto site = Parallel::Index::getIndex(x, y, z, t);
          for (unsigned int mu = 0; mu < 4; ++mu)
          {
            INFO("Mismatch at rank " << Parallel::Communicator::getProcessRank()
                                     << ", site (" << x << ", " << y << ", "
                                     << z << ", " << t << "), mu=" << mu);
            REQUIRE(before[mu][site] == after[mu][site]);
          }
        }
      }
    }
  }

  Parallel::Communicator::setBarrierActive();
}

TEST_CASE("MPI scalar field write preserves global layout", TAGS + "[scalar]")
{
  if (!Parallel::ParallelParameters::active)
  {
    SUCCEED("Inactive MPI ranks do not participate in IO tests.");
    return;
  }

  const std::string batchName = "mpi_io_scalar_layout";
  const std::string observable = "io_scalar_layout";

  configureIOEnvironment(batchName);
  prepareBatchDirectories(batchName, observable);

  auto field = allocateScalarField();
  fillScalarField(field);

  IO::FieldIO::writeDoublesFieldToFile(field, 0, observable);
  Parallel::Communicator::setBarrierActive();

  const auto outputFile = getBatchRoot(batchName) / "scalar_fields" /
                          observable / getConfigurationFilename(0);

  if (Parallel::Communicator::getProcessRank() == 0)
  {
    CAPTURE(outputFile.string());
    REQUIRE(std::filesystem::exists(outputFile));
    REQUIRE(std::filesystem::file_size(outputFile) ==
            static_cast<std::uintmax_t>(Parameters::getLatticeSize()) *
              sizeof(double));

    std::ifstream input(outputFile, std::ios::binary);
    REQUIRE(input.is_open());

    std::vector<double> rawField(Parameters::getLatticeSize(), 0.0);
    input.read(reinterpret_cast<char*>(rawField.data()),
               static_cast<std::streamsize>(rawField.size() * sizeof(double)));
    REQUIRE(input.gcount() ==
            static_cast<std::streamsize>(rawField.size() * sizeof(double)));

    for (unsigned int t = 0; t < Parameters::getNTemporal(); ++t)
    {
      for (unsigned int z = 0; z < Parameters::getNSpatial(); ++z)
      {
        for (unsigned int y = 0; y < Parameters::getNSpatial(); ++y)
        {
          for (unsigned int x = 0; x < Parameters::getNSpatial(); ++x)
          {
            const auto globalIndex = getGlobalScalarIndex(x, y, z, t);
            INFO("Mismatch at global site (" << x << ", " << y << ", " << z
                                             << ", " << t << ")");
            REQUIRE(rawField[globalIndex] == encodeScalarValue(x, y, z, t));
          }
        }
      }
    }
  }

  Parallel::Communicator::setBarrierActive();
}
