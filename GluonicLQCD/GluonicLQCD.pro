TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    actions/action.cpp \
    actions/wilsongaugeaction.cpp \
    actions/luscheraction.cpp \
    observables/correlator.cpp \
    observables/plaquette.cpp \
    observables/topologicalcharge.cpp \
    observables/energydensity.cpp \
    observables/tools/observablestorer.cpp \
    observables/mastersampler.cpp \
    observables/mastersamplertopcxyz.cpp \
    observables/latticeactionchargedensity.cpp \
    observables/supersampler.cpp \
    math/complex.cpp \
    math/matrices/su2.cpp \
    math/matrices/su3.cpp \
    math/matrices/su3matrixgenerator.cpp \
    math/exponentiation/su3exp.cpp \
    math/exponentiation/expluscher.cpp \
    math/exponentiation/taylor2exp.cpp \
    math/exponentiation/taylor4exp.cpp \
    math/exponentiation/taylorexp.cpp \
    parallelization/neighbours.cpp \
    parallelization/neighbourlist.cpp \
    parallelization/index.cpp \
    parallelization/communicator.cpp \
    parallelization/parallelparameters.cpp \
    io/observablesio.cpp \
    io/fieldio.cpp \
    tests/testsuite.cpp \
    tests/performancetests.cpp \
    config/parameters.cpp \
    config/configloader.cpp \
    config/sysprint.cpp \
    flow/flow.cpp

HEADERS += \
    system.h \
    actions/action.h \
    actions/wilsongaugeaction.h \
    actions/actions.h \
    actions/luscheraction.h \
    observables/correlator.h \
    observables/plaquette.h \
    observables/topologicalcharge.h \
    observables/energydensity.h \
    observables/tools/observablestorer.h \
    observables/mastersampler.h \
    observables/observables.h \
    observables/mastersamplertopcxyz.h \
    observables/latticeactionchargedensity.h \
    observables/supersampler.h \
    math/functions.h \
    math/complex.h \
    math/matrices/su2.h \
    math/matrices/su3.h \
    math/matrices/su3matrixgenerator.h \
    math/exponentiation/su3exp.h \
    math/exponentiation/expluscher.h \
    math/exponentiation/taylor2exp.h \
    math/exponentiation/taylor4exp.h \
    math/latticemath.h \
    math/flowexpfunctions.h \
    math/lattice.h \
    math/exponentiation/taylorexp.h \
    parallelization/neighbours.h \
    parallelization/neighbourlist.h \
    parallelization/index.h \
    parallelization/communicator.h \
    parallelization/parallel.h \
    parallelization/parallelparameters.h \
    io/observablesio.h \
    io/fieldio.h \
    tests/testsuite.h \
    tests/performancetests.h \
    tests/test.h \
    config/parameters.h \
    config/configloader.h \
    config/sysprint.h \
    flow/flow.h \
    lib/json.hpp

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# Adds processor specific optimizations
QMAKE_CFLAGS += -march=native
QMAKE_CXXFLAGS += -march=native
QMAKE_CXXFLAGS_RELEASE += -march=native

# Adds O3 optimizations
QMAKE_CFLAGS += -O3
QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS_RELEASE += -O3

# Forces C++11
QMAKE_CFLAGS += -std=c++11
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS_RELEASE += -std=c++11

# Removes flags
QMAKE_CFLAGS -= -O2
QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O2
