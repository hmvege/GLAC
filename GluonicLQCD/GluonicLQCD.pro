TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    actions/action.cpp \
    actions/wilsongaugeaction.cpp \
    observables/correlator.cpp \
    observables/plaquette.cpp \
    observables/topologicalcharge.cpp \
    observables/energydensity.cpp \
    observables/tools/observablestorer.cpp \
    math/functions.cpp \
    math/links.cpp \
    math/complex.cpp \
    math/matrices/su2.cpp \
    math/matrices/su3.cpp \
    math/matrices/su3matrixgenerator.cpp \
    math/exponentiation/su3exp.cpp \
    math/exponentiation/expluscher.cpp \
    math/exponentiation/taylor2exp.cpp \
    math/exponentiation/taylor4exp.cpp \
    parallelization/neighbours.cpp \
    parallelization/neighbourlist.cpp \
    parallelization/index.cpp \
    parallelization/communicator.cpp \
    tests/testsuite.cpp \
    tests/unittests.cpp \
    tests/performancetests.cpp \
    io/observablesio.cpp \
    io/fieldio.cpp \
    config/parameters.cpp \
    flow/flow.cpp \
    config/configloader.cpp \
    config/sysprint.cpp \
    observables/mastersampler.cpp

HEADERS += \
    system.h \
    actions/action.h \
    actions/wilsongaugeaction.h \
    observables/correlator.h \
    observables/plaquette.h \
    observables/topologicalcharge.h \
    observables/energydensity.h \
    observables/tools/observablestorer.h \
    math/functions.h \
    math/links.h \
    math/complex.h \
    math/matrices/su2.h \
    math/matrices/su3.h \
    math/matrices/su3matrixgenerator.h \
    math/exponentiation/su3exp.h \
    math/exponentiation/expluscher.h \
    math/exponentiation/taylor2exp.h \
    math/exponentiation/taylor4exp.h \
    math/latticemath.h \
    parallelization/neighbours.h \
    parallelization/neighbourlist.h \
    parallelization/index.h \
    parallelization/communicator.h \
    tests/testsuite.h \
    tests/unittests.h \
    tests/performancetests.h \
    io/observablesio.h \
    io/fieldio.h \
    config/parameters.h \
    flow/flow.h \
    config/configloader.h \
    lib/json.hpp \
    config/sysprint.h \
    math/flowexpfunctions.h \
    actions/actions.h \
    parallelization/parallel.h \
    math/lattice.h \
    observables/mastersampler.h \
    observables/observables.h

#QMAKE_PRE_LINK = hpclink #??

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = hpclink mpicc

QMAKE_CFLAGS += -O3 -std=c++11 $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# Removes flags
QMAKE_CFLAGS -= -O2
QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O2
