TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
#CONFIG += optimize_full # OK?

SOURCES += main.cpp \
    system.cpp \
    actions/action.cpp \
    actions/wilsongaugeaction.cpp \
    observables/correlator.cpp \
    observables/plaquette.cpp \
    observables/topologicalcharge.cpp \
    observables/clover.cpp \
    observables/energydensity.cpp \
    observables/observablesampler.cpp \
    math/functions.cpp \
    math/links.cpp \
    math/complex.cpp \
    math/matrices/su2.cpp \
    math/matrices/su3.cpp \
    math/matrices/su3matrixgenerator.cpp \
    math/exponentiation/su3exp.cpp \
    parallelization/neighbours.cpp \
    parallelization/neighbourlist.cpp \
    parallelization/index.cpp \
    flow/flow.cpp \
    unittests/testsuite.cpp \
    unittests/unittests.cpp \
    math/exponentiation/taylorexp.cpp \
    math/exponentiation/expmorningstar.cpp \
    math/exponentiation/expluscher.cpp
HEADERS += \
    system.h \
    actions/action.h \
    actions/wilsongaugeaction.h \
    observables/correlator.h \
    observables/plaquette.h \
    observables/topologicalcharge.h \
    observables/clover.h \
    observables/energydensity.h \
    observables/observablesampler.h \
    math/functions.h \
    math/links.h \
    math/complex.h \
    math/matrices/su2.h \
    math/matrices/su3.h \
    math/matrices/su3matrixgenerator.h \
    math/exponentiation/su3exp.h \
    parallelization/neighbours.h \
    parallelization/neighbourlist.h \
    parallelization/index.h \
    flow/flow.h \
    unittests/testsuite.h \
    unittests/unittests.h \
    math/exponentiation/taylorexp.h \
    math/exponentiation/expmorningstar.h \
    math/exponentiation/expluscher.h

#LIBS += -llapack -lblas -larmadillo

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += -O3 -std=c++11 $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += -O3 -std=c++11 $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

# Removes flags
QMAKE_CFLAGS -= -O2
QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS_RELEASE -= -O2


# Following to make openmp usable on linux
#QMAKE_LFLAGS += -fopenmp

# Following to make openmp usable on mac
#QMAKE_LDFLAGS += -L/usr/local/opt/llvm/lib -Wl,-rpath,/usr/local/opt/llvm/lib

# Following used to make armadillo usable on mac
#LIBS += -L/usr/local/lib -larmadillo
#INCLUDEPATH += /usr/local/include

#INCLUDEPATH += -I/usr/local/include
#INCLUDEPATH += -L/usr/local/lib
#compileCommand-I/usr/local/include -L/usr/local/lib -llapack -lblas -larmadillo
