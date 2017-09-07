TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    system.cpp \
    actions/action.cpp \
    functions.cpp \
    links.cpp \
    complex.cpp \
    unittests.cpp \
    correlators/correlator.cpp \
    correlators/plaquette.cpp \
    actions/wilsongaugeaction.cpp \
    matrices/su2.cpp \
    matrices/su3.cpp \
    matrices/su3matrixgenerator.cpp

HEADERS += \
    system.h \
    actions/action.h \
    functions.h \
    links.h \
    complex.h \
    unittests.h \
    correlators/correlator.h \
    correlators/plaquette.h \
    actions/wilsongaugeaction.h \
    matrices/su2.h \
    matrices/su3.h \
    matrices/su3matrixgenerator.h

#LIBS += -llapack -lblas -larmadillo

## MPI Settings
#QMAKE_CXX = mpicxx
#QMAKE_CXX_RELEASE = $$QMAKE_CXX
#QMAKE_CXX_DEBUG = $$QMAKE_CXX
#QMAKE_LINK = $$QMAKE_CXX
#QMAKE_CC = mpicc

#QMAKE_CFLAGS += $$system(mpicc --showme:compile)
#QMAKE_LFLAGS += $$system(mpicxx --showme:link)
#QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
#QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

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
