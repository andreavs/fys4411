TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    atom.cpp \
    wavefunction.cpp \
    vmcsolver.cpp \
    vmcsetup.cpp \
    random.cpp \
    slaterdeterminant.cpp \
    orbital.cpp \
    jastrow.cpp \
    zignor.cpp \
    zigrandom.cpp \
    blocking.cpp \
    moleculeorbital.cpp

HEADERS += \
    atom.h \
    wavefunction.h \
    vmcsolver.h \
    vmcsetup.h \
    vmc.h \
    random.h \
    slaterdeterminant.h \
    orbital.h \
    jastrow.h \
    zigrandom.h \
    zignor.h \
    blocking.h \
    moleculeorbital.h

LIBS += -larmadillo \
    -fopenmp \
    -llapack \
    -lblas

QMAKE_CXXFLAGS += -fopenmp -std=c++0x
QMAKE_CXXFLAGS_RELEASE = $$QMAKE_CXXFLAGS
QMAKE_CXXFLAGS_DEBUG = $$QMAKE_CXXFLAGS

release {
DEFINES += ARMA_NO_DEBUG
QMAKE_LFLAGS -= -O1
QMAKE_LFLAGS += -O3
QMAKE_LFLAGS_RELEASE -= -O1
QMAKE_LFLAGS_RELEASE += -O3
QMAKE_CXXFLAGS -= -O2
QMAKE_CXXFLAGS += -O3
QMAKE_CXXFLAGS_RELEASE -= -O2
QMAKE_CXXFLAGS_RELEASE += -O3
}

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
