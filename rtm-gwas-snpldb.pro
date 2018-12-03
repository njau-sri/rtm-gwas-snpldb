TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

QMAKE_CXXFLAGS += -std=c++11 -fopenmp
QMAKE_LFLAGS += -static -fopenmp

SOURCES += \
        main.cpp \
    rtm_gwas_snpldb.cpp \
    cmdline.cpp \
    util.cpp \
    vcf.cpp \
    block_gabriel.cpp \
    haplotype.cpp

HEADERS += \
    cmdline.h \
    split.h \
    util.h \
    vcf.h \
    block_gabriel.h \
    haplotype.h \
    version.h
