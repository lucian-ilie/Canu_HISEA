#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET    := hisea
SOURCES  := ./src/Util/Globals.cpp \
            ./src/Util/Options.cpp \
            ./src/Util/SequenceReader.cpp \
            ./src/Util/MemoryTracker.cpp \
            ./src/ThirdParty/gzstream/gzstream.C \
            ./src/Algorithm/KmerMap.cpp \
            ./src/Algorithm/Alignment.cpp \
            ./src/Hisea.cpp

SRC_INCDIRS  := ./src/Util ./src/ThirdParty/gzstream ./src/Algorithm
SRC_DEFS  := UINT64 _GLIBCXX_PARALLEL

TGT_LDFLAGS := -lz
SRC_CXXFLAGS  :=  -std=c++0x -fopenmp
TGT_CXXFLAGS  :=  -std=c++0x -fopenmp 

SUBMAKEFILES :=

