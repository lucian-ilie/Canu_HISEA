#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)/bin
endif

TARGET   := correctOverlaps
SOURCES  := correctOverlaps.C \
            correctOverlaps-Correct_Frags.C \
            correctOverlaps-Read_Olaps.C \
            correctOverlaps-Redo_Olaps.C \
            correctOverlaps-Prefix_Edit_Distance.C

SRC_INCDIRS  := .. ../AS_UTL ../stores ../overlapInCore/liboverlap

TGT_LDFLAGS := -L${TARGET_DIR}
TGT_LDLIBS  := -lcanu
TGT_PREREQS := libcanu.a

SUBMAKEFILES :=
