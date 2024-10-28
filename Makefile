CXX = g++

CPPFLAGS = -g -std=c++11 -O3 -ffast-math -DTIMING 
BUILD_PATH = $(shell pwd)

OBJECT = $(BUILD_PATH)/tools.o $(BUILD_PATH)/read.o $(BUILD_PATH)/write.o $(BUILD_PATH)/status.o
SO_OBJECT = $(BUILD_PATH)/librabbitbamtools.so $(BUILD_PATH)/librabbitbamread.so $(BUILD_PATH)/librabbitbamwrite.so
TARGET = $(BUILD_PATH)/rabbitbam

LIB =
INCLUDE = ./htslib 
HTSLIB_INSTALL_PATH = /home/user_home/ylf/someGit/bam_deps/htslib-1.20-install
LIBDEFLATE_INSTALL_PATH = /home/user_home/ylf/someGit/bam_deps/libdeflate-1.20-install

LIB += $(HTSLIB_INSTALL_PATH)/lib
LIB += $(HTSLIB_INSTALL_PATH)/lib64
LIB += $(LIBDEFLATE_INSTALL_PATH)/lib
LIB += $(LIBDEFLATE_INSTALL_PATH)/lib64

INCLUDE += $(HTSLIB_INSTALL_PATH)/include
INCLUDE += $(LIBDEFLATE_INSTALL_PATH)/include


SHARE = -fPIC -shared
#SHARE = -static

INC_FLAGS = $(addprefix -I,$(INCLUDE)) 
LIB_FLAGS = $(addprefix -L,$(LIB)) -fopenmp -lpthread -lz -lhts -ldeflate



$(BUILD_PATH)/rabbitbam: $(OBJECT) block_mul.cpp
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) block_mul.cpp -o $(BUILD_PATH)/rabbitbam $(OBJECT) $(LIB_FLAGS)

$(BUILD_PATH)/status.o: $(BUILD_PATH)/tools.o Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h $(SHARE) -o $(BUILD_PATH)/status.o $(BUILD_PATH)/tools.o $(LIB_FLAGS)

$(BUILD_PATH)/write.o : $(BUILD_PATH)/tools.o BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp  $(BUILD_PATH)/tools.o $(SHARE) -o $(BUILD_PATH)/write.o $(LIB_FLAGS)
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp  $(BUILD_PATH)/librabbitbamtools.so $(SHARE) -o $(BUILD_PATH)/librabbitbamwrite.so $(LIB_FLAGS)

$(BUILD_PATH)/read.o : $(BUILD_PATH)/tools.o BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h $(SHARE) -o $(BUILD_PATH)/read.o  $(BUILD_PATH)/tools.o $(LIB_FLAGS)
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h $(SHARE) -o $(BUILD_PATH)/librabbitbamread.so  $(BUILD_PATH)/librabbitbamtools.so $(LIB_FLAGS)

$(BUILD_PATH)/tools.o :
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamTools.cpp BamTools.h Buffer.cpp Buffer.h config.h $(SHARE) -o $(BUILD_PATH)/tools.o $(LIB_FLAGS)
	$(CXX) $(CPPFLAGS) $(INC_FLAGS) BamTools.cpp BamTools.h Buffer.cpp Buffer.h config.h $(SHARE) -o $(BUILD_PATH)/librabbitbamtools.so $(LIB_FLAGS)

.PHONY: clean

clean:
	rm -rf $(TARGET) $(OBJECT) $(SO_OBJECT)

