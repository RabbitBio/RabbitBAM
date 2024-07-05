CXX = g++
CPPFLAGS = -g -std=c++11 -lpthread -ldeflate -lz -O3 -ffast-math -flto=full -DTIMING -lhts
BUILD_PATH = $(shell pwd)
OBJECT = $(BUILD_PATH)/tools.o $(BUILD_PATH)/read.o $(BUILD_PATH)/write.o $(BUILD_PATH)/status.o
SO_OBJECT = $(BUILD_PATH)/librabbitbamtools.so $(BUILD_PATH)/librabbitbamread.so $(BUILD_PATH)/librabbitbamwrite.so
TARGET = $(BUILD_PATH)/rabbitbam
INCLUDE = ./htslib
HTSLIB_INSTALL_PATH = /home/user_home/ylf/someGit/htslib_1.20_install/
LIB = $(HTSLIB_INSTALL_PATH)/lib
SHARE = -fPIC -shared

$(TARGET): $(OBJECT) block_mul.cpp
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) block_mul.cpp -o $@ $(OBJECT)

$(BUILD_PATH)/status.o: Overrepresent.cpp Overrepresent.h Duplicate.cpp Duplicate.h BamStatus.cpp BamStatus.h $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) Overrepresent.cpp Duplicate.cpp BamStatus.cpp $(SHARE) -o $@

$(BUILD_PATH)/write.o: BamWriter.cpp BamWriter.h BamWriteCompress.cpp BamWriteCompress.h $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamWriter.cpp BamWriteCompress.cpp $(SHARE) -o $@

$(BUILD_PATH)/read.o: BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h BamRead.cpp BamRead.h BamReader.cpp BamReader.h $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamCompleteBlock.cpp BamCompress.cpp BamBlock.cpp BamRead.cpp BamReader.cpp $(SHARE) -o $@

$(BUILD_PATH)/tools.o: BamTools.cpp BamTools.h Buffer.cpp Buffer.h config.h
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamTools.cpp Buffer.cpp $(SHARE) -o $@

# Optional: Create shared libraries
$(BUILD_PATH)/librabbitbamtools.so: $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) -o $@ $(SHARE) $(BUILD_PATH)/tools.o

$(BUILD_PATH)/librabbitbamread.so: $(BUILD_PATH)/read.o $(BUILD_PATH)/librabbitbamtools.so
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) -o $@ $(SHARE) $(BUILD_PATH)/read.o $(BUILD_PATH)/librabbitbamtools.so

$(BUILD_PATH)/librabbitbamwrite.so: $(BUILD_PATH)/write.o $(BUILD_PATH)/librabbitbamtools.so
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) -o $@ $(SHARE) $(BUILD_PATH)/write.o $(BUILD_PATH)/librabbitbamtools.so

.PHONY: clean

clean:
	rm -rf $(TARGET) $(OBJECT) $(SO_OBJECT)
	rm -rf ./output.fq
	rm -rf ./log.txt
	rm -rf ./BamStatus.html

