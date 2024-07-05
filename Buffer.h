//
// Created by 赵展 on 2021/2/22.
//
#ifndef BAM_BUFFER_H
#define BAM_BUFFER_H
#include <mutex>
#include <thread>

#include <iostream>
#include <fstream>

using namespace  std;
class BufferConfig{
public:
    BufferConfig();
    BufferConfig(int Buffer_number, int consumerpack_number, int Maxn);
    ~BufferConfig();
public:
    int Buffer_number;
    int write_number;
    int consumerpack_number;
    int Maxn;
};

class Buffer {
public:
    Buffer();
    Buffer(BufferConfig *config,ofstream *fout);
    pair<char *,int> getCap();
    void initoutput(int id,int pos);
    void output();
    bool is_complete();
    void complete_thread();

public:
    BufferConfig *config;
    ofstream *fout;
    mutex mtx;
    char **buffer;
    int *write;
    int *pos;
    int write_bg;
    int write_ed;
    int *capacity;
    int cap_bg;
    int cap_ed;
};


#endif //BAM_BUFFER_H
