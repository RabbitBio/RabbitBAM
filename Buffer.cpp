#include "Buffer.h"

using namespace std;

BufferConfig::BufferConfig() {};

BufferConfig::BufferConfig(int Buffer_number, int consumerpack_number, int Maxn) {
    this->Buffer_number = Buffer_number;
    this->write_number = Buffer_number + 50;
    this->consumerpack_number = consumerpack_number;
    this->Maxn = Maxn;
}

BufferConfig::~BufferConfig() {}

Buffer::Buffer(BufferConfig *config, ofstream *fout) {
    this->fout = fout;
    this->config = config;
    this->buffer = new char *[this->config->Buffer_number];
    this->pos = new int[this->config->Buffer_number];
    for (int i = 0; i < this->config->Buffer_number; i++)
        this->buffer[i] = new char[this->config->Maxn + 10000];
    this->write = new int[this->config->write_number];
    this->write_bg = 0;
    this->write_ed = 0;
    this->capacity = new int[this->config->write_number];
    this->cap_bg = 0;
    this->cap_ed = this->config->Buffer_number;
    for (int i = this->cap_bg; i < this->cap_ed; i++) this->capacity[i] = i;
}

pair<char *, int> Buffer::getCap() {

    mtx.lock();
    while (cap_ed == cap_bg) {
        mtx.unlock();
        this_thread::sleep_for(chrono::milliseconds(1));
        mtx.lock();
    }
    int id = capacity[cap_bg];
    cap_bg = (cap_bg + 1) % config->write_number;
    mtx.unlock();
    return pair<char *, int>(buffer[id], id);
}

void Buffer::initoutput(int id, int pos) {
    mtx.lock();
    write[write_ed] = id;
    this->pos[id] = pos;
    write_ed = (write_ed + 1) % this->config->write_number;
    mtx.unlock();
}

void Buffer::output() {
    mtx.lock();
    while (write_bg != write_ed) {
        this->fout->write(buffer[write[write_bg]], pos[write[write_bg]]);
        capacity[cap_ed] = write[write_bg];
        cap_ed = (cap_ed + 1) % config->write_number;
        write_bg = (write_bg + 1) % config->write_number;
    }
    mtx.unlock();
}

bool Buffer::is_complete() {
    return (!config->consumerpack_number && write_bg == write_ed);
}

void Buffer::complete_thread() {
    mtx.lock();
    config->consumerpack_number -= 1;
    mtx.unlock();
}