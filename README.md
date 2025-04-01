## ![logo](./docs/source/_static/RabbitBAM.jpg)

An ultra-fast I/O framework for BAM files.

## Dependancy

- gcc 8.5.0 or newer
- [htslib](https://github.com/samtools/htslib) 1.15 or newer
- [libdeflate](https://github.com/ebiggers/libdeflate) 1.12 or newer

## Installation

```bash
git clone https://github.com/RabbitBio/RabbitBAM
cd RabbitBAM
bash configure.sh <path-to-htslib-installation-directory> <path-to-libdeflate-installation-directory>
source env.sh
make clean && make
```

## Benchmark

refer to `./rabbitbam -h`

## Usage

To use RabbitBAM in other software, you need to first compile RabbitBAM according to the [Installation](#installation) instructions, then link it to your new software and include the appropriate header files. 

For detailed usage instructions, please refer to: https://rabbitbam.readthedocs.io/en/latest.

Here is an example of file reading and writing (see [RabbitBAM-QC](https://github.com/RabbitBio/RabbitBAM-QC) and [RabbitBAM-SORT](https://github.com/RabbitBio/RabbitBAM-SORT) for more examples):

### Source:

```c++
#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <thread>

#include "BamTools.h"
#include "BamReader.h"
#include "BamWriter.h"

int main(int argc, char* argv[]) {
    std::string input_file;
    std::string output_file;
    int n_thread = 1;

    int opt;
    while ((opt = getopt(argc, argv, "i:t:o:")) != -1) {
        switch (opt) {
            case 'i':
                input_file = optarg;
                break;
            case 't':
                n_thread = std::atoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " -i input_file -t thread_num -o output_file\n";
                return EXIT_FAILURE;
        }
    }

    if (input_file.empty() || output_file.empty() || n_thread <= 0) {
        std::cerr << "Invalid arguments. Usage: " << argv[0] << " -i input_file -t thread_num -o output_file\n";
        return EXIT_FAILURE;
    }

    std::cout << "Input file: " << input_file << "\n";
    std::cout << "Total thread number: " << n_thread << "\n";
    std::cout << "Output file: " << output_file << "\n";

    int compress_level = 6;
    int block_size = 256;
    bool single_parser = true;

    int n_thread_read = n_thread / 5;
    int n_thread_write = n_thread - n_thread_read;
    if(n_thread_read <= 0) n_thread_read = 1;
    if(n_thread_write <= 0) n_thread_write = 1;

    BamReader *reader = new BamReader(input_file, n_thread_read, single_parser);
    BamWriter *writer = new BamWriter(output_file, reader->getHeader(), n_thread_write, compress_level, block_size, single_parser);

    bam1_t *b;
    if ((b = bam_init1()) == NULL) {
        fprintf(stderr, "[E::%s] Out of memory allocating BAM struct.\n", __func__);
    }

    long long num = 0;
    while (reader->getBam1_t(b)) {
        num++;
        writer->write(b);
    }

    writer->over();
    printf("num %lld\n", num);

    return 0;
}
```



### Compile:

```bash
export HTSLIB_INSTALL_PATH=<path-to-htslib-installation-directory>
export RABBITBAM_INSTALL_PATH=<path-to-RabbitBAM-installation-directory>
export LIBDEFLATE_INSTALL_PATH=<path-to-libdeflate-installation-directory>
export LD_LIBRARY_PATH=$RABBITBAM_INSTALL_PATH:$LD_LIBRARY_PATH

g++ -o main main.cpp -I$HTSLIB_INSTALL_PATH/include -I$RABBITBAM_INSTALL_PATH/htslib -I$RABBITBAM_INSTALL_PATH -I$LIBDEFLATE_INSTALL_PATH/include -L$HTSLIB_INSTALL_PATH/lib -L$RABBITBAM_INSTALL_PATH -L$LIBDEFLATE_INSTALL_PATH/lib -lhts -lz -fopenmp -lpthread -lrabbitbamtools -lrabbitbamread -lrabbitbamwrite
```

### Run:

```bash
time ./main -i in.bam -o out.bam -t 40
```

Note: This is for the single parser case. If you need multiple parsers, please refer to the example `parallel_test` in `block_mul.cpp`.
