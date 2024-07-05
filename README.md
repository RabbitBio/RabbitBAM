# RabbitBAM

An ultra-fast I/O framework for BAM files.

## Dependancy

- gcc 8.5.0 or newer
- [htslib](https://github.com/samtools/htslib) 1.15 or newer
- [libdeflate](https://github.com/ebiggers/libdeflate) 1.12 or newer

## Installation

```bash
git clone https://github.com/RabbitBio/RabbitBAM
cd RabbitBAM
bash configure.sh <path-to-htslib-installation-directory>
source env.sh
make clean && make
```

## Usage

refer to `./RabbitBAM -h`