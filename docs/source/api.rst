API Reference
=============

.. image:: _static/pipeline-rbam.png
   :alt: pipeline-rbam
   :align: center
   :width: 100%

BAM Reading
-----------

**Constructors**

.. code-block:: cpp

   BamReader(std::string file_name, int n_thread = 8, bool is_tgs = false);

Creates a BAM reader for the specified file.

- `n_thread`: Number of threads used for decompression and parsing (default: 8)  
- `is_tgs`: Whether the input BAM file is from third-generation sequencing

.. code-block:: cpp

   BamReader(std::string file_name, int read_block, int compress_block, int compress_complete_block, int n_thread = 8, bool is_tgs = false);

Creates a BAM reader with fine-grained control over memory usage.

- `read_block`, `compress_block`, `compress_complete_block`: Maximum number of blocks in Queue1, Queue2, and Queue3 respectively  
- Useful for tuning memory and performance in large-scale pipelines

**Functions**

.. code-block:: cpp

   bool getBam1_t(bam1_t *b);

Retrieves the next parsed BAM record.  
Returns `true` if successful, or `false` if end-of-file is reached.

BAM Writing
-----------

**Constructors**

.. code-block:: cpp

   BamWriter(std::string file_name, int threadNumber = 1, int level = 6, int BufferSize = 200, bool is_tgs = false);

Creates a BAM writer that writes directly to a file.

- `threadNumber`: Number of compression threads  
- `level`: Compression level (default: 6)  
- `BufferSize`: Maximum number of blocks in Queue4 and Queue5
- `is_tgs`: Whether the BAM is from third-generation sequencing

.. code-block:: cpp

   BamWriter(std::string file_name, sam_hdr_t *hdr, int threadNumber = 1, int level = 6, int BufferSize = 200, bool is_tgs = false);

Initializes a writer with both file path and BAM header.

**Functions**

.. code-block:: cpp

   void hdr_write(sam_hdr_t *hdr);

Writes the BAM header to the output stream.  
Must be called before writing records if not set in the constructor.

.. code-block:: cpp

   void bam_write(bam1_t *b);

Writes a single BAM record to the output.

.. code-block:: cpp

   void over();

Finalizes and closes the output stream.  
Should be called after all records are written.


Usage Example
-------------

Below is a C++ demo program that illustrates how to use the above APIs.
**Source:**

.. code-block:: cpp

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
               case 'i': input_file = optarg; break;
               case 't': n_thread = std::atoi(optarg); break;
               case 'o': output_file = optarg; break;
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

**Compile:**

.. code-block:: bash

   export HTSLIB_INSTALL_PATH=<path-to-htslib-installation-directory>
   export RABBITBAM_INSTALL_PATH=<path-to-RabbitBAM-installation-directory>
   export LIBDEFLATE_INSTALL_PATH=<path-to-libdeflate-installation-directory>
   export LD_LIBRARY_PATH=$RABBITBAM_INSTALL_PATH:$LD_LIBRARY_PATH

   g++ -o main main.cpp \
     -I$HTSLIB_INSTALL_PATH/include \
     -I$RABBITBAM_INSTALL_PATH/htslib \
     -I$RABBITBAM_INSTALL_PATH \
     -I$LIBDEFLATE_INSTALL_PATH/include \
     -L$HTSLIB_INSTALL_PATH/lib \
     -L$RABBITBAM_INSTALL_PATH \
     -L$LIBDEFLATE_INSTALL_PATH/lib \
     -lhts -lz -fopenmp -lpthread \
     -lrabbitbamtools -lrabbitbamread -lrabbitbamwrite

**Run:**

.. code-block:: bash

   time ./main -i in.bam -o out.bam -t 40
