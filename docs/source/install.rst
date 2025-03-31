Installation
============

Dependencies
------------

- gcc 8.5.0 or newer
- `htslib 1.15 or newer <https://github.com/samtools/htslib>`_
- `libdeflate 1.12 or newer <https://github.com/ebiggers/libdeflate>`_

Compilation
------------

.. code-block:: bash

   git clone https://github.com/RabbitBio/RabbitBAM
   cd RabbitBAM
   bash configure.sh <path-to-htslib-installation-directory> <path-to-libdeflate-installation-directory>
   source env.sh
   make clean && make


Test
---------

Run the following command to check whether RabbitBAM was installed successfully:

.. code-block:: bash

   ./rabbitbam -h



