# NGS_pre-process
Tools for pre-processing of short-read NGS data

These tools were used for preprocessing of NGS reads for de novo assembly of the Oryza officinalis genome.


- Removal of duplicated read pairs was performed with an in-house C++ program.
remDupFastq.zip

- The low-quality reads were filtered out by an empirically optimized custom C++ program.
qvFiltFastq.zip

- Low frequency k-mers were probably derived from base call error and the reads containing them were filtered out with a custom C++ program.
kmerPurge_pair_unorderedmap.cpp

BoostC++ Libraries are required for remDupFastq and kmerPurge_pair_unorderedmap.cpp.
http://www.boost.org/

You can compile kmerPurge_pair_unorderedmap.cpp as below.

$ g++ kmerPurge_pair_unorderedmap.cpp -o kmerPurge_pair_unorderedmap -lpthread -O3 
