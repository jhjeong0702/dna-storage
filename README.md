# Cooperative Sequence Clustering and Decoding for DNA Storage System with Fountain Codes

Here, codes and data files used for my DNA storage study is stored.

## Encoding
To encode data, you may use LT_code.m in the encoder directory. This is a Matlab file, and we used MATLAB R2019a for this study.
The encoding strategy is from "DNA Fountain enables a robust and eficient storage architecture", Erlich and Zielinski, Science 2017.
For the constraints of GC-content and homopolymer-run, you can change them by fixing the value of Max_run_length, Min_GC_content, and MAX_GC_content in the LT_code.m file.
- For our non-constrained pool, we used Max_run_length = 152, Min_GC_content = 0, and MAX_GC_content = 1.
- For our constrained pool, we used Max_run_length = 3, Min_GC_content = 0.45, and MAX_GC_content = 0.55.
The result is stored in the original_files directory.

## Sequencing
The sequenced FASTQ files are stored in the fastq_files directory. Due to the limit of the data size of github repository, I uploaded them with partitioned compression. Forward reads and reverse reads of the both non-constrained pool and constrained pool are stored.

## Decoding
I devided the decoding stage into clustering and LT-decoding stages. For the clustering, the code clustering_stage.cpp is uploaded in the decoder directory. For the compiler, I used Microsoft Visual Studio 2015, and the compile result is uploaded as clustering.exe.
 - To run the clustering.exe, you may need PEAR algorithm from "PEAR: a fast and accurate Illumina Paired-End reAd mergeR", Zhang et al, Bioinformatics 2014. We recommend to use the max_length and min_length to be same (which we used 152), and the number of N option to be 0 when merging the forward and reverse reads.
 - Then, you may run the clustering_example.bat file for the clustering. For the non-constrained pool, each encoded oligo sequence had Hamming distances of at least 81 to every other sequence, so we discarded the reads that has minimum distance of 3-79 from other clusters. However, for the constrained pool, each encoded oligo sequence had Hamming distances of at least 82 to every other sequence, except for only 2 sequences having Hamming distance of 17. Therefore, we discarded the reads that has minimum distance of 3-15 and 19-80 from other clusters.
 - If you run the clustering_example.bat file, you may type in the input Fastq file name, and the result files (raw_cnt%d.txt, raw_read_symbol_cnt_discarded%d.txt, raw_read_symbol_cnt_ocNum%d.txt, raw_read_symbol_cnt%d.txt, and raw_read%d.txt) will come out. From them, only raw_cnt%d.txt and raw_read%d.txt are used for the LT decoding.
After running the batch file, you may run the RS_correction_decoding.m code in the decoder directory. This is also a Matlab file, and we used MATLAB R2019a for this study. You may change the start_index and end_index in the code depending on the number of Fastq files used for the decoding. The Decoding_final_result.txt will be generated as a result.

