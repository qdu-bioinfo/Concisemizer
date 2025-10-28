This is a seed extraction algorithm. It builds upon the minimizer and conducts further seed screening. This project includes the following documents:

1.minimap2-concisemizer.zip is a modified version based on minimap2. The minimizer module in the original version has been changed to concisemizer. All modifications to the original version are marked as "Concisemizer_edit". For the installation and operation methods of the tool, please refer to minimap2(https://github.com/lh3/minimap2).

2.The files "concisemizer_lexicon.cpp", "concisemizer_random.cpp", "minimizer_lexicon.cpp" and "minimizer_random.cpp" contain the C++ implementation codes for the algorithms used for minimizing and simplifying the seed extraction process, and they calculate the density factor and E-hit rate. The header comments of this source code provide explanations of the operation methods. 

Additionally: The file "Supplementary1_The results of density and E-hits.xlsx" contains the density factors and E-hit counts for six reference genomes under different k_mer lengths and window sizes. The file "Supplementary2_The results compared with minimap2.xlsx" contains the comparison results of time and space for sequence alignment between minimap2 and minimap2-concisemizer.
