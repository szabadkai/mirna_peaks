mirna_peaks
===========

little pipeline, to detect de novo miRNA-s from degradome and small RNA-seq, based on peak detection. (Not yet working properly, absolute no guarantee)

Instructions:
1. Trimm your small RNAs  with (cutadapt, reaper, ect).
2. Detect peaks with peak_setect.py.
3. Use intersectbed to discover the intersection of your two libraries.
4. with fasta_parse.py extract the sequence of your intersections.
5. To confirm your miRNAs use Blastn against known miRNAs (use -tag blastn-short, and format 6) 
6. To parse the blast result use parse_blast.py