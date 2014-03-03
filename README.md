mirna_peaks
===========

requirements:
=============
Python<br>
NuPy<br>
intersectbed<br>
blastn<br>

Little pipeline, to detect de novo miRNA-s from degradome and small RNA-seq, based on peak detection. (Not yet working properly, absolutely no guarantee)<br>

Instructions:
=============
1. Trimm your small RNAs  with (cutadapt, reaper, ect).<br>
2. Detect peaks with peak_detect.py.<br>
3. Use intersectbed to discover the intersection of your two libraries.<br>
4. with fasta_parse.py extract the sequence of your intersections.<br>
5. To confirm your miRNAs use Blastn against known miRNAs (use -tag blastn-short, and format 6)<br> 
6. To parse the blast result use parse_blast.py<br>