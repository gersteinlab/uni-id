#!/bin/sh
# ----------------------------------------------------------------------------
# This file is part of the uni-id partitioning toolkit.
#
# Portions derived from ptools (Gürsoy et al., "Data Sanitization to Reduce
# Private Information Leakage from Functional Genomics", Cell 2020), which is
# distributed under the MIT License. The original ptools copyright and MIT
# permission notice are retained below; modifications for uni-id are also
# released under the MIT License.
#
# MIT License. Permission is hereby granted, free of charge, to any person
# obtaining a copy of this software and associated documentation files, to deal
# in the Software without restriction, including the rights to use, copy,
# modify, merge, publish, distribute, sublicense, and/or sell copies. The above
# notice and this permission notice shall be included in all copies. THE
# SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
# ----------------------------------------------------------------------------

bam_path=$1
reference_fasta=$2


bam_basename=$(basename "$bam_path")
bam_prefix=${bam_basename%.bam}

samtools view "${bam_path}" | awk '{if ($6~/N/) {print $0}}' > withN.sam
samtools view "${bam_path}" | awk '{if ($6!~/N/) {print $0}}' > withoutN.sam
samtools view -H "${bam_path}" > header.txt
python3 $(which getSeq_genome_wN.py) "${reference_fasta}" header.txt withN.sam | samtools view -h -bS - > withN.p.bam
python3 $(which getSeq_genome_woN.py) "${reference_fasta}" header.txt withoutN.sam | samtools view -h -bS - > withoutN.p.bam
samtools merge "${bam_prefix}".p.bam withN.p.bam withoutN.p.bam
samtools sort "${bam_prefix}".p.bam -o "${bam_prefix}".sorted.p.bam
rm withN.sam withoutN.sam
rm header.txt
rm withN.p.bam withoutN.p.bam
rm "${bam_prefix}".p.bam
