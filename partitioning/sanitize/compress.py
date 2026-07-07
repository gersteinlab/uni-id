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
"""
comp_decomp.py

Compression on a specific file, using sys (in: raw file, out: Mycompdata.txt)

Decompression of compressed file to original file (in: Mycompdata.txt, out: Mydecompdata.txt)

"""

import zlib, sys, time, base64

# Compression of raw file
rawfile = sys.argv[1]
outfile = sys.argv[2]
fp = open(rawfile, "rb")
text = fp.read()

print("Raw size:", sys.getsizeof(text))

compressed = zlib.compress(text, 9)
print("compressed size:", sys.getsizeof(compressed))

savecomp = open(outfile, "wb")
savecomp.write(compressed)
savecomp.close()
