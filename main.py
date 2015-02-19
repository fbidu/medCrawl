__author__ = 'fbidu'

import os
import sys

try:
    from Bio import Entrez
except ImportError:
    print("Biopython not found!")
    sys.exit(1)


