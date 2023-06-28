import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO
from re import findall as refindall
from matplotlib import pyplot as plt
from tabulate import tabulate

try '-i P07327.fasta'