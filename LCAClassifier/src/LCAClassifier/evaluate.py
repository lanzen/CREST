'''
Created on Mar 7, 2017

@author: andersl
'''
from LCAClassifier.taxa import CRESTree
from LCAClassifier.classify import ClassificationTree
from LCAClassifier.config import config
import sys

from Bio import Phylo

def main():
    
    database="silvamod"
    
    mapFile = ("%s/%s.map" %
               (config.DATABASES[database], database))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[database], database))
    
    
    tree = ClassificationTree(treFile, mapFile)
    
    fastafile = sys.argv[1]
    n = tree.getNode("GQ241320")    
    p = tree.getPath(n)
    print p
    textPath = [str(taxa) for taxa in p]
    print textPath
    
    # Read fasta header
    # Check classification for accession
    # Check real taxonomy path
    # Compare and count FP, TP, FN (the tricky part)
    
    # Summarise results 
            


if __name__ == '__main__':
    main()