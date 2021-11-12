'''
Created on Dec 3, 2020

Script for filling out gaps in a table of taxa at different ranks

@author: alanzen
'''

import sys
from LCAClassifier.taxa import CRESTree
from LCAClassifier.config import config
from LCAClassifier.classify import ClassificationTree



def main():
    inTable = open(sys.argv[1],"r")
    
    db="silvamod128"
    
    mapFile = ("%s/%s.map" %
               (config.DATABASES[db], db))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[db], db))

    reftree = ClassificationTree(treFile, mapFile)
    
    
    firstLine=True
    for line in inTable:
        if firstLine:
            firstLine=False
            print(line[:-1])
        else:
            items = line.split(";")
            items[-1] = items[-1].replace("\n","")
            taxa = items[1:8]
            lowestRank = 0
            for i in range(7,0,-1): # 7 -- 1 speceis -> kingdom
                if (i>lowestRank and taxa[i]):
                    lowestRank = i
                if (taxa[i]):
                    lastTaxa = taxa[i]
                elif (i<lowestRank):
                    # blank taxa above specified
                    silvaChild = reftree.getNode(lastTaxa)
                    if not silvaChild:
                        sys.stderr.write("Warning: %s not found in SILVA" % taxa[i])                        
                    else:
                        taxa[i] = reftree.getParent(silvaChild).name
                        
            print(";".join([items[0],taxa,items[9:]]))
                        
                    