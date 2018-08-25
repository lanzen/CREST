'''
Converts a five column NDS file into four readable by nds2CREST.py by adding
extra taxa from the last coloumn into the second when agreeing with the last 
child

Created on Aug 24, 2018

@author: alanzen
'''
    
import sys
import re

if __name__ == "__main__":
    longNDS = open(sys.argv[1],"r")
    for line in longNDS:
        items = line[:-1].split("\t")
        acc = items[0]
        short = items[1]
        name = items[2]
        long = items[3]
        
        bchloro = "Bacteria_Bacteria (superkingdom)/Terrabacteria_Cyanobacteria/Oxyphotobacteria/Chloroplast"
        short = short.replace(bchloro, "Eukaryota/Chloroplast")
        bmit = "Bacteria_Bacteria (superkingdom)/Proteobacteria (superphylum)_Proteobacteria/Alphaproteobacteria/Rickettsiales/Mitochondria"
        short = short.replace(bmit, "Eukaryota/Mitochondria")
        
        silvaTaxa = re.split('[_/;]', short)
        ncbiTaxa = re.split('[_/;]', long)
        if silvaTaxa[-1] in ncbiTaxa:
            corpos = ncbiTaxa.index(silvaTaxa[-1])
            if corpos > len(ncbiTaxa):
                silvaTaxa += ncbiTaxa[corpos+1:]
        
        print "\t".join([acc,"/".join(silvaTaxa),name])
