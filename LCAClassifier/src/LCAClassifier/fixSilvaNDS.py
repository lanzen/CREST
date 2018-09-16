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
        silva = items[1]
        name = items[2]
        ncbi = items[3]
        
        bchloro = "Bacteria_Bacteria (superkingdom)/Terrabacteria_Cyanobacteria/Oxyphotobacteria/Chloroplast"
        silva = silva.replace(bchloro, "Eukaryota/Chloroplast")
        bmit = "Bacteria_Bacteria (superkingdom)/Proteobacteria (superphylum)_Proteobacteria/Alphaproteobacteria/Rickettsiales/Mitochondria"
        silva = silva.replace(bmit, "Eukaryota/Mitochondria")
        
        if ncbi[-1] == ";":
            ncbi = ncbi[:-1]
        silvaTaxa = re.split('[_/;]', silva)
        ncbiTaxa = re.split('[_/;]', ncbi)
        if silvaTaxa[-1] in ncbiTaxa:
            corpos = ncbiTaxa.index(silvaTaxa[-1])
            if ((corpos +1) < len(ncbiTaxa)) and len(silvaTaxa) < len(ncbiTaxa):       
                silvaTaxa += ncbiTaxa[corpos+1:]
            
       
        
        print "\t".join([acc,"/".join(silvaTaxa),name])
