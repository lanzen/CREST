from LCAClassifier.taxa import CRESTree
from LCAClassifier.classify import ClassificationTree
from LCAClassifier.config import config
import sys
from Bio import Phylo
import random

def parserFasta(filename,acs):
    """This function is written by LAPORTE Amelie.
    Args:
       filename:The FASTA file of the population you want to compare to the reference.
       acs:The list containing the accession numbers. 
       pth:The list containing the taxonomic pathways.
        
    Returns:
       The function return a list for the accession numbers and one for the pathways.
    """
    #open the FASTA file
    sourceText=open(filename,"r")
    #create a list which contain one line per index
    source=sourceText.read().splitlines(1)
    #loop to only keep the headers
    for i in range(0,len(source),2):
        if '>' in source[i]:
            #to store accession number in list
            acs.append(getIndex(source[i],">","\n",1,0))
    sourceText.close()
    return acs

def getIndex(txt,mot1,mot2,int1,int2):
    """This function is written by LAPORTE Amelie.

    Args:
       txt:correspond to the genbank entry stored in memory with the readFlatFile function.
       mot1, mot2: the words used as index to separate the information.
       int1,int2: the number of characters to delete in the separated information. 
        
    Returns:
       The function return the same thing than if we've had used the index method but simplified to have a clearer code.
  
    """
    txt=txt[txt.index(mot1)+int1:txt.index(mot2,txt.index(mot1))+int2]
    return txt

def main():
    #Creation of the phylogenetic tree from the database
    database="unite"
    mapFile = ("%s/%s.map" %
               (config.DATABASES[database], database))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[database], database))
    
    
    reftree = ClassificationTree(treFile, mapFile)
    
    unite_file=sys.argv[1]
    
    fungi_family_file = sys.argv[2]
    
    unitefamily=[]
    unite_lst=[]

    fungi_txt=open(fungi_family_file,'r')
    fungi=fungi_txt.read().splitlines(0)
    for i in range(0,len(fungi)):
        unitefamily.append(reftree.getNode(fungi[i]))
    for i in range(0,len(unitefamily)):
        if unitefamily[i]!=None:
            unite_lst.append(unitefamily[i])
    
    
    children=[]
    for i in range(0,len(unite_lst)):
        children.append(reftree._getAllChildren(unite_lst[i]))
    
    print children[0]
    print unite_lst[0]
    
    accs_list=[]
    
    parserFasta(unite_file, accs_list)



if __name__ == '__main__':
    main()