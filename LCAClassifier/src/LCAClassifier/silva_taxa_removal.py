'''
Created on 9 mai 2017

@author: amelie
'''

from LCAClassifier.classify import ClassificationTree
from LCAClassifier.taxa import CRESTree
from LCAClassifier.config import config
import sys
import random
from Bio import SeqIO

def outputFasta(data,filename,filename2,lst):
    """Arg:
        data: the fasta file to treat.
        filename: the name of the file which contains only the sequence of the selected family and its children's.
        filename2: the name of the file which contains all the remaining sequences.
        lst: the list of accession numbers of the family to delete and its children's.
       Returns:
        Two output files
    """
    
    fasta_sequences = SeqIO.parse(open(data),'fasta')
    
    dict_acs={}
    
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        
    
    #Create a dictionary with the accession numbers as keys.
        dict_acs[name]=sequence

    
    #Output 1 : The Fasta file containing only the sequences of the selected family.

    with open(filename,"w") as output1:
        for i in lst:
            if i in dict_acs:
                print >> output1, '>'+i,'\n', dict_acs[i]
    output1.close()
    
    #Output2: The Fasta file containing all the sequences except the one of the family.

    with open(filename2,"w") as output2:
        for i in lst:
            if i in dict_acs:
                del dict_acs[i]
        output2.writelines('>{}\n{}\n'.format(k,v) for k,v in dict_acs.items())
    output2.close()

def main():
    
   
    database="silvamod128"
    mapFile = ("%s/%s.map" %
                (config.DATABASES[database], database))
    treFile = ("%s/%s.tre" %
                (config.DATABASES[database], database))

    reftree = ClassificationTree(treFile, mapFile)
    ct = CRESTree(reftree)
    
    silva128=sys.argv[1]
    
    #Get all the node at a specific level
    all_taxa=ct.getChildrenByRank(CRESTree.FAMILY)
    print len(all_taxa)
    
    #get a list with a hundred randomly selected taxa without replacement from the previous list.
    hundred_taxa=[]
    
    nombre=int(sys.argv[2])

	
    while len(hundred_taxa)<nombre:
        random_taxa=random.choice(all_taxa)
        if random_taxa not in hundred_taxa:
            hundred_taxa.append(random.choice(all_taxa))
    
    children_ID=[]
    ch_ID=[]
    accs_fin=[]
    

    
    #Get the children belonging to the selected taxa.
    for node in hundred_taxa:
        fam_ID=[k for k,v in reftree.nodeIDs.iteritems() if v==node]
        print node,'\t',fam_ID
        children=ct._getAllChildren(node)
        
        for i in children:
            children_ID.append([k for k, v in reftree.nodeIDs.iteritems() if v is i])
            
        for j in range(0,len(children_ID)):
            for l in children_ID[j]:
                ch_ID.append(l)
        taxa_ID=fam_ID+children_ID
        print len(children_ID), len(fam_ID), len(taxa_ID)
        
        map_file=open(mapFile,'r')
        for line in map_file:
            map_text=line.split("\t")
            map_acc = map_text[1]
            map_id = map_text[0]
            if map_id in taxa_ID:
                accs_fin.append(map_acc)
        print accs_fin
        
        taxon=str(node).replace(" ", "_")
        
        filename1="sequence_deleted_"+taxon+".fasta"
        filename2="troncated_database_"+taxon+".fasta"
        outputFasta(silva128, filename1,filename2, accs_fin)
        children_ID=[]
        ch_ID=[]
        accs_fin=[]

        map_file.close()




if __name__ == '__main__':
    main()
