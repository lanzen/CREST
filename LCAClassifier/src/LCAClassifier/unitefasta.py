'''
@author: Amelie Laporte

Created date: 03-20-2017
Modified date: 03-21-2017 

Description: 
Find the accession numbers of all the fungi families, choose one randomly to take it out (+children) of the database.
Output to FASTA files: One with only the sequences of the family and its children,
one with the whole database without the family and its children.
'''

from LCAClassifier.classify import ClassificationTree
from LCAClassifier.config import config
import sys
import random
import re

def getIndex(txt,mot1,mot2,int1,int2):
    """Args:
        txt:correspond to the genbank entry stored in memory with the readFlatFile function.
        mot1, mot2: the words used as index to separate the information.
        int1,int2: the number of characters to delete in the separated information. 
        
       Returns:
        The function return the same thing than if we've had used the index method but simplified to have a clearer code.
  
    """
    txt=txt[txt.index(mot1)+int1:txt.index(mot2,txt.index(mot1))+int2]
    return txt


def outputFasta(data,filename,filename2,lst):
    """Arg:
        data: the fasta file to treat.
        filename: the name of the file which contains only the sequence of the selected family and its children's.
        filename2: the name of the file which contains all the remaining sequences.
        lst: the list of accession numbers of the family to delete and its children's.
       Returns:
        Two output files
    """
    #Open the UNITE Fasta file:
    unite_txt=open(data,'r')
    unite=unite_txt.read().splitlines(1)
    
    #Set local variables
    dict_acs={}
    acs=[] 
    
    #Parse to get all the accession numbers of the Fasta file.
    for i in range(0,len(unite),2):
        if '>' in unite[i]:
            acs.append(getIndex(unite[i],">","\n",0,0))
    
    #Create a dictionary with the accession numbers as keys.
    for i in range(0,len(acs)):
        dict_acs[acs[i]]=''
    
    #Add the corresponding sequence of an accession number.
    for i in range(0,len(unite),2):
        for j in range(0,len(acs)):
            if acs[j] in unite[i]:
                dict_acs[acs[j]]=getIndex(unite[i+1],'','\n',0,0)

                
    #Delete the values of the family to delete.
    for i in lst:
        i=">"+i
        if i in dict_acs.keys():
            del dict_acs[i]
    
    #Output 1 : The Fasta file containing only the sequences of the selected family.
    with open(filename,"w") as output1:
        for j in range(0,len(unite),2):
            for i in lst:
                if i in unite[j]: 
                    print >> output1, unite[j]+unite[j+1].replace("\n","")
   
    output1.close()
    
    #Output2: The Fasta file containing all the sequences except the one of the family.
    with open(filename2,"w") as output2:
        output2.writelines('{}\n{}\n'.format(k,v) for k,v in dict_acs.items())
    output2.close()

    unite_txt.close()

def main():
    
    #Creation of the phylogenetic tree from the database
    database="unite"
    mapFile = ("%s/%s.map" %
               (config.DATABASES[database], database))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[database], database))
    
    
    reftree = ClassificationTree(treFile, mapFile)
    
    #UNITE Fasta file of the database
    unite_file=sys.argv[1]
    
    #List of ""all"" the fungi families
    fungi_family_file = sys.argv[2]
    
    unitefamily=[]
    unite_lst=[]

    fungi_txt=open(fungi_family_file,'r')
    fungi=fungi_txt.read().splitlines(0)
    
    #To get the nodes of the fungi families
    for i in range(0,len(fungi)):
        unitefamily.append(reftree.getNode(fungi[i]))
        
    for i in range(0,len(unitefamily)):
        if unitefamily[i]!=None:
            unite_lst.append(unitefamily[i])
    
    #to get the IDs of the fungi families
    family_ID=[]
    for i in unite_lst:
        #list of list
        family_ID.append([k for k, v in reftree.nodeIDs.iteritems() if v == i])
    
    #to transform into a list:
    fam_ID=[]
    for i in range(0,len(family_ID)):
        for j in family_ID[i]:
            fam_ID.append(j)
    
    random_fam=[]
    
    #to choose randomly a family to delete of the reference
    random_fam.append(random.choice(fam_ID))

    rd_fam=reftree.getNodeByID(random_fam[0])
    print "The",rd_fam,"family will be deleted"

    
    #to get the children of this family
    children=reftree._getAllChildren(rd_fam)
    print "Searching for all the organisms belonging to the",rd_fam,"family."

    
    #get the ID of the children
    children_ID=[]
    for i in children:
        children_ID.append([k for k, v in reftree.nodeIDs.iteritems() if v == i])
    
    ch_ID=[]
    for i in range(0,len(children_ID)):
        for j in children_ID[i]:
            ch_ID.append(j)

    #get a list of the IDs of the whole family
    complete_fam=ch_ID+random_fam
    
    #search to IDs into the map file
    map_file=open(mapFile,'r')
    map_text=map_file.read().split("\t")
    
    #get all the accession numbers of the family 
    accs=[]
    for i in range(0,len(map_text),3):
        for j in complete_fam:
            if j in map_text[i]:
                accs.append(map_text[i+1])

    #get only the accession numbers that are actual numbers (to search into Fasta file)
    accPatterns = [re.compile("\d\d\d\d\d\d\d\d\Z"), re.compile("\D\D\d\d\d\d\d\d\Z"),
                   re.compile("\D\D\D\D\d\d\d\d\d\d\d\d\d\Z"),
                   re.compile("\D\D\D\D\d\d\d\d\d\d\d\d\Z")]
            
    accs_fin=[]
    for i in accs:
        if (accPatterns[0].search(i) or accPatterns[1].search(i) or accPatterns[2].search(i) or accPatterns[3].search(i)):
            accs_fin.append(i)
    
    print "The output files are in creation."
    #Output the two fasta files
    outputFasta(unite_file, "sequence_deleted.fasta","troncated_database.fasta", accs_fin)
    
    print "The output files are available"
    

if __name__ == '__main__':
    main()