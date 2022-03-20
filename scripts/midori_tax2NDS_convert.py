import sys

ranks = ["kingdom","phylum","class","order","family","genus","species"]
sm138_pr2 = "silvamod138_PR2.nds"

# Read NDS file from silvamod138_PR2 and get all unique superkingdom,kingdom,phyla combos
nds = open(sm138_pr2,"r")
smpr2_taxonomy = {}
for line in nds:
    taxa = (line[line.find("\t"):line.rfind("\t")]).split("/")
    # Iterate over what is supposed to be phylum, class, order, because
    # sometimes there are intermediate ranks in the NDS
    for i in range(3,7):
        if len(taxa)>i and not (taxa[i] in smpr2_taxonomy):
            sk_k_p = "/".join(taxa[2:i+1])
            smpr2_taxonomy[taxa[i]] = sk_k_p
            #print("DEBUG: "+sk_k_p)

nds.close()

not_found = []

# Read RDP formatted (poorly) FASTA for Midori
fasta=open(sys.argv[1],"r")

for line in fasta:
    if line.startswith(">"):        
        acc = line[1:line.find(".")]
        taxa = line[line.find("root_1")+7:-1]
        tp = taxa.split(";")
        cleanTaxa = []
        
        for taxon in tp:
            ## Remove rank names that appear at random (at least in RDP version)
            for rank in ranks:
                taxon = taxon.replace(rank+"_","")
                
            ## Remove _NCBI ID
            taxon = taxon[:taxon.find("_")]
            
            # Add taxon
            cleanTaxa.append(taxon)
            
        ## Add Superkingdom and kingdom to first taxon (phylum) using silvamod138_PR2
        ph = cleanTaxa[1]
        if ph in smpr2_taxonomy:
            ph = smpr2_taxonomy[ph]
        elif not ph in not_found:
            sys.stderr.write("Phylum %s not in Silvamod138_PR2\n" % ph)
            not_found.append(ph)
            ph = "%s (superkingdom)/%s (kingdom)/%s" %(ph,ph,ph)            
        cleanTaxa[1] = ph
    
        print ("%s\t%s\t%s" % (acc, "/".join(cleanTaxa[:-1]), cleanTaxa[-1]))

fasta.close()
    
