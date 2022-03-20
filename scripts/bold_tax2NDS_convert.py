import sys

tax=open(sys.argv[1],"r")

for line in tax:
    lp = line[:-1].split("\t")
    acc = lp[0]
    taxa = lp[1].replace("_"," ")
    if taxa.endswith(";"):
        taxa = taxa[:-1]
    tp = taxa.split(";")
    sp = tp[-1]
    parents = tp[:-1]

    print ("%s\t%s\t%s" % (acc, "/".join(parents), sp))

tax.close()
    
