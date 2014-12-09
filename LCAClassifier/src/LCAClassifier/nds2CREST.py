"""Script to convert and NDS file from Silva into a map, synonyms and taxonomy
file for MEGANs.

Assumes  accession as last item in fasta header.
"""

import sys

from LCAClassifier.arbor import ARBor, Tree
import Bio.SeqIO


if len(sys.argv) < 4:
    print "Use python nds2megan output-stub fasta-file changes-file ndsfile(s)"

else:
    stub = sys.argv[1]
    
    print "...Initiating ARB Tree and reading rank lists"
    st = ARBor(name=stub)
    
    i = 1
    for nds in sys.argv[4:]:
        if "eukaryot" in nds.lower():
            st.parseSilvaNDS(nds, eukaryotic=True)
            print "Parsing eukaryotic file %s" % nds
        else:
            st.parseSilvaNDS(nds, eukaryotic=False)
            print "Parsing prokaryotic file %s" % nds
    
        i += 1
    
    print "...Making manual changes"
    
    changes = open(sys.argv[3], "r")
    for line in changes:
        st.processChangesMetadata(line)
    
    
    print "..Writing tree"
    treeFile = open(stub + "_Tree.txt", 'w')
    st.printAsTree(popDataset=False, showLeaves=True, printFile=treeFile)
    treeFile.close()
    
    print "...Writing fasta"
    
    newFasta = open(stub + ".fasta", 'w')
    rdpFasta = open(stub + "_rdp.fasta", 'w')
    added = []
    for record in Bio.SeqIO.parse(open(sys.argv[2]), "fasta"):
    
        headerItems = record.description.split()
        acc = headerItems[-1]
        if "." in acc:
            acc = acc[:acc.find(".")]
        if not (acc in st.rejected or acc in added):
            node = st.getNode(acc)
            if node:
                phyl = node.getPhylogenyRDPStyle(root=False, limit=Tree.GENUS)
                if node.parent.getHighestRank() >= Tree.GENUS:
                    rdpFasta.write(">%s %s\n%s\n" % (acc, phyl, record.seq))
                newFasta.write(">%s\n%s\n" % (acc, record.seq))
                added.append(acc)
            else:
                sys.stderr.write("Warning: cannot find %s in taxonomy. "
                                 "Skipping!\n" % acc)
    
        #DEBUG else: print "out: %s" % acc
    
    newFasta.close()
    rdpFasta.close()
    
    print "...writing MEGAN files"
    st.writeConfigFiles(map_=stub + ".map", tree=stub + ".tre",
                        rdpTraining=stub + "_rdpTraining.txt")
    
    ##st.writeConfigFiles(stub+".map", stub+".synonyms.txt", stub+".tre")
