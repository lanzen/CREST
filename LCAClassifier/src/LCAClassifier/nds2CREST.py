"""Script to convert and NDS file from Silva into a map, synonyms and taxonomy
file for MEGANs.

Assumes  accession as last item in fasta header.
"""

from LCAClassifier.arbor import ARBor
from Bio import Phylo
from optparse import OptionParser

def main():
    
    parser = OptionParser()
    
    parser.add_option("-o", "--db_name",
                      dest="name",
                      type="string",
                      default=None,
                      help="Name (stub) of derived CREST reference database")
    
    parser.add_option("-i", "--seq_file",
                      dest="fastaFile",
                      type="string",
                      default=None,
                      help="Input reference sequence file (in FASTA format)")
    
    parser.add_option("-c", "--changes_file",
                      dest="changeFile",
                      type="string",
                      default=None,
                      help="Tab separated manual changes file for processing (optional)")
    
    
    parser.add_option("-r", "--rank_file",
                      dest="rankFile",
                      type="string",
                      default=None,
                      help="Tab-separated text file containing explicit ranks for taxa")
    
    parser.add_option("-g", "--greengenes_ranks",
                      action="store_true",dest="GGRankInfo",
                      default=False,
                      help="Greengenes style explicit rank info in NDS files")
    
    parser.add_option("-p", "--nucleus_only",
                      action = "store_false", dest="euk_rearrange",
                      default=True,                    
                      help="Do not rearrange eukaryotic 16S (plastids / mitochondria) as separate domains")       
    
    parser.add_option("-m", "--megan",
                      action="store_true",dest="megan_out",
                      default=False,
                      help="Write MEGAN formatted .map and .tre files (extra)")
    
    parser.add_option("-s", "--sintax",
                      action="store_true",dest="sintax_out",
                      default=False,
                      help="Write SINTAX formatted FASTA (extra)")
                      
    
    
    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.error("You need to specify at least one NDS file in tab-separated format")
        
    if not options.fastaFile:
        parser.error("Specify a the input reference sequence file") 
        
    if not options.name:
        parser.error("Specify a name (stub) of the output reference database")
        
    if options.GGRankInfo and options.rankFile:
        parser.error("Option --greengenes_ranks is incompatible with supplying a separate rank file")        
          
        
    print "...Initiating ARB Tree and reading rank lists"
    st = ARBor(name=options.name, rearrangeOrganelles=options.euk_rearrange,
               GGRankInfo=options.GGRankInfo )
        
    i=1
    for nds in args:
        st.parseSilvaNDS(nds)
        i += 1
        
    if options.changeFile:
        print "...Making manual changes"        
        st.processChangesMetadata(options.changeFile)
            
    if options.rankFile or options.GGRankInfo:
        print "...Fixing rank structure"
        st.enforceRankStructure(options.rankFile)
        
    print "...Writing ASCII representations"
    asciiout = open(options.name+"_ascii_tree.txt","w")
    indout = open(options.name+"_indented_tree.txt","w")
    Phylo.draw_ascii(st.tree, file=asciiout)
    indout.write(str(st.tree))
    asciiout.close()
    indout.close()
        
    print "...Writing fasta"
    st.writeFasta(inFile = options.fastaFile, outFile = options.name + ".fasta")
        
    print "...writing CREST config files"
    st.writeConfigFiles(mapFile=options.name + ".map", treeFile=options.name + ".tre")    
    
    if options.megan_out:
        print "...writing MEGAN config files"
        st.writeMEGANFiles(mapFile=options.name + "_megan.map", treeFile=options.name + "_megan.tre")
    
    if options.sintax_out:
        print "...writing SINTAX training data (FASTA format)"
        st.writeFasta(inFile = options.fastaFile, outFile = options.name + "_sintax.fasta", 
                      sintax=True)

if __name__ == "__main__":

    main()