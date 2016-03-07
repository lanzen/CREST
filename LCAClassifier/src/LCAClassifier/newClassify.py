'''
Created on Feb 26, 2016

Classifier based on assignment to LCA (Lowest Common Ancestor), using BLAST
homology search results to a reference database (currently bit-score).

CREST 3.0 major update that supports taxon-specific assignment cutoff, requires OTU-table 
(or otherwise assumes abundance of 1), and enforces strict rank hierarchy

@author: andersl
'''

#TODO Inherit Biopython tree instead
#TODO Use Biopython parser to read tree
#TODO new .map format no longer like MEGAN but with cutoffs, keep legacy for MEGAN
#TODO limit options  


import sys
import os
import gzip
from optparse import OptionParser

from Bio import SeqIO
from Bio.Blast import NCBIXML

from LCAClassifier.newTaxa import CRESTree
from LCAClassifier.config import config

import biom.parse
import biom.table

BSR_DEFAULT = 0.98
MS_DEFAULT = 155

class AssignmentInfo:
    
    """
    Vectors representing assignments (primary), accumulated abundance, and diversity
    Further, singletons and doubletons are recorded to calculate Chao1. Same dimension
    and order as datasets in classifier.
    """
    
    def __init__(self, dims):
        """Create empty AssignmentInfo holder with appropriate dimensions"""        
        self.otus = []
        self.dims = dims
        #dims = len(datasets)
        self.totalAbundance = [0 for i in range(dims)]
        self.assignedAbundance = [0 for i in range(dims)]
        self.otuRichness = [0 for i in range(dims)]
        self.singletons = [0 for i in range(dims)]
        self.doubletons = [0 for i in range(dims)]
        
    def addFromOTU(self, otu, primary = True):
        """Updates abundance, diversity etc. from OTU. Primary also adds to assigned abundances"""
        if not len(otu.abundances) == self.dims:
            sys.stderr.write("Error: Incorrect dimensions in OTU table for %s" %otu.name)
            return
        
        self.totalAbundance = [self.totalAbundance[i] + otu.abundances[i] 
                                   for i in range(self.dims)]
        
        if primary:
            self.assignedAbundance = [self.totalAbundance[i] + otu.abundances[i] 
                                   for i in range(self.dims)]
        
        for i in range(self.dims):                
            if otu.abundances[i]>0:
                self.otuRichness[i]+=1
                if otu.abundances[i]==1:
                    self.singletons[i]+=1
                elif otu.abundances[i]==2:
                    self.doubletons[i]+=1        
    
    def getAbundances(self):
        """@return total abundance list"""
        return self.totalAbundance
        
    def getAssignedAbundances(self):
        """@return assigned abundances"""
        return self.assignedAbundance
        
    def getDiversity(self):
        """@return abundance list"""
        return self.otuRichness
    
    def getChaoEstimate(self):
        """@return chao1 estimate based on number of singletons, doubletons and total div."""
        chao=[]
        for i in range(self.dims):
            if self.doubletons[i]>0:
                chao.append((float(self.singletons[i]*self.singletons[i]) / 
                             float(2*self.doubletons[i]))+ self.otuRichness[i])
            else:
                chao.append(None)
                
        return chao
        
class OTU:
    """Represents an OTU, which has a sequence, abundance across datasets, and classification"""
    def __init__(self, name, abundances, sequence=None, 
                 quality=None, classification=None):
        self.name = name #given by dict in classifier?
        self.sequence = sequence
        self.abundances = abundances #same dimension and oreder as classifier datasets
        self.quality = quality
        self.classification = classification
    
    def classifyTo(self, node):
        self.classification = node
    

class LCAClassifier():
    """A classifier instance. Requires a tree for assignments"""

    def __init__(self, name, tree, datasetNames = None, otus={}):
        self.name = name
        self.tree = tree
        if datasetNames:
            self.datasets = datasetNames
        else:
            self.datasets = [name]            
        
        self.bsr = BSR_DEFAULT
        self.ms = MS_DEFAULT
        self.otus = otus
        
    def assignOTU(self, otu, node, primary=True):
        """
        Accepts and OTU instance and classifies it to the given node in the 
        reference tree. Traverses the tree downwards to add abundances and 
        diversity to parent nodes, creating or manipulating AssignmentInfo
        instances as neccesary
        """
        
        otu.classifyTo(node)
        if hasattr(node,"assignments"):
            ai = AssignmentInfo(len(self.datasets)) 
        else:
            ai = otu.assignments
            ai.addFromOTU(otu, primary)
            #Traverse tree downwards to the root
            if (node is not self.tree.root) and (node is not self.tree.noHits):
                self.assignOTU(otu, self.tree.getParent(node), primary=False)
        
    
    def classify_records(self, records, verbose=False, minFilter=True,
               euk_filter=False, noUnknowns=False):
        """Accepts a of biopython blast iterator and carries out LCA assignments."""

        for record in records:
            qName = record.query.split(" ")[0]
            
            if qName in self.otus.keys():
                otu = self.otus[qName]
            elif len(self.datasets)==1:
                # OTUs are created from Fasta files in case missing in OTU table
                # (or in case OTU table is not given). Othwerise a new one is created
                sys.stderr.write("Warning! Cannot find %s in sequences / otu list. Setting read \
                abundance to 1" %qName)
                otu = OTU(name=qName, abundances=[1])
            else:
                # An OTU table was given but this OTU still does not show up - ignore it!
                sys.stderr.write("Warning! Cannot find %s in OTU table. Ignoring" %qName)
                continue
            
            if not (record.alignments and record.alignments[0].hsps[0].bits >= self.ms):                
                self.assignOTU(otu, self.tree.noHits)
                
            else:
                #lcaNode = self.tree.common_ancestor(#the hits inside the threshold))
                                            
                best_hsp = record.alignments[0].hsps[0]
                topScore = best_hsp.bits
                if not otu.sequence:                
                    otu.sequence = str(best_hsp.query).replace("-", "")
                
                hitname = record.alignments[0].hit_def.split()[0]
                # Find node of best hit - this works because always accession number
                bestNode = self.tree.getNode(hitname)
                
                if not bestNode:
                    sys.stderr.write("ERROR: Best-scoring node %s not found! Cannot assign \n" %
                                     hitname)
                    sys.stderr.write("read (incorrectly formatted database?) %s\n" % qName)
                    continue                
                
                # Iterate through rest of hits until falling below treshold
                allHitNodes = [bestNode]
                for a in record.alignments[1:]:
                    if a.hsps[0].bits < float(topScore) * self.bsr:
                        break
                    hitname = a.hit_def.split()[0]
                    n = self.tree.getNode(hitname)

                    if not n:
                        sys.stderr.write("Node %s not found! Ignoring. \
                        (incorrectly formatted database?\n" % n)
                    else:
                        allHitNodes.append(n)
                
                lcaNode = self.tree.common_ancestor(allHitNodes)        
                
                # Take a look at similarity, print info if verbose and
                # kick down if filter
                hsp_sim = (float(best_hsp.identities) /
                           float(best_hsp.align_length))
                if verbose and hsp_sim >= .99:
                    print ("Read %s is %s percent similar to %s" %
                           (qName, hsp_sim * 100,
                            record.alignments[0].hit_def))

                # Push down to lowest possible and create new node unless opted out
                while (minFilter and hsp_sim < lcaNode.assignmentMin and 
                       lcaNode is not self.tree.root):
                    if verbose:
                        print ("OTU %s cannot be assigned to %s (similarity=%s)" %
                                   (qName, lcaNode.name, hsp_sim))                    
                    parent = lcaNode.getParent()
                    if hsp_sim < parent.assignmentMin or noUnknowns:
                        lcaNode = parent
                    else:
                        # Create unknown node - a new one for each OTU
                        i=1
                        u_name = "Unknown %s %s %s" % (lcaNode.name, self.tree.getRank(lcaNode), i)                        
                        while self.tree.getNode(u_name):
                            i+=1
                            u_name = "Unknown %s %s %s" % (lcaNode.name, 
                                                           self.tree.getRank(lcaNode), i)
                        
                        lcaNode = self.tree.addNode(u_name, parent=parent)
                
                # Last check whether to remove because eukaryote and if not assign  
                if euk_filter and lcaNode.get_path[0].name == "Eukaryota":
                    self.assignOTU(otu, self.tree.noHits)
                else:
                    self.assignOTU(otu, lcaNode)          


    def setBitscoreRange(self, percent):
        self.bsr = 1 - float(percent) / 100

    def setMinScore(self, minScore):
        self.ms = minScore
    
    def writeBIOMTable(self, bTable,biomOut):
 
        taxonomy_dict = {}
        for obs in bTable.ObservationIds:
            try:
                aNode = self.read_node_assignments[obs]
                name_list = aNode.getPhylogenyNameList()
                taxonomy_dict[obs] = {"taxonomy":name_list}
            except:
                sys.stderr.write("Problem with assignment of",obs,": classification not found!")
                taxonomy_dict[obs] = {"taxonomy":self.noHits.getPhylogenyNameList()}
        bTable.addObservationMetadata(taxonomy_dict)
        biomOut.write(bTable.getBiomFormatJsonString("CREST"))
        
    def printAllSequences(self, fastaOut):
        for seqName, node in self.read_node_assignments.iteritems():
            try:
                sn = self.seqs[seqName]
            except:
                sn = self.seqs[seqName[:seqName.find("_")]]
                
            fastaOut.write(">%s %s\n%s\n" % (seqName, 
                                           node.getPhylogenyRDPStyle(root=False), 
                                           sn))            
            
    def writeOTUsWithAssignments(self, otusOut, abundances, datasets, sep="\t"):
        otusOut.write("OTU"+sep)
        for ds in datasets:
            otusOut.write(ds+sep)
        otusOut.write("classification\n")
        for seqName in abundances.keys():
            otusOut.write(seqName)
            for ab in abundances[seqName]:
                otusOut.write("%s%i" % (sep, ab))
            taxonomy = self.read_node_assignments[seqName].getPhylogenyRDPStyle(root=False)
            otusOut.write("%s%s\n" % (sep, taxonomy))


def main():

    parser = OptionParser()
    
    parser.add_option("-o", "--output-dir",
                      dest="directory",
                      type="string",
                      default="CREST_Results",
                      help="name of output directory for classification results ")


    parser.add_option("-d", "--dbname",
                      dest="dbname",
                      type="string",
                      default="silvamod",
                      help="the taxonomy to be used (defaut = silvamod)")

    parser.add_option("-r", "--range",
                      dest="bitScoreRange",
                      type="int",
                      default=2,
                      help=("bitscore-range (the range of blast hits to find "
                            "LCA of, given in percent drop from highest "
                            "score; default = 2)"))

    parser.add_option("-s", "--minscore",
                      dest="minScore",
                      type="int",
                      default=MS_DEFAULT,
                      help=("minimum bit-score given as an integer; "
                            "default = 155"))

    parser.add_option("-n", "--normalisetobase",
                      action="store_true", dest="normBase",
                      default=False,
                      help=("Normalise domain level and lower to those "
                            "classified at base level"))

    parser.add_option("-f", "--nofilter",
                      action="store_false", dest="minFilter",
                      default=True,
                      help=("deactivate minimum identity filter preventing "
                            "classification to higher ranks when a minimum "
                            "rank-identity is not met (3% for species, 5% for "
                            "genera, 10% for family"))

    parser.add_option("-m", "--minabundance",
                      dest="mintaxonabundance",
                      default=0.0,
                      help=("Minimum relative abundance for a taxon to be included in "
                            "the relative abundance table (measured in the dataset where " 
                            "it is most abundant; default=0."))
    

    parser.add_option("-t", "--otuin",
                      dest="otus",
                      type="str",
                      default=None,
                      help=("OTU-table in tab- or comma-separated format, " 
                           "specifying the distribution across datasets "
                           "for the sequences classified. Unless specified, "
                           "one dataset per Blast output file (XML) is assumed. "))
    
    parser.add_option("-b", "--biom",
                      dest="biom",
                      type="str",
                      default=None,
                      help=("OTU-table in BIOM-format, " 
                           "specifying the distribution across datasets "
                           "for the sequences classified. Unless specified, "
                           "one dataset per Blast output file (XML) is assumed. "
                           "Cannot be used with option --otuin."))
    
    parser.add_option("-i", "--fastain",
                      dest="fastafile",
                      type="str",
                      default=None,
                      help=("FASTA-file specifying *all* sequences "
                            "to be classified. (Unless specified, fasta-"
                            "output is written using BLAST alignments only.)"
                            "Note: Abundance data will be taken from sequence"
                            "names as given in this file unless OTU table or BIOM"
                            "file is also specified! This is not compatible with"
                            "using several XML files to distinguish datasets"))
    
    
    parser.add_option("-q", "--qualin",
                      dest="qualfile",
                      type="str",
                      default=None,
                      help=("Output sequences from submitted FASTA quality " 
                      "file with annotation"))   
    
    
    parser.add_option("-a", "--fastaout",
                      action="store_true", dest="fastaOut",
                      default=False,
                      help="Print classified sequences in FASTA format for each dataset")   
    
  
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose",
                      default=False,
                      help=("Verbose mode (output information about "
                            "interesting assignments"))

    parser.add_option('-k', '--noeukaryotes',
                      action="store_true",dest='eukfilter',
                      default=False,
                      help="Filter any eukaryotic assignemnts (for comparison purpose)")
    
    parser.add_option('-u', '--nounknowns',
                      action="store_true",dest='noUnknowns',
                      default=False,
                      help="Do not create \"Unknown\" nodes when min. similarity filter \
                      was activated, pushing assignment to lower rank")

    parser.add_option('-c', '--config',
                      dest='config',
                      default=None,
                      help=('Use custom configuration file'))                      

    (options, args) = parser.parse_args()

    #Check arguments and create output directory

    if len(args) < 1:
        parser.error("No blast-file specified!")
        
    for blastFile in args:
        if not os.path.isfile(blastFile):
            parser.error("Blast output %s does not exist!" % blastFile)
        
    if options.directory[0] == "/":
        stub = os.path.join(options.directory)
    else:
        stub = os.path.join(os.getcwd(),options.directory)
    if not os.path.exists(stub):
        os.makedirs(stub)
    else:
        parser.error("Directory %s already exists! Please enter a different output directory "
                         "(using option -o)" % options.directory)
        
    
    #Use custom configuration if config file given
    
    if options.config is not None:
        config.configure(options.config)
    
    #Read (non-BIOM) OTU Table, if provided
    if options.otus:
        #Parse first line for datasets, remaining for abundances (script)
        l=0
        otuFile=open(options.otus,"r")
        otus={} #Dictionary with sequence names and list for each row
        for line in otuFile:
            l+=1
            if l==1 or line[0]=="#":
                if "\t" in line:
                    sep="\t"
                elif "," in line:
                    sep=","
                elif ";" in line:
                    sep=";"
                else:
                    parser.error(("No suitable delimeter found in OTU table %s"
                              % options.otus))
                datasets=line[:-1].split(sep)[1:]
                
            elif sep:
                abd=line.split(sep)
                seqname=abd[0]
                seq_abundances = []
                for a in abd[1:]:
                    try: 
                        ab = int(a)
                    except:
                        ab = int(float(a))
                    seq_abundances.append(ab)
                otus[seqname] = OTU(name=seqname, abundances = seq_abundances)
            
        otuFile.close()
        
    #Read BIOM-format OTU-table, if provided
    elif options.biom:
        otus={} #Dictionary with sequence names and list for each row
        biomFile = open(options.biom, "r")
        bTable = biom.parse.parse_biom_table(biomFile)
        biomFile.close()
        datasets = list(bTable.SampleIds)
        
        for obs in bTable.ObservationIds:
            seq_abundances=[]
            for ds in datasets:
                n_reads = int(bTable.getValueByIds(obs,ds))
                seq_abundances.append(n_reads)
            otus[obs] = OTU(name=obs, abundances=seq_abundances)

    # Take care of fasta file. If not otus or biom also create new OTUs
    # Otherwise assign seq
    if options.fastafile:
        if not (options.otus or options.biom):
            otus = {}
        fstream = open(options.fastafile, 'r')
        for seq_record in SeqIO.parse(fstream, "fasta"):
            if (options.otus or options.biom):
                if otus.has_key(seq_record.id):
                    otus[seq_record.id].sequence = seq_record.seq
                else:
                    sys.stderr.write("Warning: sequence %s not found in OTUs. Skipping" 
                                     %seq_record.id)
            else:
                if "numreads=" in seq_record.id:
                    readPopulation = int(seq_record.id[seq_record.id.find("numreads=") +
                                             len("numreads="):])
                     
                elif "size=" in seq_record.id:
                    readPopulation = int(seq_record.id[seq_record.id.find("size=") +
                                             len("size="):-1])
              
                else:
                    readPopulation = 1
                otus[seq_record.id] = OTU(name=seq_record.id, sequence=seq_record.seq,
                                          abundances=[readPopulation])
   
    if options.qualfile:
        if not options.fastafile:
            parser.error("Qual file given without sequence data!")
        qstream = open(options.qualfile, 'r')
        for q_record in SeqIO.parse(qstream, "qual"):
            otus[q_record.id].quality = q_record


#           
            
            
    #-------- Chaos follows --------
                    
   
    for a in args:
        if a.endswith(".gz"):
            results = gzip.open(a,"r")
        else:
            results = open(a)

        #Parse blast (only XML)
        records = NCBIXML.parse(results)

        #Classify!
#         if options.otus or options.biom:
#             lca.assign(records, self.datasets, abundances, minFilter=options.minFilter
#                        verbose=options.verbose, euk_filter=options.eukfilter,
#                        noUnknowns=options.noUknowns)
#         else: 
#             pass
            #One blast file per dataset 
            #TODO remove this option?
            
#             name = a.replace(".xml", "").replace(".XML", "")
#             name = name.replace(".gz","")
#             self.datasets.append(name)
            #lca.assign(records, self.datasets=[name], abundances=None, 
             #          verbose=options.verbose, euk_filter=options.eukfilter)
        
        results.close()

    #Crop trees from those nodes with no reads asssigned in ANY LCAClassifier
    lca.pruneUnassigned()
    
    #Initiate classifier instance from tree
    
    lca = LCAClassifier("reftree", 
                        fastafile=options.fastafile, qualfile=options.qualfile,
                        otus=options.otus)

    mapFile = ("%s/%s.map" %
               (config.DATABASES[options.dbname], options.dbname))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[options.dbname], options.dbname))

    lca.initFrom(mapFile, treFile) #, options.synFile)

    lca.setBitscoreRange(options.bitScoreRange)
    lca.setMinScore(options.minScore)

    #Write result overviews
    assignmentsCount = open(os.path.join(stub,"All_Assignments.tsv"), 'w')
    totalCount = open(os.path.join(stub,"All_Cumulative.tsv"), 'w')
    rFile = open(os.path.join(stub,"Richness.tsv"),"w")
    raFile = open(os.path.join(stub,"Relative_Abundance.tsv"),"w")
    
    
        
    lca.printAssignedCount(assignmentsCount, self.datasets)
    lca.printTotalCount(totalCount, self.datasets)
    #lca.printCompositionAlternative(altFile, datasets, options.outputNovel)
    lca.printRichness(rFile, datasets)
    lca.printRelativeAbundances(raFile, datasets, 
                                minimalMaxAbundance=float(options.mintaxonabundance),
                                normaliseToBase=options.normBase)
    
    assignmentsCount.close()
    totalCount.close()
    rFile.close()
    raFile.close()
    
    
    #Add taxonomy to BIOM file if used
    if options.biom:
        biomOut = open(os.path.join(stub,os.path.basename(options.biom)),"w")    
        lca.writeBIOMTable(bTable,biomOut)
        biomOut.close()
    
    #Add taxonomy to OTU table if supplied
    if options.otus:
        otusOut = open(os.path.join(stub,os.path.basename(options.otus)),"w")
        lca.writeOTUsWithAssignments(otusOut, abundances=abundances, datasets=datasets, sep=sep)
        otusOut.close()
    
    #Write output file-by file
    
   
    for name in datasets:
        
        treeFile = open(os.path.join(stub,(name + "_Tree.txt")), 'w')
        lca.printAsTree(popDataset=name, showLeaves=True, printFile=treeFile)
        treeFile.close()
        
        
        if options.rdpOut:
            rdpFile = open(os.path.join(stub,(name + "_Assignments.tsv")), 'w')
            lca.printAssignmentsRDPStyle(name, rdpFile)
            rdpFile.close()       
        
        if options.fastaOut:
            faFile = open(os.path.join(stub,(name + "_Assignments.fasta")), 'w')
            #or lca.printAllSequences(allFasta) ???
            lca.printAssignmentsRDPFasta(name, faFile)
            faFile.close()
        
        if options.qualfile:
            qualFile= open(os.path.join(stub,(name + "_Assignments.fasta.qual")), 'w')
            lca.printAssignmentsRDPQual(node=lca.root, dataset=name, printFile=qualFile)
            qualFile.close()
        
        tabFile = open(os.path.join(stub,(name + "_Composition.tsv")), 'w')
        
        for level in Tree.depths:
            tabFile.write("Assingments at %s level\n\n" %
                          Tree.depths[level])
            lca.printPopulationsAtDepth(level, outFile=tabFile, dataset=name,
                                        normaliseToBase=options.normBase)
            tabFile.write("\n")
        tabFile.close()
        


if __name__ == "__main__":

    main()
