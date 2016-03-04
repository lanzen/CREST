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
        dims = len(datasets)
        self.totalAbundance = [0 for i in range(dims)]
        self.assignedAbundance = [0 for i in range(dims)]
        self.otuRichness = [0 for i in range(dims)]
        self.singletons = [0 for i in range(dims)]
        self.doubletons = [0 for i in range(dims)]
        
    def addFromOTU(self, otu, primary = False):
        """Updates abundance, diversity etc. from OTU. Primary also adds to assigned abundances"""
        if not len(otu.abundances) == self.dims:
            sys.std.err("Error: Incorrect dimensions in OTU table for %s" %otu.name)
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
        for i in range(dims):
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

    def __init__(self, name, tree, datasetNames = None,
                 minFilter=True, otus={}):
        self.name = name
        self.tree = tree
        if datasetNames:
            self.datasets = datasetNames
        else:
            self.datasets = [name]            
        
        self.bsr = BSR_DEFAULT
        self.ms = MS_DEFAULT
        self.minFilter = minFilter        
        self.otus = otus
        
    def assignOTU(self, otu, node):
        """
        Accepts and OTU instance and classifies it to the given node in the 
        reference tree. Traverses the tree downwards to add abundances and 
        diversity to parent nodes, creating or manipulating AssignmentInfo
        instances as neccesary
        """
        
        otu.classifyTo(node)
        if hasattr(node,"assignments"):
            ai = AssignmentInfo
        else:
            pass
            # get existing
        #add abunances etc.
        #...
        
    
    def classify_records(self, records, abundances=None, verbose=False, 
               euk_filter=False):
        """Accepts a of biopython blast iterator and carries out LCA
        assignments to a given dataset. If abundances given, it must have same order
        for all lists as the datasets list"""

        for record in records:
            qName = record.query.split(" ")[0]
            
            read_abundances={}
            
            #Determine read population from from otus if given
            if abundances:
                try:
                    i=0
                    if qName in abundances.keys():
                        seq_abundances=abundances[qName]
                    else:                   
                        qFix = qName[:qName.find("_")]
                        seq_abundances=abundances[qFix]
                        
                    for ds in datasets:
                        read_abundances[ds] = seq_abundances[i]
                        i+=1
                except:
                    print "Warning: Cannot find %s in OTU table!" %qName
            
            #Else determine read population from its name / annotation.
            else:
                if "_" in qName:
                    try:
                        readPopulation = int(qName.split("_")[-1].replace(".00",
                                                                          ""))
                    except:
                        readPopulation = 1
                elif "numreads=" in qName:
                    readPopulation = int(qName[qName.find("numreads=") +
                                               len("numreads="):])
                    
                elif "size=" in qName:
                    readPopulation = int(qName[qName.find("size=") +
                                               len("size="):-1])
                #Or set to 1, if we cannot find
                else:
                    readPopulation = 1
                
                if datasets:
                    for ds in datasets:
                        read_abundances[ds]=readPopulation
                else:
                    ra=readPopulation
            
            # Check for minimum score and any alignemnts
            if (record.alignments and read_abundances and
                record.alignments[0].hsps[0].bits >= self.ms):
                
                best_hsp = record.alignments[0].hsps[0]
                topScore = best_hsp.bits
                if self.seqs and qName in self.seqs.keys():
                    qSeq = self.seqs[qName]
                else:
                    qSeq = str(best_hsp.query).replace("-", "")
                hitname = record.alignments[0].hit_def.split()[0]
                node = self.getNode(hitname)
                if not node:
                    sys.stderr.write("Best-scoring node %s not found!\n" %
                                     hitname)
                    sys.stderr.write("Cannot assign read %s\n" % qName)
                else:
                    parents = node.getPhylogeny()[1:]

                    # Iterate through rest of hits until falling below treshold
                    for a in record.alignments[1:]:
                        if a.hsps[0].bits < float(topScore) * self.bsr:
                            break
                        hitname = a.hit_def.split()[0]
                        n = self.getNode(hitname)

                        if not n:
                            sys.stderr.write("Node " + hitname +
                                             " not found! Ignoring.\n")
                        else:
                            p = n.parent
                            # iterate through parents until found in the
                            # parents list
                            while p not in parents:
                                p = p.parent
                            parents = parents[parents.index(p):]

                    # Take a look at similarity, print info if verbose and
                    # kick up if filter
                    hsp_sim = (float(best_hsp.identities) /
                               float(best_hsp.align_length))
                    if verbose and hsp_sim >= .99:
                        print ("Read %s is %s percent similar to %s" %
                               (qName, hsp_sim * 100,
                                record.alignments[0].hit_def))

                    if self.minFilter:
                        maxRankLimit = Tree.SPECIES
                        maxRank = maxRankLimit
                        d = maxRankLimit
                        ranks = sfLimits.keys()
                        ranks.sort()
                        ranks.reverse()
                        for rank in ranks:
                            if hsp_sim < sfLimits[rank]:
                                maxRank = rank - 1
                            else:
                                break

                        while (maxRank < Tree.SPECIES and
                               maxRank < parents[0].getHighestRank()):
                            d = min(parents[0].getHighestRank(), maxRankLimit)
                            if verbose:
                                print ("Read %s cannot be assigned to "
                                       "rank %s (similarity=%s)" %
                                       (qName, Tree.depths[d],
                                        hsp_sim))
                            parents = parents[1:]

                        if d < maxRankLimit:
                            novelName = ("Unknown %s %s" %
                                         (parents[0].name, Tree.depths[d]))
                            nn = self.getNode(novelName)
                            if nn:
                                novelNode = nn
                            else:
                                depth = parents[0].getHighestRank() + 1
                                novelNode = Node(novelName, parent=parents[0],
                                                 depth=depth)
                                self.addNode(novelNode)
                            parents = [novelNode] + parents

                    # Handle assignment
                    read = Read(qName, seq=qSeq)
                    
                    if euk_filter and self.getNode("Eukaryota") in parents:
                        parents = [self.noHits]
                    
                    if datasets:
                        for ds in datasets:
                            ra=read_abundances[ds]
                            if ra>0:
                                parents[0].assignRead(read, dataset=ds, 
                                                  abundance=ra,
                                                  primary=True, recursive=True)
                    else: 
                        parents[0].assignRead(read, dataset=None, 
                                                  abundance=ra,
                                                  primary=True, recursive=True)
                    self.read_node_assignments[qName] = parents[0]
            
            #Below min. score
            elif read_abundances:
                #No hits
                if self.seqs and qName in self.seqs.keys():
                    qSeq = self.seqs[qName]
                elif record.alignments:
                    qSeq = record.alignments[0].hsps[0].query.replace("-", "")
                else:
                    qSeq = None
                nhr = Read(name=qName, seq=qSeq)
                if datasets:
                    for ds in datasets:
                        ra=read_abundances[ds]
                        if ra>0:
                            self.noHits.assignRead(nhr, dataset=ds, abundance=ra,
                                                primary=True)
                            self.root.assignRead(nhr, dataset=ds, abundance=ra, 
                                             primary=False)
                else:
                    self.noHits.assignRead(nhr, dataset=ds, abundance=ra,
                                                primary=True)
                    
                    self.root.assignRead(nhr, dataset=ds, abundance=ra, 
                                             primary=False)
                    
                self.read_node_assignments[qName] = self.noHits


    def setBitscoreRange(self, percent):
        self.bsr = 1 - float(percent) / 100

    def setMinScore(self, minScore):
        self.ms = minScore
        
    def printAssignmentsRDPQual(self, node, dataset=None, printFile=None,
                                 newTabStyle=False):
        assignments = node.getAssignment(dataset)
        if assignments and assignments.primReads:

            for r in assignments.primReads:
                toPrint = (">%s\t%s\n" %
                           (r.name,
                            node.getPhylogenyRDPStyle(
                                      root=False, newTabStyle=newTabStyle)))
                for line in self.qual[str(r.name)].format("qual").split("\n")[1:]:
                    toPrint+=(line + " ")
                if printFile:
                        printFile.write(toPrint + "\n")
                else:
                    print toPrint
        if assignments:
            for child in node.children:
                self.printAssignmentsRDPQual(node=child, dataset=dataset, 
                                              printFile=printFile,
                                              newTabStyle=newTabStyle)
                
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

#     parser.add_option("-u", "--includeunknown",
#                       action="store_true", dest="outputNovel",
#                       default=False,
#                       help=('When the minimum percent identity is is not met, '
#                             'sequences are classified as "unknown". By '
#                             'default this is not hidden from the composition '
#                             'table.'))

    parser.add_option("-m", "--minabundance",
                      dest="mintaxonabundance",
                      default=0.0,
                      help=("Minimum relative abundance for a taxon to be included in "
                            "the relative abundance table (measured in the dataset where " 
                            "it is most abundant; default=0."))
    
    parser.add_option("-i", "--fastain",
                      dest="fastafile",
                      type="str",
                      default=None,
                      help=("FASTA-file specifying *all* sequences "
                            "to be classified. (Unless specified, fasta-"
                            "output is written using BLAST alignments only.)"))
    
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

    parser.add_option('-c', '--config',
                      dest='config',
                      default=None,
                      help=('Use custom configuration file.'))                      

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
        
    #Initiate classifier instance from tree
    
    lca = LCAClassifier("reftree", minFilter=options.minFilter,
                        fastafile=options.fastafile, qualfile=options.qualfile,
                        otus=options.otus)

    mapFile = ("%s/%s.map" %
               (config.DATABASES[options.dbname], options.dbname))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[options.dbname], options.dbname))

    lca.initFrom(mapFile, treFile) #, options.synFile)

    lca.setBitscoreRange(options.bitScoreRange)
    lca.setMinScore(options.minScore)
    
    #Read (non-BIOM) OTU Table, if provided
    if options.otus:
        #Parse first line for datasets, remaining for abundances (script)
        l=0
        otuFile=open(options.otus,"r")
        abundances={} #Dictionary with sequence names and list for each row
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
                    sep=None                
                if sep:
                    self.datasets=line[:-1].split(sep)[1:]
                
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
                abundances[seqname]=seq_abundances
                
            else:
                parser.error(("No suitable delimeter found in OTU table %s"
                              % options.otus))
            
        otuFile.close()
        
    #Read BIOM-format OTU-table, if provided
    elif options.biom:
        biomFile = open(options.biom, "r")
        bTable = biom.parse.parse_biom_table(biomFile)
        biomFile.close()
        self.datasets = list(bTable.SampleIds)
        
        abundances = {}
        for obs in bTable.ObservationIds:
            seq_abundances=[]
            for ds in self.datasets:
                n_reads = int(bTable.getValueByIds(obs,ds))
                seq_abundances.append(n_reads)
            abundances[obs] = seq_abundances
                
   
    for a in args:
        if a.endswith(".gz"):
            results = gzip.open(a,"r")
        else:
            results = open(a)        

        #Parse blast (only XML)
        records = NCBIXML.parse(results)

        #Classify!
        if options.otus or options.biom:
            lca.assign(records, self.datasets, abundances, verbose=options.verbose, euk_filter=options.eukfilter)
        else: 
            #One blast file per dataset
            
            name = a.replace(".xml", "").replace(".XML", "")
            name = name.replace(".gz","")
            self.datasets.append(name)
            lca.assign(records, self.datasets=[name], abundances=None, 
                       verbose=options.verbose, euk_filter=options.eukfilter)
        
        results.close()

    #Crop trees from those nodes with no reads asssigned in ANY LCAClassifier
    lca.pruneUnassigned()

    #Write result overviews
    assignmentsCount = open(os.path.join(stub,"All_Assignments.tsv"), 'w')
    totalCount = open(os.path.join(stub,"All_Cumulative.tsv"), 'w')
    rFile = open(os.path.join(stub,"Richness.tsv"),"w")
    raFile = open(os.path.join(stub,"Relative_Abundance.tsv"),"w")
    
    if options.fastafile:
        allFaName = os.path.basename(options.fastafile)
        allFaName = allFaName[:allFaName.find(".")] +"_Assigned.fasta"
        allFasta = open(os.path.join(stub,allFaName),"w")
        lca.printAllSequences(allFasta)
        allFasta.close() 
        
        #         if fastafile:
#             fstream = open(fastafile, 'r')
#             for seq_record in SeqIO.parse(fstream, "fasta"):
#                 self.seqs[seq_record.id] = seq_record.seq
#                 
#         if qualfile:
#             qstream = open(qualfile, 'r')
#             for q_record in SeqIO.parse(qstream, "qual"):
#                 self.qual[q_record.id] = q_record
        
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
