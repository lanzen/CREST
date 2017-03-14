'''
Created on Feb 26, 2016

Classifier based on assignment to LCA (Lowest Common Ancestor), using BLAST
homology search results to a reference database (currently bit-score).

CREST 3.0 major update that supports taxon-specific assignment cutoff, requires OTU-table 
(or otherwise assumes abundance of 1), and enforces strict rank hierarchy

@author: andersl
'''

# TODO add support for ignoring super- and kingdom in  


import sys
import os
import gzip
from optparse import OptionParser

from Bio import SeqIO
from Bio.Blast import NCBIXML

from Bio import Phylo

import re

from LCAClassifier.taxa import CRESTree
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
        self.abundances = abundances #same dimension and order as classifier datasets
        self.sequence = sequence        
        self.quality = quality
        self.classification = classification
    
    def classifyTo(self, node):
        self.classification = node
    

class ClassificationTree(CRESTree):  
    
    """A tree for classification initated from tre- and map files"""
      
    def __init__(self, trefile, mapfile):
                #Invoke superclass constructor 
        CRESTree.__init__(self, Phylo.read(trefile,"newick", values_are_confidence=True, rooted=True))
        
        self.noHits = Phylo.BaseTree.Clade(name="No hits")
        self.addNode(self.noHits, self.root)
        
        # Fill nodeID dict
        print "Making node list..."
        allNodes = self._getAllChildren(self.root)
        for c in allNodes:
            self.nodeIDs[c.name] = c
        print "..done"
        
        # the two accession number patterns that exist
        accPatterns = [re.compile("\D\D\d\d\d\d\d\d\Z"), re.compile("\D\d\d\d\d\d\Z"),
                       re.compile("\D\D\D\D\d\d\d\d\d\d\d\d\d\Z")]
        
        # Read nodes from .map file (id\t name\t cutoff)
        print "Reading .map file"
        theMap = open(mapfile, "r")
                
        for line in theMap:
            parts = line.split("\t")
            nodeID = parts[0]
            name = parts[1]
            minPercent = float(parts[2])
            if (minPercent > 100) or (minPercent<0):
                minPercent = 0  
            
            n = self.getNodeByID(nodeID)
            if n:
                self.nodeNames[name] = n                
                # Unless this is just an accession, update node name (accessions are only mapped)
                if not (accPatterns[0].match(name) or accPatterns[1].match(name) or 
                         accPatterns[2].match(name)):
                    self.assignmentMin[name] = minPercent
                    n.name = name
            else:
                sys.stderr.write("Error: Node %s (%s) not found in tree!\n" % (nodeID,name))                
                                
        theMap.close()
        
    def pruneUnassigned(self, node):
        """
        Recursively removes all clades that do not have assignments ("assignments" attribute).
        from root up
        """
        
        if not hasattr(node,"assignments"):
            self.deleteNode(node, moveUpChildren=False)
        
        else:
            children = self.getImmediateChildren(node)
            for c in children:
                self.pruneUnassigned(c)
    

class LCAClassifier():
    """A classifier instance. Requires a tree for assignments"""

    def __init__(self, name, tree, datasets = None, otus={}):
        self.name = name
        self.tree = tree
        if datasets:
            self.datasets = datasets
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
        
        if primary:
            otu.classifyTo(node)
        if hasattr(node,"assignments"):
            ai = node.assignments            
        else:
            ai = AssignmentInfo(len(self.datasets))
            node.assignments = ai
            
        ai.addFromOTU(otu, primary)
        #Traverse tree downwards to the root
        if (node is not self.tree.tree.root) and (node is not self.tree.noHits):
            self.assignOTU(otu, self.tree.getParent(node), primary=False)
    
    
    def classify_records(self, records, verbose=False, minFilter=True,
               euk_filter=False, noUnknowns=False, uniqueDataset=None):
        """
        Accepts a of biopython blast iterator and carries out LCA assignments.
        Arguments correspond to parser options except uniqueDataset, which is 
        given if the dataset name comes from a BLAST file rather than OTU table
        (abundance set to 1 and only in the unique dataset i.e. reads rather than OTUs)
        """

        for record in records:
            qName = record.query.split(" ")[0]
            
            if qName in self.otus:
                otu = self.otus[qName]
            elif uniqueDataset:
                # OTUs are created from Fasta files in case missing in OTU table
                # (or in case OTU table is not given). Othwerise a new one is created

                abFulFix = [int(uniqueDataset==self.datasets[i]) for i in range(len(self.datasets))]
                otu = OTU(name=qName, abundances=abFulFix)
            else:
                # An OTU table was given but this OTU still does not show up - ignore it!
                sys.stderr.write("Warning! Cannot find %s in OTU table. Ignoring\n" %qName)
                continue
            
            
            if not (record.alignments and record.alignments[0].hsps[0].bits >= self.ms):                
                self.assignOTU(otu, self.tree.noHits)                
                
            else:                                  
                best_hsp = record.alignments[0].hsps[0]
                topScore = best_hsp.bits
                if not otu.sequence:                
                    otu.sequence = str(best_hsp.query).replace("-", "")
                
                besthitname = record.alignments[0].hit_def.split()[0]      
                
                # Iterate through rest of hits until falling below treshold
                allHitNodes = [besthitname]
                for a in record.alignments[1:]:
                    if a.hsps[0].bits < float(topScore) * self.bsr:
                        break
                    hitname = a.hit_def.split()[0]                   
                    allHitNodes.append(hitname)
                
                # Lowest common ancestor
                lcaNode = self.tree.getCommonAncestor(allHitNodes)
                                                
                if not lcaNode:
                    sys.stderr.write("DEBUG: No LCA. Assigning to No Hits:\n%s" % allHitNodes)
                    self.assignOTU(otu, self.tree.noHits)
                    continue
                                 
                # Take a look at similarity, print info if verbose and
                # kick down if filter
                hsp_sim = (float(best_hsp.identities) /
                           float(best_hsp.align_length))*100
                if verbose and hsp_sim >= 99:
                    print ("Read %s is %s percent similar to %s" %
                           (qName, hsp_sim,
                            record.alignments[0].hit_def))

                # Push down to lowest possible and create new node unless opted out
                while (minFilter and lcaNode.name in self.tree.assignmentMin 
                       and hsp_sim < self.tree.assignmentMin[lcaNode.name] and 
                       lcaNode is not self.tree.root):
                    if verbose:
                        print ("OTU %s cannot be assigned to %s (similarity=%s)\n" %
                                   (qName, lcaNode.name, hsp_sim))                   
                    parent = self.tree.getParent(lcaNode)
                    if hsp_sim < self.tree.assignmentMin[parent.name] or noUnknowns:
                        lcaNode = parent
                    else:
                        # Create unknown node - a new one for each OTU
                        i=1
                        u_name = "Unknown %s %s %s" % (lcaNode.name, self.tree.getRank(lcaNode), i)                        
                        while u_name in self.tree.nodes:
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
    
    def writeBIOMTable(self, bTable, biomOut):
        taxonomy_dict = {}
        for obs in bTable.ObservationIds:
            try:
                aNode = self.otus[obs]
                name_list = self.tree.getFormattedTaxonomyPath(aNode)
                taxonomy_dict[obs] = {"taxonomy":name_list}
            except:
                sys.stderr.write("Problem with assignment of",obs,": classification not found!")
                taxonomy_dict[obs] = {"taxonomy":self.tree.getFormattedTaxonomyPath(self.tree.noHits)}
        bTable.addObservationMetadata(taxonomy_dict)
        biomOut.write(bTable.getBiomFormatJsonString("CREST"))
        
    def writeAllSequences(self, fastaOut, writeQual=False):
        for key in self.otus:
            otu = self.otus[key]
            try:
                if writeQual:
                    seq=otu.quality
                else:
                    seq=otu.sequence
                fastaOut.write(">%s %s\n%s\n" % (otu.name, self.tree.getFormattedTaxonomyPath(otu.classification),                                               
                                                 seq))
            except:
                sys.stderr.write("Problem with assignment of",otu,": classification not found!")
                fastaOut.write(">%s %s\n%s\n" % (otu.name, "No classification!",                                               
                                                 otu.sequence))
            
    def writeOTUsWithAssignments(self, otusOut, sep="\t"):
        otusOut.write("OTU"+sep)
        for ds in self.datasets:
            otusOut.write(ds+sep)
        otusOut.write("classification\n")
        for key in self.otus.keys():
            otu = self.otus[key]
            otusOut.write(otu.name+sep)
            otusOut.write (sep.join(str(oa) for oa in otu.abundances)) 
            otusOut.write("%s%s\n" % (sep,self.tree.getFormattedTaxonomyPath(otu.classification)))
            
    def writeAssignedCount(self, outFile, countType="totalAbundance", depths=CRESTree.depths):
        """
        Writes all assignments, either direct trefileones (not including children), total ones (incl.), or
        richness (no. of OTUs assigned)
        """        
        datasetheaders = "\t".join(ds for ds in self.datasets)       
        outFile.write("Rank\tTaxonpath\tTaxon\t%s\n" % datasetheaders)

        for i in depths:
            nodes = self.tree.getChildrenByRank(i)
            for node in nodes:                
                if not hasattr(node, "assignments"):
                    sys.stderr.write("Warning: node %s has no assignments" % node.name)
                    continue
                aToPrint = getattr(node.assignments,countType)
                formattedAb =  "\t".join(str(x) for x in aToPrint)
                outFile.write("%s\t%s\t%s\t%s\n" % (CRESTree.depths[i],
                                        self.tree.getFormattedTaxonomyPath(node),
                                        node.name, formattedAb))
                
    def writeChaoEstimates(self, outFile, depths=CRESTree.depths):
        datasetheaders = "\t".join(ds for ds in self.datasets)       
        outFile.write("Rank\tTaxonpath\tTaxon\t%s\n" % datasetheaders)

        for i in depths:
            nodes = self.tree.getChildrenByRank(i)
            for node in nodes:                
                if not hasattr(node, "assignments"):
                    sys.stderr.write("Error: node %s has no assignments!\n" % node.name)
                    continue
                chaos = node.assignments.getChaoEstimate()
                formattedAb =  "\t".join(str(f) for f in chaos)
                outFile.write("%s\t%s\t%s\t%s\n" % (CRESTree.depths[i],
                                        self.tree.getFormattedTaxonomyPath(node),
                                        node.name, formattedAb))
                
    def writeRelativeAbundances(self, outFile,minimalMaxAbundance=0, depths=CRESTree.depths): 
        
        if not (hasattr(self.tree.tree.root, "assignments") or 
                sum(self.tree.tree.root.assignments.getAbundances())==0):
            sys.stderr.write("ERROR: No assignments done for any dataset!\n")
            return
        rootAbRaw = [float(i) for i in self.tree.tree.root.assignments.getAbundances()]
        i=0
        rootAb=[]
        keep = []
        datasetheaders = ""
        for ab in rootAbRaw:            
            if ab>0:
                rootAb.append(ab)
                keep.append(True)
                datasetheaders += self.datasets[i]+"\t"
            else:
                keep.append(False)
            i+=1
        
        outFile.write("Rank\tTaxonpath\tTaxon\t%s\n" % datasetheaders[:-1])

        for depth in depths:
            nodes = self.tree.getChildrenByRank(depth)
            dn = CRESTree.depths[depth]
            for node in nodes:                
                if not hasattr(node,"assignments"):
                    sys.stderr.write("ERROR: node %s has no assignments" % node.name)
                    continue
                taRaw = node.assignments.getAbundances()
                i=0
                ra=[]
                for a in taRaw:                    
                    if keep[i]:
                        ra.append(float(a)/rootAb[i])
                    i+=1
                
                #Check min. abundance (in at least one dataset)
                aboveMin = True in [a>=minimalMaxAbundance for a in ra]
                if aboveMin:
                    formattedAb =  "\t".join(str(f) for f in ra)
                    outFile.write("%s\t%s\t%s\t%s\n" % (dn,
                                        self.tree.getFormattedTaxonomyPath(node),
                                        node.name, formattedAb))
                    
def main():

    parser = OptionParser()
    
    parser.add_option("-o", "--output-dir",
                      dest="directory",
                      type="string",
                      default="CREST_Results",
                      help="Name of output directory for classification results ")

    parser.add_option("-d", "--dbname",
                      dest="dbname",
                      type="string",
                      default="silvamod",
                      help="Taxonomy and reference database to be used (defaut = silvamod)")

    parser.add_option("-r", "--range",
                      dest="bitScoreRange",
                      type="int",
                      default=2,
                      help=("Bitscore-range (the range of blast hits to find "
                            "LCA of, given in percent drop from highest "
                            "score; default = 2)"))

    parser.add_option("-s", "--minscore",
                      dest="minScore",
                      type="int",
                      default=MS_DEFAULT,
                      help=("Minimum bit-score given as an integer; "
                            "default = 155"))

    parser.add_option("-f", "--nofilter",
                      action="store_false", dest="minFilter",
                      default=True,
                      help=("Deactivate minimum identity filter preventing "
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
                           "for the sequences classified."))
    
    parser.add_option("-b", "--biom",
                      dest="biom",
                      type="str",
                      default=None,
                      help=("OTU-table in BIOM-format, " 
                           "specifying the distribution across datasets "
                           "for the sequences classified."
                           "Cannot be used with option --otuin."))
    
    parser.add_option("-i", "--fastain",
                      dest="fastafile",
                      type="str",
                      default=None,
                      help=("FASTA-file specifying all sequences "
                            "to be classified. "
                            "Note: Abundance data will be taken from sequence "
                            "names as given in this file unless OTU table or BIOM "
                            "file is also specified! This is not compatible with "
                            "using several XML files to distinguish datasets"))
    
    
    parser.add_option("-q", "--qualin",
                      dest="qualfile",
                      type="str",
                      default=None,
                      help=("Output sequences from submitted FASTA quality " 
                      "file with annotation"))       
    
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
            datasets = [options.fastafile[:options.fastafile.find(".")]]
        fstream = open(options.fastafile, 'r')
        for seq_record in SeqIO.parse(fstream, "fasta"):
            if (options.otus or options.biom):
                if otus.has_key(seq_record.id):
                    otus[seq_record.id].sequence = seq_record.seq
                else:
                    sys.stderr.write("Warning: sequence %s not found in OTUs. Skipping\n" 
                                     %seq_record.id)
            else:
                if "numreads=" in seq_record.id:
                    readPopulation = int(seq_record.id[seq_record.id.find("numreads=") +
                                             len("numreads="):-1])
                     
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

    # Control that not fasta file is given in combination with several blast files
    if len(args) > 1 and options.fastafile and not (options.otus or options.biom) :
        parser.error("FASTA file used to generate OTUs but several blast outputs given representing"
                     "individual datasets!")
        
    # In some cases use blast outputs to create dataset info
    elif not options.fastafile and not options.otus and not options.biom:
        datasets = []
        for a in args:
            datasets.append(a[:a.find(".")])    
    
    mapFile = ("%s/%s.map" %
               (config.DATABASES[options.dbname], options.dbname))
    treFile = ("%s/%s.tre" %
               (config.DATABASES[options.dbname], options.dbname))

    reftree = ClassificationTree(treFile, mapFile)
    if options.otus or options.biom or options.fastafile:
        lca = LCAClassifier(name = "classifier instance", tree = reftree,
                            datasets = datasets, otus=otus)
    else:
        lca = LCAClassifier(name = "classifier instance", tree = reftree,
                            datasets = datasets)
        
    lca.setBitscoreRange(options.bitScoreRange)
    lca.setMinScore(options.minScore)
    
   
    # Parse blast files! 
    for a in args:
        if a.endswith(".gz"):
            results = gzip.open(a,"r")
        else:
            results = open(a)
        #Parse blast (only XML)
        records = NCBIXML.parse(results)        
        if options.otus or options.biom or options.fastafile:
            ud = None
        else:
            ud = a[:a.find(".")]     

        lca.classify_records(records, minFilter=options.minFilter,
                        verbose=options.verbose, euk_filter=options.eukfilter,
                        noUnknowns=options.noUnknowns, uniqueDataset=ud)
        results.close()

    # Remove empty nodes
    lca.tree.pruneUnassigned(lca.tree.root)
    
    #Write result overviews
    assignmentsCount = open(os.path.join(stub,"All_Assignments.tsv"), 'w')
    totalCount = open(os.path.join(stub,"All_Cumulative.tsv"), 'w')
    rFile = open(os.path.join(stub,"Richness.tsv"),"w")
    chaoFile = open(os.path.join(stub,"Chao1.tsv"),"w")
    raFile = open(os.path.join(stub,"Relative_Abundance.tsv"),"w")
    
    lca.writeAssignedCount(assignmentsCount, "assignedAbundance")
    lca.writeAssignedCount(totalCount, "totalAbundance")
    lca.writeAssignedCount(rFile, "otuRichness")
    lca.writeChaoEstimates(chaoFile)
    lca.writeRelativeAbundances(raFile,minimalMaxAbundance=float(options.mintaxonabundance))    
    assignmentsCount.close()
    totalCount.close()
    rFile.close()
    chaoFile.close()
    raFile.close()    
    
    #Add taxonomy to BIOM file if used
    if options.biom:
        biomOut = open(os.path.join(stub,os.path.basename(options.biom)),"w")    
        lca.writeBIOMTable(bTable,biomOut)
        biomOut.close()
    
    #Add taxonomy to OTU table (new in v3: always do also if not supplied)
    if options.otus:
        otusOut = open(os.path.join(stub,os.path.basename(options.otus)),"w")        
    else:
        otusOut = open(os.path.join(stub,os.path.basename("otus.csv")),"w")
        sep="\t"
    lca.writeOTUsWithAssignments(otusOut, sep=sep)
    otusOut.close()
    
    #Write fasta file (new in v3: always do this)
    if options.fastafile:
        allFaName = os.path.basename(options.fastafile)
        allFaName = allFaName[:allFaName.find(".")] +"_Assigned.fasta"
    else:
        allFaName = "Aligned_assignments.fasta"
    allFasta = open(os.path.join(stub,allFaName),"w")
    lca.writeAllSequences(allFasta)
    allFasta.close()    
  
    if options.qualfile: 
        qualName = os.path.basename(options.qualfile)
        qualName = qualName[:qualName.find(".")] +"_Assigned.fasta.qual"
        qualFile= open(os.path.join(stub,(qualName)), 'w')
        lca.writeAllSequences(qualFile,writeQual=True) 
        qualFile.close()      

if __name__ == "__main__":

    main()
