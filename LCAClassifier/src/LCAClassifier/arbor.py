import sys
import re
import os
import Queue
import Bio.SeqIO

from LCAClassifier.taxa import CRESTree



# TODO completely revise to comply with new map format and enforce stricter family names


sfLimits = {CRESTree.SUBSPECIES: 1.0, # 97
            CRESTree.SPECIES: .99, # 97
            CRESTree.GENUS: .97, # 95
            CRESTree.FAMILY: .95, # 90
            CRESTree.ORDER: .90, # 85
            CRESTree.CLASS: .85, # 80
            CRESTree.PHYLUM: .80,   # 75
            CRESTree.KINGDOM: 0,   
            CRESTree.SUPERKINGDOM: 0,
            CRESTree.DOMAIN: 0,
            CRESTree.META: 0,
            CRESTree.ROOT: 0
            }

ggCodes = {"d": CRESTree.DOMAIN,
           "k": CRESTree.KINGDOM,
           "p": CRESTree.PHYLUM,
           "c": CRESTree.CLASS,
           "o": CRESTree.ORDER,
           "f": CRESTree.FAMILY,
           "g": CRESTree.GENUS,
           "s": CRESTree.SPECIES} 

class ARBor(CRESTree):
    """A Tree with specific functions for handling data exported from ARB"""

    nonSpeciesKeys = ["uncultured", "sp.", "unidentified", "metagenome",
                      "enrichment", "clone", "Unknown"]
    
    def __init__(self, name=None, rearrangeOrganelles=True, GGRankInfo=False):
        CRESTree.__init__(self)        
        
        self.rearrangeOrganelles = rearrangeOrganelles
        self.GGMode = GGRankInfo
        
        self.nodeIDs[1] = self.root        
        self.tree.root.nodeID = self.lastID = 1
        #self.rejected = [] # accessions not to add
        
        # Accession as key, node to which it points as value
        self.accessions = {}
        
        # Node name as key, depth as value
        self.ranks = {} 
        
        # Basic structure 
        if self.rearrangeOrganelles:
            self.mainGenome = self.addNode(nodename="Main genome", parent=self.root)
            #self.organelles = self.addNode(nodename="Organelle", parent=self.root)
            self.plastid = self.addNode(nodename="Chloroplast", parent = self.root)
            self.mito = self.addNode(nodename="Mitochondria", parent = self.root)
            
        else:        
            self.cellOrg = self.addNode(nodename="Cellular organisms", parent=self.root)

    def addNode(self, nodename, parent, assignmentMin=0):
        newNode =  CRESTree.addNode(self, nodename, parent, assignmentMin=assignmentMin)
        if newNode:
            self.lastID += 1        
            nodeID = self.lastID
            newNode.nodeID = nodeID
            self.nodeIDs[nodeID] = newNode
        return newNode

    def parseSilvaNDS(self, filename):
        """Parses an NDS Export file from Silva with accession, silva tax. and
        NCBI Tax."""
        
        plastidNames = ["chloroplast","Chloroplast","Plastid","plastid"]
        organelleFind = dict((k,self.plastid) for k in plastidNames)
        for k in ["Mitochondrion","mitochondrion","mitochondria","Mitochondria"]:
            organelleFind[k] = self.mito
        
        ndsFile = open(filename, 'r')

        for line in ndsFile:            
            
            parts = line.split("\t")
            
            accession = parts[0]
            if "." in accession:
                accession = accession[:accession.find(".")]
            
            taxonomy = parts[1].replace("; ",";")
            if self.GGMode:            
                taxa = re.split('[/;]', taxonomy)
            else:
                taxa = re.split('[_/;]', taxonomy)            
            
            # Deal with some unneccessary starting nodes:
            while taxa[0].lower()=="root" or taxa[0].lower=="cellular organisms":
                taxa = taxa[1:]
            
            # Deal with NCBI name / species name
            ncbi_name = parts[2].replace("\n", "")
            
            if self._useNCBINameInTaxonomy(ncbi_name):
                taxa.append(ncbi_name)
            else:
                while taxa and not self._useNCBINameInTaxonomy(taxa[-1]):
                    taxa.pop()
                        
            # Deal with organells structure
            # (the most sensible is to include an organelle node 
            # named "chloroplast" or "mitochondria" anywherein the taxonomy
            # as in Silva, but the possibility of an extra fourth column
            # kept for backwards compability)
            
            if self.rearrangeOrganelles:
                parent = self.mainGenome
                for k in organelleFind:
                    if ((len(parts) > 3) and k in parts[3]) or k in taxonomy:
                        parent = organelleFind[k]
                    if k in taxa:            
                        taxa.pop(taxa.index(k))               
            else:
                parent = self.cellOrg
                
            # Fix organelle names
            for taxon in taxa:               
                if parent == self.plastid:
                    taxon+=(" (Chloroplast)")
                elif parent == self.mito:
                    taxon+=(" (Mitochondrion)")
                    
                                                                               
            # Now go through taxa one by one again            
 
                depth = self.getDepth(parent)+1
                
                if depth in sfLimits:
                    alim = sfLimits[depth]
                else:
                    alim = 0
                    
                if depth > max(CRESTree.depths):
                    depth = max(CRESTree.depths)
                
                rank = CRESTree.depths[depth]
                
                #Remove explicit rank from greengenes node since not useful
                if self.GGMode and len(taxon)>3 and taxon[1:3]=="__":                    
                    taxon = taxon[3:]
                    dsym = taxon[0:1]                   
                
                #--- Some common problems before default check and adding of taxa
                #Check if taxon is numerical or empty and then replace with parent name
                if len(taxon)== 0 or (len(re.findall('[0-9]', taxon)) == len(taxon)):
                    taxon = "%s (%s)" % (parent.name, rank)
                
                #Fix parent of self ambigousity error
                elif taxon == parent.name or taxon.startswith("Unknown "):
                    taxon = "%s (%s)" % (parent.name, rank)
                
                #Fix incertae sedis only issues (still in 106)
                if taxon == "Incertae Sedis" or taxon == "Incertae_sedis":
                    taxon = "%s Incertae Sedis" % parent.name                
                
                
                while (taxon in self.nodeNames) and self.getParent(self.getNode(taxon)) is not parent:                     
                    
                        node = self.getNode(taxon)
                        # Set new parent if parent matches
                        #print "DEBUG: Parent %s, node: %s" %(parent.name, node.name)
                        taxon += " (%s)" % parent.name
                
                if taxon in self.nodeNames:
                    parent = self.getNode(taxon)
                else:
                    parent = self.addNode(nodename=taxon,parent=parent,
                                      assignmentMin=alim)
                 
                if self.GGMode and dsym:
                    self.ranks[taxon] = ggCodes[dsym]
                 
            # Add accession to parent
            if not accession in self.accessions:
                self.accessions[accession] = parent
            else:
                #self.rejected += accession
                sys.stderr.write("Warning: Accession %s occurs more than once.\n" % accession)
        ndsFile.close()    

    def processChangesMetadata(self, changeFile):
        """Process metadata update in tab-sep format"""
        # Changing of rank can no longer be supported due to the explicit rank criteria
        changes = open(changeFile, "r")
        for line in changes:
            try:
                parts = line.split(";")
                #process parts
                oldName = parts[0]
                newName = parts[1]
                newParent = parts[2]
                #rank = parts[3]
            except:
                sys.stderr.write("Warning: Incorrect manual change line: %s\n" %
                                 line)
                return
            
            if not newParent and not oldName:
                sys.stderr.write("Error: Incomplete information for taxon %s" % newName)
                return
            
            if newParent:
                parent = self.getNode(newParent)
                if not parent:
                    sys.stderr.write("Cannot find parent node %s\n" % newParent)
                    return
                
                if newName and not oldName:
                    # Add new node
                    if self.getNode(newName):
                        sys.stderr.write("Warning: taxon %s is already present - not added" % newName)                    
                    
                else:
                    newNode = self.addNode(nodeName=newName, parent=parent,
                                           assignmentMin=sfLimits(parent.getDepth+1))
                    self.addNode(newNode)
                    print "Added new taxon: %s" % newName            
            
            elif not newName and not newParent:
                # Remove node
                if not self.getNode(oldName):
                    sys.stderr.write("Warning: Taxon %s not found and cannot be deleted" % oldName)
                else:
                    self.deleteNode(self.getNode(oldName), False)
                    print "Deleted taxon: %s" % oldName    
            
            elif not ".." in newName:
                # Move or rename taxon
                n = self.getNode(oldName)
                if not n:
                    if self.getNode(newName):
                        sys.stderr.write("Warning: Already moved / renamed: %s" % oldName)
                    else:
                        sys.stderr.write("Error: Cannot find taxon %s\n" % oldName)
                        return
                else:
                    if newName:
                        self.renameNode(n, newName)
                        print "Renamed %s to %s" % (oldName, newName)
                    if newParent and not (newParent == n.parent.name):
                        print ("Moving %s from %s to %s" %
                               (oldName, n.parent.name, newParent))
                        self.moveNode(n, self.getNode(newParent))
    
            #Control shorthand annotation with ..
            else:
                twoNodes = oldName.split("..")
                firstParentName = self.getParent(self.getNode(twoNodes[0])).name
                secParentName = self.getParent(self.getNode(twoNodes[1])).name
                if not (firstParentName == parent and secParentName == parent):
                    sys.stderr.write("Warning: Taxons %s not moved properly in "
                                     "NDS file!!\n" % oldName)
                
    
    def enforceRankStructure(self, rankFile=None):
        """
        Read rank info and fix all nodes, removing or adding intermediate nodes
        @rankFile external tab-separated file with explicit ranks.
        """
        
        depthNames = dict((v,k) for k, v in CRESTree.depths.iteritems())        
        
        if rankFile:
            rr = open(rankFile,"r")
            for line in rr:
                lp = line[:-1].split("\t")
                if (lp[0]) in self.nodeNames and lp[1] in depthNames:
                    self.ranks[lp[0]] = depthNames[lp[1]]
                    
            rr.close()            
               
        # Do breadth first traversal and delete intermediate taxa that are not
        # part of available ranks to enforce rank structure. Detected when 
        # child node is deeper than it should be. Remember parents for
        # later reassignment of accessions
        
        parentsToDeleted = {} # name as key, parent node as value
        print "DEBUG: Root=",self.root
        delQ = Queue.Queue()
        for c in self.tree.root.clades:
            delQ.put(c)
        while not delQ.empty():
            current = delQ.get()
            print "DEBUG:",current
            properRank = self._getProperRank(current)
             
            if properRank:   
                while self.getDepth(current) > properRank:
                    parent = self.getParent(current)
                    grandparent = self.getParent(parent)                    
                    self.deleteNode(parent, moveUpChildren=True)
                    parentsToDeleted[parent.name] = grandparent
                    if properRank>=max(sfLimits):
                        self.assignmentMin[current.name] = sfLimits[properRank]
                    else:
                        self.assignmentMin[current.name] = 0
                    
            for child in self.getImmediateChildren(current):
                delQ.put(child)                
        
        # Go through accessions and check if that node is if in deletedNodes, if so re-assign to parent
        # (until parent is no longer in deletedNodes itself)        
        for acc in self.accessions:                        
            while self.accessions[acc].name in parentsToDeleted:
                self.accessions[acc] = parentsToDeleted[self.accessions[acc].name]
        
        
        self._addIntermediateRanks(self.tree.root)         

    def _addIntermediateRanks(self, node):
        
        # Do recursive depth-first and to intermediate nodes if too shallow
        
        properRank = self._getProperRank(node)
        if properRank:
            while self.getDepth(node) < properRank:
                parentName = "%s %s incertae sedis" %(node.name, CRESTree.depths[self.getDepth(node)])
                intermediate = self.addNode(parentName, self.getParent(node), sfLimits[properRank-1])
                self.moveNode(node, intermediate)                                
                
        for child in self.getImmediateChildren(node):
            self._addIntermediateRanks(child)
            
    def _getProperRank(self, node):
        """ Find depth of node from name or self.ranks dict"""            
        name = node.name
        if " (" in name:
            name = name[:name.find(" (")]                        
        if name in self.ranks: # other tests incl. "ales" "aceae"
            return self.ranks[name]
        elif name[-4:] == "ales":
            return CRESTree.ORDER
        elif name[-4:] == "ceae":
            return CRESTree.FAMILY
        else: 
            return None        
    
    def writeFasta(self, inFile, outFile):    
        newFasta = open(outFile , 'w')
    
        for record in Bio.SeqIO.parse(open(inFile), "fasta"):
        
            headerItems = record.description.split()
            acc = headerItems[-1]
            if "." in acc:
                acc = acc[:acc.find(".")]
            
            #if not (acc in self.rejected):
            if acc in self.accessions:
                #id = self.accessions[acc].nodeID
                newFasta.write(">%s\n%s\n" % (acc, record.seq))              
            else:
                sys.stderr.write("Warning: cannot find %s in taxonomy. "
                                 "Skipping!\n" % acc)
        newFasta.close()
    
    def writeConfigFiles(self, mapFile=None, treeFile=None):
        """Writes metadata files for MEGAN in tabular format and (special)        
        @mapFile the map file name
        @treeFile Newick treeFile format file name
        @limitTo only include the following accessions in map file
        """
        mapFile = open(mapFile, "w")        
        
        treeList = [self.root.nodeID, ";"]
        
        self._writeMEGANData(node=self.root, mapFile=mapFile, treeList=treeList)
        
        # also write accessions in map file mapped to their nodes
        for acc in self.accessions:
            anode = self.accessions[acc]
            assMin = self.assignmentMin[anode.name]
            mapFile.write("%s\t%s\t-1\t%s\n" % (anode.nodeID, acc, assMin))
       
        mapFile.close()

        treeFile = open(treeFile, "w")
        treeFile.write(''.join([str(item) for item in treeList]))
        treeFile.close()


    def _writeMEGANData(self, node, mapFile, treeList):
        # Self recursive, depth-first type
        
        if node.name in  self.assignmentMin:
            assMin = self.assignmentMin[node.name]
        else:
            assMin = 0
        mapFile.write("%s\t%s\t-1\t%s\n" % (node.nodeID, node.name, assMin))
      
        # Example:
        # 1. Root;
        # 2. (Bacteria)Root
        # 3. ((Firmicutes)Bacteria)Root;
        # 4. ((Firmicutes,Proteobacteria)Bacteria)Root;

        if node is self.root:
            pass
        else:
            parent = self.getParent(node)
            i = treeList.index(parent.nodeID)
            # If "," or "(" or nothing before
            # parent -> insert "(New)" after it
            if i == 0 or treeList[i - 1] == '(' or treeList[i - 1] == ',':
                treeList.insert(i, ')')
                treeList.insert(i, node.nodeID)
                treeList.insert(i, '(')
            #If ")" before parent -> insert ",New" before it
            elif  treeList[i - 1] == ')':
                treeList.insert(i - 1, node.nodeID)
                treeList.insert(i - 1, ',')
            else:
                sys.stderr.write("Problem with tree structure! Breaking.")

        for child in self.getImmediateChildren(node):
            self._writeMEGANData(node=child, mapFile=mapFile, treeList=treeList)

    # NDS helper methods
    def _useNCBINameInTaxonomy(self, ncbi_name):
        for nsKey in ARBor.nonSpeciesKeys:
            if nsKey in ncbi_name:
                return False
        return True
    
        