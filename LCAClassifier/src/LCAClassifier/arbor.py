import sys
import re
import Queue
import Bio.SeqIO

from LCAClassifier.taxa import CRESTree
from posix import access


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

def isInt(str):
    try:
        int(str)
        return True
    except ValueError:
        return False


class ARBor(CRESTree):
    """A Tree with specific functions for handling data exported from ARB"""

    nonSpeciesKeys = ["uncultured", "sp.", "unidentified", "metagenome",
                      "environmental","enrichment",
                       "clone", "Unknown", "sp.", "spp.", "bacteri"]
    
    def __init__(self, name=None, rearrangeOrganelles=True, GGRankInfo=False, NCBIColumn=True):
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
        
        if self.rearrangeOrganelles:
            plastidNames = ["chloroplast","Chloroplast","Plastid"]
            organelleFind = dict((k,self.plastid) for k in plastidNames)
            for k in ["Mitochondrion","mitochondrion","mitochondria","Mitochondria"]:
                organelleFind[k] = self.mito
        
        ndsFile = open(filename, 'r')

        for line in ndsFile:            
            
            parts = line[:-1].split("\t")
            
            accession = parts[0]
            if "." in accession:
                accession = accession[:accession.find(".")]
            
            #print "DEBUG %s" % accession
            
            taxonomy = parts[1].replace("; ",";")
            if self.GGMode:            
                taxa = re.split('[/;,]', taxonomy)
            else:
                taxa = re.split('[_/;]', taxonomy)
            
            # Deal with some unneccessary starting nodes:
            while taxa[0].lower()=="root" or taxa[0].lower=="cellular organisms":
                taxa = taxa[1:]
            
            # Deal with NCBI name / species name            
            if self.NCBIColumn and self._useNameInTaxonomy(taxa=taxa, ncbi_name=parts[2]):
                taxa.append(parts[2])
            else:
                while taxa and not self._useNameInTaxonomy(taxa=taxa):
                    taxa.pop()
                    
            # Remove all empty or only numerical taxa
            remove = []
            for i in range(0,len(taxa)):
                if len(taxa[i])== 0 or (len(re.findall('[0-9]', taxa[i])) == len(taxa[i])):
                    remove.append(taxa[i])
            for r in remove:
                taxa.pop(taxa.index(r))

                        
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
                        print "DEBUG: %s : %s" % (taxa, organelleFind[k])
                    if k in taxa:            
                        taxa.pop(taxa.index(k))   
                        
                # Fix organelle names
                for i in range(0,len(taxa)):
                    if parent == self.plastid:
                        taxa[i]+=(" (Chloroplast)")
                    elif parent == self.mito:
                        taxa[i]+=(" (Mitochondrion)")
                       
            else:
                parent = self.cellOrg
            
                                          
            # Now go through taxa one by one again            
            for i in range(0,len(taxa)):
                depth = self.getDepth(parent)+1
                
                if depth in sfLimits:
                    alim = sfLimits[depth]
                else:
                    alim = 0
                    
                if depth > max(CRESTree.depths):
                    depth = max(CRESTree.depths)
                
                rank = CRESTree.depths[depth]
                
                #Remove explicit rank from greengenes node since not useful
                dsym = ""
                if self.GGMode and len(taxa[i])>3:                    
                    dsym = taxa[i][0:1]
                    if taxa[i][1:3]=="__": # Original GG
                        taxa[i] = taxa[i][3:]
                    elif taxa[i][1]==":":
                        taxa[i] = taxa[i][2:]
                                                                 
                
                #Fix parent of self ambigousity error
                if taxa[i] == parent.name or taxa[i].startswith("Unknown "):
                    taxa[i] = "%s (%s)" % (parent.name, rank)
                    print "DEBUG: -----------"
                    print "DEBUG", line[:-1]
                    print "DEBUG: same name as parent - adding %s" %taxa[i]                    
                    print "DEBUG", taxa
                    print "DEBUG: -----------"
                
                #Fix incertae sedis only issues (still in 106)
                if taxa[i].lower() == "incertae sedis" or taxa[i] == "Incertae_sedis":
                    taxa[i] = "%s incertae sedis" % parent.name                
                
                
                while (taxa[i] in self.nodeNames) and self.getParent(self.getNode(taxa[i])) is not parent:                     
                    
                        node = self.getNode(taxa[i])
                        # Set new parent if parent matches
                        #print "DEBUG: Parent %s, node: %s" %(parent.name, node.name)
                        taxa[i] += " (%s)" % parent.name
                
                if taxa[i] in self.nodeNames:
                    parent = self.getNode(taxa[i])
                else:
                    parent = self.addNode(nodename=taxa[i],parent=parent,
                                      assignmentMin=alim)
                 
                if self.GGMode and dsym:
                    self.ranks[taxa[i]] = ggCodes[dsym]
                 
            # Add accession to parent
            if not accession in self.accessions:
                self.accessions[accession] = parent
            else:
                #self.rejected += accession
                sys.stderr.write("Warning: Accession %s occurs more than once!\n" % accession)
                
        ndsFile.close()    

    def processChangesMetadata(self, changeFile):
        """Process metadata update in semicolon-sep format"""
        # Changing of rank can no longer be supported due to the explicit rank criteria
        changes = open(changeFile, "r")
        for line in changes:
            try:
                parts = line[:-1].split(";")
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
                
                else:          
                    # Move node
                    n = self.getNode(oldName)
                    if not n:
                        if self.getNode(newName):
                            sys.stderr.write("Warning: Already moved / renamed: %s" % oldName)
                        else:
                            sys.stderr.write("Error: Cannot find taxon %s\n" % oldName)
                            return                       
                    if newName:                        
                        self.renameNode(n, newName)
                        print "Renamed %s to %s" % (oldName, newName)
                    if newParent and not (newParent == self.getParent(n).name):
                        print ("Moving %s from %s to %s" %
                               (oldName, self.getParent(n).name, newParent))
                        self.moveNode(n, self.getNode(newParent))
                
            elif not newName and not newParent:
                # Remove node
                if not self.getNode(oldName):
                    sys.stderr.write("Warning: Taxon %s not found and cannot be deleted" % oldName)
                else:
                    self.deleteNode(self.getNode(oldName), False)
                    print "Deleted taxon: %s" % oldName   
            

    
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
        # child node is deeper or shallower than it should be. Add or remove 
        # Parents. Remember parents of deleted for later reassignment of accessions
        
        parentsToDeleted = {} # name as key, parent node as value
        #print "DEBUG: Root=",self.root
        delQ = Queue.Queue()
        delQOrganelle = Queue.Queue()
        if self.rearrangeOrganelles:
            delQ.put(self.mainGenome)
            delQOrganelle.put(self.mito)
            delQOrganelle.put(self.plastid)
        else: 
            for c in self.tree.root.clades:
                delQ.put(c)
        while not delQ.empty() or not delQOrganelle.empty():
            if not delQ.empty(): 
                current = delQ.get()
            else:
                current = delQOrganelle.get()            
            
            #print "DEBUG:",current
            properRank = self._getProperRank(current)
             
            while properRank and self.getDepth(current) > properRank:
                tp = [str(n) for n in self.getPath(current)] #DEBUG
                print "DEBUG: Too low: %s" % tp 
                parent = self.getParent(current)
                grandparent = self.getParent(parent)                    
                self.deleteNode(parent, moveUpChildren=True)
                parentsToDeleted[parent.name] = grandparent
                if properRank>=max(sfLimits):
                    self.assignmentMin[current.name] = sfLimits[max(sfLimits)]
                elif properRank<min(sfLimits):
                    self.assignmentMin[current.name] = 0
                else:
                    self.assignmentMin[current.name] = sfLimits[properRank]
                
                tp = [str(n) for n in self.getPath(current)] #DEBUG
                print "DEBUG: --> %s" % tp
            
            while properRank and self.getDepth(current) < properRank:
                tp = [str(n) for n in self.getPath(current)] #DEBUG
                print "DEBUG: Too high: %s" % tp
                
                parent = self.getParent(current)
                actualRank = self.getDepth(current)
                actualRankName = CRESTree.depths[actualRank]
                if parent.name.startswith("Bacteria") and actualRank<CRESTree.PHYLUM:
                    parentName = "Bacteria (%s)" % actualRankName
                elif parent.name.startswith("Archaea") and actualRank<CRESTree.PHYLUM:
                    parentName = "Archaea (%s)" % actualRankName
                elif properRank == CRESTree.SPECIES and " " in current.name and not "(" in current.name:
                    genusName = current.name[:current.name.find(" ")]                    
                    if actualRank == CRESTree.GENUS and genusName not in self.nodeNames:
                        parentName = genusName
                    else:
                        parentName = "%s %s incertae sedis" % (current.name, actualRankName)
                elif " (Chloroplast)" in current.name:
                    parentName = self._adaptFromMainGenome(node=current, rank=properRank, 
                                                           add=" (Chloroplast)",genomeType=self.plastid)
                elif " (Mitochondrion)" in current.name:
                    parentName = self._adaptFromMainGenome(node=current, rank=properRank, 
                                                           add=" (Mitochondrion)",genomeType=self.mito)
                        
                else:
                    parentName = "%s %s incertae sedis" %(current.name, actualRankName)
                
                if parentName == self.getParent(current).name:
                    break
                
                elif parentName in self.nodeNames:
                    intermediate = self.getNode(parentName)
                else:
                    intermediate = self.addNode(parentName, parent, sfLimits[properRank-1])
                self.moveNode(current, intermediate)
                self.assignmentMin[current.name] = sfLimits[properRank]
                
                tp = [str(n) for n in self.getPath(current)] #DEBUG                
                print "DEBUG: --> %s" % tp      
                    
            for child in self.getImmediateChildren(current):
                delQ.put(child) 
        
        # Go through accessions and check if that node is if in deletedNodes, if so re-assign to parent
        # (until parent is no longer in deletedNodes itself)        
        for acc in self.accessions:                      
            while self.accessions[acc].name in parentsToDeleted:
                self.accessions[acc] = parentsToDeleted[self.accessions[acc].name]
        
        
         #self._addIntermediateRanks(self.tree.root)         #Must be done in parallel! :(
         
    def _adaptFromMainGenome(self, node, rank, add, genomeType):        
        
        newParent = genomeType
        realName = node.name.replace(add,"")
        if not realName in self.nodeNames:
            print "DEBUG: Not in 18S taxonomy"
            return self.getParent(node).name
        realPath = self.getPath(self.getNode(realName))
        nodePath = self.getPath(node)
        for i in range(1,min(rank,len(realPath))-1):
            pName = realPath[i].name+add
            if len(nodePath)> i and pName == nodePath[i].name:
                newParent= nodePath[i]
                print "DEBUG: Agreement: %s" % (nodePath[i].name)
            else:
                if len(nodePath)>i:
                    print "DEBUG: 18S: %s, 16S: %s" %(realPath[i].name, nodePath[i].name)
                if pName in self.nodeNames:
                    newParent = self.getNode(pName)
                    print "DEBUG: New Parent set to existing: %s" % pName
                else:
                    newParent = self.addNode(pName, newParent, 0)
                    print "DEBUG: New Parent added: %s" % pName
                
        for i in range(len(realPath)-1,rank-1):
            nn = "%s %s incertae sedis" % (node.name, CRESTree.depths[i])
            if nn in self.nodeNames:
                newParent = self.getNode(nn)
                print "DEBUG: New Parent set to existing: %s" % nn
            else:
                newParent = self.addNode(nn,newParent,0)
                print "DEBUG: New Parent added: %s" % nn
        
        return newParent.name

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
            inBrackets = name[name.rfind("(")+1:name.rfind(")")]
            if not inBrackets in CRESTree.depths.values()+["superphylum","Chloroplast"+"Mitochondrion"]:
                name = name[:name.find(" (")]                        
        if name in self.ranks: # other tests incl. "ales" "aceae"
            #print "DEBUG: %s in rank file as %s" % (name, self.ranks[name])
            return self.ranks[name]
        elif self.isParent(node) and name[-4:] == "ales":
            return CRESTree.ORDER
        elif self.isParent(node) and name[-4:] == "ceae" and not name[-7:] == "phyceae":
            return CRESTree.FAMILY
        else: 
            return None
        
    def writeFasta(self, inFile, outFile, sintax=False):
        added = []    
        newFasta = open(outFile , 'w')
    
        for record in Bio.SeqIO.parse(open(inFile), "fasta"):            
        
            headerItems = record.description.split()
            acc = headerItems[0]
            if "." in acc:
                acc = acc[:acc.find(".")]
            
            #if not (acc in self.rejected):
            if  acc in self.accessions:
                if acc not in added:                    
                    if not sintax:
                        newFasta.write(">%s\n%s\n" % (acc, record.seq))
                    else:
                        sintaxpath = self.getSintaxFormattedTaxPath(self.accessions[acc])
                        newFasta.write(">%s;%s;\n%s\n" % (acc, sintaxpath, record.seq))
                    added.append(acc)
            else:
                sys.stderr.write("Warning: cannot find %s in taxonomy. "
                                 "Skipping!\n" % acc)
        newFasta.close()
        newFasta.close()
    
    def writeConfigFiles(self, mapFile, treeFile):
        """Writes metadata files for CREST in tabular format and (special)        
        @mapFile the map file name
        @treeFile Newick treeFile format file name
        """
        mapFile = open(mapFile, "w")
        treeList = [self.root.nodeID, ";"]
        
        self._writeMapAndTre(node=self.root, mapFile=mapFile, treeList=treeList)
                        
        # also write accessions in map file mapped to their nodes
        for acc in self.accessions:
            anode = self.accessions[acc]
            #assMin = self.assignmentMin[anode.name]
            mapFile.write("%s\t%s\t-1\t%s\n" % (anode.nodeID, acc, -1))
       
        mapFile.close()
        
        treeFile = open(treeFile, "w")
        treeFile.write(''.join([str(item) for item in treeList]))
        treeFile.close()
        
    def writeMEGANFiles(self, mapFile, treeFile):
        """Writes metadata files for MEGAN same format as CREST config 
        except all accessions are included in the tree and each .map file
        entry maps uniquely to the tree        
        @mapFile the map file name
        @treeFile Newick treeFile format file name        
        """
        mapFile = open(mapFile, "w")
        treeList = [self.root.nodeID, ";"]
        
        # Add all accessions as nodes
        for acc in self.accessions:
            self.addNode(acc, self.accessions[acc], -1)
        
        self._writeMapAndTre(node=self.root, mapFile=mapFile, treeList=treeList)
       
        mapFile.close()
        
        treeFile = open(treeFile, "w")
        treeFile.write(''.join([str(item) for item in treeList]))
        treeFile.close()


    def _writeMapAndTre(self, node, mapFile, treeList):
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
            self._writeMapAndTre(node=child, mapFile=mapFile, treeList=treeList)            
    

    # NDS helper methods
    def _useNameInTaxonomy(self, taxa, ncbi_name=None):               
        
        i=1
        if ncbi_name:
            if not " " in ncbi_name:
                print "DEBUG: discarding incomplete sp. name: %s, parent: %s" % (ncbi_name, taxa[-1])
                return False
            
            
            spParts = ncbi_name.split(" ")
            if taxa[-1] in ["Eukaryota", "Chloroplast", "Mitochondrion", "Mitochondria"]:
                # Silva v128 organelle
                for nsKey in ARBor.nonSpeciesKeys:
                    if nsKey in ncbi_name:
                        return False              
                taxa.append(spParts[0])
                 
            else:                
                while isInt(taxa[-i]) and len(taxa)>i:
                    i+=1 
                                       
                if (not taxa[-i].startswith(spParts[0])):
                    print "DEBUG: discarding sp. name with wrong genus: %s, parent: %s" % (ncbi_name, taxa[-i])
                    return False
        
            name = ncbi_name        
        else:
            name = taxa[-i]
        
        for nsKey in ARBor.nonSpeciesKeys:
            if nsKey in name:
                return False
        
        #print "DEBUG: keeping sp. name %s" % (name)
        return True
    
        
