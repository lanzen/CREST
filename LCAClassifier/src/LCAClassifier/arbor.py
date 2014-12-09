import sys
import re
import os

from LCAClassifier.taxa import Tree, Node, Read


class ARBor(Tree):
    """A Tree with specific functions for handling data exported from ARB"""

    nonSpeciesKeys = ["uncultured", "sp.", "unidentified", "metagenome",
                      "enrichment", "clone"]

    def __init__(self, name=None):
        Tree.__init__(self, Node(name="root", nodeID=1), name)
        self.cellOrg = Node(name="Cellular organisms",
                            parent=self.root, nodeID=2)
        self.addNode(self.cellOrg)
        self.idCount = 10000000
        self.problemNodes = []
        self.rejected = []
        self.phylumNames = []
        self.classNames = []
        self.orderNames = []
        self.familyNames = []
        self.genusNames = []
        
        try:
            rankDir = (os.environ['LCATaxonomyDir'] + "/ranks/")
    
            ncbiRanks = {rankDir + "all_ncbi_phylum": self.phylumNames,
                         rankDir + "all_ncbi_class": self.classNames,
                         rankDir + "all_ncbi_order": self.orderNames,
                         rankDir + "all_ncbi_family": self.familyNames,
                         rankDir + "all_ncbi_genus": self.genusNames}
    
            for filename in ncbiRanks.keys():
                readFile = open(filename, 'r')
                for line in readFile:
                    ncbiRanks[filename].append(line[:-1])
                readFile.close()
        except: 
            pass

    def deleteNode(self, node, moveUpChildren=False):
        if isinstance(node, str):
            nodename = node
            node = self.getNode(nodename)

        if (not node) or node not in self.nodes:
            if node:
                nodename = node.name
            sys.stderr.write("Node \"%s\" not found\n" % nodename)
            return

        if moveUpChildren:
            for child in node.children:
                self.moveNode(child, node.parent)
            if node.getAssignment():
                for read in node.getAssignment().reads:
                    node.parent.assignRead(read)
        else:
            self.removeReads(node)

        i = node.parent.children.index(node)
        del node.parent.children[i]

    def removeReads(self, node):
        """Remove reads from read list recursively."""
        if node.getAssignment():
                for read in node.getAssignment().reads:
                    self.rejected.append(read.name)
                    for node in node.children:
                        self.removeReads(node)

    def parseSilvaNDS(self, filename, eukaryotic=True):
        """Parses an NDS Export file from Silva with accession, silva tax. and
        NCBI Tax."""
        ndsFile = open(filename, 'r')

        for line in ndsFile:
            self._readNDSLine(line, eukaryotic)
        ndsFile.close()

        for problemNode in self.problemNodes:
            if not problemNode.children:
                self.deleteNode(problemNode)
            else:
                self.renameNode(problemNode,
                                "%s (%s)" % (problemNode.name,
                                             problemNode.parent.name))
        self.problemNodes = []

    def _readNDSLine(self, line, eukaryotic=True, altTax=True, GGMode=True):
        parts = line.split("\t")
        accession = parts[0]
        if "." in accession:
            accession = accession[:accession.find(".")]
        taxonomy = parts[1]
        ncbi_name = parts[2]
        plast = False
        mito = False
        if (eukaryotic and len(parts) > 3):
            # or (GGMode and "Chloroplast" in taxonomy)
            if "Chloroplast" in parts[3]:
                plast = True
            elif "mitochondria" in parts[3].lower():
                mito = True


        parent = self.cellOrg
        depth = Tree.META

        if len(parts) < 1:
            print "Problem:\n%s" % line
            return
        if GGMode:
            taxonomy = taxonomy.replace("; ", ";")
            taxa = re.split('[/;]', taxonomy)
        else:
            taxa = re.split('[_/;]', taxonomy)
        ncbi_name = ncbi_name.replace("\n", "")

        # Do not use eukaryotic reads in non-eukaryotic mode
        if not eukaryotic and not GGMode:
            if taxa[0] == "Eukarya" or ("Chloroplast" in taxa) \
                or ("mitochondria" in taxa):
                self.rejected.append(accession)
                return

        # Ignore uncultured groups and handle like clustered to parent.
        if taxa[-1] == "uncultured":
            taxa = taxa[:-1]

        # Alt tax. fix
        if eukaryotic and altTax:
            if accession in self.rejected:
                del self.rejected[self.rejected.index(accession)]
            alt = ["Eukaryota"]
            if plast or (" plastid" in ncbi_name.lower()) or \
            ("chloroplast" in ncbi_name.lower()):
                alt.append("Plastid")
                plast = True
            elif  mito or ("mitochondrion" in ncbi_name.lower()):
                alt.append("Mitochondrion")
                mito = True
            else:
                alt.append("Nucleus")

            #Put extra labels on all childs of these new groups
            extra = ""
            if mito:
                extra = " (Mitochondrion)"
            elif plast:
                extra = " (Plastid)"
            for taxon in taxa[1:]:
                if len(taxon) > 0:
                    alt.append(taxon + extra)
            taxa = alt

        #Use NCBI Taxonomy species name if meaningful
        species = True
        for key in ARBor.nonSpeciesKeys:
            if key in ncbi_name:
                species = False
                break

        if species:
            if eukaryotic:
                if (plast or mito or
                    (" plastid" in ncbi_name.lower()) or
                    (" mitochondrion" in ncbi_name.lower())):
                    name_only = ncbi_name.replace(" Plastid", "")
                    name_only = name_only.replace(" plastid", "")
                    name_only = name_only.replace(" Mitochondrion", "")
                    name_only = name_only.replace(" mitochondrion", "")
                else:
                    name_only = ncbi_name
                    ncbi_name = name_only + " nucleus"
                if not altTax:
                    taxa.append(name_only)

            taxa.append(ncbi_name)
            # --- Remove line to use synonyms file instead
            taxa.append(accession)  # Makes synoynms file redundant.
            # ---
        else:
            #Not species
            taxa.append(accession)  # Remove line to use synonyms file instead
        
        i = 0
        parents = []
        for taxon in taxa:
            # If we are in GGMode then there will be a species name and accession 
            #that should not be touched:
            if (GGMode and len(taxon) > 3 and not
                ((taxon is taxa[-1]) or (taxon is taxa[-2] and species))):
                
                GGfix = True
                taxon = taxon[3:]
                dsym = taxon[0]
            else:
                GGfix = False

            #Check if taxon is numerical or empty
            if len(re.findall('[0-9]', taxon)) == len(taxon) and not GGMode:
                pass
            else:
                depth += 1

                #Fix parent of self ambigousity error
                if taxon in parents:
                    if depth > Tree.SPECIES or taxa[0] == 'Eukarya':
                        taxon += " (group)"
                    else:
                        taxon += " (%s)" % Tree.depths[depth]

                #Fix incertae sedis only issues (still in 106)
                if taxon == "Incertae Sedis" or taxon == "Incertae_sedis":
                    taxon = "%s Incertae Sedis" % parent.name

                #Find or create node
                node = self.getNode(taxon)
                if node:
                    # Fix parent conflicts
                    if node.parent is not parent:
                        if node not in self.problemNodes:
                            sys.stderr.write("Warning: Conflicting placement "
                                             "of node %s" % taxon)
                            sys.stderr.write("( parent: %s) Parent of "
                                             "existing: %s\n\n" %
                                             (node.parent, parent))
                            self.problemNodes.append(node)

                        # Two accession number nodes with different taxonomy
                        # not allowed. Skip these.
                        if taxon == accession:
                            return

                        # Otherwise try to find different name
                        altname = "%s (%s)" % (taxon, parent.name)
                        node = self.getNode(altname)
                        alt = 2

                        # Backup plan if we are still not ok
                        while node and not node.parent is parent:
                            altname = "%s (%s)" % (taxon, alt)
                            alt += 1
                            node = self.getNode(altname)

                        # Add new node with alternative name.
                        if not node:
                            node = Node(name=altname, parent=parent,
                                        nodeID=self.newID(), depth=depth)
                            self.addNode(node)
                    else:
                        depth = node.getHighestRank()
                else:
                    # Fix depth issues
                    rank = depth

                    if (taxon == "Plastid" or
                        taxon == "Nucleus" or
                        taxon == "Mitochondrion"):
                        rank = Tree.NORANK
                        depth = Tree.DOMAIN

                    #Find depth
                    elif taxon is taxa[-1] and depth < Tree.SPECIES:
                        if not species:
                            rank = Tree.SPECIES
                        else:
                            rank = Tree.SUBSPECIES
                    elif (taxon is taxa[-2] and
                          species and
                          depth < Tree.SPECIES):
                        rank = Tree.SPECIES
                    elif taxon is taxa[0]:
                        depth = rank = Tree.DOMAIN
                    elif GGfix:
                        if dsym == "k":
                            rank = Tree.DOMAIN
                        elif dsym == "p":
                            rank = Tree.PHYLUM
                        elif dsym == "c":
                            rank = Tree.CLASS
                        elif dsym == "o":
                            rank = Tree.ORDER
                        elif dsym == "f":
                            rank = Tree.FAMILY
                        elif dsym == "g":
                            rank = Tree.GENUS
                        elif dsym =="s":
                            rank = Tree.SPECIES
                    elif taxon in self.phylumNames:
                        depth = rank = Tree.PHYLUM
                    elif taxon in self.classNames:
                        depth = rank = Tree.CLASS
                    elif taxon in self.orderNames or taxon[-4:] == "ales":
                        depth = rank = Tree.ORDER
                    elif taxon in self.familyNames or taxon[-4:] == "ceae":
                        depth = rank = Tree.FAMILY
                    elif taxon in self.genusNames:
                        depth = rank = Tree.GENUS

                    p = parent
                    if not GGMode:
                        while (rank > Tree.META and
                               p.depth and
                               p.depth >= rank):
                            p.depth = Tree.NORANK
                            if p.parent:
                                p = p.parent
                                

                    #Add node
                    node = Node(name=taxon, parent=parent,
                                nodeID=self.newID(), depth=rank)
                    self.addNode(node)

                parent = node
                parents.append(node.name)

                if taxon is taxa[-1]:
                    #Associate accession with taxa
                    node.assignRead(Read(accession))
                i += 1

    def processChangesMetadata(self, line):
        """Process one line of metadata update in tab-sep format"""
        try:
            parts = line.split(";")
            #process parts
            oldName = parts[0]
            newName = parts[1]
            newParent = parts[2]
            rank = parts[3]
        except:
            sys.stderr.write("Warning: Incorrect manual change line: %s\n" %
                             line)
            return

        parent = self.getNode(newParent)
        if newParent and not parent:
            sys.stderr.write("Cannot find parent node %s\n" % newParent)
            return

        #___Fix depth___
        depth = 0
        #if rank is specified: use this
        if rank:
            for d in ARBor.depths.keys():
                if ARBor.depths[d] == rank:
                    depth = d
                    break

        elif parent:
            depth = parent.depth + 1

        #____

        # Add new node
        if not oldName:
            if self.getNode(newName):
                print "Already present: %s" % newName
            else:
                if not parent:
                    print ("Warning: Cannot add taxon %s as parent %s "
                           "does not exist" % (newName, newParent))
                else:
                    newNode = Node(name=newName, parent=parent,
                                   nodeID=self.newID(), depth=depth)
                    self.addNode(newNode)
                    print "Added new taxon: %s" % newName

        # Remove node
        elif not newName and not newParent:
            if not self.getNode(oldName):
                print "Already deleted: %s" % oldName
            else:
                self.deleteNode(self.getNode(oldName), False)
                print "Deleted taxon: %s" % oldName

        # Move or rename taxon
        elif not ".." in newName:
            n = self.getNode(oldName)
            if not self.getNode(oldName):
                if self.getNode(newName):
                    print "Already moved / renamed: %s" % oldName
                else:
                    sys.stderr.write("Cannot find node: %s\n" % oldName)
            else:
                if newName:
                    self.renameNode(n, newName)
                    print "Renamed %s to %s" % (oldName, newName)
                if newParent and not (newParent == n.parent.name):
                    print ("Moving %s from %s to %s" %
                           (oldName, n.parent.name, newParent))
                    self.moveNode(n, self.getNode(newParent))
                    n.depth = depth

        #Control shorthand annotation with ..
        else:
            twoNodes = oldName.split("..")
            firstParentName = self.getNode(twoNodes[0]).parent.name
            secParentName = self.getNode(twoNodes[1]).parent.name
            if not (firstParentName == parent and secParentName == parent):
                sys.stderr.write("Warning: Taxons %s not moved properly in "
                                 "NDS file!!\n" % oldName)

    def writeConfigFiles(self, map_=None, synonyms=None, tree=None,
                         rdpTraining=None):
        """Writes metadata files for MEGAN in tabular format and (special)
        Newick tree format.
        """
        if map_:
            map_ = open(map_, "w")
        if synonyms:
            synonyms = open(synonyms, "w")
        if rdpTraining:
            rdp = open(rdpTraining, "w")
        else:
            rdp = None

        # Need to find a list with good insertion options here.
        # One item = separator or nodeID.
        if tree:
            treeList = [self.root.nodeID, ";"]
        else:
            treeList = None

        self._writeMEGANData(node=self.root, map_=map_,
                             syn=synonyms, treeList=treeList, rdp=rdp)

        if map_:
            map_.close()
        if synonyms:
            synonyms.close()

        if tree:
            tree = open(tree, "w")
            tree.write(''.join([str(item) for item in treeList]))
            tree.close()

        if rdpTraining:
            rdp.close()

    def _writeMEGANData(self, node, map_, syn, treeList, rdp):
        # Self recursive, depth-first type
        if map_:
            if node.depth > Tree.SUBSPECIES:
                rankCode = Tree.mapCodes[Tree.SUBSPECIES]
            else:
                try:
                    rankCode = Tree.mapCodes[node.depth]
                except:
                    sys.stderr.write("Error: Depth key not found: %s\n"
                                     % node.depth)
                    rankCode = Tree.mapCodes[Tree.NORANK]

            map_.write("%s\t%s\t-1\t%s\n" % (node.nodeID, node.name, rankCode))

        if rdp and node.children and node.parent and node.depth < Tree.SPECIES:
            depth = len(node.getPhylogeny()) - 1
            if node.parent.isRoot():
                pID = 0
            else:
                pID = node.parent.nodeID
            if (not node.depth) or node.depth == Tree.NORANK:
                rankName = "norank"
            elif node.depth < Tree.GENUS and not node.children[0].children:
                node.depth = Tree.GENUS
                rankName = "genus"
            else:
                rankName = Tree.depths[node.depth]
            rdp.write("%s*%s*%s*%s*%s\n" %
                      (node.nodeID, node.name, pID, depth - 1, rankName))

        if syn and node.getAssignment():
            for read in node.getAssignment().reads:
                syn.write("%s\t%s\n" % (read.name, node.nodeID))

        if treeList:
            # Example:
            # 1. Root;
            # 2. (Bacteria)Root
            # 3. ((Firmicutes)Bacteria)Root;
            # 4. ((Firmicutes,Proteobacteria)Bacteria)Root;

            if node is self.root:
                pass
            else:
                i = treeList.index(node.parent.nodeID)
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

        for child in node.children:
            self._writeMEGANData(node=child, map_=map_, syn=syn,
                                 treeList=treeList, rdp=rdp)
