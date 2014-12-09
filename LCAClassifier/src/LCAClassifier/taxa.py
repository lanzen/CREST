import sys
import re
from types import StringType, IntType

from numpy import longfloat

import biom.table


class Tree:

    NORANK = 0
    META = 1
    DOMAIN = 2
    PHYLUM = 3
    CLASS = 4
    ORDER = 5
    FAMILY = 6
    GENUS = 7
    SPECIES = 8
    SUBSPECIES = 9
    SUB2 = 10

    depths = {META: 'base',
              DOMAIN: 'domain',
              PHYLUM: 'phylum',
              CLASS: 'class',
              ORDER: 'order',
              FAMILY: 'family',
              GENUS: 'genus',
              SPECIES: 'species',
              #SUBSPECIES: 'subspecies',
              #SUB2: 'sub2',
              }

    mapCodes = {NORANK: 0,
                META: 0,
            DOMAIN: 0,
            PHYLUM: 2,
            CLASS: 3,
            ORDER: 4,
            FAMILY: 5,
            GENUS: 98,
            SPECIES: 100,
            SUBSPECIES: 101}

    """A node holder for retrieving and manipulating nodes more easily"""

    def __init__(self, rootNode=None, name=None):
        if not rootNode:
            rootNode = Node(name="root", nodeID=1)
        self.root = rootNode
        self.name = name
        self.nodes = {self.root.name: self.root}
        self.nodesNoBrackets = {}
        self.idCount = 10000000
        

    def addNode(self, node):
        self.nodes[node.name] = node

    def getNode(self, name):
        if name in self.nodes:
            return self.nodes[name]
        return None

    def getCloseNode(self, searchName):
        """Return the node with a name close to searchName, first by ignoring
        words in brackets in Silva, then by cutting words off in the end of
        the searchName"""

        if not self.nodesNoBrackets:
            for nn in self.nodes.keys():
                nbName = nn
                if "(" in nn:
                    nbName = nn[0:nn.find("(") - 1]
                    self.nodesNoBrackets[nbName] = self.getNode(nn)
                else:
                    self.nodesNoBrackets[nn] = self.getNode(nn)

        while (" " in searchName and not searchName in self.nodesNoBrackets):
            space = searchName.rfind(" ")
            searchName = searchName[0:space]

        if searchName in self.nodesNoBrackets:
            return self.getNode(self.nodesNoBrackets[searchName])

        return None

    def newID(self):
        self.idCount += 1
        return self.idCount

    def printAsTree(self, popDataset=False, showLeaves=True, printFile=None):
        """Print as tree.

        If popDataste is None or a dataset name, the no. assignments will be
        printed.
        """
        _printAsTree(self.root, popDataset, showLeaves, printFile=printFile)

    def moveNode(self, node, newParent):
        if type(node) is StringType:
            node = self.getNode(node)
        i = node.parent.children.index(node)
        del node.parent.children[i]
        if type(newParent) is StringType:
            newParent = self.getNode(newParent)

        node.parent = newParent
        newParent.addChild(node)

    def renameNode(self, node, newName):
        """Rename node to newName. Updates self.nodes dictionary"""
        if type(node) is StringType:
            nodename = node
            node = self.getNode(nodename)
            if not node:
                sys.stderr.write("Node \"%s\" not found\n" % nodename)
                return
        del self.nodes[node.name]
        node.name = newName
        self.addNode(node)

    def deleteNode(self, node, moveUpChildren=False):
        """Removes node from tree structure. Still retains the disconnected
        node in self.nodes dict.

        @param moveUpChildren: if true, move up the children of the node to
                               its parent before deleting
        """
        #TODO test without retaining node

        if type(node) is StringType:
            nodename = node
            node = self.getNode(nodename)
            if not node:
                sys.stderr.write("Node \"%s\" not found\n" % nodename)
                return

        if moveUpChildren:
            for child in node.children:
                self.moveNode(child, node.parent)

        i = node.parent.children.index(node)
        del node.parent.children[i]
        del self.nodes[node.name]

    def mergeNodes(self, nodeList, newName, newID=None):
        """Creates a new node with newName, then moves all child nodes of the
        nodes in the nodeList to the new node.

        Checks that all nodes in nodeList has the same parent.
        """
        kids = []
        reads = []

        first = nodeList[0]
        parent = first.parent
        depth = first.depth
        if first.population is None:
            pop = None
        else:
            pop = 0
        if first.singles is None:
            singles = None
        else:
            singles = 0
        if first.doubles is None:
            doubles = None
        else:
            doubles = 0

        for node in nodeList:
            if node.parent != parent:
                return "Nodes could not be merged due to different parents"

            for read in node.reads:
                reads.append(read)

            for child in node.children:
                kids.append(child)

            if pop is not None:
                pop += node.population
            if singles is not None:
                singles += node.population
            if doubles is not None:
                doubles += node.population

        if len(reads) == 0:
            reads = None

        newNode = Node(name=newName, parent=parent, population=pop,
                       reads=reads, singles=singles, doubles=doubles,
                       nodeID=newID, depth=depth)
        self.addNode(newNode)

        for kid in kids:
            self.moveNode(kid, newNode)
        for node in nodeList:
            self.deleteNode(node)

    def pruneUnassigned(self):
        """Removes all nodes except those with names matching nodesAssigned"""
        toRemove = []
        self.root._pruneUnassigned(toRemove)
        for node in toRemove:
            try:
                self.deleteNode(node, moveUpChildren=False)
            except:
                sys.stderr.write("Weird error: Could not remove node %s\n" %
                                 node.name)

    def initFrom(self, mapFile, treFile, synonymsFile=None):
        """Reads a MEGAN .tre, .theMap and (optionally) synonyms file and adds
        to self from this under root node"""

        theMap = open(mapFile, "r")

        codeMaps = {}
        for key in Tree.mapCodes.keys():
            codeMaps[Tree.mapCodes[key]] = key
        # Beacause of unfortunate use of NORANK code for MEGAN compability:
        codeMaps[0] = Tree.NORANK

        self.nodeIDs = {}
        self.nodeIDs[self.root.nodeID] = self.root.name
        self.ranks = {}
        for line in theMap:
            parts = line.split("\t")
            nodeID = int(parts[0])
            self.nodeIDs[nodeID] = parts[1]
            try:
                self.ranks[nodeID] = codeMaps[int(parts[-1])]
            except:
                sys.stderr.write("Error: Depth key not found: %s\n" %
                                 int(parts[-1]))
                self.ranks[nodeID] = None
        theMap.close()

        tre = open(treFile, 'r')

        tree = ''
        for line in tre:
            tree += line
        tre.close()

        parent = self.root
        this = ''
        for i in range(1, len(line) + 1):
            char = line[-i]
            if char == ';':
                pass
            elif char == ',' or char == '(' or char == ')':
                if not this == '':
                    nodeID = int(this)
                    this = ''
                    if nodeID == 1:
                        n = self.root
                    else:
                        d = self.ranks[nodeID]
                        # Fix for MEGAN compability forcing all ranks < phylum
                        # to be 2
                        if d == 0:
                            pl = len(parent.getPhylogeny())
                            if pl < Tree.PHYLUM:
                                d = pl
                        #--Fix--
                        n = Node(name=self.nodeIDs[nodeID], depth=d,
                                 nodeID=nodeID, parent=parent)
                        self.addNode(n)

                if char == '(':
                    parent = parent.parent
                elif char == ')':
                    parent = n
            else:
                this = char + this

        if synonymsFile:
            syn = open(synonymsFile, "r")
            for line in syn:
                parts = line.split("\t")
                sName = parts[0]
                nodeID = int(parts[1])
                name = self.nodeIDs[nodeID]
                if name != sName:
                    synNode = Node(name=sName, depth=self.ranks[nodeID],
                                   nodeID=self.newID(),
                                   parent=self.getNode(name))
                    self.addNode(synNode)

    def initFromIndentedText(self, indFile, ncbiMode=False):
        treefile = open(indFile, 'r')
        lastNode = self.root
        lastIndent = -1
        i = 0
        sub = 0

        for line in treefile:
            i += 1

            if len(line) > 1:
                ind = 0
                while line[ind] == ' ':
                    ind += 1
                name = line[ind:].replace("\n", "")
                if ":" in name:
                    name = name[:name.find(":")]

                if sub > 0:
                    if ind <= sub:
                        sub = 0
                    else:
                        ind -= 2

                if  (ncbiMode and
                     (sub == 0 and
                      ('group' in name or
                       'subdivisions' in name))):
                    sub = ind
                elif ind == 0:
                    self.root.name = name
                else:
                    ind = ind / 2
                    for d in range(lastIndent + 1 - ind):
                        lastNode = lastNode.parent

                    thisNode = Node(name=name, parent=lastNode)
                    self.addNode(thisNode)

                    lastIndent = ind
                    lastNode = thisNode

    def getChildrenByRank(self, rank):
        """Recursively add taxa to rankChildren if:
        1) the rank is right, or
        2) the rank is too high, but the correct depth is missing (i e genera
        directly under domain)
        (only if there is no taxa in between, ow back one)
        3) the taxa are terminal and have no rank, but the correct rank is
        missing (only if the highest rank is one above that seaarched for)
        """
        rankChildren = []
        rankChildren = self.root._getChildrenByRank(rankChildren, rank)
        return rankChildren

    def printPopulationsAtDepth(self, depth=PHYLUM, outFile=None,
                                ignoreAcc=True, dataset=None,
                                normaliseToBase=True,
                                useRankData=True):
        """Print number assigned, unique and Chao estimate for each taxa in
        rank, tab-separated text format
        """

        if not type(depth) is IntType:
            #Find depth
            for i in self.dephts:
                if self.depths[i] == depth:
                    depth = i
                    break

        if not outFile:
            print "Parent taxa\tTaxa\tAbundance\tShare\tUnique\tChao"
        else:
            outFile.write("Parent taxa\tTaxa\tAbundance"
                          "\tShare\tUnique\tChao\n")

        total = 0
        totalU = 0
        if useRankData:
            nodes = self.getChildrenByRank(depth)
        else:
            nodes = self.root.getChildrenByDepth(depth=depth,
                                                 ignoreAcc=ignoreAcc)
            nodes = self._fixRanks(nodes, depth)

        if not self.root.getAssignment(dataset):
            sys.stderr.write("No assignments in %s\n" % dataset)
            return

        if normaliseToBase and depth > self.META:
            rootNode = self.getNode("Cellular organisms")
        else:
            rootNode = self.root

        rpop = rootNode.getAssignment(dataset).population

        for node in nodes:
            if "Unknown" in node.name:
                pass
            else:
                a = node.getAssignment(dataset)
                if a:
                    npop = a.population
                    nu = a.numberAssigned()
                    if not "Unknown" in node.name:
                        total += npop
                        totalU += nu
                else:
                    npop = 0
                    nu = 0

                percent = longfloat(npop) / longfloat(rpop)
                if not outFile:
                    print("%s\t%s\t%s\t%s\t%s\t%s" %
                          (node.parentPrintName, node.name,
                           npop, percent, nu, node.chaoEstimate(dataset)))
                else:
                    outFile.write("%s\t%s\t%s\t%s\t%s\t%s\n" %
                                  (node.parentPrintName, node.name, npop,
                                   percent, nu, node.chaoEstimate(dataset)))

        percent = longfloat(total) / longfloat(rpop)
        if not outFile:
            print "\tTotal\t%s\t%s\t%s" % (total, percent, totalU)
        else:
            outFile.write("\tTotal\t%s\t%s\t%s\n" % (total, percent, totalU))
        ut = rpop - total
        percent = longfloat(ut) / longfloat(rpop)
        utU = rootNode.getAssignment(dataset).numberAssigned() - totalU
        if not outFile:
            print ("\tUnclassified at %s level\t%s\t%s\t%s" %
                   (self.depths[depth], ut, percent, utU))
        else:
            outFile.write("\tUnclassified at %s level\t%s\t%s\t%s\n" %
                          (self.depths[depth], ut, percent, utU))

#     def printCompositionAlternative(self, outFile=None, datasets=[],
#                               outputNovel=False):
#         """Print a number of datasets in alternative format
#         (from Hakon 22-2-2012)       
#         """
#         
#         ## bugfix 17-1-2013: taxa with missing ranks repeated by assigned col.
#          
#         datasetheaders = ""
#         for dataset in datasets:
#             datasetheaders += ("Sum_" + dataset + "\t" +
#                                "Assigned_" + dataset + "\t")
#         datasetheaders = datasetheaders[:-1]
#         if not outFile:
#             print "Level\tTaxonpath\tTaxon\t%s" % datasetheaders
#         else:
#             outFile.write("Level\tTaxonpath\tTaxon\t%s\n" % datasetheaders)
#             
#         listedNodes = []
# 
#         for i in self.depths:
#             nodes = self.getChildrenByRank(i)
#             for node in nodes:
#                 if ("Unknown" in node.name and not outputNovel) or (node in listedNodes):
#                     pass
#                 else:
#                     listedNodes.append(node)
#                     row = "%s\t%s\t%s" % (self.depths[i],
#                                         node.getPhylogenyRDPStyle(),
#                                         node.name)
#                     for dataset in datasets:
#                         a = node.getAssignment(dataset)
#                         if a:
#                             setSum = a.population
#                             setAssigned = a.primaryPopulation()
#                         else:
#                             setSum = 0
#                             setAssigned = 0
#                         row += "\t%s\t%s" % (setSum, setAssigned)
# 
#                     if not outFile:
#                         print(row)
#                     else:
#                         outFile.write(row + "\n")
                        
    def printAssignedCount(self, outFile=None, datasets=[]):
                
        ## 2014-12-08 - due to problem with unknowns never counted in assigned
         
        datasetheaders = ""
        for dataset in datasets:
            datasetheaders += (dataset + "\t")
        datasetheaders = datasetheaders[:-1]
        if not outFile:
            print "Level\tTaxonpath\tTaxon\t%s" % datasetheaders
        else:
            outFile.write("Level\tTaxonpath\tTaxon\t%s\n" % datasetheaders)
            
        listedNodes = []

        for i in self.depths:
            nodes = self.getChildrenByRank(i)
            for node in nodes:
                if ("Unknown" in node.name) or (node in listedNodes):
                    pass
                else:
                    listedNodes.append(node)
                    row = "%s\t%s\t%s" % (self.depths[i],
                                        node.getPhylogenyRDPStyle(),
                                        node.name)
                    for dataset in datasets:
                        a = node.getAssignment(dataset)
                        if a:
                            setAssigned = a.primaryPopulation()
                            if i<Tree.SPECIES and not node.name == "No hits" :
                                uc_name = "Unknown %s %s" % (node.name, Tree.depths[node.depth+1])
                                unknown_children = node.getChildByName(1, uc_name)
                                if unknown_children:
                                    u_assigned = unknown_children.getAssignment(dataset).primaryPopulation()
                                    setAssigned += u_assigned                                    
                        else:
                            setAssigned = 0
                        row += "\t%s" % (setAssigned)

                    if not outFile:
                        print(row)
                    else:
                        outFile.write(row + "\n")
                        
    def printTotalCount(self, outFile=None, datasets=[]):
                
        ## 2014-12-08 - due to problem with unknowns never counted in assigned
         
        datasetheaders = ""
        for dataset in datasets:
            datasetheaders += (dataset + "\t")
        datasetheaders = datasetheaders[:-1]
        if not outFile:
            print "Level\tTaxonpath\tTaxon\t%s" % datasetheaders
        else:
            outFile.write("Level\tTaxonpath\tTaxon\t%s\n" % datasetheaders)
            
        listedNodes = []

        for i in self.depths:
            nodes = self.getChildrenByRank(i)
            for node in nodes:
                if ("Unknown" in node.name) or (node in listedNodes):
                    pass
                else:
                    listedNodes.append(node)
                    row = "%s\t%s\t%s" % (self.depths[i],
                                        node.getPhylogenyRDPStyle(),
                                        node.name)
                    for dataset in datasets:
                        a = node.getAssignment(dataset)
                        if a:
                            setSum = a.population                            
                        else:
                            setSum = 0                            
                        row += "\t%s" % (setSum)

                    if not outFile:
                        print(row)
                    else:
                        outFile.write(row + "\n")
                        
    def printRichness(self, outFile=None, datasets=[]):
        datasetheaders = ""
        for dataset in datasets:
            datasetheaders += (dataset + "\t")
        datasetheaders = datasetheaders[:-1]
        if not outFile:
            print "Level\tTaxonpath\tTaxon\t%s" % datasetheaders
        else:
            outFile.write("Level\tTaxonpath\tTaxon\t%s\n" % datasetheaders)
            
        listedNodes = []

        for i in self.depths:
            nodes = self.getChildrenByRank(i)
            for node in nodes:
                if ("Unknown" in node.name) or (node in listedNodes):
                    pass
                else:
                    listedNodes.append(node)
                    row = "%s\t%s\t%s" % (self.depths[i],
                                        node.getPhylogenyRDPStyle(),
                                        node.name)
                    for dataset in datasets:
                        a = node.getAssignment(dataset)
                        if a:
                            r = a.numberAssigned()
                        else:
                            r = 0
                        row += "\t%s" % r

                    if not outFile:
                        print(row)
                    else:
                        outFile.write(row + "\n")
    
    def printRelativeAbundances(self, outFile=None, datasets=[], 
                                minimalMaxAbundance=0.0, normaliseToBase=False):
        datasetheaders = ""
        for dataset in datasets:
            datasetheaders += (dataset + "\t")
        datasetheaders = datasetheaders[:-1]
        if not outFile:
            print "Level\tTaxonpath\tTaxon\t%s" % datasetheaders
        else:
            outFile.write("Level\tTaxonpath\tTaxon\t%s\n" % datasetheaders)
            
        listedNodes = []

        for i in self.depths:
            nodes = self.getChildrenByRank(i)
            for node in nodes:
                if ("Unknown" in node.name) or (node in listedNodes):
                    pass
                else:
                    listedNodes.append(node)
                    row = "%s\t%s\t%s" % (self.depths[i],
                                        node.getPhylogenyRDPStyle(),
                                        node.name)
                    
                    maxPop=0.0
                    for dataset in datasets:
                        if normaliseToBase:
                            rc = self.getNode("Cellular organisms")
                        else:
                            rc = self.root
                        all_assignments = rc.getAssignment(dataset)
                        if all_assignments:
                            
                            rpop = all_assignments.population
                        else:
                            rpop = 0
                            sys.stderr.write("Warning: Dataset %s has no assignments" %dataset)
                        a = node.getAssignment(dataset)
                        if a:
                            npop = a.population
                        else:
                            npop = 0
                        try:
                            ra = (longfloat(npop) / longfloat(rpop))
                        except:
                            sys.stderr.write("Problem with longfloat dividing %s with %s" % (npop, rpop))
                            ra = float ( npop / rpop )
                            sys.stderr.write("Using", ra)
                        if ra > maxPop:
                            maxPop = ra
                        row += "\t%f" % ra
                    if (maxPop >= minimalMaxAbundance):
                        if not outFile:
                            print(row)
                        else:
                            outFile.write(row + "\n")


    def printAssignmentsRDPStyle(self, dataset=None, printFile=None,
                                 newTabStyle=False):
        self.root._printAssignmentsRDPStyle(dataset, printFile, seq=False,
                                            newTabStyle=newTabStyle)

    def printAssignmentsRDPFasta(self, dataset=None, printFile=None,
                                 newTabStyle=False):
        self.root._printAssignmentsRDPStyle(dataset, printFile, seq=True,
                                            newTabStyle=newTabStyle)
        
        
    

    def _fixRanks(self, nodes, depth):
        # Old method to check depth when not known from map file
        nodesToRemove = []
        nodesAppend = []
        speciesPattern = re.compile('^[A-Z][a-z]+ [a-z]+')
        for node in nodes:

            i = 1
            newKids = node.children
            if (((depth == Tree.ORDER and node.name[-4:] != "ales") or
                 (depth == Tree.FAMILY and node.name[-4:] != "ceae")) and not
                ("incertae sedis" in node.name.lower() or
                 "division" in node.name.lower() or
                 "group" in node.name.lower() or
                 re.findall(speciesPattern, node.name))):

                while newKids:
                    if ((depth == Tree.ORDER and
                         newKids[0].name[-4:] == "ales") or
                        (depth == Tree.FAMILY and
                         newKids[0].name[-4:] == "ceae")):
                        nodesAppend += newKids
                        nodesToRemove.append(node)
                        newKids = None
                    elif re.findall(speciesPattern, newKids[0].name):
                        newKids = None
                    else:
                        kidsLeft = False
                        for newKid in newKids:
                            if newKid.children:
                                kidsLeft = True
                                break
                        if kidsLeft:
#                            print ("DEBUG: Checking for %s in first "
#                                   "child %s" %
#                                   (Tree.depths[depth], newKids[0].name)
                            i += 1
                            newKids = node.getChildrenByDepth(i)
                        else:
                            newKids = None

            elif (depth == Tree.GENUS and
                  node.parent.name[-4:] != "ceae" and not
                  ("incertae sedis" in node.name.lower() or
                   "division" in node.name.lower() or
                   "group" in node.name.lower() or
                   re.findall(speciesPattern, node.name) or
                   (node.children and
                    re.findall(speciesPattern, node.children[0].name)))):
                while newKids:
                    if (newKids[0].parent.name[-4:] == "ceae" or
                        (newKids[0].children and
                         re.findall(speciesPattern,
                                    newKids[0].children[0].name))):

                        nodesAppend += newKids
                        nodesToRemove.append(node)
                        newKids = None
                    elif re.findall(speciesPattern, newKids[0].name):
                        newKids = None
                    else:
                        kidsLeft = False
                        for newKid in newKids:
                            if newKid.children:
                                kidsLeft = True
                                break
                        if kidsLeft:
#                            print ("DEBUG: Checking for genus in first "
#                                   "child %s under %s" %
#                                   (newKids[0].name, newKids[0].parent.name))
                            i += 1
                            newKids = node.getChildrenByDepth(i)
                        else:
                            newKids = None

            elif (depth == Tree.SPECIES and not
                  ("incertae sedis" in node.name.lower() or
                   "division" in node.name.lower() or
                   "group" in node.name.lower()) and not
                  re.findall(speciesPattern, node.name)):
                while newKids:
                    if re.findall(speciesPattern, newKids[0].name):
                        nodesAppend += newKids
                        nodesToRemove.append(node)
                        newKids = None
                    else:
                        kidsLeft = False
                        for newKid in newKids:
                            if newKid.children:
                                kidsLeft = True
                                break
                        if kidsLeft:
                            i += 1
                            newKids = node.getChildrenByDepth(i)
                        else:
                            newKids = None

            if node.parent:
                node.parentPrintName = node.parent.name

        for node in nodesToRemove:
            del nodes[nodes.index(node)]
        for node in nodesAppend:
            pName = node.parent.name
            for level in node.getPhylogeny()[2:1 - depth]:
                pName = "%s : %s" % (level.name, pName)
            node.parentPrintName = pName
            nodes.append(node)
        return nodes


class Node:
    """Assigns new node. For root node parent is None."""

    def __init__(self, name=None, parent=None, nodeID=None, depth=None):

        self.name = name
        self.parent = parent
        self.assignments = {}
        self.children = []
        self.nodeID = nodeID
        self.depth = depth
        if parent:
            self.parent.addChild(childNode=self)

    def getAssignment(self, dataset=None):
        if dataset in self.assignments.keys():
            return self.assignments[dataset]
        else:
            return None
        
    def getPrimaryAssignedReads(self, dataset=None):
        """Return all reads assigned to dataset provided, or all reads
        with duplicates removes if dataset is None"""
        if dataset in self.assignments.keys():
            return self.assignments[dataset].primReads
        else:
            aReads = []
            for ds in self.assignments.keys():
                for read in self.assignments[ds].primReads:
                    if not read in aReads:
                        aReads.append(read)
            return aReads

    def assignRead(self, read, dataset=None, abundance=None,
                   primary=False, recursive=False):
        assignment = self.getAssignment(dataset)
        if assignment:
            assignment.addRead(read, readPopulation=abundance, primary=primary)
        else:
            assignment = Assignment(read, readPopulation=abundance, primary=primary)
            self.assignments[dataset] = assignment
        if recursive:
            if self.parent:
                self.parent.assignRead(read, dataset, abundance,
                                       primary=False, recursive=True)

    def addChild(self, childNode):
        self.children += ([childNode])

    def getChildrenByDepth(self, depth=1, ignoreAcc=True):
        accNumberPattern = re.compile('[A-Z]+\d\d\d\d\d')
        if depth == 0:
            return [self]
        tempChildren = self.children
        for i in range(1, depth):
            dChildren = []
            for child in tempChildren:
                dChildren += (child.children)
            tempChildren = dChildren

        if ignoreAcc:
            temp2 = []
            for node in tempChildren:
                if not re.findall(accNumberPattern, node.name):
                    temp2.append(node)
            tempChildren = temp2

        return tempChildren

    def getChildByName(self, depth, name):
        tc = self.getChildrenByDepth(depth)
        for c in tc:
            if c.name == name:
                return c

    def addRead(self, read):
        self.assignRead(read)

    def chaoEstimate(self, dataset=None):
        if self.getAssignment(dataset):
            return self.getAssignment(dataset).chaoEstimate()
        else:
            return None

    def isRoot(self):
        return not self.parent

    def getPhylogeny(self, limit=False):
        """Return all nodes until root in a list, in order of decreasing depth
        (ACC,Species,Genus..Root).
        """
        phyl = [self]
        while not phyl[-1].isRoot():
            phyl.append(phyl[-1].parent)

        if limit:
            while phyl[0].depth > limit:
                phyl = phyl[1:]

        return phyl
    
    def getPhylogenyNameList(self):
        "Return all node names until root in a list"
        phyl = self
        p_list=[]
        while not phyl.isRoot():
            p_list.append(phyl.name)
            phyl = phyl.parent
        p_list.reverse()
        return p_list

    def getPhylogenyRDPStyle(self, root=True, limit=False, newTabStyle=False):
        """Return all nodes until root in a semicolon-separated list"""

        parentNameList = ""
        phyl = self.getPhylogeny(limit)
        if not root:
            phyl = phyl[:-1]
        #if not last: phyl=phyl[1:]
        for p in phyl:
            parentNameList = "%s;%s" % (p.name, parentNameList)
        return parentNameList[:-1]
    
    def getHighestRank(self):
        """Return the highest rank to which this node is classified, assuming
        intermediate parents without rank.
        """
        maxRank = 0
        for node in self.getPhylogeny():
            if node.depth > maxRank:
                maxRank = node.depth
        return maxRank

    def _getChildrenByRank(self, rankChildren, rank):
        """Recursively add to rankChildren if the depth is right, too high
        or if it has no more kids and the depth is. Otherwise iterate with
        all all kids.
        """
        # If the rank is right, or no more children and the rank remains
        # unknown and not found (once)

        if ((self.depth and self.depth >= rank) or
            (((not self.depth) or
              self.depth == Tree.NORANK) and
             (not self.children) and
             self.getHighestRank() == rank - 1)):
                if self.parent:
                    p = self.parent
                    self.parentPrintName = p.name
                    while p.parent and (not p.depth or p.depth == Tree.NORANK):
                        p = p.parent
                        self.parentPrintName = (p.name + " : " +
                                                self.parentPrintName)
                else:
                    self.parentPrintName = ""
                rankChildren.append(self)
        else:
            "Not matching!"
            for child in self.children:
                rankChildren = child._getChildrenByRank(rankChildren, rank)
        return rankChildren

    def _pruneUnassigned(self, toRemove):
        #Exhaustive recursion child by child stopping when node removed
        if not self.isRoot() and not self.assignments:
            toRemove.append(self)
        else:
            for child in self.children:
                child._pruneUnassigned(toRemove)

    def _printAssignmentsRDPStyle(self, dataset, printFile,
                                  seq=None, newTabStyle=False):
        assignments = self.getAssignment(dataset)
        if assignments and assignments.primReads:

            for r in assignments.primReads:
                toPrint = ("%s\t%s" %
                           (r.name,
                            self.getPhylogenyRDPStyle(
                                      root=False, newTabStyle=newTabStyle)))
                if seq and r.seq:
                    toPrint = ">%s\n%s" % (toPrint, r.seq)
                if (seq and r.seq) or not seq:
                    if printFile:
                        printFile.write(toPrint + "\n")
                    else:
                        print toPrint
                elif not r.seq:
                    pass  # DEBUG: print "No seq for %s" % r.name
                else:
                    print "Undefined problem for seq %s" % r.seq
        if assignments:
            for child in self.children:
                child._printAssignmentsRDPStyle(dataset, printFile, seq=seq)
        

    def __repr__(self):
        if self.isRoot():
            return "Root Node: %s" % (self.name)
        else:
            return "Node: %s, parent: %s" % (self.name, self.parent.name)


class Read:
    """Represents a sequence read"""
    def __init__(self, name, seq=None):
        self.name = name
        self.seq = seq

    def __repr__(self):
        if not self.seq:
            return '>' + self.name
        else:
            return '>%s\n%s' % (self.name, self.seq)

class Assignment:
    """An assignment of a read to a node. The population is assumed to be
    indicated after the last "_" in the read name as for AmpliconNoise.
    Support for a specific dataset when several datasets are assigned to one
    tree.
    """

    def __init__(self, read, readPopulation=None, primary=False):
        self.singletons = None
        self.doubletons = None
        self.population = 0
        self.reads = []
        self.primReads = []
        self.addRead(read, readPopulation, primary)

    def addRead(self, read, readPopulation, primary=False):
        self.reads += [read]
        if primary:
            self.primReads += [read]


        if readPopulation:
            if not self.singletons:
                self.singletons = 0
            if not self.doubletons:
                self.doubletons = 0

            self.population += readPopulation
            if readPopulation == 1:
                self.singletons += 1
            elif readPopulation == 2:
                self.doubletons += 1
        else:
            self.population += 1
            self.singletons = None
            self.doubletons = None

    def numberAssigned(self):
        return len(self.reads)

    def primaryPopulation(self):
        primPop = 0
        for pr in self.primReads:
            if "_" in pr.name:
                try:
                    primPop += int(pr.name.split("_")[-1].replace(".00", ""))
                except:
                    primPop += 1
            elif "numreads=" in pr.name:
                primPop += int(pr.name[pr.name.find("numreads=") +
                                       len("numreads="):])
            else:
                primPop += 1
        return primPop

    def __len__(self):
        return len(self.reads)

    def chaoEstimate(self):
        if (self.singletons != None and self.doubletons):
            return ((self.singletons * self.singletons) /
                    (2 * self.doubletons)) + self.numberAssigned()
        else:
            return None


def _printAsTree(node, pop, showLeaves, space="", printFile=None):
    if showLeaves or node.children:
        if pop or pop is None:
            a = node.getAssignment(pop)
            if a:
                s = "%s%s: %s" % (space, node.name, a.population)
                if printFile:
                    printFile.write(s + "\n")
                else:
                    print s
        else:
            s = "%s%s" % (space, node.name)
            if printFile:
                printFile.write(s + "\n")
            else:
                print s

    if node.children:
        for c in node.children:
            _printAsTree(c, pop, showLeaves, space=space + "  ",
                         printFile=printFile)
