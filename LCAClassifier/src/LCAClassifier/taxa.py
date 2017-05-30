import sys
from types import StringType, IntType
from Bio import Phylo


class CRESTree:
    """ Instances have specific methods and attributes needed for classification and 
    stores the actual CachingTree instance as self.tree"""

    ROOT = 0
    META = 1
    DOMAIN = 2
    SUPERKINGDOM = 3
    KINGDOM = 4
    PHYLUM = 5
    CLASS = 6
    ORDER = 7
    FAMILY = 8
    GENUS = 9
    SPECIES = 10
    SUBSPECIES = 11

    depths = {ROOT: 'root',
              META: 'meta',
              DOMAIN: 'domain',
              SUPERKINGDOM: 'superkingdom',
              KINGDOM: 'kingdom',
              PHYLUM: 'phylum',
              CLASS: 'class',
              ORDER: 'order',
              FAMILY: 'family',
              GENUS: 'genus',
              SPECIES: 'species',
              SUBSPECIES: 'strain',
              }



    """Reference tree initiated from newick and map file"""

    def __init__(self, rootedTree=None):
        
        if rootedTree:
            self.tree = rootedTree                        
            self.root = self.tree.root     
        else:
            self.tree = Phylo.BaseTree.Tree(name="Tree")
            self.root = Phylo.BaseTree.Clade(name="root")
            self.tree.root = self.root
                
        self.nodeNames = {} # name as key, node as value
        self.nodeIDs = {} # id as key, node as value
        self.assignmentMin = {}
        self.parents = {} # dict for caching parents
        
        
               
    def getNode(self, name):
        if not name in self.nodeNames:
            sys.stderr.write("Error: Node %s not found!\n" % name)
            return None
        else:
            return self.nodeNames[name]
        
    def getNodeByID(self, id):
        if id in self.nodeIDs:
            return self.nodeIDs[id]
        else:
            return None        
        
    def verifyNode(self, node):
        """Returns proper node if instead StringType is given as argument"""
        if isinstance(node, Phylo.BaseTree.Clade):
            return node
        elif type(node) is StringType:
            nn = node
            node = self.getNode(nn)
            if not node:
                sys.stderr.write("Verification Error: Node \"%s\" not found\n" % nn)
                return
            else:
                return node
        else: 
            sys.stderr.write("Verification Error: Node %s is a strange instance \n" %node)
            return        
        

    def addNode(self, nodename, parent, assignmentMin=0):
        #insert instead of append?
        if nodename in self.nodeNames:
            sys.stderr.write("Node name \"%s\" is not unique - not added\n" % nodename)
            return None
        node = Phylo.Newick.Clade(name=nodename)
        parent = self.verifyNode(parent)
        if not parent:
            return None
        parent.clades.append(node)
        self.nodeNames[nodename] = node
        self.assignmentMin[node.name] = assignmentMin
        self.parents[node] = parent
        return node
        
    def getImmediateChildren(self, node):
        node = self.verifyNode(node)
        children = []
        for c in node.clades:
            children.append(c)
        return children
    
    def getParent(self, node):
        node = self.verifyNode(node)
        if not node:
            return None
        
        if node in self.parents:
            return self.parents[node] 
        
        else:        
            p = self.tree.get_path(node)
            if p and len(p) > 1:
                parent = p[-2]
                self.parents[node] = parent
                return parent
            else:
                return self.tree.root
    
    def getRank(self, node):
        """Return the name of the rank of the node"""
        depth = self.getDepth(node)        
        if depth > CRESTree.SUBSPECIES:
            depth = CRESTree.SUBSPECIES
        return CRESTree.depths[depth]
    
    def getDepth(self, node):
        """Return the depth of the rank of the node"""
        node = self.verifyNode(node)
        pth = self.getPath(node)
        if len(pth)==1 and pth[0] is self.tree.root:
            return 0
        else:
            return len(pth)
    
    def getPath(self, node):
        """ Uses cached parent info if available to look up bath backwards instead of the 
        slow implementation of Phylo"""
        plist = [node]
        if node is self.tree.root:
            return plist
        parent = self.getParent(node)
        while parent and parent is not self.tree.root:
            plist = [parent] + plist
            parent = self.getParent(parent)
        if not parent:
            sys.stderr.write("Cannot find parent beyond %s\n" % plist)
        else:
            return plist    
          
    def getFormattedTaxonomyPath(self, node):
        """Return string formatted, separated by semicolon"""        
        p = self.getPath(node)
        return ";".join(clade.name for clade in p)
    
    def getSintaxFormattedTaxPath(self, node):
        """Return string formatted, separated by semicolon with nice letters before"""
        tp = self.getPath(node)
        taxString = "tax="
        for clade in tp:
            depth = self.getDepth(clade)
            if depth > CRESTree.META  and depth < CRESTree.SUBSPECIES and depth != CRESTree.SUPERKINGDOM:
                letter = CRESTree.depths[depth][0]
                taxString+=("%s:%s," % (letter, clade.name.replace(" ","_")))
        
        return taxString[:-1]

    def deleteNode(self, node, moveUpChildren=False):
        """Removes node from tree structure.

        @param moveUpChildren: if true, move up the children of the node to
                               its parent before deleting
        """
        node = self.verifyNode(node)
        if not node:            
            return
                
        parent = self.getParent(node)
        
        #self.tree.clearAllPaths()        
        if not node.is_terminal():
            children = self.getImmediateChildren(node)
            if moveUpChildren:
                parent.clades.extend(children)
            for child in children:
                if moveUpChildren:
                    self.parents[child] = parent
                else:
                    self.deleteNode(child, moveUpChildren=False)
        
        parent.clades.pop(parent.clades.index(node))
            
        if node.name in self.nodeNames:
            del self.nodeNames[node.name]
            
        if node in self.parents:
            del self.parents[node]

            
    def moveNode(self, node, newParent):
        """Moves node to new parent. Does not touch the nodeMapping"""
        node = self.verifyNode(node)
        if not node:
            return
        newParent = self.verifyNode(newParent)
        if not newParent:
            return
        oldParent = self.getParent(node)
        oldParent.clades.pop(oldParent.clades.index(node))
        newParent.clades.append(node)
        self.parents[node] = newParent
        #self.tree.clearAllPaths()
    
    def renameNode(self, node, newName):
        """Renames node by changing local name and key in nodeMapping"""
        node = self.verifyNode(node)
        if not node:
            return
        self.nodeNames[node.name] = newName
        node.name = newName
            
    def mergeNodes(self, nodeList, newName, newID=None):
        """Creates a new node with newName, then moves all child nodes of the
        lowest common ancestor identified in nodeList, deleting the nodes in the list
        """
        parent = self.common_ancestor(nodeList)        
        mergeNode = self.addNode(newName, parent)
        
        for node in nodeList:
            self.moveNode(node, mergeNode)
            self.deleteNode(node, moveUpChildren=True)    
    

    def getChildrenByRank(self, rank):
        """Recursively add taxa to rankChildren"""
        if rank==0:
            return [self.tree.root]
        children = self.getImmediateChildren(self.tree.root)
        for i in range(1,rank):
            grandchildren = []
            for child in children:
                grandkids = self.getImmediateChildren(child)
                for grandkid in grandkids:
                    grandchildren.append(grandkid)
            children = grandchildren
        return children
    
    def getCommonAncestor(self, nodeNames):
        """Accepts a list of accessions (corresponding to tre file) incl. duplicates and returns LCA"""
        nodes = []
        if not nodeNames:
            sys.stderr.write("ERROR: Empty list of node names for LCA received!\n")
            return None
        for nn in nodeNames:
            n = self.getNode(nn)
            if not n:
                sys.stderr.write("ERROR: Node %s not found in reference database\n" % nn)
                return None
            if not n in nodes:
                nodes.append(n)
                
        #return self.tree.common_ancestor(nodes)
        # Replacing the above with own implementation due to being extremely slow
        
        if len(nodes) == 1:
            return nodes[0]
        
        paths = []
        # Get paths for all nodes
        for n in nodes:
            p = self.getPath(n)
            paths.append(p)
            
        # Sort so that shortest path is first, start with its terminal node
        paths = sorted(paths,key=len)
        lcaPath = paths[0]
        others = paths[1:]
        
        # Go up in list until the node is in all of the paths 
        for o in others:
            while not lcaPath[-1] in o:
                if len(lcaPath) == 1:
                    return self.tree.root
                else:
                    lcaPath.pop()
        return lcaPath[-1]
    
    
    def _getAllChildren(self, node, children=None, parent = None):
        if children is None: children = []
        for c in node.clades:
            self._getAllChildren(c,children, parent = node)
        children.append(node)
        self.parents[node] = parent
        return children