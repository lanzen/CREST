import sys
import re
from types import StringType, IntType
from Bio import Phylo
from numpy import longfloat

import biom.table

#TODO - completely revise code to inherit from Phylo.Tree and also
# assign OTUs rather than reads, using a smarter structure

class CRESTree(Phylo.BaseTree.Tree):

    NORANK = 0
    DOMAIN = 1
    SUPERKINGDOM = 2
    KINGDOM = 3
    PHYLUM = 4
    CLASS = 5
    ORDER = 6
    FAMILY = 7
    GENUS = 8
    SPECIES = 9
    SUBSPECIES = 10

    depths = {NORANK: 'root',
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

    meganCodes = {NORANK: 0,
            DOMAIN: 0,
            PHYLUM: 2,
            CLASS: 3,
            ORDER: 4,
            FAMILY: 5,
            GENUS: 98,
            SPECIES: 100,
            SUBSPECIES: 101}

    """Reference tree initiated from newick and map file"""

    def __init__(self, trefile, mapfile, name=None):
        
        self = Phylo.read(trefile,"newick")
        self.name = name
        self.rooted = True        
        self.nodeNames = {} #dict of real name to each node
        self.assignmentMin = {} #dict of min. similarity in percent (int) for assignment
        
        self.noHits = Phylo.BaseTree.Clade(name="No hits")
        #self.addNode(self.noHits, self.root)
        
        # Read nodeNames from .map file (id\t name\t cutoff)
        theMap = open(mapfile, "r")
        for line in theMap:
            parts = line.split("\t")
            nodeID = parts[0]
            name = parts[1]
            minPercent = float(parts[2])
            targetNode = self.find_clades(nodeID).next()
            if targetNode:
                self.nodeNames[name] = targetNode
                self.assignmentMin[name] = minPercent
                targetNode.name = name
            else:
                sys.stderr.write("Error: Node %s not found in tree!\n" %
                                 nodeID)                
        theMap.close()
        
        
    def verifyNode(self, node):
        """Returns proper node if instead StringType is given as argument"""
        if isinstance(node, Phylo.BaseTree.Clade):
            return node
        elif type(node) is StringType:
            nn = node
            node = self.getNode(nn)
            if not node:
                sys.stderr.write("Node \"%s\" not found\n" % nn)
                return
            else:
                return node
        else: 
            sys.stderr.write("Node %s is a strange instance (should never happen)")
            return        
        

    def addNode(self, nodename, parent):
        #insert instead of append?
        node = Phylo.Newick.Clade(name=nodename)
        parent = self.verifyNode(parent)
        if not parent:
            return
        parent.clades.append(node)
        self.nodeNames[node.name] = node
        self.assignmentMin[node.name] = 0
        return node

    def getNode(self, name):
        if name in self.nodeNames:
            return self.nodes[name]
        else:
            sys.stderr.write("Error: Node %s not found!" % name)
            return None    
        
    def getImmediateChildren(self, node):
        node = self.verifyNode(node)
        return node.clades
    
    def getParent(self, node):
        node = self.verifyNode(node)
        if not node:
            return None
        path = self.get_path(node)
        if len(path) > 1:
            return path[2]
        else:
            return self.root            

    
    def deleteNode(self, node, moveUpChildren=False):
        """Removes node from tree structure. Still retains the disconnected
        node in self.nodes dict.

        @param moveUpChildren: if true, move up the children of the node to
                               its parent before deleting
        """
        node = self.verifyNode(node)
        if not node:                
            return

        if moveUpChildren or node.is_terminal():
            self.collapse(node)
            if node.name in self.nodeNames:
                del self.nodeNames[node.name]        
        else:
            for child in self.getImmediateChildren(node):
                self.deleteNode(child, moveUpChildren=False)
            
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
        newParent.append(node)        
            
            
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
            return [self.root, self.noHits]
        children = self.getImmediateChildren(self.root)
        for i in range(1,rank):
            grandchildren = []
            for child in children:
                grandchildren.append(self.getImmediateChildren(child))
            children = grandchildren
        return children    


