import sys

from LCAClassifier.rdpclassifierParse import RDPTree, RDPRead
from LCAClassifier.taxa import Node


class RDPQIIMETree(RDPTree):

    def parseQIIME(self, outFile, confidence=0.5, weighFlows=True):

        outFile = open(outFile, 'r')
        i = 0

        #Format:
        #IO95TX02.787F-MID-1_s60_c01_T400_s30_c08_2_3
        #Bacteria;Proteobacteria;Deltaproteobacteria;Desulfobulbaceae    0.82

        for line in outFile:
            #Add read
            #print line
            sl = line.split("\t")
            rName = sl[0]
            taxa = sl[1].split(";")
            conf = sl[2]

            read = RDPRead(name=rName, reverse=None)
            node = parent = self.root

            levels = len(taxa)
            for i in range(levels):
                nodeName = taxa[i].replace("\"", "")
                p = parent
                while p is not self.root:
                    if nodeName == p.name:
                        nodeName += " (%s)" % p.name
                        break
                    else:
                        p = p.parent

                if conf >= confidence:
                    node = self.getNode(nodeName)
                    if not node:
                        node = Node(parent=parent, name=nodeName, depth=i)
                        self.addNode(node)
                    parent = node
            node.assignRead(read, primary=True, recursive=True)


if __name__ == "__main__":
        if len(sys.argv) < 3:
            print ("Use qiimeAssignmentParse.py rdp_out_file confidence "
                   "(0-1.0) ['tree'] ")
        else:
            tree = False
            mt = RDPQIIMETree(name=sys.argv[1])
            for i in range(3, len(sys.argv)):
                if sys.argv[i] == 'tree':
                    tree = True
            mt.parseQIIME(sys.argv[1], confidence=(sys.argv[2]))  # (not tree))
            mt.pruneUnassigned()

        if tree:
            dn = []
            for n in mt.root.children:
                dn.append(n)

            co = Node(parent=mt.root, name="Cellular organisms", population=0,
                                        reads=[], singles=0, doubles=0)

            for d in dn:
                mt.moveNode(d, co)

            mt.printAsTree()
        else:
            for level in RDPTree.depths.keys():
                print("Assingments at %s level" % RDPTree.depths[level])
                print
                mt.printPopulationsAtDepth(level, normaliseToBase=False)
                print
