import sys

from taxa import Tree, Node, Read


class RDPRead(Read):

    def __init__(self, name, reverse):
        Read.__init__(self, name)
        self.reverse = reverse


class RDPTree(Tree):

    META = 0
    DOMAIN = 1
    PHYLUM = 2
    CLASS = 3
    ORDER = 4
    FAMILY = 5
    GENUS = 6
    SPECIES = 7
    SUB = 8

    depths = {  # META: 'meta',
              DOMAIN: 'domain',
              PHYLUM: 'phylum',
              CLASS: 'class',
              ORDER: 'order',
              FAMILY: 'family',
              GENUS: 'genus',
              SPECIES: 'species',
              SUB: 'subspecies'}

    def __init__(self, name=None):
        #print "Confidence level: %s" % confidence
        Tree.__init__(self, rootNode=Node(name="Root", depth=0), name=name)

    def parseRDPClassifier(self, outFile, confidence=0.5, weighFlows=True,
                           ignoreSubs=True):
        outFile = open(outFile, 'r')
        i = 0

        #Format:

        #IO95TX02.787F-MID-1_s60_c01_T400_s30_c08_2_3  Bacteria  domain   1.0 \
        #"Proteobacteria"phylum    0.77    Deltaproteobacteria    class   0.61\
        # Desulfobacterales  order 0.55    Desulfobulbaceae    family    0.44 \
        # Desulfopila    genus    0.43

        for line in outFile:
            if ";" in line or "\t" in line:
                line = line.replace(" - ", "")
                line = line.replace("\t-\t", "\t")
                line = line.replace("\t\t", "\t")
                sl = line.split("\t")
                reverse = None
                rName = sl[0]

                read = RDPRead(name=rName, reverse=reverse)
                node = parent = self.root
                levels = (len(sl) - 1) / 3
                #Problems at genus level that can have same name as phylum
                #if levels > 6: levels = 6
                for i in range(levels):
                    nodeName = sl[i * 3 + 1].replace("\"", "")
                    level = sl[i * 3 + 2]
                    conf = sl[i * 3 + 3]
                    if "%" in conf:
                        conf = conf[:-1] / 100
                    if ignoreSubs and level[0:3] == "sub":
                        pass
                    else:
                        p = parent
                        while p is not self.root:
                            if nodeName == p.name:
                                nodeName += " (%s)" % level
                                break
                            else:
                                p = p.parent

                        if conf >= confidence:

                            node = self.getNode(nodeName)
                            if not node:
                                node = Node(parent=parent, name=nodeName,
                                            depth=i)
                                self.addNode(node)
                            parent = node
                        else:
                            break
                node.assignRead(read, primary=True, recursive=True)

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "Use rdptree.py rdp_out_file confidence (0-1.0) ['tree'] "
    else:
        tree = False
        mt = RDPTree(name=sys.argv[1])
        for i in range(3, len(sys.argv)):
            if sys.argv[i] == 'tree':
                tree = True

        mt.parseRDPClassifier(sys.argv[1], confidence=(sys.argv[2]),
                              ignoreSubs=True)  # (not tree))
        mt.pruneUnassigned()

    if tree:
        dn = []
        for n in mt.root.children:
            dn.append(n)

        co = Node(parent=mt.root, name="Cellular organisms")

        for d in dn:
            mt.moveNode(d, co)

        mt.printAsTree(popDataset=None)
    else:
        for level in RDPTree.depths.keys():
            print("Assingments at %s level" % RDPTree.depths[level])
            print
            mt.printPopulationsAtDepth(level, normaliseToBase=False)
            print
