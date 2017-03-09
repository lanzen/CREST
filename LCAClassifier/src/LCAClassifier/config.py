import ConfigParser
import os.path
import sys


def buildout_path():
    return os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))


class LCAClassifierConfig(object):

    def __init__(self):
        self.DATABASES = None
        #run configure
        self.configure()

    def configure(self, config=None):

        if config is None:
            # Assuming script being run from buildout path
            config_path = os.path.join('parts', 'etc', 'lcaclassifier.conf')
            bp = buildout_path()
            if "/src" in bp:
                bp = bp.replace("/src","")
            config = os.path.join(bp, config_path)            

        if os.path.isfile(config):
            parser = ConfigParser.ConfigParser()
            parser.read(config)
            self.DATABASES = dict(parser.items('Databases'))
        else:
            raise Exception("Unable to find %s" % config)

config = LCAClassifierConfig()
