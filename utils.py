from configparser import SafeConfigParser
import os

def parse_config(filename):
    """ Parse config file, exit if it can't be found """
    if not os.path.exists(filename):
        sys.exit("Cannot find configuration file '{}'".format(filename))
    config = SafeConfigParser(allow_no_value=True)
    config.optionxform = str # Make sure keys are case sensitive
    config.read(filename)
    return config

