import unittest
import os

def test_suite():
    thisdir = os.path.dirname(__file__)
    return unittest.defaultTestLoader.discover(thisdir,top_level_dir=thisdir)
