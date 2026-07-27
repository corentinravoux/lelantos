"""lelantos test suite."""

import unittest
import os

def test_suite():
    """Discover and return the lelantos unittest test suite."""
    thisdir = os.path.dirname(__file__)
    return unittest.defaultTestLoader.discover(thisdir,top_level_dir=thisdir)
