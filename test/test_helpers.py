import fitsio
import os
import numpy as np
import sys
import unittest
import shutil, tempfile
from pkg_resources import resource_filename


class AbstractTest(unittest.TestCase):
    """
        Class with Helper functions for the picca unit tests
    """



    def compare_fits(self, path1, path2, nameRun=""):
        """
            Compares all fits files in 2 directories against each other

        Args:
            path1 (str): path where first set of fits files lies
            path2 (str): path where second set of fits files lies
            nameRun (str, optional): A name of the current run for identification. Defaults to "".
        """
        m = fitsio.FITS(path1)
        self.assertTrue(os.path.isfile(path2), "{}".format(nameRun))
        b = fitsio.FITS(path2)

        self.assertEqual(len(m), len(b), "{}".format(nameRun))

        for i, _ in enumerate(m):

            ###
            r_m = m[i].read_header().records()
            ld_m = []
            for el in r_m:
                name = el['name']
                if len(name) > 5 and name[:5] == "TTYPE":
                    ld_m += [el['value'].replace(" ", "")]
            ###
            r_b = b[i].read_header().records()
            ld_b = []
            for el in r_b:
                name = el['name']
                if len(name) > 5 and name[:5] == "TTYPE":
                    ld_b += [el['value'].replace(" ", "")]

            self.assertListEqual(ld_m, ld_b, "{}".format(nameRun))

            for k in ld_m:
                d_m = m[i][k][:]
                d_b = b[i][k][:]
                if d_m.dtype in ['<U23',
                                'S23']:  # for fitsio old version compatibility
                    d_m = np.char.strip(d_m)
                if d_b.dtype in ['<U23',
                                'S23']:  # for fitsio old version compatibility
                    d_b = np.char.strip(d_b)
                self.assertEqual(d_m.size, d_b.size,
                                "{}: Header key is {}".format(nameRun, k))
                if not np.array_equal(d_m, d_b):
                    diff = d_m - d_b
                    diff_abs = np.absolute(diff)
                    w = d_m != 0.
                    diff[w] = np.absolute(diff[w] / d_m[w])
                    allclose = np.allclose(d_m, d_b)
                    self.assertTrue(
                        allclose,
                        "{}: Header key is {}, maximum relative difference is {}, maximum absolute difference is {}".
                        format(nameRun, k, diff.max(), diff_abs.max()))

        m.close()
        b.close()

        return


    @classmethod
    def load_requirements(cls, picca_base):
        """
            Loads reqirements file from picca_base
        """
        req = {}

        if sys.version_info > (3, 0):
            path = picca_base + '/requirements.txt'
        else:
            path = picca_base + '/requirements-python2.txt'
        with open(path, 'r') as f:
            for l in f:
                l = l.replace('\n', '').replace('==', ' ').replace('>=',
                                                                ' ').split()
                assert len(
                    l) == 2, "requirements.txt attribute is not valid: {}".format(
                        str(l))
                req[l[0]] = l[1]
        return req

    @classmethod
    def send_requirements(cls, req):
        """
            Compares requirements in req to currently loaded modules
        """
        for req_lib, req_ver in req.items():
            try:
                local_ver = __import__(req_lib).__version__
                if local_ver != req_ver:
                    print(
                        "WARNING: The local version of {}: {} is different from the required version: {}"
                        .format(req_lib, local_ver, req_ver))
            except ImportError:
                print("WARNING: Module {} can't be found".format(req_lib))

        return

    @classmethod
    def setUpClass(cls):
        """
            sets up directory structure in tmp
        """
        cls._branchFiles = tempfile.mkdtemp() + "/"
        cls.produce_folder(cls)
        cls.picca_base = resource_filename('picca',
                                           './').replace('py/picca/./', '')
        cls.send_requirements(cls.load_requirements(cls.picca_base))
        cls._masterFiles = cls.picca_base + '/py/picca/test/data/'
        cls._test=True


    @classmethod
    def tearDownClass(cls):
        """
            removes directory structure in tmp
        """
        if os.path.isdir(cls._branchFiles):
            shutil.rmtree(cls._branchFiles, ignore_errors=True)
