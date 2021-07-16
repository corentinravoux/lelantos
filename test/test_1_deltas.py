import unittest
import os
import numpy as np
import fitsio
import healpy

from .test_helpers import AbstractTest


class TestDelta(AbstractTest):
    """
        Test case for picca_deltas.py
    """
    def produce_folder(self):
        """
            Create the necessary folders
        """

        lst_fold = [
            "/Products/", "/Products/Spectra/", "/Products/Spectra_MiniSV/",
            "/Products/Delta_Pk1D/", "/Products/Delta_Pk1D/Delta/",
            "/Products/Delta_Pk1D/Log/", "/Products/Delta_Pk1D_MiniSV/",
            "/Products/Delta_Pk1D_MiniSV/Delta/",
            "/Products/Delta_Pk1D_MiniSV/Log/", "/Products/Delta_LYA/",
            "/Products/Delta_LYA/Delta/", "/Products/Delta_LYA/Log/"
        ]

        for fold in lst_fold:
            if not os.path.isdir(self._branchFiles + fold):
                os.mkdir(self._branchFiles + fold)

        return

    def produce_cat(self, nObj, name="cat", thidoffset=0):
        """
            produces a fake catalog for testing
        """


        ### Create random catalog
        ra = 10. * np.random.random_sample(nObj)
        dec = 10. * np.random.random_sample(nObj)
        plate = np.random.randint(266, high=10001, size=nObj)
        mjd = np.random.randint(51608, high=57521, size=nObj)
        fiberid = np.random.randint(1, high=1001, size=nObj)
        thid = np.arange(thidoffset + 1, thidoffset + nObj + 1)
        z_qso = (3.6 - 2.0) * np.random.random_sample(nObj) + 2.0

        ### Save
        out = fitsio.FITS(self._branchFiles + "/Products/" + name + ".fits",
                          'rw',
                          clobber=True)
        cols = [ra, dec, thid, plate, mjd, fiberid, z_qso]
        names = ['RA', 'DEC', 'THING_ID', 'PLATE', 'MJD', 'FIBERID', 'Z']
        out.write(cols, names=names, extname='CAT')
        out.close()

        if self._test:
            path1 = self._masterFiles + "/test_delta/" + name + ".fits"
            path2 = self._branchFiles + "/Products/" + name + ".fits"
            self.compare_fits(path1, path2, "produce cat")

        return

    def produce_forests(self):
        """
            randomly creates Lya forests for testing
        """
        nside = 8

        ### Load DRQ
        vac = fitsio.FITS(self._branchFiles + "/Products/cat.fits")
        ra = vac[1]["RA"][:] * np.pi / 180.
        dec = vac[1]["DEC"][:] * np.pi / 180.
        thid = vac[1]["THING_ID"][:]
        plate = vac[1]["PLATE"][:]
        mjd = vac[1]["MJD"][:]
        fiberid = vac[1]["FIBERID"][:]
        vac.close()

        ### Get Healpy pixels
        pixs = healpy.ang2pix(nside, np.pi / 2. - dec, ra)

        ### Save master file
        path = self._branchFiles + "/Products/Spectra/master.fits"
        head = {}
        head['NSIDE'] = nside
        cols = [thid, pixs, plate, mjd, fiberid]
        names = ['THING_ID', 'PIX', 'PLATE', 'MJD', 'FIBER']
        out = fitsio.FITS(path, 'rw', clobber=True)
        out.write(cols, names=names, header=head, extname="MASTER TABLE")
        out.close()

        ### Log lambda grid
        logl_min = 3.550
        logl_max = 4.025
        logl_step = 1.e-4
        log_lambda = np.arange(logl_min, logl_max, logl_step)

        ### Loop over healpix
        for p in np.unique(pixs):

            ### Retrieve objects from catalog and produce fake spectra
            p_thid = thid[(pixs == p)]
            p_fl = np.random.normal(loc=1.,
                                    scale=1.,
                                    size=(log_lambda.size, p_thid.size))
            p_iv = np.random.lognormal(mean=0.1,
                                       sigma=0.1,
                                       size=(log_lambda.size, p_thid.size))
            p_am = np.zeros((log_lambda.size, p_thid.size)).astype(int)
            p_am[np.random.random_sample(size=(log_lambda.size,
                                               p_thid.size)) > 0.90] = 1
            p_om = np.zeros((log_lambda.size, p_thid.size)).astype(int)

            ### Save to file
            p_path = self._branchFiles + "/Products/Spectra/pix_" + str(
                p) + ".fits"
            out = fitsio.FITS(p_path, 'rw', clobber=True)
            out.write(p_thid, header={}, extname="THING_ID_MAP")
            out.write(log_lambda, header={}, extname="LOGLAM_MAP")
            out.write(p_fl, header={}, extname="FLUX")
            out.write(p_iv, header={}, extname="IVAR")
            out.write(p_am, header={}, extname="ANDMASK")
            out.write(p_om, header={}, extname="ORMASK")
            out.close()

        return


    def send_delta(self):
        """
            Test the continuum fitting routines on randomly generated eBOSS mock data
        """
        import picca.bin.picca_deltas as picca_deltas

        ### Send
        cmd = "picca_deltas.py"
        cmd += " --in-dir " + self._branchFiles + "/Products/Spectra/"
        cmd += " --drq " + self._branchFiles + "/Products/cat.fits"
        cmd += " --out-dir " + self._branchFiles + "/Products/Delta_LYA/Delta/"
        cmd += " --iter-out-prefix " + self._branchFiles + \
            "/Products/Delta_LYA/Log/delta_attributes"
        cmd += " --log " + self._branchFiles + "/Products/Delta_LYA/Log/input.log"
        cmd += " --nproc 1"
        picca_deltas.main(cmd.split()[1:])

        ### Test
        if self._test:
            path1 = self._masterFiles + "/test_delta/delta_attributes.fits.gz"
            path2 = self._branchFiles + "/Products/Delta_LYA/Log/delta_attributes.fits.gz"
            self.compare_fits(path1, path2, "picca_deltas.py")

        return


    def test_delta(self):
        """
            wrapper around send_delta test to produce the mock datasets
        """
        np.random.seed(42)
        self.produce_cat(nObj=1000)
        self.produce_forests()
        self.produce_cat(nObj=1000, name="random", thidoffset=1000)
        self.send_delta()
        return



if __name__ == '__main__':
    unittest.main()
