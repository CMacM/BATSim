import fitsio
import numpy as np


class WCS(object):
    """
    Read SIP distortion coefficients from a FITS image header.

    The parsed coefficients are stored in ``ab`` with shape
    ``(2, order + 1, order + 1)``.  The first plane contains the A polynomial
    coefficients and the second plane contains the B polynomial coefficients.
    """

    def __init__(self, file_name=None):
        """
        Load SIP WCS coefficients from a FITS file.

        Parameters
        ----------
        file_name : str
            Path to a FITS file containing ``A_ORDER``/``B_ORDER`` and
            ``A_i_j``/``B_i_j`` header cards.
        """
        if file_name is None:
            print("Provide the filename")

        _, header = fitsio.read(file_name, header=True)

        self._set_ab(header)

        return

    def _set_ab(self, header):
        """
        Populate the SIP A/B coefficient array from a FITS header.

        Parameters
        ----------
        header : mapping
            FITS header object containing SIP polynomial metadata.
        """
        a_order = int(header.get("A_ORDER"))
        b_order = int(header.get("B_ORDER"))
        order = max(a_order, b_order)
        a = [float(header.get(f"A_{i}_{j}", 0)) for i in range(order + 1) for j in range(order + 1)]
        a = np.array(a).reshape((order + 1, order + 1))
        b = [float(header.get(f"B_{i}_{j}", 0)) for i in range(order + 1) for j in range(order + 1)]
        b = np.array(b).reshape((order + 1, order + 1))
        a[1, 0] += 1.0
        b[0, 1] += 1.0
        self.ab = np.array([a, b])

        return
