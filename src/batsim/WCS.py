import fitsio
import numpy as np


class WCS(object):
    def __init__(self, file_name=None):
        if file_name is None:
            print("Provide the filename")

        _, header = fitsio.read(file_name, header=True)

        self._set_ab(header)

        return

    def _set_ab(self, header):
        a_order = int(header.get("A_ORDER"))
        b_order = int(header.get("B_ORDER"))
        order = max(a_order, b_order)
        a = [
            float(header.get(f"A_{i}_{j}", 0))
            for i in range(order + 1)
            for j in range(order + 1)
        ]
        a = np.array(a).reshape((order + 1, order + 1))
        b = [
            float(header.get(f"B_{i}_{j}", 0))
            for i in range(order + 1)
            for j in range(order + 1)
        ]
        b = np.array(b).reshape((order + 1, order + 1))
        a[1, 0] += 1.0
        b[0, 1] += 1.0
        self.ab = np.array([a, b])

        return
