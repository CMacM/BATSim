import batsim
import galsim
import numpy as np

def test_get_flux_vec():
    sb_obj = galsim.Sersic(n=4, half_light_radius=0.5)
    trans_obj = sb_obj.shear(g1=0.1, g2=0.2).shift(0.5, 0.5)

    try:
        sb_flux = batsim._gsinterface.getFluxVec(
            scale=0.2, gsobj=sb_obj._sbp, xy_coords=np.array([[0.1, 0.1], [0.2, 0.2]])
        )
    except Exception as e:
        print("Error in getFluxVec for Sersic profile:", e)

    try:
        trans_flux = batsim._gsinterface.getFluxVec(
            scale=0.2, gsobj=trans_obj._sbp, xy_coords=np.array([[0.1, 0.1], [0.2, 0.2]])
        )
    except Exception as e:
        print("Error in getFluxVec for Transform profile:", e)

if __name__ == "__main__":
    test_get_flux_vec()

