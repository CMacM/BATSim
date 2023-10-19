"""Two sample benchmarks to compute runtime and memory usage.

For more information on writing benchmarks:
https://asv.readthedocs.io/en/stable/writing_benchmarks.html."""

import batsim
import galsim
import time

def compare_shear_speed():
    
    # create galaxy
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for batsim
    start = time.time()
    
    # generate stamp, shear and sample
    gal_stamp = batsim.stamp.Stamp(nn=64, scale=0.2, centering='galsim')
    lens = batsim.transforms.LensTransform(gamma1=0.2, center=[-0.5*scale, -0.5*scale])
    gal_stamp.transform_grids(lens)
    bat_array = gal_stamp.sample_galaxy(sersic_gal)
    
    end = time.time()
    
    bat_time = end-start
    
    # start timing for galsim
    start = time.time()
    
    gal_shear = gal.shear(g1=0.2)
    gal_array = gal_shear.drawImage(nx=64, ny=64, scale=0.2).array()
    
    end = time.time()
    
    gal_time = end-start
    
    return bat_time/gal_time

if __name__ == "__main__":
    test_affine_galsim()
    
    
