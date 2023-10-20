"""Benchmarks to check computation time and memory usage for batsim."""
import batsim.stamp as batstamp
import batsim.transforms as batforms
import galsim
import time

def time_shear_speed():
    
    # create galaxy
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for batsim
    start = time.time()
    
    # generate stamp, shear and sample
    gal_stamp = batstamp.Stamp(nn=640, scale=0.2, centering='galsim')
    lens = batforms.LensTransform(gamma1=0.2, gamma2=0, kappa=0,
                                  center=[-0.5*0.2, -0.5*0.2])
    gal_stamp.transform_grids(lens)
    bat_array = gal_stamp.sample_galaxy(gal)
    
    end = time.time()
    
    bat_time = end-start
    
    # start timing for galsim
    start = time.time()
    
    gal_shear = gal.shear(g1=0.2)
    gal_array = gal_shear.drawImage(nx=640, ny=640, scale=0.2).array
    
    end = time.time()
    
    gal_time = end-start
    
    return {'batsim time': bat_time, 'galsim time' : gal_time}

if __name__ == "__main__":
    test_affine_galsim()
    
    
