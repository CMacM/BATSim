"""Benchmarks to check computation time and memory usage for batsim."""
import batsim.stamp as batstamp
import batsim.transforms as batforms
import galsim
import time

def time_shear_speed(nn=64, scale=0.2):
    
    # create galaxy
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for batsim
    start = time.time()
    
    # generate stamp, shear and sample
    gal_stamp = batstamp.Stamp(nn=nn, scale=scale, centering='galsim')
    lens = batforms.LensTransform(gamma1=0.2, gamma2=0, kappa=0,
                                  center=[-0.5*0.2, -0.5*0.2])
    gal_stamp.transform_grids(lens)
    bat_array = gal_stamp.sample_galaxy(gal)
    
    end = time.time()
    
    bat_time = end-start
    
    # start timing for galsim
    start = time.time()
    
    gal_shear = gal.shear(g1=0.2)
    gal_array = gal_shear.drawImage(nx=nn, ny=nn, scale=scale).array
    
    end = time.time()
    
    gal_time = end-start
    
    return {'batsim time': bat_time, 'galsim time' : gal_time}

def time_ia_speed(nn=128, scale=0.1):
    '''
        Times the time taken to apply an IA shear to a 128x128 pixel
        image and compares it with the time required to apply an
        affine shear.
    '''
    # initialise a galaxy object
    gal = galsim.Sersic(n=1.5, half_light_radius=1.5, flux=40)
    
    # start timing for ia transform
    ia_start = time.time()
    
    # crete stamp, apply transform, and sample gal obj
    ia_stamp = batstamp.Stamp(nn=nn, scale=scale)
    ia = batforms.IaTransform(scale=scale, hlr=1.5, phi=0.2)
    ia_stamp.transform_grids(ia)
    ia_gal = ia_stamp.sample_galaxy(gal)
    
    # stop timing
    ia_end = time.time()
    ia_time = ia_end - ia_start
    
    # get equivalent shear at hlr to pass to lens transform
    g1, g2 = ia.get_g1g2(1.5,0)
    
    # start timing for lens
    aff_start = time.time()
    
    # crete stamp, apply transform, and sample gal obj
    lens_stamp = batstamp.Stamp(nn=nn, scale=scale)
    lens = batforms.LensTransform(gamma1=g1, gamma2=g2, kappa=0)
    lens_stamp.transform_grids(lens)
    lens_gal = lens_stamp.sample_galaxy(gal)
    
    # stop timing
    aff_end = time.time()
    aff_time = aff_end - aff_start
    
    return {'IA time' : ia_time, 'Lens time' : aff_time}
    

if __name__ == "__main__":
    time_shear_speed()
    time_ia_speed()
    
    
