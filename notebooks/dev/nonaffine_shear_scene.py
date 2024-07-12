import galsim
import batsim
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from parallelbar import progress_starmap
from multiprocessing import cpu_count

# Function to shift the image and fill the edges with zeros
def shift_image(image, shift_x, shift_y):
    shifted_image = np.roll(image, shift=(int(shift_y), int(shift_x)), axis=(0, 1))

    if shift_x > 0:
        shifted_image[:, :int(shift_x)] = 0
    elif shift_x < 0:
        shifted_image[:, int(shift_x):] = 0

    if shift_y > 0:
        shifted_image[:int(shift_y), :] = 0
    elif shift_y < 0:
        shifted_image[int(shift_y):, :] = 0

    return shifted_image

# define function to be parallelized
def shear_and_shift(n_pix, pixel_scale, gal, psf_obj, shift_x, shift_y, gamma1, gamma2, kappa, F1, F2, G1, G2):

    # Flexion transform object to be passed to BATSim
    flexion_transform = batsim.FlexionTransform(
        gamma1=gamma1,
        gamma2=gamma2,
        kappa=kappa,
        F1=F1,
        F2=F2,
        G1=G1,
        G2=G2,
    )

    # simulate galaxy and get image with galaxy centered
    img = batsim.simulate_galaxy(
        ngrid=n_pix, # number of pixels in each dimension
        pix_scale=pixel_scale, # arcsec/pixel
        gal_obj=gal, # galsim galaxy object
        transform_obj=flexion_transform, # transform class telling BATSim what coordinate to request flux at
        psf_obj=psf_obj, # PSF to convolve with image after shearing
        draw_method="auto", # Include pixel response
        force_ngrid=False, # Force the number of pixels in each dimension for simulation if calculated size is smaller
        maximum_num_grids=8000 # Controls maximum possible sim size, 
                                # may lead to clipping of flux at edges of light profile for large galaxies if set too low 
    )
    
    # shift the image so the galaxy is not at the center
    img = shift_image(img, shift_x, shift_y)

    return img

# define parameters to be passed to function
n_pix = 512
pixel_scale = 0.1
psf_obj = galsim.Moffat(beta=3.5, fwhm=0.8, trunc=4*0.8)

# load in cosmos galaxy catalogue
cosmos = galsim.COSMOSCatalog()

# for an example we just randomise everything
rng = np.random.RandomState(1234)
inds = rng.choice(len(cosmos), 100)
gal_list = cosmos.makeGalaxy(index=inds, gal_type='parametric')

# position of galaxies
shift_x = rng.uniform(-n_pix*0.5, n_pix*0.5, len(gal_list))
shift_y = rng.uniform(-n_pix*0.5, n_pix*0.5, len(gal_list))

# shear parameters
gamma1 = rng.uniform(-0.01, 0.01, len(gal_list))
gamma2 = rng.uniform(-0.01, 0.01, len(gal_list))
kappa = rng.uniform(-0.01, 0.01, len(gal_list))
F1 = rng.uniform(-0.01, 0.01, len(gal_list))
F2 = rng.uniform(-0.01, 0.01, len(gal_list))
G1 = rng.uniform(-0.01, 0.01, len(gal_list))
G2 = rng.uniform(-0.01, 0.01, len(gal_list))

# Force galsim into single threading so we can batch process galaxies
original_omp_num_threads = os.environ.get("OMP_NUM_THREADS", None)
os.environ["OMP_NUM_THREADS"] = "1"

# construct list of arguments to be passed to paralellised function
args_list = [
    (
        n_pix, 
        pixel_scale, 
        gal, 
        psf_obj, 
        shift_x[i], 
        shift_y[i], 
        gamma1[i], 
        gamma2[i], 
        kappa[i], 
        F1[i], 
        F2[i], 
        G1[i], 
        G2[i]
        ) for i, gal in enumerate(gal_list)
]

# collect resulting images in a list
outcome = progress_starmap(shear_and_shift, args_list, n_cpu=cpu_count()//2) # use half the cores

# return OMP env variable to default state
if original_omp_num_threads is None:
    del os.environ['OMP_NUM_THREADS']
else:
    os.environ['OMP_NUM_THREADS'] = original_omp_num_threads

# sum all images to get final image
img = np.sum(outcome, axis=0)

# plot the final image
norm = matplotlib.colors.Normalize(vmin=0, vmax=0.2) # random norm values I find work well
plt.imshow(img, origin='lower', norm=norm)
plt.savefig("galaxy_scene_nonaffine.png")