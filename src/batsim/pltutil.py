import numpy as np
import galsim
import matplotlib.pyplot as plt
from astropy.visualization import simple_norm

colors = [
    "#000000",
    "#1976D2",
    "#E53935",
    "#43A047",
    "#673AB7",
    "#4DD0E1",
    "#E91E63",
    "#F2D026",
    "#333333",
    "#9E9E9E",
    "#FB8C00",
    "#FFB300",
    "#795548",
]

cblues = ["#004c6d", "#346888", "#5886a5", "#7aa6c2", "#9dc6e0", "#c1e7ff"]
creds = ["#DC1C13", "#EA4C46", "#F07470", "#F1959B", "#F6BDC0", "#F8D8E3"]


def make_figure_axes(ny=1, nx=1, square=True):
    """Makes figure and axes

    Args:
        ny (int):       number of subplots in y direction
        nx (int):       number of subplots in y direction
        square (bool):  whether using square plot
    """
    if not isinstance(ny, int):
        raise TypeError("ny should be integer")
    if not isinstance(nx, int):
        raise TypeError("nx should be integer")
    axes = []
    if ny == 1 and nx == 1:
        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(ny, nx, 1)
        axes.append(ax)
    elif ny == 2 and nx == 1:
        if square:
            fig = plt.figure(figsize=(6, 11))
        else:
            fig = plt.figure(figsize=(6, 7))
        ax = fig.add_subplot(ny, nx, 1)
        axes.append(ax)
        ax = fig.add_subplot(ny, nx, 2)
        axes.append(ax)
    elif ny == 1 and nx == 2:
        fig = plt.figure(figsize=(11, 6))
        for i in range(1, 3):
            ax = fig.add_subplot(ny, nx, i)
            axes.append(ax)
    elif ny == 1 and nx == 3:
        fig = plt.figure(figsize=(18, 6))
        for i in range(1, 4):
            ax = fig.add_subplot(ny, nx, i)
            axes.append(ax)
    elif ny == 1 and nx == 4:
        fig = plt.figure(figsize=(20, 5))
        for i in range(1, 5):
            ax = fig.add_subplot(ny, nx, i)
            axes.append(ax)
    elif ny == 2 and nx == 3:
        fig = plt.figure(figsize=(15, 8))
        for i in range(1, 7):
            ax = fig.add_subplot(ny, nx, i)
            axes.append(ax)
    elif ny == 2 and nx == 4:
        fig = plt.figure(figsize=(20, 8))
        for i in range(1, 9):
            ax = fig.add_subplot(ny, nx, i)
            axes.append(ax)
    else:
        raise ValueError("Do not have option: ny=%s, nx=%s" % (ny, nx))
    return fig, axes

def determine_cuts(data):
    """
    Determine min_cut and max_cut for the data using median and standard deviation.

    Parameters:
        data (ndarray): 2D numpy array containing the image data
        sigma (int): Number of standard deviations to use for max_cut

    Returns:
        min_cut, max_cut: Calculated cuts
    """
    min_cut = np.percentile(np.ravel(data), 5)
    max_cut = np.percentile(np.ravel(data), 98)
    return min_cut, max_cut


def make_plot_image(data):
    min_cut, max_cut = determine_cuts(data)
    sn = simple_norm(data, "asinh", asinh_a=0.1, min_cut=min_cut, max_cut=max_cut)
    fig = plt.imshow(data, aspect="equal", cmap="RdYlBu_r", origin="lower", norm=sn)
    return fig


def stitch_images(images, direction='horizontal', spacing=None):
    """
    Stitch multiple images together to create a single composite image.

    Args:
        images (list): A list of images to be stitched together.
        direction (str, optional): The direction of stitching. Can be 'horizontal' or 'vertical'. Defaults to 'horizontal'.
        spacing (int, optional): The spacing between images. If None, images will be stitched with no gap between them. Defaults to None.

    Returns:
        galsim.ImageF: The stitched composite image.

    Raises:
        None

    """
    # read in sizes of the individual images
    # MUST BE SAME FOR ALL IMAGES RIGHT NOW
    nx = images[0].xmax
    ny = images[0].ymax
    scale = images[0].scale

    # Check is images should be stitched with no gap between them
    if spacing is None:
        # Check for direction of stitching
        if direction=='horizontal':
            # Determine required size of final 'super' image
            # And create an empty canvas to draw to
            Nx = len(images) * nx
            super_image = galsim.ImageF(Nx, ny, scale=scale)
            # Counter to ensure correct bounds are applied
            i = 0
            for image in images:
                # Determine the bounds within which image should be
                # placed in super_image and then place
                bounds = galsim.BoundsI(xmin = 1 + (i*nx), 
                                        xmax = nx + (i*nx),
                                        ymin = 1,
                                        ymax = ny)
                                       
                super_image.setSubImage(bounds, image)
                
                i = i + 1 # update for next iteration
            return super_image
        
        # Same as above but for vertical stitching
        elif direction=='vertical': 
            # Determine required size of final 'super' image
            # And create an empty canvas to draw to
            Ny = len(images) * ny
            super_image = galsim.ImageF(nx, Ny, scale=scale)
            # Counter to ensure correct bounds are applied
            i = 0
            for image in images:
                # Determine the bounds within which image should be
                # placed in super_image and then place
                bounds = galsim.BoundsI(xmin = 1, 
                                        xmax = nx,
                                        ymin = 1 + (i*ny),
                                        ymax = ny + (i*ny))
                                       
                super_image.setSubImage(bounds, image)
            
                i = i + 1 # update for next iteration
            return super_image    
    # TODO: Allow for empty space to be inserted between images
    
def split_image(image, nsplit, direction='horizontal', spacing=None):
    ''' Utility function which can be used to split galsim images
    into smaller individual stamps.'''
    # read in sizes of the image
    Nx = image.xmax
    Ny = image.ymax
    scale = image.scale
    
    # Check if images have gaps between them
    if spacing is None:
        # Check for direction of spliting
        if direction=='horizontal':
            # Determine required size of each split image
            nx = Nx / nsplit
            ny = Ny
            # Counter to ensure correct bounds are applied
            i = 0
            # list to contain split images
            split_images = []
            for i in range(nsplit):
                split_image = galsim.ImageF(nx, ny, scale=scale)
                # Determine bounds within which to get the sub image
                bounds = galsim.BoundsI(xmin = 1 + (i*nx), 
                                        xmax = nx + (i*nx),
                                        ymin = 1,
                                        ymax = ny)  
                
                sub_image = image.subImage(bounds)
                split_image.copyFrom(sub_image)
                # Determine the bounds within which image should be
                split_images.append(split_image)
                
                i = i + 1 # update for next iteration 
            return split_images
        
        # Same as above but for vertical stitching
        elif direction=='vertical':
            # Determine required size of each split image
            ny = Ny / nsplit
            nx = Nx
            # Counter to ensure correct bounds are applied
            i = 0
            # list to contain split images
            split_images = []
            for i in range(nsplit):
                split_image = galsim.ImageF(nx, ny, scale=scale)
                # Determine bounds within which to get the sub image
                bounds = galsim.BoundsI(xmin = 1, 
                                        xmax = nx,
                                        ymin = 1 + (i*ny),
                                        ymax = ny + (i*ny))   
                
                sub_image = image.subImage(bounds)
                split_image.copyFrom(sub_image)
                # Determine the bounds within which image should be
                split_images.append(split_image)
                
                i = i + 1 # update for next iteration
            return split_images
    # TODO: Allow for empty space to be inserted between images