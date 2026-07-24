"""Legacy plotting utilities used by BATSim examples and validation notebooks."""

import galsim
import matplotlib.pyplot as plt
import numpy as np
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
    """
    Create a Matplotlib figure with one of BATSim's preset subplot layouts.

    Parameters
    ----------
    ny : int, optional
        Number of subplot rows.
    nx : int, optional
        Number of subplot columns.
    square : bool, optional
        Select the taller preset size for the supported ``2 x 1`` layout.

    Returns
    -------
    tuple
        ``(fig, axes)`` where ``fig`` is a Matplotlib figure and ``axes`` is a
        flat list of axes in creation order.

    Raises
    ------
    TypeError
        Raised if ``ny`` or ``nx`` is not an integer.
    ValueError
        Raised if the requested layout is not one of the supported presets.
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
    Determine display cuts for image data using fixed percentiles.

    Parameters
    ----------
    data : ndarray
        Image array used to calculate display cuts.

    Returns
    -------
    tuple of float
        ``(min_cut, max_cut)`` from the 5th and 98th percentiles of the
        flattened input data.
    """
    min_cut = np.percentile(np.ravel(data), 5)
    max_cut = np.percentile(np.ravel(data), 98)
    return min_cut, max_cut


def make_plot_image(data):
    """
    Display image data with BATSim's default asinh plotting normalisation.

    Parameters
    ----------
    data : ndarray
        Two-dimensional image array to display.

    Returns
    -------
    matplotlib.image.AxesImage
        Image artist returned by ``matplotlib.pyplot.imshow``.
    """
    min_cut, max_cut = determine_cuts(data)
    sn = simple_norm(data, "asinh", asinh_a=0.1, vmin=min_cut, vmax=max_cut)
    fig = plt.imshow(data, aspect="equal", cmap="RdYlBu_r", origin="lower", norm=sn)
    return fig


def stitch_images(images, direction="horizontal", spacing=None):
    """
    Stitch equal-sized Galsim image objects into a single composite image.

    Parameters
    ----------
    images : sequence of galsim.Image
        Images to stitch together. All images are expected to share the same
        dimensions and pixel scale.
    direction : {"horizontal", "vertical"}, optional
        Direction in which images are concatenated.
    spacing : None, optional
        Placeholder for future gap support. Only ``None`` is currently
        implemented.

    Returns
    -------
    galsim.ImageF or None
        Composite image when ``spacing`` is ``None`` and ``direction`` is
        supported; otherwise ``None``.
    """
    # read in sizes of the individual images
    # MUST BE SAME FOR ALL IMAGES RIGHT NOW
    nx = images[0].xmax
    ny = images[0].ymax
    scale = images[0].scale

    # Check is images should be stitched with no gap between them
    if spacing is None:
        # Check for direction of stitching
        if direction == "horizontal":
            # Determine required size of final 'super' image
            # And create an empty canvas to draw to
            Nx = len(images) * nx
            super_image = galsim.ImageF(Nx, ny, scale=scale)
            # Counter to ensure correct bounds are applied
            i = 0
            for image in images:
                # Determine the bounds within which image should be
                # placed in super_image and then place
                bounds = galsim.BoundsI(xmin=1 + (i * nx), xmax=nx + (i * nx), ymin=1, ymax=ny)

                super_image.setSubImage(bounds, image)

                i = i + 1  # update for next iteration
            return super_image

        # Same as above but for vertical stitching
        elif direction == "vertical":
            # Determine required size of final 'super' image
            # And create an empty canvas to draw to
            Ny = len(images) * ny
            super_image = galsim.ImageF(nx, Ny, scale=scale)
            # Counter to ensure correct bounds are applied
            i = 0
            for image in images:
                # Determine the bounds within which image should be
                # placed in super_image and then place
                bounds = galsim.BoundsI(xmin=1, xmax=nx, ymin=1 + (i * ny), ymax=ny + (i * ny))

                super_image.setSubImage(bounds, image)

                i = i + 1  # update for next iteration
            return super_image


def split_image(image, nsplit, direction="horizontal", spacing=None):
    """
    Split a stitched GalSim image into equal-sized image stamps.

    Parameters
    ----------
    image : galsim.Image
        Input image to split.
    nsplit : int
        Number of equal pieces to extract along the split direction.
    direction : {"horizontal", "vertical"}, optional
        Axis along which the input image is split.
    spacing : None, optional
        Placeholder for future gap support. Only ``None`` is currently
        implemented.

    Returns
    -------
    list of galsim.ImageF or None
        List of split images when ``spacing`` is ``None`` and ``direction`` is
        supported; otherwise ``None``.
    """
    # read in sizes of the image
    Nx = image.xmax
    Ny = image.ymax
    scale = image.scale

    # Check if images have gaps between them
    if spacing is None:
        # Check for direction of spliting
        if direction == "horizontal":
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
                bounds = galsim.BoundsI(xmin=1 + (i * nx), xmax=nx + (i * nx), ymin=1, ymax=ny)

                sub_image = image.subImage(bounds)
                split_image.copyFrom(sub_image)
                # Determine the bounds within which image should be
                split_images.append(split_image)

                i = i + 1  # update for next iteration
            return split_images

        # Same as above but for vertical stitching
        elif direction == "vertical":
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
                bounds = galsim.BoundsI(xmin=1, xmax=nx, ymin=1 + (i * ny), ymax=ny + (i * ny))

                sub_image = image.subImage(bounds)
                split_image.copyFrom(sub_image)
                # Determine the bounds within which image should be
                split_images.append(split_image)

                i = i + 1  # update for next iteration
            return split_images
