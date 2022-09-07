#!/usr/bin/env python
import argparse
from os import fspath, path

import numpy as np
from matplotlib.widgets import CheckButtons


def plot(
    nmap_filename,
    slc_stack_filename=None,
    slc_stack_bands=(1, 2, 3),
    slc_filename=None,
    ps_filename=None,
    alpha_nmap=0.8,
    cmap_nmap="Reds_r",
    block=False,
):
    """Interactively plot the neighborhood map over an SLC image.

    Click on a pixel to see the neighborhood map.
    """
    import matplotlib.pyplot as plt
    from matplotlib import patches

    slc = _load_slc(slc_stack_filename, slc_stack_bands, slc_filename)

    fig, ax = plt.subplots()
    ax.imshow(_scale_mag(slc), cmap="gray")

    # Also overlay the PS pixels, if passed
    ps = _load_ps(ps_filename)
    if ps is not None:
        ax2 = fig.add_axes((0.9, 0.9, 0.1, 0.1))
        ps_img = np.ma.masked_where(ps == 0, ps)
        axim_ps = ax.imshow(ps_img, cmap="winter_r", alpha=0.95, interpolation="nearest")

        labels = ["PS"]
        check = CheckButtons(ax2, labels, [True])

        def show_ps(label):
            # Flip the state of the checkbox
            # showing = not showing
            axim_ps.set_visible(not axim_ps.get_visible())
            fig.canvas.draw()

        check.on_clicked(show_ps)
        num_base_images = 2
    else:
        num_base_images = 1

    ny, nx = _get_windows(nmap_filename)

    def onclick(event):
        # Ignore right/middle click, clicks off image
        if event.button != 1 or not event.inaxes:
            return
        if event.inaxes != ax:
            return
        # Check if the toolbar has zoom or pan active
        # https://stackoverflow.com/a/20712813
        # MPL version 3.3: https://stackoverflow.com/a/63447351
        # if mpl.__version__ >= "3.3":
        state = fig.canvas.manager.toolbar.mode
        if state != "":  # Zoom/other tool is active
            return

        # Save limits to restore after adding neighborhoods
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        row, col = int(event.ydata), int(event.xdata)
        # Somehow clicked outside image, but in axis
        if row >= slc.shape[0] or col >= slc.shape[1]:
            return

        n = load_neighborhood(nmap_filename, row, col)
        n_img = np.ma.masked_where(n == 0, n)
        extent = _get_extent(row, col, ny, nx)

        # Remove old neighborhood images
        if len(event.inaxes.get_images()) > num_base_images:
            # event.inaxes.get_images()[1].remove()
            event.inaxes.get_images()[-1].remove()

        ax.imshow(
            n_img,
            cmap=cmap_nmap,
            alpha=alpha_nmap,
            extent=extent,
            origin="lower",
        )
        # Remove old neighborhood patches
        for p in ax.patches:
            p.remove()
        # add a rectangle around the neighborhood
        rect = patches.Rectangle(
            (extent[0], extent[2]),
            1 + 2 * nx,
            1 + 2 * ny,
            linewidth=1,
            edgecolor="r",
            facecolor="none",
        )
        ax.add_patch(rect)

        # Restore original viewing bounds
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        fig.canvas.draw()

    fig.canvas.mpl_connect("button_press_event", onclick)


    plt.show(block=block)


def load_neighborhood(filename, row, col):
    """Get the neighborhood of a pixel as a numpy array

    Parameters
    ----------
    row : int
        Row of the pixel
    col: int
        column of pixel

    Returns
    -------
    neighborhood : numpy array (dtype = np.bool)
    """
    from osgeo import gdal

    ny, nx = _get_windows(filename)
    # nbands = ds.RasterCount
    # Shape: (nbands, 1, 1)
    ds = gdal.Open(fspath(filename), gdal.GA_ReadOnly)
    pixel = ds.ReadAsArray(xoff=col, yoff=row, xsize=1, ysize=1)
    ds = None

    pixel = pixel.ravel()
    assert pixel.view("uint8").shape[0] == 4 * len(pixel)
    # 1D version of neighborhood, padded with zeros
    neighborhood = np.unpackbits(pixel.view("uint8"), bitorder="little")
    wx = 2 * nx + 1
    wy = 2 * ny + 1
    ntotal = wx * wy
    # assert np.all(neighborhood[ntotal:] == 0)
    return neighborhood[:ntotal].reshape((wy, wx))


def _get_windows(filename):
    from osgeo import gdal

    ds = gdal.Open(fspath(filename), gdal.GA_ReadOnly)
    if "ENVI" in ds.GetMetadataDomainList():
        meta = ds.GetMetadata("ENVI")
    else:
        meta = ds.GetMetadata()
    nx = int(meta["HALFWINDOWX"])
    ny = int(meta["HALFWINDOWY"])
    ds = None
    return ny, nx


def _get_extent(row, col, ny, nx):
    """Get the row/col extent of the window surrounding a pixel."""
    # Matplotlib extent is (left, right, bottom, top)
    # Also the extent for normal `imshow` is shifted by -0.5
    return col - nx - 0.5, col + nx + 1 - 0.5, row - ny - 0.5, row + ny + 1 - 0.5


def _load_ps(ps_filename):
    from osgeo import gdal

    if ps_filename is not None and path.exists(ps_filename):
        ds = gdal.Open(fspath(ps_filename), gdal.GA_ReadOnly)
        ps = ds.ReadAsArray()
        ds = None
    else:
        ps = None
    return ps


def _load_slc(slc_stack_filename, stack_bands, slc_filename):
    from osgeo import gdal

    if slc_stack_filename is not None:
        ds = gdal.Open(fspath(slc_stack_filename), gdal.GA_ReadOnly)
        # Average the power of the complex bands
        slc = np.zeros((ds.RasterYSize, ds.RasterXSize), dtype=np.float32)
        for b in stack_bands:
            slc += np.abs(ds.GetRasterBand(b).ReadAsArray()) ** 2
        slc = np.sqrt(slc / len(stack_bands))
        ds = None
    else:
        # If one SLC is provided, use that
        ds = gdal.Open(fspath(slc_filename), gdal.GA_ReadOnly)
        slc = np.abs(ds.ReadAsArray())
        ds = None
    return slc


def _scale_mag(img, exponent=0.3, max_pct=99.95):
    """Scale the magnitude of complex radar image for better display"""
    out = np.abs(img) ** exponent
    max_val = np.nanpercentile(out, max_pct)
    return np.clip(out, None, max_val)


def get_cli_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(
        description="Interactively view neighborhood maps over an SLC."
    )
    parser.add_argument(
        "-n",
        "--nmap-file",
        help="Neighborhood map file",
        required=True,
    )
    parser.add_argument(
        "--slc-stack-file",
        help="SLC stack file",
    )
    parser.add_argument(
        "--slc-stack-bands",
        nargs="+",
        type=int,
        default=[1, 2, 3, 4, 5],
        help="Bands to use from SLC stack file",
    )
    parser.add_argument(
        "--slc-file",
        help="Alternative background: a single SLC filename.",
    )
    parser.add_argument(
        "--ps-file",
        default="PS/ps_pixels",
        help="PS file to overlay",
    )
    parser.add_argument(
        "--cmap-nmap",
        default="Reds_r",
        help="Colormap for neighborhood map",
    )
    parser.add_argument(
        "--alpha-nmap",
        type=float,
        default=0.8,
        help="Alpha value for neighborhood map",
    )
    parser.add_argument(
        "--no-block",
        action="store_true",
        help="Don't block the matplotlib window (default is blocking)",
    )
    return parser.parse_args()


def run_cli():
    args = get_cli_args()
    plot(
        nmap_filename=args.nmap_file,
        slc_stack_filename=args.slc_stack_file,
        slc_stack_bands=args.slc_stack_bands,
        slc_filename=args.slc_file,
        ps_filename=args.ps_file,
        cmap_nmap=args.cmap_nmap,
        alpha_nmap=args.alpha_nmap,
        block=not args.no_block,
    )


if __name__ == "__main__":
    run_cli()
