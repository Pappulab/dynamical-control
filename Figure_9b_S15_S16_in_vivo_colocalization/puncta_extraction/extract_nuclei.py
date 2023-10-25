#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import bigfish.stack as stack
from cellpose import models
from cellpose import io
from cellpose.io import imread
from pprint import pprint
from argparse import ArgumentParser
from skimage import measure
from scipy import ndimage
from tabulate import tabulate
from skimage.filters import try_all_threshold
from skimage.measure import regionprops, regionprops_table
from skimage.filters import threshold_otsu, threshold_local, try_all_threshold
from skimage.restoration import (denoise_tv_chambolle, denoise_bilateral,
                                 denoise_wavelet, estimate_sigma)
from cellpose import plot as cellplot



DEFAULT_IMAGE_EXTENSION = '.tif'

# Supremely useful answer to convert a DataFrame to a TSV file: 
# https://stackoverflow.com/a/35974742/866930
def to_fwf(df, fname, header_lines=None):
    content = tabulate(df.values.tolist(), 
                       list(df.columns), 
                       tablefmt="plain", 
                       floatfmt='.5f',
                       numalign='left')
    with open(fname, 'w') as ofile:
        if header_lines is not None:
            ofile.write('\n'.join(header_lines))
            ofile.write('\n')
        ofile.write(content)

# Patch `pandas` to use our custom method
pd.DataFrame.to_fwf = to_fwf


def load_images(directory, extension):
    d = directory
    e = extension
    files = [os.path.join(d, f) for f in os.listdir(directory) if f.endswith(e)]
    return list(sorted(files))


def main():
    parser = ArgumentParser()
    parser.add_argument('-d',
                        '--directory',
                        help='The directory containing the nuclei images (TIF).',
                        type=str)
    parser.add_argument('-e',
                        '--extension',
                        help='The extension of the nuclei images.',
                        type=str,
                        default=DEFAULT_IMAGE_EXTENSION)
    args = parser.parse_args()

    if args.directory is None or len(args.directory) == 0:
        raise RuntimeError('No directory provided.')

    tif_files = load_images(args.directory, args.extension)
    if len(tif_files) == 0:
        raise RuntimeError('Please provide a directory that contains TIF files.')

    arrays_path = os.path.join(os.getcwd(), 'arrays_nuclei')
    plots_path = os.path.join(os.getcwd(), 'plots_nuclei')
    results_path = os.path.join(os.getcwd(), 'results_nuclei')
    if not os.path.exists(arrays_path):
        os.makedirs(arrays_path)

    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

    if not os.path.exists(results_path):
        os.makedirs(results_path)

    model = models.Cellpose(gpu=False, model_type='cyto')

    block_size = 3
    channels = [[0,0]]
    ndim = 2
    voxel_size_nm = (108, 108, 200)
    puncta_radius_nm = (500, 500)
    puncta_radius_px = (puncta_radius_nm[0]/voxel_size_nm[0], puncta_radius_nm[1]/voxel_size_nm[1])

    for index, tfile in enumerate(tif_files, start=1):
        rna = stack.read_image(tfile)

        fig = plt.figure(figsize=(12, 5))
        masks, flows, styles, diams = model.eval(rna, diameter=10, channels=channels)

        # determine the locations of the centroids
        regions = regionprops(masks)

        ofilename = os.path.basename(tfile).replace('.tif', '_cellpose.txt')
        header = [ f'# Filename: {ofilename}',
                   f'# Voxel size (nm): {voxel_size_nm[:ndim]}',
                   '# Area quantity: voxels',
                   f'# Num nuclei detected: {len(regions)}']
        cellpose_savename = os.path.join(results_path, ofilename)
        print(cellpose_savename)
        
        d = dict()
        d['Y'] = list()
        d['X'] = list()
        d['AREA'] = list()
        d['MAJOR_AXIS_LENGTH'] = list()
        d['MINOR_AXIS_LENGTH'] = list()

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111)
        for region in regions:
            y, x = region.centroid
            d['Y'].append(y)
            d['X'].append(x)
            d['AREA'].append(region.area)
            d['MAJOR_AXIS_LENGTH'].append(region.axis_major_length)
            d['MINOR_AXIS_LENGTH'].append(region.axis_minor_length)
            ax.plot(x, 2048 - y, marker='o', markersize=4, color='red', alpha=0.5)
        df = pd.DataFrame(d)
        df.to_fwf(cellpose_savename, header_lines=header)

        fname, ext = os.path.splitext(os.path.basename(tfile))
        savename = os.path.join(arrays_path, f'{fname}.npz')
        np.savez(savename, masks=masks, flows=flows, styles=styles, diams=diams)
        cellplot.show_segmentation(fig, rna, masks, flows[0], channels=channels)
        fig.tight_layout()
        fig.savefig(os.path.join(plots_path, f'{fname}.png'))



if __name__ == '__main__':
    main()
