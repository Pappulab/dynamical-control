#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from tabulate import tabulate


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


pd.DataFrame.to_fwf = to_fwf


FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    directory = os.path.join(os.path.dirname(FILE_DIR), 'results')
    replicate_range = [1,2,3]
    names = ['Control Rep %d', 'Mutant Rep %d']
    filename = '100___BNI1_CLN3.tsv'
    puncta_radius_px = 150/108
    puncta_radius_nm = (150, 150)
    cluster_radius_nm = 300
    min_puncta_per_cluster = 2
    voxel_size_nm = (108, 108, 200)

    for name in names:
        for rep in replicate_range:
            path = os.path.join(directory, name % rep, '100', filename)
            print(path)
            df = pd.read_csv(path, sep='\s+')
            for _, row in df.iterrows():
                fname = row['Filename']
                p = os.path.dirname(path).replace(directory, 'arrays')
                p = os.path.join(p, fname.replace('.tif', '_puncta.npy'))
                a = np.load(p, allow_pickle=True)
                print(path)
                print(fname)
                if 'BNI1' in fname:
                    assert row['Num-puncta-1'] == len(a)
                else:
                    assert row['Num-puncta-2'] == len(a)

                header = [ f'# Filename: {fname}',
                           f'# Puncta radius (pixels): {puncta_radius_px}',
                           f'# Puncta radius (nm): {puncta_radius_nm}',
                           f'# Cluster radius (nm): {cluster_radius_nm}',
                           f'# Min puncta per cluster: {min_puncta_per_cluster}',
                           f'# Voxel size (nm): {voxel_size_nm}']

                fpath, ext = os.path.splitext(os.path.basename(fname))
                clustered_puncta_savename = os.path.join(os.path.dirname(path).replace(directory, 'results'), f'{fpath}_puncta.txt')
                d = dict()
                d['Y']          = a[:,0]
                d['X']          = a[:,1]
                df = pd.DataFrame(d)
                df.to_fwf(clustered_puncta_savename, header_lines=header)
                #print(clustered_puncta_savename, 111)


if __name__ == '__main__':
    main()
