#!/usr/bin/env python
import os
import numpy as np
import pandas as pd


FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    directory = os.path.join(os.path.dirname(FILE_DIR), 'results')
    replicate_range = [1,2,3]
    names = ['Control Rep %d', 'Mutant Rep %d']
    filename = '100___BNI1_CLN3.tsv'

    for name in names:
        for rep in replicate_range:
            path = os.path.join(directory, name % rep, '100', filename)
            df = pd.read_csv(path, sep='\s+')
            for _, row in df.iterrows():
                fname = row['Filename']
                p = os.path.dirname(path).replace(directory, 'arrays')
                p = os.path.join(p, fname.replace('.tif', '_puncta.npy'))
                a = np.load(p, allow_pickle=True)
                f = os.path.join(os.path.dirname(path).replace(directory, 'results'), fname.replace('.tif', '_puncta.txt'))
                d = pd.read_csv(f, sep='\s+', comment='#')
                print(path)
                print(fname)
                if 'BNI1' in fname:
                    assert row['Num-puncta-1'] == len(a)
                    assert row['Num-puncta-1'] == len(d)
                else:
                    assert row['Num-puncta-2'] == len(a)
                    assert row['Num-puncta-2'] == len(d)
                print(d.tonumpy() == a)


if __name__ == '__main__':
    main()
