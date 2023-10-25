#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    directory = os.path.join(os.path.dirname(FILE_DIR),'results')
    result_file = '100___BNI1_CLN3.tsv'
    radius_range = range(100, 260, 10)
    all_results = dict()
    for radius in radius_range:
        p = os.path.join(directory, f'Test_{radius}', result_file)
        df = pd.read_csv(p, sep='\s+')
        all_results[radius] = df

    filenames = '001 002 007'.split()
    proteins = 'BNI1 CLN3'.split()

    protein = proteins[0]
    colors = 'red blue green purple'.split()
    for suffix in filenames:
        all_frac1 = list()
        all_frac2 = list()
        all_spots1 = list()
        all_spots2 = list()
        for radius in radius_range:
            df = all_results[radius]
            sel = df[df['Filename'].str.startswith(protein)]
            row = sel[sel['Filename'].str.endswith(f'{suffix}.tif')]
            frac1 = float(row['Colocalization-1'])
            frac2 = float(row['Colocalization-2'])
            spots1 = int(row['Num-puncta-1'])
            spots2 = int(row['Num-puncta-2'])
            
            all_frac1.append(frac1)
            all_frac2.append(frac2)
            all_spots1.append(spots1)
            all_spots2.append(spots2)

            print(suffix, protein, frac1, frac2, spots1, spots2)

        width = 0.4
        x = np.arange(len(radius_range))
        fig = plt.figure(figsize=(18, 8))
        ax = fig.add_subplot(111)
        rects1 = ax.bar(x - width/2, all_frac1, width, label=proteins[0], color='red')
        rects2 = ax.bar(x + width/2, all_frac2, width, label=proteins[1], color='blue')
        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        labels = ['%d' % r for r in radius_range]
        ax.set_xticks(x, labels)
        ax.set_title(suffix)
        ax.set_xlabel('Titration Radius (nm)')
        ax.set_ylabel('Fraction colocalized')
        ax.legend()
        fig.tight_layout()
        fig.savefig(f'frac_colocalized_{suffix}.png')
            
        fig = plt.figure(figsize=(18, 8))
        ax = fig.add_subplot(111)
        rects1 = ax.bar(x - width/2, all_spots1, width, label=proteins[0], color='green')
        rects2 = ax.bar(x + width/2, all_spots2, width, label=proteins[1], color='purple')
        ax.bar_label(rects1, padding=3)
        ax.bar_label(rects2, padding=3)
        labels = ['%d' % r for r in radius_range]
        ax.set_xticks(x, labels)
        ax.set_title(suffix)
        ax.set_xlabel('Titration Radius (nm)')
        ax.set_ylabel('Num Puncta')
        ax.legend()
        fig.tight_layout()
        fig.savefig(f'puncta_{suffix}.png')


if __name__ == '__main__':
    main()
