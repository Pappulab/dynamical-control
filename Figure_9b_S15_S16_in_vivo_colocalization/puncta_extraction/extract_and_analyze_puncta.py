#!/usr/bin/env python
import matplotlib
matplotlib.use('Agg')

import os
import sys
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

import bigfish
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.multistack as multistack
import bigfish.plot as plot

from pprint import pprint
from tabulate import tabulate
from collections import namedtuple
from configparser import ConfigParser
from argparse import ArgumentParser


# =====================================================================================================================
# CONSTANTS


PROJECT_ROOT                        = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PLOTS_DIR                   = os.path.join(PROJECT_ROOT, 'plots')
DEFAULT_RESULTS_DIR                 = os.path.join(PROJECT_ROOT, 'results')
DEFAULT_ARRAYS_DIR                  = os.path.join(PROJECT_ROOT, 'arrays')
DEFAULT_CONFIG_FILENAME             = 'config.ini'
DEFAULT_CROP_PERCENTAGE             = 100
DEFAULT_PROTEINS                    = 'BNI1 CLN3'.split()
DEFAULT_PROTEIN_TYPES               = 'bidi ctrl'.split()
DEFAULT_VOXEL_SIZE                  = (108, 108, 200)
DEFAULT_IMAGE_FILETYPE              = '.tif'
DEFAULT_MIN_PUNCTA                  = 2
DEFAULT_PLOT_TYPE                   = 'png'
DEFAULT_PUNCTA_RADIUS               = (150, 150)
DEFAULT_COLOCALIZATION_THRESHOLD    = 300
DEFAULT_CLUSTER_RADIUS              = 300
DEFAULT_NUM_CORES                   = 2


# =====================================================================================================================
# UTILITY FUNCTIONS


def read_config(config_filename):
    if not os.path.exists(os.path.expanduser(config_filename)):
        print(f'[WARNING]: Config file does not exist at: "{config_filename}".')
        return

    cp = ConfigParser(allow_no_value=True,
                      converters={'list': 
                                  lambda x: [e.strip() for e in x.split(',')]})
    cp.read(config_filename)
    return cp


# ---------------------------------------------------------------------------------------------------------------------    


def validate_config(config):
    expected_keys = ('crop_percentage proteins protein_types '
                     'voxel_size_nm data_directory min_puncta_per_cluster '
                     'plot_type puncta_radius_nm colocalization_threshold_nm '
                     'cluster_radius_nm interactive_mode num_cores').split()
    keys_with_type = {'voxel_size_nm': int, 
                      'crop_percentage': int,
                      'min_puncta_per_cluster': int,
                      'puncta_radius_nm': int,
                      'colocalization_threshold_nm': int,
                      'plot_type': str,
                      'data_directory': str,
                      'cluster_radius_nm': int,
                      'interactive_mode': bool,
                      'num_cores': int}
    keys_with_required_lengths = {'puncta_radius_nm': {'min': 2, 'max': 3},
                                  'voxel_size_nm': {'min': 2, 'max': 3}}
    if 'default' not in config.sections():
        raise RuntimeError('Default section not found in config. Exiting.')

    parsed_config = dict()
    for key in expected_keys:
        values = config['default'].getlist(key)
        if key in keys_with_type:
            key_type = keys_with_type[key]
            if len(values) == 1:
                parsed_config[key] = key_type(values[0])
            else:
                parsed_config[key] = [key_type(v) for v in values]
        else:
            parsed_config[key] = values
            if len(values) == 1 and len(values[0]) == 0:
                print(f'Warning: key "{key}" does not have a value assigned.')
                parsed_config[key] = None            

    # Check lengths
    for key in keys_with_required_lengths:
        values = parsed_config[key]
        min_length = keys_with_required_lengths[key]['min']
        max_length = keys_with_required_lengths[key]['max']
        key_as_flag = '--%s' % key.replace('_', '-')
        if len(values) < min_length:
            raise RuntimeError(f'The passed argument for "{key_as_flag}" '
                               f'contains less than {min_length} arguments.')
        if len(values) > max_length:
            raise RuntimeError(f'The passed argument for "{key_as_flag}" '
                               f'contains more than {max_length} arguments.')
    return parsed_config


# ---------------------------------------------------------------------------------------------------------------------    


def validate_command_line_args(args):
    cannot_be_none = ('crop_percentage data_directory proteins protein_types '
                      'voxel_size min_puncta_per_cluster plot_type '
                      'puncta_radius_nm colocalization_threshold_nm '
                      'cluster_radius_nm interactive_mode num_cores').split()
    parsed = vars(args)
    for key in parsed:
        value = parsed[key]
        if key in cannot_be_none and value is None:
            key_as_flag = '--%s' % key.replace('_', '-')
            raise RuntimeError(f'The passed argument for "{key_as_flag}" '
                                'cannot be None.')
    return parsed        


# ---------------------------------------------------------------------------------------------------------------------


def create_directories(config_dict):
    basedir = os.path.basename(config_dict['data_directory'])
    percent = config_dict['crop_percentage']
    config_dict['ARRAYS_DIR'] = os.path.join(DEFAULT_ARRAYS_DIR, 
                                             basedir, 
                                             f'{percent}')
    config_dict['PLOTS_DIR'] = os.path.join(DEFAULT_PLOTS_DIR, 
                                            basedir, 
                                            f'{percent}')
    config_dict['RESULTS_DIR'] = os.path.join(DEFAULT_RESULTS_DIR, 
                                              basedir, 
                                              f'{percent}')
    dirs = [config_dict['ARRAYS_DIR'], 
            config_dict['PLOTS_DIR'], 
            config_dict['RESULTS_DIR']]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
    return config_dict


# ---------------------------------------------------------------------------------------------------------------------


def assign_config(args):
    config = dict()

    if args.config_file is not None:
        raw_config = read_config(args.config_file)
        if raw_config is not None:
            config = validate_config(raw_config)
        else:
            raise RuntimeError('Unable to read config file.')
    else:
        config = validate_command_line_args(args)
    return config


# ---------------------------------------------------------------------------------------------------------------------


def parse_config(config_dict):
    # For easier use - allows dot notation
    Config = namedtuple('Config', config_dict.keys())
    config = Config(**config_dict)
    return config


# ---------------------------------------------------------------------------------------------------------------------


# Calculate the centroid from an array of 2D points:
# https://stackoverflow.com/a/23021198/866930
def centeroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    return sum_x/length, sum_y/length


# ---------------------------------------------------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------------------------------------------------


def load_files(directory, names_of_interest, filetype):
    files_of_interest = list()
    name_table = {name:list() for name in names_of_interest}
    names = set(names_of_interest)
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        name_parts = filename.replace('-', '_').split('_')
        intersection = names.intersection(name_parts)
        if len(intersection) > 0 and filename.endswith(filetype):
            files_of_interest.append(filepath)
            name = list(intersection)[0]
            name_table[name].append(filepath)
    for name in name_table:
        values = list(sorted(name_table[name]))
        name_table[name] = values
    return name_table, files_of_interest


# ---------------------------------------------------------------------------------------------------------------------


def validate_mutants_and_WT(name_table, protein_types):
    type_table = dict()
    num_files = list()
    for protein_name in name_table:
        files = name_table[protein_name]
        num_files.append(len(files))
    
    if len(list(set(num_files))) > 1:
        raise RuntimeError("There number of files between proteins differs.")


# =====================================================================================================================
# ANALYSIS FUNCTIONS AND CLASSES


class DetectAndExportPuncta(object):
    def __init__(self, image_filename, config):
        self.image_filename = image_filename
        self.config = config

        self.filename = os.path.basename(image_filename)
        self.fpath = os.path.splitext(self.filename)[0]
        self.ARRAYS_DIR = config.ARRAYS_DIR
        self.PLOTS_DIR = config.PLOTS_DIR
        self.RESULTS_DIR = config.RESULTS_DIR
        self.rna = None

        self.ndim = len(config.puncta_radius_nm)
        self.puncta_radius_nm = config.puncta_radius_nm
        self.puncta_radius_px = None
        self.voxel_size = self.config.voxel_size_nm[:self.ndim]
        self.puncta = None
        self.clusters = None
        self.sorted_puncta_post_clustering = None
    

    def export_puncta_and_clusters(self):
        header = [ f'# Filename: {self.filename}',
                   f'# Puncta radius (pixels): {self.puncta_radius_px}', 
                   f'# Puncta radius (nm): {self.config.puncta_radius_nm}', 
                   f'# Cluster threshold (nm): {self.config.cluster_radius_nm}', 
                   f'# Min puncta per cluster: {self.config.min_puncta_per_cluster}', 
                   f'# Voxel size (nm): {self.config.voxel_size_nm[:self.ndim]}']

        puncta_savename = os.path.join(self.RESULTS_DIR, f'{self.fpath}_puncta.txt')
        d = dict()
        d['Y'] = self.puncta[:,0]
        d['X'] = self.puncta[:,1]
        df = pd.DataFrame(d)
        df.to_fwf(puncta_savename, header_lines=[l for l in header if 'Min' not in l])

        clusters_savename = os.path.join(self.RESULTS_DIR, f'{self.fpath}_clusters.txt')
        d = dict()
        d['Y']          = self.clusters[:,0]
        d['X']          = self.clusters[:,1]
        d['NUM_PUNCTA'] = self.clusters[:,2]
        d['CLUSTER_ID'] = self.clusters[:,3]
        df = pd.DataFrame(d)
        df.to_fwf(clusters_savename, header_lines=header)

        header = [ f'# Filename: {self.filename}',
                   f'# Puncta radius (pixels): {self.puncta_radius_px}', 
                   f'# Puncta radius (nm): {self.config.puncta_radius_nm}', 
                   f'# Cluster threshold (nm): {self.config.cluster_radius_nm}', 
                   f'# Min puncta per cluster: {self.config.min_puncta_per_cluster}', 
                   f'# Voxel size (nm): {self.config.voxel_size_nm}']

        clustered_puncta_savename = os.path.join(self.RESULTS_DIR, f'{self.fpath}_clustered_puncta.txt')
        d = dict()
        d['Y']          = self.sorted_puncta_post_clustering[:,0]
        d['X']          = self.sorted_puncta_post_clustering[:,1]
        d['CLUSTER_ID'] = self.sorted_puncta_post_clustering[:,2]
        df = pd.DataFrame(d)
        df.to_fwf(clustered_puncta_savename, header_lines=header)


    def stage1_detect_puncta(self):
        # Read image and begin analysis
        rna = stack.read_image(self.image_filename)
        self.rna = rna
        print("smfish channel")
        print("\r shape: {0}".format(rna.shape))
        print("\r dtype: {0}".format(rna.dtype))

        # puncta radius
        puncta_radius_px = detection.get_object_radius_pixel(voxel_size_nm=self.voxel_size, 
                                                             object_radius_nm=self.puncta_radius_nm, 
                                                             ndim=self.ndim)
        self.puncta_radius_px = puncta_radius_px
        print('Puncta radius in pixels: ', puncta_radius_px)
        
        # LoG filter
        rna_log = stack.log_filter(rna, sigma=puncta_radius_px)

        # local maximum detection
        mask = detection.local_maximum_detection(rna_log, min_distance=puncta_radius_px)

        # thresholding
        threshold = detection.automated_threshold_setting(rna_log, mask)
        puncta, _ = detection.spots_thresholding(rna_log, mask, threshold)
        self.puncta = puncta

        print("Detected puncta:")
        print("\r shape: {0}".format(puncta.shape))
        print("\r dtype: {0}".format(puncta.dtype))
        print("\r threshold: {0}".format(threshold))

        puncta_savepath = os.path.join(self.PLOTS_DIR, f'{self.fpath}_puncta')
        elbow_savepath = os.path.join(self.PLOTS_DIR, f'{self.fpath}_elbow')
        reference_puncta_savepath = os.path.join(self.PLOTS_DIR, f'{self.fpath}_reference')

        plot.plot_detection(self.rna,
                            puncta, 
                            contrast=True, 
                            show=self.config.interactive_mode, 
                            path_output=puncta_savepath,
                            ext=self.config.plot_type)

        plot.plot_elbow(images=self.rna, 
                        voxel_size=self.config.voxel_size_nm[:self.ndim], 
                        spot_radius=self.puncta_radius_nm,
                        show=self.config.interactive_mode,
                        path_output=elbow_savepath,
                        ext=self.config.plot_type)

        # save results
        puncta_npy_savename = os.path.join(self.ARRAYS_DIR, f'{self.fpath}_puncta.npy')
        stack.save_array(puncta, puncta_npy_savename)


    def stage2_dense_region_decomposition(self):
        reference_puncta_savepath = os.path.join(self.PLOTS_DIR, f'{self.fpath}_reference')

        # gaussian kernel
        kernel_size = detection.get_object_radius_pixel(voxel_size_nm=self.voxel_size, 
                                                        object_radius_nm=self.puncta_radius_nm, 
                                                        ndim=self.ndim)
        large_kernel_size = tuple([kernel_size_ * 5 for kernel_size_ in kernel_size])

        # denoising
        rna_denoised = stack.remove_background_gaussian(self.rna, sigma=large_kernel_size)

        # reference puncta
        reference_puncta = detection.build_reference_spot(image=rna_denoised,
                                                          spots=self.puncta,
                                                          voxel_size=self.voxel_size, 
                                                          spot_radius=self.puncta_radius_nm,
                                                          alpha=0.7)

        plot.plot_reference_spot(reference_puncta,
                                 show=self.config.interactive_mode,
                                 path_output=reference_puncta_savepath)

        # fit a gaussian function on the reference puncta
        (sigma_yx, 
        amplitude, 
        background) = detection.modelize_spot(reference_spot=reference_puncta, 
                                              voxel_size=self.voxel_size,
                                              spot_radius=self.puncta_radius_nm)

        # detect dense regions
        (regions_to_decompose, 
        puncta_out_regions, 
        region_size) = detection.get_dense_region(image=rna_denoised, 
                                                  spots=self.puncta,
                                                  voxel_size=self.voxel_size,
                                                  spot_radius=self.puncta_radius_nm,
                                                  beta=1)

        # precompute gaussian function values
        max_grid = max(200, region_size + 1)
        precomputed_gaussian = detection.precompute_erf(ndim=self.ndim,
                                                        voxel_size=self.voxel_size,
                                                        sigma=(sigma_yx, sigma_yx),
                                                        max_grid=max_grid)

        # simulate gaussian mixtures
        puncta_in_regions, _ = detection.simulate_gaussian_mixture(image=rna_denoised,
                                                                   candidate_regions=regions_to_decompose,
                                                                   voxel_size=self.voxel_size,
                                                                   sigma=(sigma_yx, sigma_yx),
                                                                   amplitude=amplitude,
                                                                   background=background,
                                                                   precomputed_gaussian=precomputed_gaussian)

        puncta_post_decomposition = np.concatenate((puncta_out_regions, 
                                                    puncta_in_regions[:, :self.ndim]), 
                                                    axis=0)
        self.puncta_post_decomposition = puncta_post_decomposition
        print("Detected puncta before decomposition")
        print("\r shape: {0}".format(self.puncta.shape))
        print("\r dtype: {0}".format(self.puncta.dtype), "\n")
        print("Detected puncta after decomposition")
        print("\r shape: {0}".format(puncta_post_decomposition.shape))
        print("\r dtype: {0}".format(puncta_post_decomposition.dtype))

        reference_puncta_npy_savename = os.path.join(self.ARRAYS_DIR, f'{self.fpath}_reference.npy')
        stack.save_array(reference_puncta, reference_puncta_npy_savename)


    def stage3_detect_clusters(self):
        clusters_savepath = os.path.join(self.PLOTS_DIR, f'{self.fpath}_clusters')

        puncta_post_clustering, clusters = detection.detect_clusters(spots=self.puncta_post_decomposition, 
                                                                     voxel_size=self.voxel_size, 
                                                                     radius=self.config.cluster_radius_nm, 
                                                                     nb_min_spots=self.config.min_puncta_per_cluster)
        self.clusters = clusters
        
        print("Detected puncta after clustering")
        print("\r shape: {0}".format(puncta_post_clustering.shape))
        print("\r dtype: {0}".format(puncta_post_clustering.dtype), "\n")
        print("Detected clusters")
        print("\r shape: {0}".format(clusters.shape))
        print("\r dtype: {0}".format(clusters.dtype))

        # data appears to be of the form y, x, num, id
        sorted_puncta_indices = puncta_post_clustering[:,self.ndim].argsort(kind='mergesort')
        sorted_puncta_post_clustering = puncta_post_clustering[sorted_puncta_indices]
        self.sorted_puncta_post_clustering = sorted_puncta_post_clustering

        # export our progress thus far
        self.export_puncta_and_clusters()
        
        # plot
        plot.plot_detection(self.rna, 
                            spots=[self.puncta_post_decomposition, clusters[:,:self.ndim]], 
                            shape=["circle", "polygon"], 
                            radius=[3, 6], 
                            color=["red", "blue"],
                            linewidth=[1, 2], 
                            fill=[False, True], 
                            contrast=True,
                            show=self.config.interactive_mode,
                            path_output=clusters_savepath)

        clusters_npy_savename = os.path.join(self.ARRAYS_DIR, f'{self.fpath}_clusters.npy')
        stack.save_array(clusters, clusters_npy_savename)


    def detect_and_export_puncta_and_clusters(self):
        self.stage1_detect_puncta()
        self.stage2_dense_region_decomposition()
        self.stage3_detect_clusters()


# ---------------------------------------------------------------------------------------------------------------------


def colocalize(image_filepath,
               config_dict,
               protein1, 
               protein2, 
               protein1_puncta, 
               protein2_puncta, 
               suffix):
    config = parse_config(config_dict)
    basename = os.path.basename(image_filepath).replace('.tif', '')
    plots_savepath = config.PLOTS_DIR
    arrays_savepath = config.ARRAYS_DIR
    results_savepath = config.RESULTS_DIR
    ndim = len(config.puncta_radius_nm)

    (puncta_1_colocalized, 
    puncta_2_colocalized, 
    distances, 
    indices_1, 
    indices_2, 
    threshold) = multistack.detect_spots_colocalization(spots_1=protein1_puncta, 
                                                        spots_2=protein2_puncta,
                                                        voxel_size=config.voxel_size_nm[:ndim],
                                                        threshold=config.colocalization_threshold_nm,
                                                        return_indices=True,
                                                        return_threshold=True)

    print("Colocalized puncta")
    print("\r shape 1: {0}".format(puncta_1_colocalized.shape))
    print("\r shape 2: {0}".format(puncta_2_colocalized.shape))
    print("\r distances: {0}".format(distances.shape))
    print("\r indices 1: {0}".format(indices_1.shape))
    print("\r indices 2: {0}".format(indices_2.shape))
    print("\r threshold: {0:0.2f} nm".format(threshold))
    
    print('Puncta 1:', len(puncta_1_colocalized)/(len(protein1_puncta)))
    print('Puncta 2: ', len(puncta_2_colocalized)/len(protein2_puncta))
    print(len(protein1_puncta), len(protein2_puncta))

    d = dict()
    d['Filename']           = [os.path.basename(image_filepath)]
    d['Protein-1']          = [protein1]
    d['Protein-2']          = [protein2]
    d['Num-puncta-1']       = [len(protein1_puncta)]
    d['Num-puncta-2']       = [len(protein2_puncta)]
    d['Num-Colocalized']    = [len(puncta_1_colocalized)]
    d['Colocalization-1']   = [len(puncta_1_colocalized)/len(protein1_puncta)]
    d['Colocalization-2']   = [len(puncta_2_colocalized)/len(protein2_puncta)]
    d['Threshold']          = [threshold]
    d_savepath = os.path.join(arrays_savepath, f'{basename}___{protein1}_{protein2}___{suffix}.npy')
    np.save(d_savepath, d, allow_pickle=True)

    # plot
    image = stack.read_image(image_filepath)
    detection_savename = os.path.join(plots_savepath, f'{basename}___{protein1}_{protein2}___{suffix}')
    plot.plot_detection(image, 
                        spots=[puncta_1_colocalized, puncta_2_colocalized], 
                        radius=2, 
                        color=["magenta", "red"],
                        title=f"Magenta= {protein1} | Red = {protein2}",
                        linewidth=2, 
                        contrast=True,
                        show=config.interactive_mode,
                        path_output=detection_savename)

    colocalization_elbow_savename = os.path.join(plots_savepath, f'{basename}___{protein1}_{protein2}___elbow_{suffix}')
    plot.plot_elbow_colocalized(spots_1=protein1_puncta,
                                spots_2=protein2_puncta, 
                                voxel_size=config.voxel_size_nm[:ndim],
                                show=config.interactive_mode,
                                path_output=colocalization_elbow_savename)


# ---------------------------------------------------------------------------------------------------------------------


def detect_and_export_puncta_from_files(filepaths, config, num_processes=None):
    num_cpus = mp.cpu_count() * 2
    if num_processes is not None and type(num_processes) is int:
        num_cpus = num_processes
    pool = mp.Pool(processes=num_cpus)
    for filepath in filepaths:
        detector_and_exporter = DetectAndExportPuncta(filepath, config)
        pool.apply_async(detector_and_exporter.detect_and_export_puncta_and_clusters(), args=())
    pool.close()
    pool.join()


# ---------------------------------------------------------------------------------------------------------------------


def calculate_colocalization(protein_files_table, config_dict, num_processes=None):
    num_cpus = mp.cpu_count() * 2
    if num_processes is not None and type(num_processes) is int:
        num_cpus = num_processes
    pool = mp.Pool(processes=num_cpus)
    config = parse_config(config_dict)

    protein1, protein2 = config.proteins
    protein1_files = protein_files_table[protein1]
    protein2_files = protein_files_table[protein2]
    protein_arrays = dict()
    for protein in config.proteins:
        protein_arrays[protein] = list()

    for protein1_file, protein2_file in zip(protein1_files, protein2_files):
        basename1 = os.path.basename(protein1_file)
        basename2 = os.path.basename(protein2_file)
        path1 = os.path.join(config.ARRAYS_DIR, basename1.replace('.tif', '_puncta.npy'))
        path2 = os.path.join(config.ARRAYS_DIR, basename2.replace('.tif', '_puncta.npy'))
        protein1_puncta = np.load(path1)
        protein2_puncta = np.load(path2)
        
        bpath = os.path.basename(basename1).replace('.tif', '')
        p1 = os.path.join(config.ARRAYS_DIR, f'{bpath}___{protein1}_{protein2}___{protein1}.npy')
        cpath = os.path.basename(basename2).replace('.tif', '')
        p2 = os.path.join(config.ARRAYS_DIR, f'{cpath}___{protein1}_{protein2}___{protein2}.npy')
        protein_arrays[protein1].append(p1)
        protein_arrays[protein2].append(p2)

        # `mp.pool`` uses `mp.SimpleQueue` to pass tasks to workers. Everything passed through `mp.SimpleQueue`
        # must be pickleable. Hence, while `config` is the natural choice, we instead pass the `config_dict`
        # and convert it to a Config object. The alternative is to ensure that `Config` is global, however
        # given its dynamic initialization, that's not the best approach. This is a reasonable workaround.
        args1 = (protein1_file, config_dict, protein1, protein2, protein1_puncta, protein2_puncta, protein1,)
        args2 = (protein2_file, config_dict, protein1, protein2, protein1_puncta, protein2_puncta, protein2,)
        pool.apply_async(colocalize, args=args1)
        pool.apply_async(colocalize, args=args2)
    pool.close()
    pool.join()

    # join results
    all_dfs = list()
    for puncta_path1, puncta_path2 in zip(protein_arrays[protein1], protein_arrays[protein2]):
        p1 = np.load(puncta_path1, allow_pickle=True).item()
        p2 = np.load(puncta_path2, allow_pickle=True).item()
        all_dfs += [pd.DataFrame(p1), pd.DataFrame(p2)]
    df = pd.concat(all_dfs, ignore_index=True)
    df.to_fwf(os.path.join(config.RESULTS_DIR, f'{config.crop_percentage}___{protein1}_{protein2}.tsv'))


'''
Want to titrate the colocalization of BNI1 and CLN3 for each figure (WT / Ctrl vs Mut).
puncta is in distance of pixels not nm
'''
def titrate_colocalization(image_filepath,
                           config_dict,
                           protein1,
                           protein2,
                           protein1_puncta,
                           protein2_puncta,
                           suffix,
                           pixel_shifts):
    config = parse_config(config_dict)
    ndim = len(config.puncta_radius_nm)
    colocalizations = {p:list() for p in [protein1, protein2]}
    voxel_size = config.voxel_size_nm[:ndim]
    for pixel_shift in pixel_shifts:
        puncta1_nm = detection.convert_spot_coordinates(spots=protein1_puncta, voxel_size=voxel_size) + pixel_shift
        puncta2_nm = detection.convert_spot_coordinates(spots=protein2_puncta, voxel_size=voxel_size)

        distance_matrix = cdist(puncta1_nm, puncta2_nm)
        
        # assign spots based on their euclidean distance
        indices_1, indices_2 = linear_sum_assignment(distance_matrix)

        # get distance between colocalized spots
        distances = distance_matrix[indices_1, indices_2]

        # keep colocalized spots within a specific distance
        mask = distances <= config.colocalization_threshold_nm
        indices_1 = indices_1[mask]
        indices_2 = indices_2[mask]
        distances = distances[mask]

        # get colocalized spots
        puncta_1_colocalized = protein1_puncta[indices_1, ...]
        puncta_2_colocalized = protein2_puncta[indices_2, ...]
    
        cl_protein1 = len(puncta_1_colocalized)/len(protein1_puncta)
        cl_protein2 = len(puncta_2_colocalized)/len(protein2_puncta)
        colocalizations[protein1].append(cl_protein1)
        colocalizations[protein2].append(cl_protein2)
        print(image_filepath, pixel_shift, cl_protein1, cl_protein2)

    print(colocalizations[protein1])
    print(colocalizations[protein2])
    print()
    
    width = 0.4
    x = np.arange(len(pixel_shifts))
    fig = plt.figure(figsize=(30, 8))
    ax = fig.add_subplot(111)
    rects1 = ax.bar(x - width/2, colocalizations[protein1], width, label=protein1)
    rects2 = ax.bar(x + width/2, colocalizations[protein2], width, label=protein2)
    ax.bar_label(rects1, padding=3)
    ax.bar_label(rects2, padding=3)
    labels = ['%d' % r for r in pixel_shifts]
    ax.set_xticks(x, labels)
    ax.set_title(suffix)
    ax.set_xlabel('Pixel Shift (nm)')
    ax.set_ylabel('Fraction colocalized')
    ax.legend()
    fig.tight_layout()

    bpath = os.path.basename(image_filepath).replace('.tif', '')
    fig.savefig(os.path.join(config.PLOTS_DIR, f'pixelshifts_{bpath}_{protein1}_{protein2}.png'))


def titrate_all_colocalizations(protein_files_table, config_dict, num_processes=None):
    num_cpus = mp.cpu_count() * 2
    if num_processes is not None and type(num_processes) is int:
        num_cpus = num_processes
    pool = mp.Pool(processes=num_cpus)
    config = parse_config(config_dict)

    protein1, protein2 = config.proteins
    protein1_files = protein_files_table[protein1]
    protein2_files = protein_files_table[protein2]
    protein_arrays = dict()
    for protein in config.proteins:
        protein_arrays[protein] = list()

    pixel_shifts = range(100, 3100, 100)
    for protein1_file, protein2_file in zip(protein1_files, protein2_files):
        basename1 = os.path.basename(protein1_file)
        basename2 = os.path.basename(protein2_file)
        path1 = os.path.join(config.ARRAYS_DIR, basename1.replace('.tif', '_puncta.npy'))
        path2 = os.path.join(config.ARRAYS_DIR, basename2.replace('.tif', '_puncta.npy'))
        protein1_puncta = np.load(path1)
        protein2_puncta = np.load(path2)

        bpath = os.path.basename(basename1).replace('.tif', '')
        p1 = os.path.join(config.ARRAYS_DIR, f'{bpath}___{protein1}_{protein2}___{protein1}.npy')
        cpath = os.path.basename(basename2).replace('.tif', '')
        p2 = os.path.join(config.ARRAYS_DIR, f'{cpath}___{protein1}_{protein2}___{protein2}.npy')
        protein_arrays[protein1].append(p1)
        protein_arrays[protein2].append(p2)

        # `mp.pool`` uses `mp.SimpleQueue` to pass tasks to workers. Everything passed through `mp.SimpleQueue`
        # must be pickleable. Hence, while `config` is the natural choice, we instead pass the `config_dict`
        # and convert it to a Config object. The alternative is to ensure that `Config` is global, however
        # given its dynamic initialization, that's not the best approach. This is a reasonable workaround.
        args1 = (protein1_file, config_dict, protein1, protein2, protein1_puncta, protein2_puncta, protein1, pixel_shifts,)
        args2 = (protein2_file, config_dict, protein1, protein2, protein1_puncta, protein2_puncta, protein2, pixel_shifts,)
        #titrate_colocalization(*args1)
        pool.apply_async(titrate_colocalization, args=args1)
        pool.apply_async(titrate_colocalization, args=args2)
    pool.close()
    pool.join()
    

# =====================================================================================================================
# DOCUMENTATION DESCRIPTORS THAT REQUIRE ELABORATION


CROP_PERCENTAGE_HELP   = ('The amount to crop the protein images by for '
                          'analysis. If values less than 100%% are provided, a '
                          'full analysis will be performed, and the centroid '
                          'will be used for subsequent analysis.')

COLOCALIZATION_HELP     = ('The interpuntca distance (nm) to use for '
                           'determining colocalization.')

MIN_PUNCTA_HELP         = ('The minimum number of puncta to consider as part '
                           'of a cluster.')

CLUSTER_RADIUS_HELP     = ('The radius (nm) by which to determine whether '
                           'puncta belong to a cluster. Used in conjunction '
                           'with "--min-puncta".')

INTERACTIVE_MODE_HELP   = ('Whether or not to show plots for interactive '
                           'introspection.')


# =====================================================================================================================
# MAIN FUNCTION


def main():
    parser = ArgumentParser()
    parser.add_argument('-cf',
                        '--config-file',
                        help='The INI file containing the configuration.',
                        type=str)

    parser.add_argument('-cp',
                        '--crop-percentage',
                        help=CROP_PERCENTAGE_HELP,
                        type=int,
                        nargs='+',
                        default=DEFAULT_CROP_PERCENTAGE)

    parser.add_argument('-cr',
                        '--cluster-radius-nm',
                        help=CLUSTER_RADIUS_HELP,
                        type=int,
                        default=DEFAULT_CLUSTER_RADIUS)

    parser.add_argument('-ct',
                        '--colocalization-threshold-nm',
                        help=COLOCALIZATION_HELP,
                        type=int,
                        default=DEFAULT_COLOCALIZATION_THRESHOLD)

    parser.add_argument('-d',
                        '--data-directory',
                        help='The directory containing the data / TIFF files.',
                        type=str)

    parser.add_argument('-i',
                        '--interactive-mode',
                        help=INTERACTIVE_MODE_HELP,
                        action='store_true')

    parser.add_argument('-mp',
                        '--min-puncta-per-cluster',
                        help=MIN_PUNCTA_HELP,
                        type=int,
                        default=DEFAULT_MIN_PUNCTA)

    parser.add_argument('-nc',
                        '--num-cores',
                        help='The number of cores to use for analysis and exporting.',
                        type=int,
                        default=DEFAULT_NUM_CORES)

    parser.add_argument('-p',
                        '--proteins',
                        help='The RNA protein TIFF images to analyze.',
                        type=str,
                        nargs='+',
                        default=DEFAULT_PROTEINS)
    
    parser.add_argument('-pr',
                        '--puncta-radius-nm',
                        help='The radius (nm) of the puncta.',
                        type=int,
                        nargs='+',
                        default=DEFAULT_PUNCTA_RADIUS)

    parser.add_argument('-pt',
                        '--plot-type',
                        help='The extension to use for saving the plots.',
                        type=str,
                        default=DEFAULT_PLOT_TYPE)

    parser.add_argument('-t',
                        '--protein-types',
                        help='The types of proteins (ctrl / bidi)',
                        type=str,
                        default=DEFAULT_PROTEIN_TYPES)

    parser.add_argument('-v',
                        '--voxel-size-nm',
                        help='The voxel size (nm) of the TIF images.',
                        type=int,
                        nargs='+',
                        default=DEFAULT_VOXEL_SIZE)

    args = parser.parse_args()
    config_dict = assign_config(args)
    config_dict = create_directories(config_dict)
    config = parse_config(config_dict)
    pprint(config)

    files_by_protein, files = load_files(config.data_directory, 
                                         config.proteins, 
                                         DEFAULT_IMAGE_FILETYPE)

    files = [os.path.join(config.data_directory, f) for f in os.listdir(config.data_directory) if '.tif' in f and 'Nuclei' not in f]
    
    detect_and_export_puncta_from_files(files, config)
    calculate_colocalization(files_by_protein, config_dict)


if __name__ == '__main__':
    main()
