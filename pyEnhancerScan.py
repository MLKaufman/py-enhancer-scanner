import pyBigWig
import os
from os import path
import pandas as pd
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

from lxml import html
import requests
import statistics
import gzip
import shutil

from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from primer3plus import Design

VERSION = 0.6

# Version 1.0 Feature TODO list:
# Compare CHIP to motifs
# self.load_ChIP('track', tfactor) # plot as overlay on motif data? plot as overlay
# confirm accurate coordinates for motifs
# confirm position / other attriutes of jaspar scanner

class EnhancerScan:
    """ Class to handle scanning of bigwig files. """
    def __init__(self):
        pd.options.display.max_rows = 4000
        pd.options.display.max_columns = 4000

        self.track1_name = None
        self.track1_bw = None
        self.track1_values = None
        self.track1_min_value = None
        self.track1_max_value = None

        self.track2_name = None
        self.track2_bw = None
        self.track2_values = None
        self.track2_min_value = None
        self.track2_max_value = None

        self.track_comparison_type = None

        self.track_plot_header = None

        self.genome = None
        self.chromosome = None

        self.region_start = None
        self.region_stop = None

        self.region_values = None
        self.region_max_value = None
        self.region_mean_value = None
        self.region_median_value = None
        
        self.all_detected_peaks = None

        self.auto_threshold = None

        self.df_results = None
        self.df_motifs = None
        self.reset_results() # also generates a new dataframe for self.df_results and self.df_motifs

        ### Streamlit App functionality
        self.last_figure = None # conforms to streamlit depreciating global pyplot functionality

        print("pyEnhancerScanner version " + str(VERSION) + " beta\n")
        print("Important changes: all write functions have been renamed to save. write_stats is now save_results\n")
        print("The following tracks are locally available:")
        for track in self.list_tracks(): print(track)
        print("")
        print("To download additional tracks, use the download_tracks() function or list them with list_external_tracks().")

    def list_external_tracks(self):
        """ Alias for download"""
        print('External Track List:')
        df_quicklist = pd.read_csv('external_tracks.db', delimiter=',', header=0)
        #pd.set_option('display.max_colwidth', None) # so it doesnt truncate columns
        #print(df_quicklist.loc[:, :'Size'])
        #print("")
        print("To download one of these tracks, use download_tracks(track_num=X) where X is the track number / index.")
        print("You can specify your own download url by download_tracks(url='X').")
        return df_quicklist

    def load_track(self, genome, track1, track2=None, track_comparison_type=None, reset=True): #TODO add multitrack here?
        if reset is True:
            self.reset_results()

        if track1 is not None and genome is not None:
            
            self.genome = genome
            self.track1_name = track1
            self.track1_bw = pyBigWig.open(track1)

            self.track_plot_header = track1

            header = self.track1_bw.header()
            self.track1_min_value = header['minVal']
            self.track1_max_value = header['maxVal']

        if track2 is not None and track_comparison_type is not None:
            self.track2_bw = pyBigWig.open(track2)
            self.track2_name = track2

            track2_header = self.track2_bw.header()
            self.track2_min_value = track2_header['minVal']
            self.track2_max_value = track2_header['maxVal']


            self.track_comparison_type = track_comparison_type
            self.track_plot_header = track1 + ' ' + track_comparison_type + ' ' + track2

        else:
            return(print("Error: You must load a track by specifyinhg a bigwig track and a genome."))

    def download_tracks(self, url = '', track_num = 0):
        
        if url == '' and track_num == 0:
            self.list_external_tracks()
        
        elif url != '':
            if type(url) is int:
                df_quicklist = pd.read_csv('external_tracks.db', delimiter=',', header=0)

                url = df_quicklist.loc[url, 'URL_Path']
                self.download_url(url)
            else:
                self.download_url(url)

        elif track_num !=0:
            df_quicklist = pd.read_csv('external_tracks.db', delimiter=',', header=0)

            url = df_quicklist.loc[track_num, 'URL_Path']
            self.download_url(url)

    def download_jaspar(self, url = ''):
        if url == '':
            url = 'http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt'

        self.download_url(url)

    def list_tracks(self):
        tracks = []
        for bigwig in os.listdir(os.getcwd()):
            if bigwig.endswith(".bw"):
                if bigwig !='None':
                    tracks.append(bigwig)
        return tracks
                    
    def list_bed(self):
        list_bedfiles = []
        for bedfile in os.listdir(os.getcwd()):
            if bedfile.endswith(".bed"):
                if bedfile !='None':
                    list_bedfiles.append(bedfile)
        return list_bedfiles
    
    def load_bed(self, genome, bed_file):
        """ Load an existing bed file for analyzing regions enhancers. """
        # load a bed file
        # for each entry in bed file, create an entry in df_results.
        # calculate 
        self.reset_results()

        self.genome = genome
        ## read in bed file and convert to dataframe
        df_bed = pd.read_csv(bed_file, sep='\t', comment='t', header=None)
        header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes', 'blockStarts']
        df_bed.columns = header[:len(df_bed.columns)]
        print(df_bed)

        list_region_start = []
        list_region_stop = []

        #convert to df_results
        for index, enhancer in df_bed.iterrows():
            chromosome = enhancer['chrom']
            self.chromosome = chromosome #from bed files
            start = enhancer['chromStart']
            stop = enhancer['chromEnd']
            name = enhancer['name']
            full_name = enhancer['name']
            size = enhancer['chromEnd'] - enhancer['chromStart']
            primerf = ''
            primerftemp = ''
            primerr = ''
            primerrtemp = ''
            mean_peak_values = self.get_mean_range_values(chromosome, start, stop) #should be bqsed on track
            max_peak_values = self.get_max_range_values(chromosome, start, stop)
            sequence = self.get_sequence(self.genome, chromosome, start, stop)
            
            list_region_start.append(start)
            list_region_stop.append(stop)

            self.df_results.loc[len(self.df_results)]=([chromosome, start, stop, name, full_name, mean_peak_values, max_peak_values, size, primerf, primerftemp, primerr, primerrtemp, sequence])

        self.region_start = sorted(list_region_start)[0]
        self.region_stop = sorted(list_region_stop)[-1]

        if self.track1_bw != None:
            self.region_values = self.track1_bw.values(self.chromosome, self.region_start, self.region_stop) # but we need a track loaded....

    def enhancer_scanner(self, chromosome, region_start=0, region_stop=0, peak_width=50, peak_height='auto', merge_distance=200, final_size=None, final_mean_peak_values=None, reset=True):
        if reset is True:
            self.reset_results()

        #allow direct copy from ucsc genome browser
        if ':' in chromosome:
            temp = chromosome
            chromosome, temp = temp.split(':')
            region_start, region_stop = temp.replace(',','').split('-')
            
        self.chromosome = chromosome
        self.region_start = int(region_start)
        self.region_stop = int(region_stop)

        ### GET REGION VALUES universal for single or multiple tracks

        self.get_region_values()

        if peak_height == 'auto':
            self.auto_threshold = self.region_max_value * .25
            peak_height = self.auto_threshold
            print("")
            print("Using auto detection of peak height: ", self.auto_threshold)

        elif peak_height == 'mean':
            peak_height = float(self.region_mean_value)

        elif peak_height == 'median':
            peak_height = float(self.region_median_value)
        else:
            peak_height = float(peak_height)

        # PEAK DETECTION
        # detect peaks, returns a tuple of an array of peaks and a dictionary or properties
        peaks = signal.find_peaks(self.region_values, width = float(peak_width), height = peak_height)

        # make list of unqiue peak widths as tuples and sort
        list_widths = sorted(list(set(zip(list(peaks[1]['left_bases'] + self.region_start), list(peaks[1]['right_bases'] + self.region_start))))) # its fixed for final location from widths now
        print('Total Peaks Detected:', len(list_widths))
        #print(list_widths)
        
        #TODO: clean up:
        # merge overlapping tuples and tuples within a certain distance
        list_merged = []

        for higher in sorted(list_widths, key=lambda tup: tup[0]): #sorted by lower bound
            if not list_merged:
                list_merged.append(higher)
            else:
                lower = list_merged[-1]
                # test for intersection between lower and higher:
                # we know via sorting that lower[0] <= higher[0]
                if higher[0] <= lower[1]:
                    upper_bound = max(lower[1], higher[1])
                    list_merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                elif higher[0] - lower[1] < merge_distance: #seems to work #TODO: confirm this works and whether it needs to be run multiple times
                    upper_bound = max(lower[1], higher[1])
                    list_merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                else:
                    list_merged.append(higher)

        self.all_detected_peaks = list_merged

        # update results dataframe
        # 'chrom', 'chromStart', 'chromEnd', 'name', 'full_name', 'mean_peak_values', 'max_peak_values', 'size_bp', 'sequence'
        for peak_num, peak in enumerate(self.all_detected_peaks):
            chromosome = self.chromosome
            start = int(peak[0])
            stop = int(peak[1])
            name = 'P' + str(peak_num+1)
            full_name = 'PEAK_' + str(peak_num+1) + '_' + str(self.chromosome) + ':' + str(self.region_start) + '-' + str(self.region_stop) #TODO change the region to local coords
            size = peak[1] - peak[0]
            primerf = ''
            primerftemp = ''
            primerr = ''
            primerrtemp = ''
            mean_peak_values = self.get_mean_range_values(start-self.region_start, stop-self.region_start)
            max_peak_values = self.get_max_range_values(start-self.region_start, stop-self.region_start)
            sequence = self.get_sequence(self.genome, chromosome, start, stop)

            #columns=['chrom', 'chromStart', 'chromEnd', 'name', 'full_name', 'mean_peak_values', 'max_peak_values', 'size_bp', 'primerF', 'primerFtemp', 'primerR' 'primerRtemp','sequence']
            self.df_results.loc[len(self.df_results)]=([chromosome, start, stop, name, full_name, mean_peak_values, max_peak_values, size, primerf, primerftemp, primerr, primerrtemp, sequence])

        #filter size and mean peak values
        if final_size is not None:
            self.df_results = self.df_results[self.df_results.size_bp >= final_size]
        if final_mean_peak_values is not None:
            self.df_results = self.df_results[self.df_results.mean_peak_values >= final_mean_peak_values]

        print('Final Passed Merge/Filters:', len(self.df_results))

    def motif_scanner(self, tfactor=None, score_threshold=8, plot=True, fig_width=8, fig_height=4, jaspar_data='JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt'):
        """ Function to scan detected sequences for transcription factor motifs using JASPAR. Multiple transcription factors can be specificied with a plus sign.
        example: 'OTX2+VSX2'
        """

        if len(self.df_results) < 1:
            raise RuntimeError("You run the ecr_scanner prior to running the motif_scanner!")

        self.df_motifs = pd.DataFrame(columns=['name', 'full_name', 'tfactor', 'motif_coords', 'position', 'strand', 'score', 'rel_score'])

        file_handle = open(jaspar_data)
        tf_dict = {}
        list_tfs = []

        coords = ''
        strand = ''

        # biopython motif parser using jaspar format
        for motif in motifs.parse(file_handle, fmt="jaspar"):
            tf_dict[motif.name] = motif
        
        if tfactor is None:
            print("You must select a transcription factor to scan for.")
            print(tf_dict.keys())
        
        for tf in list(str(tfactor).split('+')):
            if tf not in tf_dict.keys():
                print(tf + " was not found in the database.\n")
                print("You must select an available transcription factor to scan for. Both upper and titlecase are valid posibilities.\n")
                print(tf_dict.keys())
        
        else:
            list_tfs = list(str(tfactor).split('+'))
            for index, enhancer in self.df_results.iterrows():

                test_seq=Seq(enhancer['sequence'])

                for tfactor in list_tfs:
                    self.df_motifs.loc[len(self.df_motifs)]=(enhancer['name'], enhancer['full_name'], tfactor, coords, 0, strand, -1, 0) #dummy value to preserve order for ploting

                    pssm = tf_dict[tfactor].pssm

                    footprint = len(pssm[0])
                    start_pos = 0
                    max_score = pssm.max
                    min_score = pssm.min
                    abs_score_threshold = (max_score - min_score) * 0.8 + min_score

                    for position, score in pssm.search(test_seq):
                        rel_score = (score - min_score) / (max_score - min_score)

                        if int(position) >= 0:
                            strand = '+'
                            start_pos = enhancer['chromStart'] + position 
                        if int(position) < 0:
                            strand = '-' 
                            start_pos = enhancer['chromEnd'] + position 

                        coords = enhancer['chrom'] + ':' + str(start_pos) + '-' + str(start_pos + footprint)

                        if score > score_threshold:
                            self.df_motifs.loc[len(self.df_motifs)]=(enhancer['name'], enhancer['full_name'], tfactor, coords, position, strand, float(score), rel_score)
        if plot is True:
            plt.figure(figsize=(fig_width, fig_height))
            for tfactor in list_tfs:
                df_plot = self.df_motifs.loc[self.df_motifs['tfactor'] == tfactor]
                plt.scatter(x=df_plot.name, y=df_plot.score, alpha=0.75)
            plt.title('JASPAR detected ' + str(list_tfs)+' motif(s)')
            plt.ylim(bottom=score_threshold-1)
            plt.legend(list_tfs, loc='center left', bbox_to_anchor=(1, 0.5))
            plt.xticks(rotation=90)
            plt.ylabel('JASPAR score')
            self.df_motifs = self.df_motifs[self.df_motifs['score'] >= 0] # drop negative scores
        else:
            self.df_motifs = self.df_motifs[self.df_motifs['score'] >= 0] # drop negative scores
            print(self.df_motifs)
    
    def plot_detected_enhancers(self, fig_width=20, fig_height=4, ymax='auto'):
        #calculate the x indexes for plotting
        x_index = [self.region_start + i for i in range(self.region_stop-self.region_start)]

        if ymax is 'auto':
            ymax = self.region_max_value

        plt.figure(figsize=(fig_width, fig_height))
        plt.plot(x_index, self.region_values)
        plt.title(self.track_plot_header + ' ' + self.genome + ' ' + self.chromosome + ': ' + str(self.region_start) + ' -> ' + str(self.region_stop-1))
        plt.xlabel('Coordinates')
        plt.ylabel('Peak Values')
        plt.ylim(0, ymax)

        for index, enhancer in self.df_results.iterrows():
            plt.annotate(enhancer['name'], ((enhancer['chromStart'] + enhancer['chromEnd'])/2, ymax*.90), fontsize=9, color='red')
            plt.axvspan(enhancer['chromStart'],enhancer['chromEnd'], 0, ymax, alpha=0.25, color='red')

        plt.xticks([self.region_start, (self.region_start+self.region_stop)/2, self.region_stop])
        plt.tight_layout()

    def plot_detected_mean_peak_values(self, sort=False):
        self.df_results.plot.bar(x='name', y='mean_peak_values', title=self.track, ylabel='Mean Peak Values')

    def plot_detected_max_peak_values(self, sort=False):
        self.df_results.plot.bar(x='name', y='max_peak_values', title=self.track, ylabel='Max Peak Values')

    def plot_detected_size(self, sort=False):
        self.df_results.plot.bar(x='name', y='size_bp', title=self.track, ylabel='Size (bp)')
    
    def plot_detected_motifs(self, fig_width=6, fig_height=4):
        pass
        # plt.figure(figsize=(fig_width, fig_height))
        # for tfactor in list_tfs:
        #     df_plot = self.df_motifs.loc[self.df_motifs['tfactor'] == tfactor]
        #     plt.scatter(x=df_plot.name, y=df_plot.score, alpha=0.75)
        # plt.title('JASPAR detected ' + str(list_tfs)+' motif(s)')
        # plt.ylim(bottom=score_threshold-1)
        # plt.legend(list_tfs, loc='center left', bbox_to_anchor=(1, 0.5))
    
    def plot_custom(self, x, y):
        fig, ax = plt.subplots()
        ax.bar(self.df_results[x], self.df_results[y])
        ax.set_title(self.track)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        self.last_figure = fig
    
    def save_bed(self, file_name):
        df_bed = pd.DataFrame(columns=['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

        #convert df_results into df_bed
        for index, row in self.df_results.iterrows():
            df_bed.loc[len(df_bed)]=([row['chrom'], row['chromStart'], row['chromEnd'], row['name'], row['max_peak_values'], 0])
        df_bed.to_csv(file_name + '.bed', sep='\t', index=None)
        print(file_name + '.bed' + " Saved!")

    def save_genbank(self, file_name):
        sequence_string = self.get_sequence(self.genome, self.chromosome, self.region_start, self.region_stop)

        sequence_object = Seq(sequence_string)

        #create a record
        genbank_record = SeqRecord(sequence_object,
                        id='123456789', # random accession number
                        name=self.chromosome + str(self.region_start) + str(self.region_stop),
                        description='Autogenerated Genbank file annotated with enhancer scanner data')

        genbank_record.annotations['molecule_type']='DNA'

        #annotate peaks
        for index, peak in self.df_results.iterrows():
            feature = SeqFeature(FeatureLocation(start=peak['chromStart']-self.region_start, end=peak['chromEnd']-self.region_start), type='Peak', strand=0)
            feature.qualifiers['label'] = peak['name']
            #feature.qualifiers['ApEinfo_revcolor'] = '#000000' #black
            #feature.qualifiers['ApEinfo_fwdcolor'] = '#000000' #black

            genbank_record.features.append(feature)

        #annotate motifs
        for index, motif in self.df_motifs.iterrows():

            coords = motif['motif_coords']
            coords = coords.split(':')[1]
            start_pos, stop_pos = coords.split('-')

            feature = SeqFeature(FeatureLocation(start=int(start_pos)-self.region_start, end=int(stop_pos)-self.region_start), type='TF', strand=0)
            feature.qualifiers['label'] = motif['tfactor']
            genbank_record.features.append(feature)


        output_file = open(file_name + '.gb', 'w')
        SeqIO.write(genbank_record, output_file, 'genbank')
        print(file_name + '.gb' + " Saved!")

    def save_results(self, file_name):
        self.df_results.to_csv(file_name + '.csv', sep=',', index=None)
        print(file_name + '.csv' + " Saved!")
    
    def save_motifs(self, file_name):
        self.df_motifs.to_csv(file_name + '.csv', sep=',', index=None)
        print(file_name + '.csv' + " Saved!")

    def reset_results(self):
        """ Internal function to reset the results and motif dataframes. """

        self.df_results = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'full_name', 'mean_peak_values', 'max_peak_values', 'size_bp', 'primerF', 'primerFtemp', 'primerR', 'primerRtemp','sequence'])
        self.df_motifs = pd.DataFrame(columns=['name', 'full_name', 'tfactor', 'motif_coords', 'position', 'strand', 'score', 'rel_score'])

    def get_mean_range_values(self, start, stop):
        return self.region_values[start:stop+1].mean()

    def get_max_range_values(self, start, stop):
        print(start, stop)
        print(self.region_values)
        print(self.region_values[start:stop])
        return self.region_values[start:stop+1].max()

        #self.bw.stats(chromosome, region_start, region_stop, type='max')[0]

    #def get_median_peak_values(self, chromosome, region_start, region_stop):
    #    return statistics.median(self.bw.values(chromosome, region_start, region_stop))

    def get_region_values(self):
        region_start = int(self.region_start)
        region_stop = int(self.region_stop)
        chromosome = self.chromosome

        # If multiple track comparison
        if self.track2_bw != None:
            region_values = None

            if self.track_comparison_type == 'subtract':
                region_values = np.subtract(np.array(self.track1_bw.values(chromosome, region_start, region_stop)), np.array(self.track2_bw.values(chromosome, region_start, region_stop)))
            if self.track_comparison_type == 'add':
                region_values = np.add(np.array(self.track1_bw.values(chromosome, region_start, region_stop)), np.array(self.track2_bw.values(chromosome, region_start, region_stop)))
            if self.track_comparison_type == 'divide':
                region_values = np.true_divide(np.array(self.track1_bw.values(chromosome, region_start, region_stop)), np.array(self.track2_bw.values(chromosome, region_start, region_stop)))
            if self.track_comparison_type == 'multiply':
                region_values = np.multiply(np.array(self.track1_bw.values(chromosome, region_start, region_stop)), np.array(self.track2_bw.values(chromosome, region_start, region_stop)))

            self.region_values = np.nan_to_num(region_values, nan=0.0, posinf=100, neginf=0.0).clip(min=0)
            
            self.region_max_value = self.region_values.max()
            self.region_mean_value = self.region_values.mean()


        else: #only one track
        # grab region values to detect

            self.region_values = np.array(self.track1_bw.values(chromosome, region_start, region_stop))
            self.region_max_value = self.region_values.max()
            self.region_mean_value = self.region_values.mean()

            #self.region_max_value = self.track1_bw(self.chromosome, self.region_start, self.region_stop, type='max')[0]
            #self.region_mean_value = self.get_mean_peak_values(self.chromosome, self.region_start, self.region_stop)
            #self.region_median_value = self.get_median_peak_values(self.chromosome, self.region_start, self.region_stop)

        print("Max peak height in this range: ",self.region_max_value)
        print("Mean peak height in this range: ", self.region_mean_value)
        print("Median peak height in this range: ", self.region_median_value)
        print("")
        print("Max peak height across whole track: ", self.track1_max_value, self.track2_max_value)
        print("Minumum peak height across whole track1: ", self.track1_min_value, self.track2_min_value)

    def get_sequence(self, genome, chromosome, start, stop, padding_start=0, padding_stop=0):
        #page = requests.get('http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chr14:48645000,48660000')
        # be careful: the DAS server uses an index of (+1) for the first base.
        
        padding_start = int(padding_start)
        padding_stop = int(padding_stop)
        start = int(start)
        stop = int(stop)
        page = requests.get('http://genome.ucsc.edu/cgi-bin/das/'+genome+'/dna?segment='+chromosome+':'+str(start-padding_start)+','+str(stop+padding_stop))
        return html.fromstring(page.content).text_content().strip().replace('\n', '')

    def ungzip(self, filepath):
        new_filepath = filepath.split('.')[-1]

        with gzip.open(filepath, 'rb') as f_in:
            with open(new_filepath, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    def download_url(self, url):
        filename = url.split('/')[-1]

        if '?' in filename:
            filename = filename.split('?')[0]
        if '=' in filename:
            filename = filename.split('=')[-1]

        filename = filename.replace('%2E', '.')
        filename = filename.replace('%2D', '-')
        filename = filename.replace('%5F', '_')

        #check if file exists already
        if path.exists(filename):
            print("This track already exists!")

        else:
            print('Downloading ', filename, '...')
            r = requests.get(url)

            with open(filename, 'wb') as f:
                f.write(r.content)
            
            if int(r.status_code) == 200:
                print('Download successful!')
            else:
                print(r.status_code)
        
            if filename.endswith('gz'):
                print('Unzipping...')
                self.ungzip(filename)
                print('Done.')
            
# phase these out
    def write_stats(self, file_name):
        print("Warning this is being renamed to .save_results()")
        self.save_results(file_name)

    def write_genbank(self, file_name):
        print("Warning this is being renamed to .save_results()")
        self.save_genbank(file_name)

    def write_bed(self, file_name):
        print("Warning this is being renamed to .save_bed()")
        self.save_bed(file_name)