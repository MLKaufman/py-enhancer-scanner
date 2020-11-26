import pyBigWig
import os
from os import path
import pandas as pd
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

VERSION = 0.4

# Version 1.0 Feature TODO list:
# TODO: Compare tracks option?
# Run via streamlit

class EnhancerScan:
    """ Class to handle scanning of bigwig files. """
    def __init__(self):
        self.chromosome = None
        self.region_start = None
        self.region_stop = None

        self.region_values = None
        self.all_detected_peaks = None

        self.reset_results()

        print("pyEnhancerScanner version " + str(VERSION) + " beta\n")
        print("The following tracks are locally available:")
        for track in self.list_tracks(): print(track)
        print("")
        print("To download additional tracks, use the download_tracks() function.")

    def list_external_tracks(self):
        """ Alias for download"""
        print('External Track List:')
        df_quicklist = pd.read_csv('external_tracks.csv', delimiter=',', header=0)
        pd.set_option('display.max_colwidth', None) # so it doesnt truncate columns
        print(df_quicklist.loc[:, :'Size'])
        print("")
        print("To download one of these tracks, use download_tracks(track_num=X) where X is the track number / index.")
        print("You can specify your own download url by download_tracks(url='X').")

    def load_track(self, track, genome, reset=True):
        if reset is True:
            self.reset_results()

        if track is not None and genome is not None:
            self.bw = pyBigWig.open(track)
            self.genome = genome
            self.track = track

            header = self.bw.header()
            self.track_min_value = header['minVal']
            self.track_max_value = header['maxVal']

        else:
            return(print("Error: You must load a track by specifyinhg a bigwig track and a genome."))

    def download_tracks(self, url = '', track_num = 0):
        
        if url == '' and track_num == 0:
            self.list_external_tracks()
        
        elif url != '':
            if type(url) is int:
                df_quicklist = pd.read_csv('external_tracks.csv', delimiter=',', header=0)

                url = df_quicklist.loc[url, 'URL_Path']
                self.download_url(url)
            else:
                self.download_url(url)

        elif track_num !=0:
            df_quicklist = pd.read_csv('external_tracks.csv', delimiter=',', header=0)

            url = df_quicklist.loc[track_num, 'URL_Path']
            self.download_url(url)

    def download_jaspar(self, url = ''):
        print(type(url))
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
        bedfiles = []
        for bedfiles in os.listdir(os.getcwd()):
            if bedfiles.endswith(".bed"):
                if bedfiles !='None':
                    bedfiles.append(bedfiles)
        return bedfiles
    
    def load_bed(self, bed_file, chromosome):
        """ Load an existing bed file for analyzing regions enhancers. """
        # load a bed file
        # for each entry in bed file, create an entry in df_results.
        # calculate 
        self.reset_results()
        self.chromosome = chromosome

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
            start = enhancer['chromStart']
            stop = enhancer['chromEnd']
            name = enhancer['name']
            full_name = enhancer['name']
            size = enhancer['chromEnd'] - enhancer['chromStart']
            mean_peak_values = self.get_mean_peak_values(chromosome, start, stop)
            max_peak_values = self.get_max_peak_values(chromosome, start, stop)
            sequence = self.get_sequence(self.genome, chromosome, start, stop)
            
            list_region_start.append(start)
            list_region_stop.append(stop)

            self.df_results.loc[len(self.df_results)]=([chromosome, start, stop, name, full_name, mean_peak_values, max_peak_values, size, sequence])

        #need to set self.chromosome, self.region_start, self.region_stop?
        self.region_start = sorted(list_region_start)[0]
        self.region_stop = sorted(list_region_stop)[-1]
        self.region_values = self.bw.values(chromosome, self.region_start, self.region_stop)

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

        self.region_max_value = self.bw.stats(self.chromosome, self.region_start, self.region_stop, type='max')[0]
        self.region_mean_value = self.get_mean_peak_values(self.chromosome, self.region_start, self.region_stop)
        self.region_median_value = self.get_median_peak_values(self.chromosome, self.region_start, self.region_stop)

        print("Max peak height in this range: ",self.region_max_value)
        print("Mean peak height in this range: ", self.region_mean_value)
        print("Median peak height in this range: ", self.region_median_value)
        print("")
        print("Max peak height across whole track: ", self.track_max_value)
        print("Minumum peak height across whole track: ", self.track_min_value)

        if peak_height is 'auto':
            peak_height = self.region_max_value * .25
            print("")
            print("Using auto detection of peak height: ", peak_height)

        if peak_height is 'mean':
            peak_height = self.region_mean_value

        if peak_height is 'median':
            peak_height = self.region_median_value

        print("")

        # PEAK DETECTION

        # grab region values to detect
        self.region_values = self.bw.values(self.chromosome, self.region_start, self.region_stop)

        # detect peaks, returns a tuple of an array of peaks and a dictionary or properties
        peaks = signal.find_peaks(self.region_values, width = peak_width, height = peak_height)

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
            full_name = 'PEAK_' + str(peak_num+1) + '_' + str(self.chromosome) + str(self.region_start) + '-' + str(self.region_stop)
            size = peak[1] - peak[0]
            mean_peak_values = self.get_mean_peak_values(chromosome, start, stop)
            max_peak_values = self.get_max_peak_values(chromosome, start, stop)
            sequence = self.get_sequence(self.genome, chromosome, start, stop)

            self.df_results.loc[len(self.df_results)]=([chromosome, start, stop, name, full_name, mean_peak_values, max_peak_values, size, sequence])

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

        df_motifs = pd.DataFrame(columns=['name', 'tf_factor', 'position', 'score', 'rel_score'])

        file_handle = open(jaspar_data)
        tf_dict = {}
        list_tfs = []

        # biopython motif parser using jaspar format
        for motif in motifs.parse(file_handle, fmt="jaspar"):
            tf_dict[motif.name] = motif
        
        if tfactor is None:
            print("You must select a transcription factor to scan for.")
            print(tf_dict.keys())
        
        else:
            list_tfs = list(str(tfactor).split('+'))
            for index, enhancer in self.df_results.iterrows():
                df_motifs.loc[len(df_motifs)]=(enhancer['name'], tfactor, 0, -1, 0)

                test_seq=Seq(enhancer['sequence'])

                for tfactor in list_tfs:
                    pssm = tf_dict[tfactor].pssm
                    max_score = pssm.max
                    min_score = pssm.min
                    abs_score_threshold = (max_score - min_score) * 0.8 + min_score
                    #print("abs score threshold: ", abs_score_threshold)

                    for position, score in pssm.search(test_seq):
                        rel_score = (score - min_score) / (max_score - min_score)
                        #print(enhancer['name'], position, score, rel_score)
                        
                        if score > score_threshold:
                            df_motifs.loc[len(df_motifs)]=(enhancer['name'], tfactor, position, float(score), rel_score)
        if plot is True:
            plt.figure(figsize=(fig_width, fig_height))
            for tfactor in list_tfs:
                df_plot = df_motifs.loc[df_motifs['tf_factor'] == tfactor]
                plt.scatter(x=df_plot.name, y=df_plot.score, alpha=0.75)
            plt.title('JASPAR detected ' + str(list_tfs)+' motif(s)')
            plt.ylim(bottom=score_threshold-1)
            plt.legend(list_tfs, loc='center left', bbox_to_anchor=(1, 0.5))        
        else:
            df_motifs = df_motifs[df_motifs['score'] >= 0] # drop negative scores
            print(df_motifs)
    
    def plot_detected_enhancers(self, fig_width=20, fig_height=4, ymax='auto'):
        #calculate the x indexes for plotting
        x_index = [self.region_start + i for i in range(self.region_stop-self.region_start)]

        if ymax is 'auto':
            ymax = self.region_max_value

        plt.figure(figsize=(fig_width, fig_height))
        plt.plot(x_index, self.region_values)
        plt.title(self.track + ' ' + self.genome + ' ' + self.chromosome + ': ' + str(self.region_start) + ' -> ' + str(self.region_stop-1))
        plt.xlabel('Coordinates')
        plt.ylabel('Peak Values')
        plt.ylim(0, ymax)

        for index, enhancer in self.df_results.iterrows():
            plt.annotate(enhancer['name'], ((enhancer['chromStart'] + enhancer['chromEnd'])/2, ymax*.90), fontsize=9, color='red')
            plt.axvspan(enhancer['chromStart'],enhancer['chromEnd'], 0, ymax, alpha=0.25, color='red')

        plt.xticks([self.region_start, (self.region_start+self.region_stop)/2, self.region_stop])
        plt.tight_layout()

    def plot_detected_mean_peak_values(self, sort=False):
        self.df_results.plot.bar(x='name', y='mean_peak_values')

    def plot_detected_max_peak_values(self, sort=False):
        self.df_results.plot.bar(x='name', y='max_peak_values')

    def plot_detected_size(self, sort=False):
        self.df_results.plot.bar(x='name', y='size_bp')
    
    def plot_detected_motifs(self):
        #TODO: plotting of detected motifs
        pass
    
    def plot_custom(self, x, y):
        self.df_results.plot.bar(x, y)

    def write_bed(self, file_name, ecr_padding=0):
        #TODO: make this use df_results instead of self.all_detected_peaks
        df_bed = pd.DataFrame(columns=['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand'])

        for peak_num, peak in enumerate(self.all_detected_peaks):
            df_bed.loc[len(df_bed)]=([self.chromosome, self.region_start+peak[0]-ecr_padding, self.region_start+peak[1]+ecr_padding, 'PEAK_' + str(peak_num+1), 0, '.'])

            df_bed.to_csv(file_name + '.bed', sep='\t', index=None)
    
    def write_stats(self, file_name):
        self.df_results.to_csv(file_name + '.csv', sep=',', index=None)

    def write_genbank(self, file_name):
        #TODO: create annotations for each primer
        #TODO: create annotations for TFBS

        sequence_string = self.get_sequence(self.genome, self.chromosome, self.region_start, self.region_stop)

        sequence_object = Seq(sequence_string)

        #create a record
        genbank_record = SeqRecord(sequence_object,
                        id='123456789', # random accession number
                        name=self.chromosome + str(self.region_start) + str(self.region_stop),
                        description='Autogenerated Genbank file annotated with enhancer scanner data')

        genbank_record.annotations['molecule_type']='DNA'

        for index, peak in self.df_results.iterrows():
            feature = SeqFeature(FeatureLocation(start=peak['chromStart']-self.region_start, end=peak['chromEnd']-self.region_start), type='Peak', strand=0)
            feature.qualifiers['label'] = peak['name']
            #feature.qualifiers['ApEinfo_revcolor'] = '#000000' #black
            #feature.qualifiers['ApEinfo_fwdcolor'] = '#000000' #black

            genbank_record.features.append(feature)

        output_file = open(file_name + '.gb', 'w')
        SeqIO.write(genbank_record, output_file, 'genbank')

    def reset_results(self):
        self.df_results = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'full_name', 'mean_peak_values', 'max_peak_values', 'size_bp', 'sequence'])

    def get_mean_peak_values(self, chromosome, region_start, region_stop):
        return self.bw.stats(chromosome, region_start, region_stop)[0]

    def get_max_peak_values(self, chromosome, region_start, region_stop):
        return self.bw.stats(chromosome, region_start, region_stop, type='max')[0]

    def get_median_peak_values(self, chromosome, region_start, region_stop):
        return statistics.median(self.bw.values(chromosome, region_start, region_stop))

    def get_sequence(self, genome, chromosome, start, stop, padding_start=0, padding_stop=0):
        #page = requests.get('http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment=chr14:48645000,48660000')
        # be careful: the DAS server uses an index of (+1) for the first base.
        
        padding_start = int(padding_start)
        padding_stop = int(padding_stop)
        start = int(start)
        stop = int(stop)
        page = requests.get('http://genome.ucsc.edu/cgi-bin/das/'+genome+'/dna?segment='+chromosome+':'+str(start-padding_start)+','+str(stop+padding_stop))
        return html.fromstring(page.content).text_content().strip().replace('\n', '')

    def get_primers(self, file_name, padding=0, primer_walk=False, infusionf=None, infusionr=None, verbose=False):
        df_primers = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd', 'name', 'full_name', 'mean_peak_values', 'max_peak_values', 'size_bp', 'primer_fwd', 'fwd_tm', 'primer_rev', 'rev_tm', 'product_size','sequence', 'padded_sequence']
        )

        success = 0
        total_count = 0

        for index, enhancer in self.df_results.iterrows():
            total_count+=1
            padded_sequence = 'NA'

            chromosome = enhancer['chrom']
            start = enhancer['chromStart']
            stop = enhancer['chromEnd']
            name = enhancer['name']
            full_name = enhancer['full_name']
            mean_peak_values = enhancer['mean_peak_values']
            max_peak_values = enhancer['max_peak_values']
            size = enhancer['size_bp']

            #primer magic
            primer_design = Design()
            template = self.get_sequence(self.genome, chromosome, start, stop, padding)

            primer_design.settings.template(template)
            primer_design.settings.as_cloning_task()
            primer_design.settings.primer_num_return(1)
            primer_design.settings.product_size((size, 10000))
            results, explain = primer_design.run()

            if verbose is True:
                print(results)
                print(explain)
                #print(dict(primer_design.params))

            if len(results.values()) > 0:
                success+=1
                print(name + " = success!")
                primer_fwd = results[0]['LEFT']['SEQUENCE']
                fwd_tm = results[0]['LEFT']['TM']
                primer_rev = results[0]['RIGHT']['SEQUENCE']
                rev_tm = results[0]['RIGHT']['TM']
                product_size = results[0]['PAIR']['PRODUCT_SIZE']

            else: #TODO: padding walking to get primers made
                primer_fwd = "not found"
                fwd_tm = "not found"
                primer_rev = "not found"
                rev_tm = "not found"
                product_size = "none"

            if infusionf is not None and infusionr is not None and primer_fwd != "not found" and primer_rev != "not found":
                primer_fwd = str(infusionf) + primer_fwd
                primer_rev = str(infusionr) + primer_rev

            sequence = enhancer['sequence']
            padded_sequence = template

            df_primers.loc[len(df_primers)]=([chromosome, start, stop, name, full_name, mean_peak_values, max_peak_values, size, primer_fwd, fwd_tm, primer_rev, rev_tm, product_size, sequence, padded_sequence])
        print("Success rate= "+str(success)+" out of "+str(total_count))
        df_primers.to_csv(file_name+'.csv', sep=',', index=None)

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
            
