import pyBigWig
import os
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt

from lxml import html
import requests

from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from Bio.Alphabet import IUPAC


from primer3plus import Design

VERSION = 0.3

class EnhancerScan:
    """ Class to handle scanning of bigwig files. """
    def __init__(self):
        self.chromosome = None
        self.region_start = None
        self.region_stop = None

        self.region_values = None
        self.all_detected_peaks = None

        self.reset_results()

        print("pyEnhancerScanner version " + str(VERSION) + " beta")
        print("The following tracks are available:\n")
        print(self.list_tracks(), "\n")
        print("Ability to load remote files: ", pyBigWig.remote)

    def load_track(self, track, genome):
        if track is not None and genome is not None:
            self.bw = pyBigWig.open(track)
            self.genome = genome

            header = self.bw.header()
            self.track_min_value = header['minVal']
            self.track_max_value = header['maxVal']

        else:
            return(print("Error: You must load a track by specifyinhg a bigwig track and a genome."))

    def list_tracks(self):
        for bigwig in os.listdir(os.getcwd()):
            if bigwig.endswith(".bw"):
                if bigwig !='None':
                    print(bigwig)

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

    def enhancer_scanner(self, chromosome, region_start, region_stop, peak_width=50, peak_height=.5, merge_distance=200, final_size=None, final_mean_peak_values=None):

        self.chromosome = chromosome
        self.region_start = region_start
        self.region_stop = region_stop

        print("Max peak height in this range: ",self.bw.stats(chromosome, region_start, region_stop, type='max')[0])
        self.region_max_value = self.bw.stats(chromosome, region_start, region_stop, type='max')[0]
        print("Max peak height in whole track: ", self.track_max_value)
        print("Minumum peak height in dataset: ", self.track_min_value)
        print("")

        # grab region values to detect
        self.region_values = self.bw.values(chromosome, region_start, region_stop)

        # detect peaks, returns a tuple of an array of peaks and a dictionary or properties

        peaks = signal.find_peaks(self.region_values, width = peak_width, height = peak_height)

        # make list of unqiue peak widths as tuples and sort
        list_widths = sorted(list(set(zip(list(peaks[1]['left_bases'] + region_start), list(peaks[1]['right_bases'] + region_start))))) # its fixed for final location from widths now
        print('Total Peaks Detected:', len(list_widths))
        print(list_widths)
        
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
        if len(self.df_results) < 1:
            raise RuntimeError("You run the ecr_scanner prior to running the motif_scanner!")

        df_motifs = pd.DataFrame(columns=['name', tfactor, 'position', 'score', 'rel_score'])

        file_handle = open(jaspar_data)
        tf_dict = {}
        # biopython motif parser using jaspar format
        for motif in motifs.parse(file_handle, "jaspar"):
            tf_dict[motif.name] = motif
        
        if tfactor is None:
            print("You must select a transcription factor to scan for.")
            print(tf_dict.keys())
        
        else:
            for index, enhancer in self.df_results.iterrows():
                df_motifs.loc[len(df_motifs)]=(enhancer['name'], tfactor, 0, -1, 0)

                test_seq=Seq(enhancer['sequence'])

                pssm = tf_dict[tfactor].pssm
                max_score = pssm.max
                min_score = pssm.min
                abs_score_threshold = (max_score - min_score) * 0.8 + min_score
                print("abs score threshold: ", abs_score_threshold)

                for position, score in pssm.search(test_seq):
                    rel_score = (score - min_score) / (max_score - min_score)
                    #print(enhancer['name'], position, score, rel_score)
                    
                    if score > score_threshold:
                        df_motifs.loc[len(df_motifs)]=(enhancer['name'], tfactor, position, float(score), rel_score)

        if plot is True:
            plt.figure(figsize=(fig_width, fig_height))
            df_motifs.plot.scatter(x='name', y='score', ylim=score_threshold-1, title='JASPAR detected ' + tfactor+' motif(s)', rot=90)
        
        else:
            df_motifs = df_motifs[df_motifs['score'] >= 0]
            print(df_motifs)
        
        #drop -1 scores before printing?
    
    def plot_detected_enhancers(self, fig_width=20, fig_height=4, ymax='auto'):
        #calculate the x indexes for plotting
        x_index = [self.region_start + i for i in range(self.region_stop-self.region_start)]

        if ymax is 'auto':
            ymax = self.region_max_value

        plt.figure(figsize=(fig_width, fig_height))
        plt.plot(x_index, self.region_values)
        plt.title(self.chromosome + ': ' + str(self.region_start) + ' -> ' + str(self.region_stop-1))
        plt.xlabel('Coordinates')
        plt.ylabel('Peak Values')
        plt.ylim(0, ymax)

        for index, enhancer in self.df_results.iterrows():
            plt.annotate(enhancer['name'], ((enhancer['chromStart'] + enhancer['chromEnd'])/2, ymax), fontsize=14, color='red')

            #plt.axhline(y=.5, xmin=self.region_start+enhancer['chromStart'], xmax=self.region_start + enhancer['chromEnd'], linestyle='solid', linewidth=5, color='red')
            #plt.plot((self.region_start+enhancer['chromStart'], 0.5), (self.region_start + enhancer['chromEnd'], 0.5), '-')

            #TODO: better annotate with line matched to actual coords and not text
            plt.annotate('[', (enhancer['chromStart'], 0), fontsize = 18, color='red')
            plt.annotate(']', (enhancer['chromEnd'], 0), fontsize = 18, color='red')

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
