# py-enhancer-scanner
python implementation of an API to facilitate the scanning, integrating, and visualizing genomic track data for regions of interest in gene regulation
---
![Demo Peaks](/images/demopeaks.png)

---

## Requirements
- python >= 3.5
- pyBigWig==0.3.17
- biopython==0.0.6
- matplotlib==3.3.2
- requests==2.24.0
- pandas==1.1.3
- lxml==4.6.1
- scipy==1.5.3
- numpy==1.19.2
- Bio==0.0.6
- primer3plus==1.0.8
- jupyter (optional)


---

## Features
- Analyze genome tracks for potential enhancer sequences
- Uses bigwig formated tracks of Conservation, ATAC, etc
- The `scipy.signal` suite of functions are used to perform singal analysis and peak calling
- Motif scanning for enrichment of transcription factors binding motifs
- Auto creates primers for amplification of identified region products with optional Gibson Assembly overhangs
- Generates plots describing the identified sequences
- Output of bed file for anotations
- Output all results to CSV file  
- TODO: Output of sequences and annotations to genbank file for import into sequence analysis tools
- TODO: Remote access of bigwig files
---
## Installation

1. Download or clone this git repository.  
2. Use your desired virtual environment solution. i.e. conda or pipx etc..  
For conda environment I currently recommend installing with python=3.7.  
3. Activte your evironment, navigate to the cloned directory and run:  
`pip install -r requirements.txt`

---
## Basic Usage
The prefered method to interact with pyEnhancerScanner is through a Jupyter notebook interface. An example of this:  
[ExampleEnhancerScanNotebook.ipynb](https://github.com/MLKaufman/py-enhancer-scanner/blob/main/ExampleEnhancerScanNotebook.ipynb)

Functionality of pyEnhancerScanner is through the creation of a scanobject from the EnhancerScanner class and choosing a data track and genome to analyze:  
`from pyEnhancerScanner import EnhancerScanner`

`scanobject = EnhancerScanner()`  

Next you load a track you want to use in the analysis:  
`scanobject.load_track('bigwigfile.bw', 'mm10)`  

The track must be in the format of bigwig, but can be any sort of conservation, ATAC, etc type data. You will need to supply your own bigwig file which can be obtained from the UCSC Genome browser or from various NCBI Geo datasets. You also supply the reference genome in order to properly grab the sequences.  

Following the creation of the scanobject. Calling the `enhancerscanner()` function will initiate the scan and peak finding.  

Example:  
`scanobject.enhancer_scanner('chr14', 48641226, 48657194, peak_width=10, peak_height=0.5, merge_distance=500, final_size=400, final_mean_conservation=0)`  

Note: `peak_width` and `peak_height` are the most basic variables to adjust to acquire peaks / regions of interest that fit your requirements. These must be set according to the height of the bigwig data track loaded in and minimum width of the desired peak at this height. See the function calls to see the definitions for the advanced settings.  

You can then plot these detected peaks using:  
`scanobject.plot_detected_enhancers()`  

Additional plots:  
`scanobject.plot_detected_size()  scanobject.plot_detected_mean_peak_values()  scanobject.plot_detected_max_peak_values()`  

Advanced:  
`scanobject.motif_scanner()` # detect transcription factor binding site predicitions  
`scanobject.get_primers()` # automatically generate primers for amplification of detected regions  
`scanobject.write_bed()` # write bed file of detected peak locations  
`scanobject.write_stats() ` # write results to CSV file
