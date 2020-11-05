# py-enhancer-scanner
python implementation of an API to facilitate the scanning, integrating, and visualizing genomic track data for regions of interest in gene regulation
---

## Requirements
- python >= 3.5
- pyBigWig
- pandas
- matplotlib
- scipy
- primer3plus
- lxml
- jupyter (optional)

---

## Features
- Analyze genome tracks for potential enhancer sequences
- Uses bigwig formated tracks of Conservation, ATAC, etc
- The `scipy.signal` suite of functions are used to perform singal analysis and peak calling
- TODO: Motif scanning for enrichment of transcription factors binding motifs
- Auto creates primers for amplification of identified region products with optional Gibson Assembly overhangs
- Generates plots describing the identified sequences
- Output of bed file for anotations
- Output all results to CSV file
- Output of sequences and annotations to genbank file for import into sequence analysis tools
---
## Installation

Download or clone this git repository.  
Use your desired virtual environment solution. i.e. conda or pipx etc..  
For conda environment I currently recommend installing with python=3.7.  
Navigate to the cloned directory and run:  
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

The track must be in the format of bigwig, but can be any sort of conservation, ATAC, etc type data. You also supply the reference genome in order to properly grab the sequences.

Following the creation of the scanobject. Calling the `enhancerscanner()` function will initiate the scan and peak finding.  

Example:  
`scanobject.enhancer_scanner('chr14', 48641226, 48657194, peak_width=10, peak_height=0.5, merge_distance=500, final_size=400, final_mean_conservation=0)`  

Note: `peak_width` and `peak_height` are the most basic variables to adjust to acquire peaks / regions of interest that fit your requirements. These must be set according to the height of the bigwig data track loaded in and minimum width of the desired peak at this height. See the function calls to see the definitions for the advanced settings.  