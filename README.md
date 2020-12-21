# py-enhancer-scanner

Run with Binder:  
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/MLKaufman/py-enhancer-scanner/HEAD)  
[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/mlkaufman/py-enhancer-scanner/main)  

python implementation of an API to facilitate the scanning, integrating, and visualizing genomic track data for regions of interest in gene regulation
---
![Demo Peaks](/demopeaks.png)

---

## Features
- Analyze genome tracks for potential enhancer sequences
- Uses bigwig formated tracks of Conservation, ATAC, etc
- Download bigwig files from an external database
- The `scipy.signal` suite of functions are used to perform singal analysis and peak calling
- Motif scanning for enrichment of transcription factors binding motifs
- Auto creates primers for amplification of identified region products with optional Gibson Assembly overhangs
- Generates plots describing the identified sequences
- Output of bed file for anotations
- Output all results to CSV file  
- Output of sequences and annotations to genbank file for import into sequence analysis tools
  
---

## Installation
 
1. Download or clone this git repository.  
2. Use your desired virtual environment solution. i.e. conda or pipx etc..  
For conda environment I currently recommend installing with python=3.7. As bioconda is not supporting 3.8 yet.  
3. Activte your evironment, navigate to the cloned directory and run:  
`pip install -r requirements.txt`  
If primer3-py fails on your system, you can try `conda install primer3-py -c bioconda`  
4. (Optional but prefered) `pip install jupyterlab`

---

## Requirements
- python >= 3.5 and < 3.8
- pyBigWig==0.3.17
- biopython==1.78
- matplotlib==3.3.2
- requests==2.24.0
- pandas==1.1.3
- lxml==4.6.1
- scipy==1.5.3
- numpy==1.19.2
- primer3plus==1.0.8
- jupyterlab (optional)

---
## Basic Usage
The prefered method to interact with pyEnhancerScanner is through a Jupyter notebook interface. An example of this:  
[ExampleEnhancerScanNotebook.ipynb](https://github.com/MLKaufman/py-enhancer-scanner/blob/main/ExampleEnhancerScanNotebook.ipynb)

Functionality of pyEnhancerScanner is through the creation of a scanobject from the EnhancerScanner class and choosing a data track and genome to analyze:  
`from pyEnhancerScan import EnhancerScan`

`scanner = EnhancerScan()`  

Next you load a track you want to use in the analysis:  
`scanner.load_track('bigwigfile.bw', 'mm10)`  
If you require downloading of new bigwig tracks. A list of external tracks can be displayed by using `scanner.download_tracks()`.  
You can select one from this list and use `scanner.download_tracks(X)`, where X is a track number from the list or a url to automatically download.

The track must be in the format of bigwig, but can be any sort of conservation, ATAC, etc type data. You will need to supply your own bigwig file which can be obtained from the UCSC Genome browser or from various NCBI Geo datasets. You must also supply the reference genome name in order to properly grab the DNA sequences.  

Following the creation of the scanobject. Calling the `enhancer_scanner()` function will initiate the scan and peak finding. This requires a chromosome and genome coordinates based on the genome and track loaded above.  

Minimal Example:  
`scanner.enhancer_scanner('chr14:48641226-48657194')`  
or  
`scanner.enhancer_scanner('chr14', '48641226', '48657194')`  
Full Example:  
`scanner.enhancer_scanner('chr14', 48641226, 48657194, peak_width=10, peak_height=0.5, merge_distance=500, final_size=400, final_mean_conservation=0)`  

Note: `peak_width` and `peak_height` are the most basic variables to adjust to acquire peaks / regions of interest that fit your requirements. These must be set according to the height of the bigwig data track loaded in and minimum width of the desired peak at this height. See the function calls to see the definitions for the advanced settings.  

You can then plot these detected peaks using:  
`scanner.plot_detected_enhancers()`  

Additional plots:  
`scanner.plot_detected_size()  scanner.plot_detected_mean_peak_values()  scanner.plot_detected_max_peak_values()`  `scanner.plot_custom()`

PyEnhancerScanner can also use existing genome annotations and detected peaks from BED files to use in downstream analysis instead of de novo scanning of a region.  

Load BED files:  
`scanner.list_bed()` # prints list of locally available bed files for loading  
`scanner.load_bed()` # loads a local bed file  


Transcription Factor Motif Scanner:  
`scanner.download_jaspar()` # downloads the 2020 version of vertebrate position weight matrices if needed  
Note that when selecting a transcription factor motif that the JASPAR database has named them according to available data. i.e. All CAPS for human and Title for mouse, etc.  
`scanner.motif_scanner()` # detect transcription factor binding site predicitions  
Use `scanner.motif_scanner('OTX2+VSX2', plot=False)` for viewing the dataframe containing the list of motifs and scores.

Primer Design:  
`scanner.get_primers()` # automatically generate primers for amplification of detected regions  

Saved Output:  
`scanner.save_bed()` # save bed file of detected peak locations  
`scanner.save_motifs()` # save csv file with motifs dataframe  
`scanner.save_results() ` # save results to CSV file  
`scanner.save_genbank()` # outputs an annotated genbank file for importing into DNA sequence programs (Benchling, Snapgene, APe, etc.)  