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

## Usage
The prefered method to interact with pyEnhancerScanner is through a Jupyter notebook interface. Functionality of pyEnhancerScanner is through the creation of a scanobject from the EnhancerScanner class and choosing a data track and genome to analyze: 
<&nbsp>
`from pyEnhancerScanner import EnhancerScanner`

`scanobject = EnhancerScanner('mm10.60way.phastCons60wayPlacental.bw', 'mm10')`  
The track must be in the format of bigwig, but can be any sort of conservation, ATAC, etc type data.

Following the creation of the scanobject. Calling the `enhancerscanner()` function will initiate the scan and peak finding.
`scanobject.enhancer_scanner('chr14', 48641226, 48657194, peak_width=10, peak_height=0.5, merge_distance=500, final_size=400, final_mean_conservation=0)`  
