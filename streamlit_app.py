import streamlit as st
import streamlit.components.v1 as components

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pyEnhancerScan import EnhancerScan
from Bio import motifs
import base64
from io import BytesIO

# TODO:
# custom track upload 

st.set_option('deprecation.showPyplotGlobalUse', False)

def download_link(object_to_download, download_filename, download_link_text):
    """
    Generates a link to download the given object_to_download.

    object_to_download (str, pd.DataFrame):  The object to be downloaded.
    download_filename (str): filename and extension of file. e.g. mydata.csv, some_txt_output.txt
    download_link_text (str): Text to display for download link.

    Examples:
    download_link(YOUR_DF, 'YOUR_DF.csv', 'Click here to download data!')
    download_link(YOUR_STRING, 'YOUR_STRING.txt', 'Click here to download your text!')

    """
    if isinstance(object_to_download,pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    # some strings <-> bytes conversions necessary here
    b64 = base64.b64encode(object_to_download.encode()).decode()

    return f'<a href="data:file/txt;base64,{b64}" download="{download_filename}">{download_link_text}</a>'

st.title('PyEnhancerScanner Streamlit App')
scanner = EnhancerScan()
st.write('interactive implementation of a Python API to facilitate the scanning, integrating, and visualizing genomic track data for regions of interest in gene regulation')
st.write('Source code and documentation: https://github.com/MLKaufman/py-enhancer-scanner')
st.write('Author: Michael Kaufman - 2021')

st.header('Track and genome to analyze:')

with st.beta_expander('View built in browser:'):
    #components.iframe('https://genome.ucsc.edu/cgi-bin/hgTracks?db=mm10', height=400, width=800, scrolling=True)
    components.iframe('https://igv.org/app/', height=500, width=700, scrolling=True)

col1, col2 = st.beta_columns([3, 1])
with col1:
    track = st.selectbox('Which track would you like to analyze?', scanner.list_external_tracks()['Track_Name'])
with col2:
    genome = st.selectbox( 'Pick a genome', ['mm10'])

if track not in scanner.list_tracks():
    index = list(scanner.list_external_tracks().index[scanner.list_external_tracks()['Track_Name'] == track])[0]
    print(track, index)
    with st.spinner('Downloading track, please wait...'):
        scanner.download_tracks(index)
    st.success('Done!')

st.header('Genomic coordinates for peak analysis:')

col1, col2 = st.beta_columns(2)
with col1:
    coords = st.text_input("Please Enter Coordinates", 'chr14:48,605,306-48,630,475')
with col2:
    user_peak_height = st.text_input('Peak height threshold [ auto, mean, median, or a value]:','auto')

scanner.load_track(track, genome)
print(user_peak_height)
scanner.enhancer_scanner(coords, peak_height=user_peak_height)

scanner.plot_detected_enhancers(fig_height=4, fig_width=10)
st.pyplot()
st.write('Track Max Height:', scanner.track_max_value, 'Track Min Height:', scanner.track_min_value, 'Auto Peak Threshold:', scanner.auto_threshold)

scanner.df_results
if st.button('Download Dataframe as CSV', 1):
    tmp_download_link = download_link(scanner.df_results, 'Results.csv', 'CSV generated! Click here to download your data!')
    st.markdown(tmp_download_link, unsafe_allow_html=True)

st.header('Additional downstream analysis:')
if st.checkbox('Plot Sizes'):
    scanner.plot_detected_size()
    st.pyplot()

if st.checkbox('Plot Max Peak Values'):
    scanner.plot_detected_max_peak_values()
    st.pyplot()

if st.checkbox('Load BED file'):
    col1, col2 = st.beta_columns(2)
    with col1:
        genome = st.selectbox( 'Pick a genome', ['mm10'], key='bedgenome')
    with col2:
        chrom = st.text_input('Chromosome:', 'chr14')

    bedfile = st.file_uploader('BED File')
    try:
        scanner.load_bed(bedfile, 'chr14')
    except ValueError:
        pass

if st.checkbox('Transcription Factor Motif Analysis'):
    tf_dict = {}
    file_handle = open('JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt')
    for motif in motifs.parse(file_handle, fmt="jaspar"):
        tf_dict[motif.name] = motif

    option = st.multiselect('Which motif would you like to analyze?', list(tf_dict.keys()), default='OTX2')
    multi = ''
    for each in option:
        if multi == '':
            multi = each
        else:
            multi = multi + '+' + each

    scanner.motif_scanner(multi, plot=True)
    st.pyplot()

    st.write(scanner.df_motifs)
    if st.button('Download Dataframe as CSV', 2):
        tmp_download_link = download_link(scanner.df_motifs, 'Motifs.csv', 'CSV generated! Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)