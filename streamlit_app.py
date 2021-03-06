import streamlit as st
import streamlit.components.v1 as components

import numpy as np
import pandas as pd
from pyEnhancerScan import EnhancerScan
from Bio import motifs
import base64
import os
from io import BytesIO

GENOME_LIST = ['mm10', 'mm9', 'hg38', 'hg19']

# TODO:
# custom track upload
# descriptions
# convert to new pyplot method

st.set_option('deprecation.showPyplotGlobalUse', False)
#st.set_option('server.maxUploadSize', 500)

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

def get_binary_file_downloader_html(bin_file, file_label='File'):
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a>'
    return href

##### SIDEBAR
st.sidebar.title('PyEnhancerScanner Streamlit App')
st.sidebar.write('interactive implementation of a Python API to facilitate the scanning, integrating, and visualizing genomic track data for regions of interest in gene regulation')

st.sidebar.header("Analysis Modes:")
sidebar_result = st.sidebar.radio("Mode", ['Single Track', 'Compare Tracks', 'Analyze BEDfile', 'Genome Browser'])

st.sidebar.write('Source code and documentation: https://github.com/MLKaufman/py-enhancer-scanner')
st.sidebar.write('Author: Michael Kaufman - 2021')

##### SINGLE TRACK ANALYSIS
if sidebar_result == 'Single Track':
    st.title('Single Track Analysis')
    st.write('Select a track and genome to analyze for potential enhancer peaks.')
    st.write('Choose your coordinates to scan a specific region.')
    st.write('Additionally these peaks can be scanned for transcription factor binding sites.')
    st.markdown('<HR>', unsafe_allow_html=True)

    scanner = EnhancerScan()
    st.header('Track and genome to analyze:')

    with st.beta_expander('View built in browser:'):
        components.iframe('https://igv.org/app/', height=500, width=700, scrolling=True)

    col1, col2 = st.beta_columns([3, 1])
    with col1:
        track = st.selectbox('Which track would you like to analyze?', scanner.list_external_tracks()['Track_Name'])
    with col2:
        genome = st.selectbox( 'Pick a genome', GENOME_LIST)

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

    scanner.load_track(genome, track)
    print(user_peak_height)
    scanner.enhancer_scanner(coords, peak_height=user_peak_height)

    scanner.plot_detected_enhancers(fig_height=4, fig_width=10)
    st.pyplot()
    st.write('Track Max Height:', scanner.track1_max_value, 'Track Min Height:', scanner.track1_min_value, 'Auto Peak Threshold:', scanner.auto_threshold)

    scanner.df_results
    if st.button('Download Dataframe as CSV', 1):
        tmp_download_link = download_link(scanner.df_results, 'Results.csv', 'CSV generated! Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

    if st.button('Download Bedfile of Peaks'):
        scanner.save_bed('Bedfile')
        st.markdown(get_binary_file_downloader_html('Bedfile.bed', 'Bedfile generated! Click here to download your data!'), unsafe_allow_html=True)

    st.header('Additional downstream analysis:')

    if st.checkbox('Plot Peaks'):
        y_axis = st.selectbox('Y Axis', ['size_bp', 'mean_peak_values', 'max_peak_values'])
        scanner.plot_custom('name', y_axis)
        st.pyplot()

    if st.checkbox('Transcription Factor Motif Analysis'):
        tf_dict = {}
        file_handle = open('JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt')
        for motif in motifs.parse(file_handle, fmt="jaspar"):
            tf_dict[motif.name] = motif

        col10, col11 = st.beta_columns([3,1])
        with col10:
            option = st.multiselect('Which motif would you like to analyze?', list(tf_dict.keys()), default='OTX2')
        
        with col11:
            threshold = st.number_input('Score Treshold', min_value=0, max_value=20, value=8, step=1)

        multi = ''
        for each in option:
            if multi == '':
                multi = each
            else:
                multi = multi + '+' + each

        scanner.motif_scanner(multi, score_threshold=threshold, plot=True)
        st.pyplot()

        st.write(scanner.df_motifs)
        if st.button('Download Dataframe as CSV', 2):
            tmp_download_link = download_link(scanner.df_motifs, 'Motifs.csv', 'CSV generated! Click here to download your data!')
            st.markdown(tmp_download_link, unsafe_allow_html=True)

        if st.button('Download Genbank'):
            scanner.save_genbank('Genbank')
            st.markdown(get_binary_file_downloader_html('Genbank.gb', 'Genbank generated! Click here to download your data!'), unsafe_allow_html=True)
        
        with st.beta_expander('View motif locations:'):
            list_motifs = scanner.df_motifs.tfactor.unique()
            motif = st.selectbox('Motif:', list_motifs)

            scanner.plot_detected_motifs(motif, fig_width=10, fig_height=4)
            st.pyplot()

##### USER TRACK ANALYSIS
if sidebar_result == 'User Track':
    st.title('User Track Analysis')
    st.write('Description/Instructions')
    st.markdown('<HR>', unsafe_allow_html=True)    

    scanner = EnhancerScan()
    st.header('Track and genome to analyze:')

    with st.beta_expander('View built in browser:'):
        components.iframe('https://igv.org/app/', height=500, width=700, scrolling=True)

    col1, col2 = st.beta_columns([3, 1])
    with col1:
        uploadedfile = st.file_uploader('Upload a track in bigwig (.bw) format', type='bw')
        #print(uploadedfile.name)
        #track = base64.b64encode(uploadedfile.read()).decode()
        track = st.selectbox('Which track would you like to analyze?', scanner.list_external_tracks()['Track_Name'])
    with col2:
        genome = st.selectbox( 'Pick a genome', GENOME_LIST)

    st.header('Genomic coordinates for peak analysis:')

    col1, col2 = st.beta_columns(2)
    with col1:
        coords = st.text_input("Please Enter Coordinates", 'chr14:48,605,306-48,630,475')
    with col2:
        user_peak_height = st.text_input('Peak height threshold [ auto, mean, median, or a value]:','auto')

    scanner.load_track(genome, track)
    print(user_peak_height)
    scanner.enhancer_scanner(coords, peak_height=user_peak_height)

    scanner.plot_detected_enhancers(fig_height=4, fig_width=10)
    st.pyplot()
    st.write('Track Max Height:', scanner.track1_max_value, 'Track Min Height:', scanner.track1_min_value, 'Auto Peak Threshold:', scanner.auto_threshold)

    scanner.df_results
    if st.button('Download Dataframe as CSV', 1):
        tmp_download_link = download_link(scanner.df_results, 'Results.csv', 'CSV generated! Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

    if st.button('Download Bedfile of Peaks'):
        scanner.save_bed('Bedfile')
        st.markdown(get_binary_file_downloader_html('Bedfile.bed', 'Bedfile generated! Click here to download your data!'), unsafe_allow_html=True)

    st.header('Additional downstream analysis:')

    if st.checkbox('Plot Peaks'):
        y_axis = st.selectbox('Y Axis', ['size_bp', 'mean_peak_values', 'max_peak_values'])
        scanner.plot_custom('name', y_axis)
        st.pyplot()

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

        if st.button('Download Genbank'):
            scanner.save_genbank('Genbank')
            st.markdown(get_binary_file_downloader_html('Genbank.gb', 'Genbank generated! Click here to download your data!'), unsafe_allow_html=True)

##### MULTI TRACK ANALYSIS
if sidebar_result == 'Compare Tracks':
    st.title('Track Comparison Analysis')
    st.write('Two tracks and be compared to generate new region values to analyze for enhancer peaks.')
    st.write('Operations include: add, subtract, divide, and multiply')
    st.markdown('<HR>', unsafe_allow_html=True)

    st.header('Track and genome to analyze:')
    scanner = EnhancerScan()

    with st.beta_expander('View built in browser:'):
        components.iframe('https://igv.org/app/', height=500, width=700, scrolling=True)

    col1, col2 = st.beta_columns([3, 1])
    with col1:
        track1 = st.selectbox('Track1:', scanner.list_external_tracks()['Track_Name'])
        operation = st.selectbox('Comparison Operation', ['subtract', 'add', 'divide', 'multiply'])
        track2 = st.selectbox('Track2:', scanner.list_external_tracks()['Track_Name'])
    with col2:
        genome = st.selectbox( 'Pick a genome', GENOME_LIST)

    if track1 not in scanner.list_tracks():
        index = list(scanner.list_external_tracks().index[scanner.list_external_tracks()['Track_Name'] == track1])[0]
        print(track1, index)
        with st.spinner('Downloading track, please wait...'):
            scanner.download_tracks(index)
        st.success('Done!')

    if track2 not in scanner.list_tracks():
        index = list(scanner.list_external_tracks().index[scanner.list_external_tracks()['Track_Name'] == track2])[0]
        print(track2, index)
        with st.spinner('Downloading track, please wait...'):
            scanner.download_tracks(index)
        st.success('Done!')

    st.header('Genomic coordinates for peak analysis:')

    col1, col2 = st.beta_columns(2)
    with col1:
        coords = st.text_input("Please Enter Coordinates", 'chr14:48,605,306-48,630,475')
    with col2:
        user_peak_height = st.text_input('Peak height threshold [ auto, mean, median, or a value]:','auto')

    scanner.load_track(genome, track1, track2, operation)

    scanner.enhancer_scanner(coords, peak_height=user_peak_height)

    scanner.plot_detected_enhancers(fig_height=4, fig_width=10)
    st.pyplot()
    #st.write('Track Max Height:', scanner.track_max_value, 'Track Min Height:', scanner.track_min_value, 'Auto Peak Threshold:', scanner.auto_threshold)

    scanner.df_results
    if st.button('Download Dataframe as CSV', 1):
        tmp_download_link = download_link(scanner.df_results, 'Results.csv', 'CSV generated! Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

    if st.button('Download Bedfile of Peaks'):
        scanner.save_bed('Bedfile')
        st.markdown(get_binary_file_downloader_html('Bedfile.bed', 'Bedfile generated! Click here to download your data!'), unsafe_allow_html=True)

##### BED FILE ANALYSIS
if sidebar_result == 'Analyze BEDfile':
    st.title('BEDfile Analysis')
    st.write('Load a BED file of previously saved results or selected regiosn to scan.')
    st.markdown('<HR>', unsafe_allow_html=True)

    scanner = EnhancerScan()

    col1, col2 = st.beta_columns(2)
    with col1:
        genome = st.selectbox( 'Pick a genome', GENOME_LIST, key='bedgenome')

    # dummy data to allow bedfile
    track = 'GSE102092_Retina_E14.5_ATAC.mm10.bw'
    coords = 'chr14:48,605,306-48,630,475'
    user_peak_height = 'auto'

    #scanner.load_track(track, genome)
    #scanner.enhancer_scanner(coords, peak_height=user_peak_height)

    bedfile = st.file_uploader('BED File')
    try:
        scanner.load_bed(genome, bedfile)
    except ValueError:
        pass

    tf_dict = {}
    file_handle = open('JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt')
    for motif in motifs.parse(file_handle, fmt="jaspar"):
        tf_dict[motif.name] = motif


    col10, col11 = st.beta_columns([3,1])
    with col10:
        option = st.multiselect('Which motif would you like to analyze?', list(tf_dict.keys()), default='OTX2')
    
    with col11:
        threshold = st.number_input('Score Treshold', min_value=0, max_value=20, value=8, step=1)


    multi = ''
    for each in option:
        if multi == '':
            multi = each
        else:
            multi = multi + '+' + each
    try:
        scanner.motif_scanner(multi, score_threshold=threshold, plot=True)
        st.pyplot()
    except RuntimeError:
        pass
        
    st.write(scanner.df_motifs)

    if st.button('Download Dataframe as CSV', 2):
        tmp_download_link = download_link(scanner.df_motifs, 'Motifs.csv', 'CSV generated! Click here to download your data!')
        st.markdown(tmp_download_link, unsafe_allow_html=True)

    if st.button('Download Genbank'):
        scanner.save_genbank('Genbank')
        st.markdown(get_binary_file_downloader_html('Genbank.gb', 'Genbank generated! Click here to download your data!'), unsafe_allow_html=True)

if sidebar_result == 'Genome Browser':
    st.title('Genome Browser')
    
    scanner = EnhancerScan()

    with st.beta_expander('View built in tracks:'):
        st.write(scanner.list_external_tracks())

    components.iframe('https://igv.org/app/', height=500, width=700, scrolling=True)