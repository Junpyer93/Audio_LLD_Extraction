# Audio_LLD_Extraction
Matlab script for LLDs extraction from an audio file

The repository contains two .m scripts: one for extracting low-level and psychoacoustic descriptors from an audio signal (main.m), and the second for extracting the same descriptors from a dataset of audio signals (mainDataset.m). The scripts work for both mono and stereo audio signals, regardless of sampling frequency or audio encoding. The extracted descriptors comply with the MPEG-7 standard, but also include non-standard descriptors and psychoacoustic descriptors.

# main.m
Upon startup, the main.m file asks for the input file name, including the path, the time interval to analyze, and the option to display partial and total plots.

# mainDataset.m
Upon startup, the mainDataset.m file asks for the folder path, the seconds number to analyze, and the option to display partial and total plots.

# LLDsfeatureStereo.m & LLDsfeatureMono.m
The script allows the analysis of the following MPEG-7 standardized low-level descriptors: AudioWaveform, AudioPower, Silence Segment, Audio Spectrum Envelope, Audio Spectrum Centroid, Audio Spectrum Spread, Audio Spectrum Flatness, Harmonic Ratio, Upper Limit of Harmonicity, Audio Fundamental Frequency, Log Attack Time, Temporal Centroid, Harmonic Spectral Centroid, Harmonic Spectral Deviation, Harmonic Spectral Spread, Harmonic Spectral Variation, Spectral Centroid, Audio Spectrum Basis, Audio Spectrum Projection.

# FeaturePlotLowLevelStereo.m & FeaturePlotLowLevelMono.m
The following scripts allow visualizing the plots of the corresponding descriptors.

# PsychoFeature.m
The script allows the analysis of some psychoacustical descriptors: Loudness, Brightness, Roughness, Sharpness, Fluctuation Strength.

# OtherFeature.m
The script allows the analysis of the following no-MPEG-7 standardized low-level descriptors: Zero Crossing Rate, Spectral Rolloff Frequency, Spectral Flux, MFCC

# License and References 

This project uses tools from the MPEG-7 Audio Reference Software for extracting low-level audio descriptors. 
The software is developed by Michael Casey and maintained by Goldsmiths College, University of London, UK. 
For more details, refer to the official MPEG-7 standard ISO/IEC 15938-4:2001 and to (http://mpeg7.doc.gold.ac.uk/mirror/index.html). I have modified some lines of the MPEG-7 scripts to better adapt them to the code of this repository.

This project uses functions from the MIRToolbox (licensed under MIT License).
The software is developed by Olivier Lartillot, Petri Toivanien, Pasi Saari and Tuomas Eerola, members of the Finnish Centre of Excellence in Interdisciplinary Music Research, University of Jyväskylä, Finland.
For more details, refer to  (https://github.com/olivierlar/mirtoolbox)

To run this project, it is necessary to install the MATLAB Audio Toolbox.



