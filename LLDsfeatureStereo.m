   function [LLDs] = LLDsfeatureStereo (folder, fs, Startsec ,Endsec, flagp)

% Extraction of low-level features according to the conventions adopted by
% MPEG7

%% Initialization

filename = folder;       % Audio signal path
Fs = fs;                 % Sample rate
ssec = Startsec;             % Seconds
esec = Endsec;
s = esec - ssec;
flagpartials = flagp;    % Flag that enables partial visualization

% Initializing Struct

LLDs = struct(); 


%% DONE

% 1. Elimination XML file

%% PARAMETER WINDOWING

% They will be used to create a uniform window with the different
% descriptors

% No Standard

% Windowsize =              1024 samples ------------ 23.21 ms
% Overlaplength = Hopsize = 341 samples ------------- 7.73 ms (1/3 Windowsize)

prompt = {'Length Window', 'Hopsize'};
dlgtitle = 'Parameter for Windowed Descriptors';
dims = [1 100];
definput = {'1024', '341'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);   

Lw = str2double(answer(1));         %Frame's length 
Hopsize = str2double(answer(2));     %Lenght between consecutive frames 


%% 1) AudioWaveform (AWF) 

% Provide signal estimation in the time domain

% For each frame it returns the maximum value and the value 
% minimum audio amplitude, taking into account that the frames are not overlapped.

% Input: 
% - Filename;
% - s = seconds;
% - Fs = Sample Frequency

% Output:
% Struct with:
% - Audio Input (samplesXchannels);
% - Sample Frequency
% - Time Vector;
% - Max Values;
% - Min Values;
% - Root First Values;
% - Variance Scalewise Values;

AWF = struct();

% Try-Catch used to account for the case where the length of the audio
% signal is inexistent

try
 [AWF.Input, AWF.Fs] = audioread(filename, [ssec*Fs esec*Fs]);                      
 flag = 0;
 N = length(AWF.Input(:,1));                                                
 AWF.t = (0:N-1)/Fs;                                                       
catch 
 [AWF.Input, AWF.Fs] = audioread(filename);                                
 flag = 1;
 N = length(AWF.Input(:,1));                                                
 AWF.t = (0:N-1)/Fs;   % Vettore dei tempi
end

% Prompt to get the waveform in compact form 

prompt = {'Scaling Ratio', 'Weight Flag', 'Weight', 'Root First'};
dlgtitle = 'Audio Waveform SeriesOfScalars';
dims = [1 100];
definput = {'256', '0', '0', '1'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);                               
scalingRatio = str2double(answer(1));
weight_flag = str2double(answer(2));
weight = str2double(answer(3));
rootFirst = str2double(answer(4));

totalSampleNum = N;
elementNum = floor(totalSampleNum/scalingRatio);
[Raw, AWF.maxValues, AWF.minValues, AWF.rootFirstValues,AWF.varianceScalewiseValues] = AudioWaveformD(AWF.Input,totalSampleNum,scalingRatio, elementNum, weight_flag, weight, 0, rootFirst);

LLDs.AWF = AWF;


%% (2)Audio Power  (AP)                     

% Allows you to measure the evolution of the signal's amplitude 
% as a function of time. Helps represent a rapid 
% representation of the spectrogram.

% For each frame it returns the instantaneous signal strength 
% calculated as the quadratic mean of the s(n) values

% SP: Lw = Hopsize 
%     l  = 0

prompt = {'Global AP', 'Instantaneous AP', 'Windowed AP', 'Series of Scalars AP'};
dlgtitle = 'Audio Power';
dims = [1 100];
definput = {'1', '1', '1', '1'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);                              
globalAP = str2double(answer(1));
instaAP = str2double(answer(2));
windowAP = str2double(answer(3));
frameAP = str2double(answer(4));


AP = struct();
 
% 1 Global AP

if globalAP == 1

    APchannel1 = rms((AWF.Input(:,1)).^2);
    APchannel2 = rms((AWF.Input(:,2)).^2);
    AP.APchannel1 = APchannel1;
    AP.APchannel2 = APchannel2;
end

% 2. Instantaneous AP

if instaAP == 1

% SP: Lw = Hopsize = 1/44100
% l = all samples    

    audioPowerData1 = zeros(N,1);
    audioPowerData2 = zeros(N,1);

    for i = 1:N
        signal1 = AWF.Input(i,1);
        signal2 = AWF.Input(i,2);
        audioPowerData1(i) = rms(signal1^2); 
        audioPowerData2(i) = rms(signal2^2);
    end

    AP.audioPowerData1pp = audioPowerData1;
    AP.audioPowerData2pp = audioPowerData2;

end

% 3. Windowed AP

if windowAP == 1
    if flag == 0 % length s
        iter = floor((s*Fs)/(Lw));
    else         % length less than s
        s = (length(AWF.Input)/Fs);
        iter = floor((s*Fs)/(Lw));
    end


    for i=0:1:(iter-1)
        InputLw1 = AWF.Input((1+(i*(Lw))):(Lw)+(i*(Lw)),1); 
        InputLw2 = AWF.Input((1+(i*(Lw))):(Lw)+(i*(Lw)),2);  
        audioPowerData1L(i+1) = rms(InputLw1.^2);
        audioPowerData2L(i+1) = rms(InputLw2.^2);
    end


    AP.audioPowerData1L = audioPowerData1L;
    AP.audioPowerData2L = audioPowerData2L;

end

% 4 AudioPower Series of Scalar

if frameAP == 1
    
    prompt = {'scalingRatio', 'elementNum', 'weight'};
    dlgtitle = 'Audio Power Series of Scalars';
    dims = [1 100];
    definput = {'','',''};
    opts.Interpreter = 'tex';
    answers = inputdlg(prompt,dlgtitle,dims,definput,opts);
    answers(1) = strrep(answers(1), '[', '');
    answers(1) = strrep(answers(1), ']', '');
    scalingRatio = cell2mat(answers(1));
    scalingRatio = str2num(scalingRatio);
    answers(2) = strrep(answers(2), '[', '');
    answers(2) = strrep(answers(2), ']', ''); 
    elementNum = cell2mat(answers(2));
    elementNum = str2num(elementNum);
    answers(3) = strrep(answers(3), '[', '');
    answers(3) = strrep(answers(3), ']', '');
    weight = cell2mat(answers(3));
    weight = str2num(weight);
    
    AP.frameAP = AudioPowerD(AWF.Input, N, Fs, scalingRatio, elementNum, weight);

end

LLDs.AP = AP;

%%   SILENCE SEGMENT             

% Define an audio segment with a "silence" label
% when there is no sound. 

% - Confidence = Confidence measure contained in the range
% Audiosegment Type [0 1] indicating the degree of certainty in 
% to which the identified silence corresponds 
% certainly in silence.

% - minDurationRef = Used to communicate a minimum
% time threshold within 
% a portion of the signal that is 
% identified as a silent segment. 
% It is usually used uniformly 
% to a decomposition segment such as 
% a parameter for extracting
% algorithm. 

% Can be used to distinguish different regions in the sound
% or to segment audio stream into different classes such as silence, speech,
% music and noise.

% Our choice ---> detectSpeech

% Identify the intervals in which the sound is absent

% We preferred to use the detectSpeech function
% identifies the intervals in which the sound is present 
% and then deduce which are the intervals in which 
% the sound is absent. 

% - STFT calculation for each frame
% - Spectral Spread and Short - Term Energy Spread calculation for
% each frame (Doc.)
% - Calculation of the relevant histograms.
% - Calculation of thresholds for histograms
% - Noise purification for Energy and Spread
% - Creation of comparison masks for Energy and Spectral Spread
% with the respective thresholds. To declare a frame containing 
% speech, a function must be above its threshold. 
% - Masks are combined. For a frame that is
% declared as speech, both the energy and the specral spread 
% must be above the threshold.
% - Regions declared as speech join if the distance between them is
% less than MergeDistance

% SP:   - Lw = 1024, Hamming, NO overlapping.
%       - mergedist = numel(window)*5
%       - resolution = hop length = numel (Window) - Overlaplength

% Input:       
% - Audio Input
% - Sample Frequency
% - Window
% - Overlap
% - MergeDistance
% - Threshold (If not specified it is calculated by histogram)

% Output:       - idxSpeech (Interval in which there is spoken)
%               - thresholds 

Silent = struct();

SampleSong = [1, s * Fs];

[Silent.idxSpeech1, Silent.thresholds] = detectSpeech(AWF.Input(:,1), Fs, 'Window', hamming(Lw));
Nspeech1 = length(Silent.idxSpeech1);

Silent.idxSilent1 = [];
if ~isempty(Silent.idxSpeech1)
    % Find silence before first speak
    if Silent.idxSpeech1(1) > SampleSong(1)
        Silent.idxSilent1 = [SampleSong(1), Silent.idxSpeech1(1) - 1];
    end
    
    % Find silence between consecutive speaks 
    for i = 1:(Nspeech1 - 1)
        startSilent = Silent.idxSpeech1(i) + 1;
        endSilent = Silent.idxSpeech1(i + 1) - 1;
        if startSilent <= endSilent
            Silent.idxSilent1 = [Silent.idxSilent1; startSilent, endSilent];
        end
    end
    
    % Find silence after speak
    if Silent.idxSpeech1(end) < SampleSong(end)
        Silent.idxSilent1 = [Silent.idxSilent1; Silent.idxSpeech1(end) + 1, SampleSong(end)];
    end
else
    % Se non c'è parlato, tutto è silenzio
    Silent.idxSilent1 = SampleSong;
end

% Channel 2
[Silent.idxSpeech2, Silent.thresholds] = detectSpeech(AWF.Input(:,2), Fs, 'Window', hamming(Lw));
Nspeech2 = length(Silent.idxSpeech2);

Silent.idxSilent2 = [];
if ~isempty(Silent.idxSpeech2)
    if Silent.idxSpeech2(1) > SampleSong(1)
        Silent.idxSilent2 = [SampleSong(1), Silent.idxSpeech2(1) - 1];
    end
    
    for i = 1:(Nspeech2 - 1)
        startSilent = Silent.idxSpeech2(i) + 1;
        endSilent = Silent.idxSpeech2(i + 1) - 1;
        if startSilent <= endSilent
            Silent.idxSilent2 = [Silent.idxSilent2; startSilent, endSilent];
        end
    end
    
    if Silent.idxSpeech2(end) < SampleSong(end)
        Silent.idxSilent2 = [Silent.idxSilent2; Silent.idxSpeech2(end) + 1, SampleSong(end)];
    end
else
    % All silent
    Silent.idxSilent2 = SampleSong;
end

LLDs.Silent = Silent;

%% BASIC SPECTRAL DESCRIPTORS

% Description of the signal power spectrum in the logarithmic frequency domain 

%% (3)Audio Spectrum Envelope 

% Generation of a reduced spectrogram of the audio signal
% original in the log frequency domain.

% The energy of the original power spectrum is added 
% within the frequency band series.

% SP: - Frequency band resolution: Expressed in octaves. If it were lower we would have 
% more frequency bands. 
% - Frequency band distribution: /Bands are logarithmically
%                               distributed between two edges 
%                               frequency of size
%                               2^rn * 1KHz with n integer value
%
%                               /By default the spectrum is distributed 
%                               between 62.5 Hz and 16 kHz taking into account
%                               also of the out of band (31.25-62.5 and 16k-Fs/2).
%                               The range is divided into octaves as the scale is
%                               Logarithmic. 
%                               /The ratio between two frequencies 
%                               at the end of an octave is
%                               2:1. Note that 16KHz = 62.5 Hz * 2^8.
% - Number of log bands:        Inside the Range [LoEdge, hiEdge]
%                                the number of log bands
%                               that correspond to r is Bin=8/r.
% - Output:                     In output I have a matrix formed
%                               by 10 columns. Each column indicates the
%                               corresponding value of the ASE
%                               in the given band (10 bands in total).

% - Window parameters: - Hopsize: 341
% - Overlap: 341
% - Window: 1024
% - Window type: Hamming

prompt = {'loEdge', 'hiEdge', 'octaveResolution'};
dlgtitle = 'Audio Power Frame';
definput = {'62.5','16000','1'};
answers = inputdlg(prompt,dlgtitle,dims,definput,opts);
attributegrp = struct();
attributegrp.loEdge = str2double(answers(1));
attributegrp.hiEdge = str2double(answers(2));
attributegrp.octaveResolution = cell2mat(answers(3));
hopSizeASE = floor(Hopsize*1000/Fs);
hopSizeASE = sprintf('PT%02dN1000F', hopSizeASE);
ASE = struct();

% Input: - audiofile: file path;
%        - hopSize: time interval between two successive frames (def. 7 ms)
%        - attributegrp: structure (HiEdge,LoEdge,Octaveresolution) (def. above)

% Output: - AudioSpectrumEnvelope: matrix ( ASE = rows, bands = columns)
%         - attributegrp: (HiEdge,LoEdge,Octaveresolution)
%         - map: frequency band mapping

[ASE.AudioSpectrumEnvelope, ASE.attributegrp, ASE.map] = AudioSpectrumEnvelopeD(filename,hopSizeASE,attributegrp);

LLDs.ASE = ASE;


%%   (4)Audio Spectrum Centroid  (ASC)    

% Provide center of gravity of a log-frequency power spectrum

% Indicates when a power spectrum is dominated by low or high frequencies
% and can be seen as an approximation of the sharpness of the signal.
% The log-frequency scale approximates the perception of frequencies
% in the human ear.

% - Window parameters: - Hopsize: 341
% - Overlap: 341
% - Window: 1024
% - Window type: Hamming

[ASC] =AudioSpectrumCentroidD(filename, hopSizeASE);

LLDs.ASC = ASC;


%%  (5)Audio Spectrum Spread    (ASS)       

% Define the second central moment of the power spectrum log-frequency.

% - Provides information about how much the spectrum is distributed
% around its centroid.
% - A low value indicates how much the spectrum is concentrated
% around the centroid, while a high value indicates
% a distribution over a large frequency range.
% - It is used to differentiate loud sounds from more tonal ones.

% - Windowing parameters: - Hopsize: 341
% - Overlap: 341
% - Window: 1024
% - Window type: Hamming

% N.B.: The centroid is also calculated within the function.

[ASS] = AudioSpectrumSpreadD(filename, hopSizeASE);

LLDs.ASS = ASS;


%%   (6)Audio Spectrum Flatness        (ASF)  

% Expresses the flatness properties of the spectral power.

% For a given frame it consists of a series of values,
% each expressing the deviation of the spectral power
% of the signal from a flat shape within a
% predefined frequency band.

% Usefulness: - It is a measure of how similar an audio signal is to white noise or how
% correlated a signal is to it.
% - High values ​​indicate noisiness. Low values ​​indicate
% a harmonic structure of the spectrum.
% - A strong deviation from the flat shape indicates tonal sounds.
% - It is used to match pairs of audio signals.

% - Windowing parameters: - Hopsize: 341
% - NO Overlap
% - Window: 1024
% - Window type: Hamming

% Input: - folder
% - hopsize (specified format)
% - Loedge
% - HiEdge
% - Flag

% Output: Rows = bands, columns = coefficients

[ASF] = AudioSpectrumFlatnessD(filename, hopSizeASE ,attributegrp.loEdge,attributegrp.hiEdge,0);

LLDs.ASF = ASF;


%% BASIC SIGNAL PARAMETERS

% - Provide a representation of the spectral power.
% - The following descriptors provide some complementary information
% describing the degree of harmonicity of the audio signals since
% the resolution of the ASE is too coarse to have
% a detailed representation of the harmonic peaks.

%%    (7)(Audio Harmonicity)     

% Definition of two measures of the harmonic properties of the
% spectrum.

% They rely on the standard fundamental frequency estimation method based on the normalized
% autocorrelation function of the signal.

% - This approach is widely used for local pitch estimation and is independent of the extraction of the audio fundamental
% frequency described later.
% - They are useful for defining the harmonic properties
% of the sound and for distinguishing harmonic from non-harmonic sounds.

% 7.1 Harmonic Ratio (HR)

% - Measurement of the proportion of harmonic components
% in the power spectrum.
% - An HR coefficient is found for every N sample frame of the
% original signal with a predefined hopsize.

% HR -> 0 White Noise;
% HR -> 1 Purely periodic;

% 7.2 Upper Limit Of Harmonicity (ULH)

% Estimation of the frequency above which the spectrum has no harmonic
% structure.

% Based on the power ratio between output and input
% of a comb filter tuned to the fundamental period
% of the signal estimated by the prec. 

% Lw = samplingRate*0.032       
% Hopsize = 1 (campione) 
% Win = Hanning

AH = struct();
ch1 = AWF.Input(:,1);
audioFilteredch1 = [];
for i = 1:size(Silent.idxSpeech1, 1) 
    startIdx = Silent.idxSpeech1(i, 1); 
    endIdx = Silent.idxSpeech1(i, 2);    
    audioFilteredch1 = [audioFilteredch1; ch1(startIdx:endIdx)]; 
end

audiowrite('output_audio_ch1.wav', audioFilteredch1, Fs);
[Input_filteredch1, AWF.Fs] = audioread(filename);
Nch1 = length(Input_filteredch1(:,1));

ch2 = AWF.Input(:,2);
audioFilteredch2 = [];
for i = 1:size(Silent.idxSpeech2, 1) 
    startIdx = Silent.idxSpeech2(i, 1); 
    endIdx = Silent.idxSpeech2(i, 2);    
    audioFilteredch2 = [audioFilteredch2; ch2(startIdx:endIdx)];
end

audiowrite('output_audio_ch2.wav', audioFilteredch2, Fs);
[Input_filteredch2, AWF.Fs] = audioread(filename);
Nch2 = length(Input_filteredch2(:,2));

[AH.harmonicRatioch1, AH.upperLimitOfHarmonicitych1] = AudioHarmonicityD(Input_filteredch1(:,1),Nch1,Fs);
[AH.harmonicRatioch2, AH.upperLimitOfHarmonicitych2] = AudioHarmonicityD(Input_filteredch2(:,2),Nch2,Fs);



% Input:                               - Input
%                                      - Number samples
%                                      - Sample Frequency
% Output:                              


LLDs.AH = AH;


%%   (8)(Audio Fundamental Frequency)           

% Estimation of the fundamental frequency in segments where it is
% assumed to be periodic.

% Pitch approximation in audio signals to obtain
% information on the harmonic structure of a periodic sound.

% - Varies according to the approach used. (Normalized
% Cross Correlation)
% - The fundamental must be searched within a specific range of
% frequencies.
% - The time window in which to identify the
% fundamental and any overlapping is also specified.
% - In some cases a confidence measure is used for the presence of
% periodicity contained between 0 and 1. For example the HR.

prompt = {'Method'};
dlgtitle = 'Audio Fundamental Frequency';
definput = {'NCF'};
answers = inputdlg(prompt,dlgtitle,dims,definput,opts);

% SP: (No Standard) -Lw = 1024, OverlapLength = 341.

method = string(answers(1));
f0.ch1 = pitch(AWF.Input(:,1),Fs,'Method',method,'WindowLength',round(Lw),'OverlapLength',round(Hopsize)); 
f0.ch2 = pitch(AWF.Input(:,2),Fs,'Method',method,'WindowLength',round(Lw),'OverlapLength',round(Hopsize));

LLDs.pitch = f0;

%% Timbral Descriptors    

% Description of the perceptual features of instrument sounds.

% Timbre: Characteristic that allows to distinguish two sounds
% that are equal in pitch, loudness and duration.
% Timbre is considered as a multidimensional feature.

%% Temporal Timbral          

% - Extracts from the signal envelope in the time domain.
% The signal envelope describes the energy change
% of the signal and is equivalent to the ADSR.

LAT = struct();

% According to the regulation, it is preferable to find the RMS envelope of the original signal
% by evaluating the RMS value frame by frame.

[LAT.Yu1, LAT.Yl1] = envelope(AWF.Input(:,1),Lw,'rms');
[LAT.Yu2, LAT.Yl2] = envelope(AWF.Input(:,2),Lw,'rms');

t = (AWF.t)';


%%    (9)Log Attack Time                         

% Calculate the time required to reach the maximum amplitude
% of the signal from a minimum time threshold.

% The main motivation is the description of individual sound samples
% from different musical instruments.

prompt = {'Threshold'};
dlgtitle = 'Log Attack Time';
definput = {'2'};
answers = inputdlg(prompt,dlgtitle,dims,definput,opts);
threshold = str2double(answers(1));

env1 = [t LAT.Yu1];
env2 = [t LAT.Yu2];
thr = threshold;

% Input: - envelope_bp (Energy Envelope (first column: time [second] | second column: value))
% - threshold_percent : percentage of maximum applied energy.

% Output: - LAT evaluated on the two channels (in Db)
%         - Start Attack
%         - Stop Attack


[LAT.lat1,LAT.StAtt1,LAT.StPAtt1] = LogAttackTimeD(env1,thr);
[LAT.lat2,LAT.StAtt2,LAT.StPAtt2] = LogAttackTimeD(env2,thr);

LLDs.LAT = LAT;

%%    (10)Temporal Centroid                 

% Find the time average of the signal energy envelope

TC = struct();

[TC.ch1] = TemporalCentroidD(env1);
[TC.ch2] = TemporalCentroidD(env2);

% Input: - envelope_bp (Energy envelope (first column: time [second] | second column: value))
% Output: - TC evaluated on the two channels

LLDs.TC = TC;


%% SPECTRAL TIMBRAL FEATURES 

% Description of the harmonic structure of the spectrum,
% in the linear space of frequency 

% Lw = 1024; Hopsize = 341 = Overlap; Window = Hamm.

% The extraction of these descriptors requires a priori
% the estimation of the fundamental and the identification
% of the harmonic components of the signal following an algorithm

% 0) Clean the signal from noise
% 1) Find FFT of the windowed signal
% 2) Consider the amplitude spectrum.
% 3) Estimate the peak frequency.
% 4) Find the peaks of the spectrum.
% 5) It is necessary to observe which of the peaks are harmonics.

% METHOD 1

% - Find the maximum values ​​of the spectrum amplitude around frequencies
%   that are multiples of f0. Harmonic peaks are found within
%   intervals centered at each multiple of f0.
% - The FFT bin kh corresponding to the h-th harmonic is found
%   by finding the maximum of the signal FFT amplitude within a
%   search interval (??) The intervals are specified in the book

% METHOD 2

% - Due to many noisy components in the signal it is not easy to find the
%   harmonic peaks, so other alternative methods are adopted that involve
%   cleaning the signal from noise.

% Values ​​are calculated locally (for the single frame) and globally
% for the entire sample

%% Extraction Harmonic Spectral Features 

% Initalization vector Harmonical Spectral Features

HSCp = 0;
HSDp = 0;
HSSp = 0;
HSVp = 0;

% Initialization of Harmonic Spectral Features vectors

if flag == 0 % length equal to s
   k = floor((s*Fs)/(Lw+Hopsize));
else         % length less than s
   s = (length(AWF.Input)/Fs);
   k = floor((s*Fs)/(Lw+Hopsize));
end
              
% Pre-allocation of local arrays

LHSC = zeros(k);
LHSD = zeros(k);
LHSS = zeros(k);
LHSV = zeros(k);

prompt = {'Method Smooth'};
dlgtitle = 'STFT';
definput = {'movmean'};
answers = inputdlg(prompt,dlgtitle,dims,definput,opts);

disp(answers(1))

% Calculating local expressions

% The STFT is calculated on first channel

for i=0:1:k-1 

InputLw1 = AWF.Input((1+(i*(Lw+Hopsize))):(Lw)+(i*(Lw+Hopsize)),1).*hamming(Lw); 
Nw1 = length(InputLw1);                               % Signal Lenght (first channel!)
disp(Nw1);
Inputsmooth1 = smoothdata(InputLw1, char(answers(1)));      % default:movmean
T = 1/Fs;                                            
tw = (0:Nw1-1)*T;                                     

% Plot: Frame width of the audio signal if enabled from Main

if strcmp(flagpartials, 'y')
   figure(3);
   plot(tw,InputLw1); 
   grid on;
   title ("Amplitude of Time Frame Audio Signal");
   xlabel("Time (t)");
   ylabel("Amplitude"');

% Plot: Frame width of the noise-free audio signal

   figure(4);
   plot(tw,Inputsmooth1);
   grid on;
   title ("Amplitude of Time Frame Audio Signal Smoothed");
   xlabel("Time (t)");
   ylabel("Amplitude");
end 

% Amplitude of the windowed signal in the frequency domain (DFT)

NFT = 2^nextpow2(Nw1); % The nextpow function defines the next power of 2 after Nw.
% In this way I fill the signal with zeros to improve the % performance of the fft.

ft = fft(Inputsmooth1,NFT);                               
a = abs(ft/NFT);                                          % Two-Sided Spectrum
a1 = a(1:NFT/2+1); 
a1(2:end-1) = 2*a1(2:end-1);
f = Fs*(0:(NFT/2))/NFT;                                   

% Plot: DFT of the windowed signal cleaned of noise

if strcmp(flagpartials, 'y')
figure(5);
plot(f,a1);
grid on;
title ("Single-Sided Amplitude Spectrum of Frame-Audio Signal Smooth");
xlabel("f(Hz)");
ylabel("|a1(f)|"');
end 

% Evaluation of spectral peaks 

[pks,locs] = findpeaks(a1); % I find positions and amplitudes of harmonics 
H = length(pks);            % number of spectral peaks 

% - Since the number of peaks, positions and amplitudes of the harmonics is
% variable by frame, I use cell-type variables to have
% peaks, positions and number of peaks outside the For. 

pksCell(i+1) = {pks}; % max number of harmonical peaks
locsCell(i+1) = {locs};
HCell(i+1) = {H}; 

%%     (11)Harmonic Spectral Centroid  (HSC) 

% Calculate the average, over the duration of the signal, weighted
% of the amplitude (on a linear scale) of the harmonic peaks
% of the spectrum.

% LHSC

LHSC(i+1) = HarmonicSpectralCentroidD(locs, pks, H); 

% HSCp Sum (outside the for, I calculate the total)

HSCp = + HarmonicSpectralCentroidD(locs, pks, H);

% Input: - locs = positions
%        - pks = peak values
%        - H = number of spectral peaks

%%    (12)Harmonic Spectral Deviation  (HSD)  


% Measurement of the deviation of harmonic peaks from the
% envelopes of the local spectrum.

% The curve clearly deals with the spectral modulation
% within the vibrato note.

% SE_v vector containing the spectral envelope corresponding to the
% positions of the harmonic peaks.

SE_v = 1:1:H;

% Identifying the spectral envelope
for j=1:1:H
    if H == 1  
        SE_v(j) = pks(j);  
    elseif j == 1
        SE_v(j) = (1/2)*(pks(j)+pks(j+1));
    elseif j == H
        SE_v(j) = (1/2)*(pks(H-1)+pks(H));
    else
        SE_v(j) = (1/3)*(pks(j-1)+pks(j)+pks(j+1));
    end
end


% LHSD

LHSD(i+1) = HarmonicSpectralDeviationD(pks, SE_v, H);

% Input: - pks = harmonic peaks
%        - SE_v = spectral envelope

% HSDp

HSDp = + HarmonicSpectralDeviationD(pks, SE_v, H);

%%    (13)Harmonic Spectral Spread      (HSS) 

% Measure of the average spectrum spread in relation to the HSC

% Vibrato modulation measurement

% LHSS

LHSS(i+1) = HarmonicSpectralSpreadD(locs, pks, H);

% HSSp

HSSp = + HarmonicSpectralSpreadD(locs, pks, H);


%%    (14)Harmonic Spectral Variation         (HSV)   

% Highlight the spectral variation between adjacent frames

% x1_v and x2_v contain the harmonic amplitude values ​​for consecutive
% frames.

% Note that the number of harmonics to be taken into account must be the same
% for x1 and x2.

if i == 0
    LHSV(i+1) = 1;
    pk1 = pks;
    HSVp = LHSV(i+1);
else 
    pk2 = pks;
    H1 = min(length(pk1),length(pk2));
    LHSV(i+1) = HarmonicSpectralVariationD(pk1, pk2, H1); %BE CAREFUL, H MUST BE THE SAME
    HSVp = + HarmonicSpectralVariationD(pk1, pk2, H1);
    pk1 = pk2;
end

% Input: - x1_v, x2_v = values ​​of the amplitudes of the harmonics of two consecutive
% frames.
%        - H = number of peaks taken into consideration

end

% Calculating final values

HSCentroid = (HSCp)/k;

HSC.HSCentroid = HSCentroid;
HSC.LHSC = LHSC;
LLDs.HSC = HSC;            % Harmonic Spectral Centroid

HSDeviation = (HSDp)/k;

HSD.HSDeviation = HSDeviation;
HSD.LHSD = LHSD;
LLDs.HSD = HSD;            % Harmonic Spectral Deviation

HSSpread = (HSSp)/k;

HSS.HSSpread = HSSpread;
HSS.LHSS = LHSS;
LLDs.HSS = HSS;            % Harmonic Spectral Spread

HSVariation = (HSVp)/k;

HSV.HSVariation = HSVariation;
HSV.LHSV = LHSV;
LLDs.HSV = HSV;            % Harmonic Spectral Variation


%%    (15)Spectral Centroid                 

% Provide the weighted average of the power of the discrete frequencies
% of the estimated spectrum along the sound segment. 

% It is used to distinguish the timbres of musical instruments
% Connected to the perceptual sharpness of the sound.
% Associated with the brightness of the sound.
                           
SC.ch1 = spectralCentroid(AWF.Input(:,1), Fs, ...
                            'Window', hamming(Lw), ...
                            'OverlapLength', Hopsize);
                        
SC.ch2 = spectralCentroid(AWF.Input(:,2), Fs, ...
                            'Window', hamming(Lw), ...
                            'OverlapLength', Hopsize);
                        
LLDs.SC = SC;

%Input:   - Audio Input
%         - Sample Frequency
%         - Window
%         - Overlap


%% SOUND CLASSIFICATION

% Project the audio signal spectrum into a
% low-dimensional representation, obtaining a
% classification in a more compact and efficient way.
%
% Some tasks in audio analysis consist in
% identifying similarities between sounds and classifying them.

% A classification system is based on the following steps:

% 1) Segmentation: Isolation of sound segments from the background in order to
% distinguish sounds from each other.

% 2) Feature Extraction: Properties of the sound useful for
% classification, it is necessary to have content-rich information (ASP,
% ASB, MFCC). The feature vector must not be too
% large, algorithms are used to reduce it.

% 3) Classification: A classifier uses the feature vector of
% reduced dimension to assign a sound to a category.

%%    (16)Audio Spectrum Basis              (ASB)   

% Reduction of statistical dependencies of observations and
% reduction of the size of features by capturing as much useful
% information as possible.

% Basis functions are extracted to best describe the
% spectral shape. Different models can be used.

% The function uses SVD and optionally ICA.

% Singular Value Decomposition (SVD): - Model used to transform a
% matrix X into three matrices U, D, V correlated
% with each other and with X. V is the
% matrix formed by the basis functions
% of size FxF (F maximum value
% of the range of logarithmic frequencies).

% The SVD transformation produces low-dimensional
% bases for the data,
% uncorrelated and also reduced in such a way that few
% basis functions can
% represent the entire spectrum

% ICA model: - Used with JADE algorithm 

prompt = {'Num_IC', 'JADE', 'hopSize', 'loEdge','hiEdge', 'octaveResolution'};
dlgtitle = 'Audio Spectrum Basis and Projection';
definput = {'10', '1', 'PT10N100F', '62.5', '16000', '1/8'};
answers = inputdlg(prompt,dlgtitle,dims,definput,opts);

Num_IC = str2double(answers(1));
flag = str2double(answers(2));
hopSize = char(answers(3));
loEdge = str2double(answers(4));
hiEdge = str2double(answers(5));
octaveResolution = str2double(answers(6));

ASB = struct();
 
[ASB.V,ASB.envBasis] = AudioSpectrumBasisD(ASE.AudioSpectrumEnvelope, Num_IC, 'Jade',flag, 'hopSize', hopSize, 'loEdge', loEdge, 'hiEdge', hiEdge, 'octaveResolution', octaveResolution);

LLDs.ASB = ASB;

% INPUT: 
% - X = AudioSpectrumEnvelope; you can directly enter the name of the
% file and the descriptor will automatically find the ASE.

% - NUM_IC = Number of components to extract, i.e. basis
% functions

% - varargin = Series of optional arguments to insert inside. For
% example Jade algorithm, used to find the basis functions (It activates if the value is 1).

% OUTPUT: 

% - ASB = Audio Basis Functions (nxNUM_IC) with n = number of frequency bins
% - env ​​= L2-norm envelope of log Spectrogram data (Euclidean Norm)


%%    (17)Audio Spectrum Projection          (ASP) 

ASP = struct();

[ASP.P,ASP.maxenv] = AudioSpectrumProjectionD(ASE.AudioSpectrumEnvelope,ASB.V, 'hopSize', hopSize, 'loEdge', loEdge, 'hiEdge', hiEdge, 'octaveResolution', octaveResolution);

% INPUT: - X = AudioSpectrumEnvelope; you can directly enter the name of the
% file and the descriptor will automatically find the ASE.

% - V = Radial basis matrix

% OUTPUT: - ASP = t x (1 + NUM_IC) matrix where each row contains 1 x L2-norm envelope
% coefficient and NUM_IC x spectral projection coefficients.

% - maxenv = maximum value of L2-norm envelope (used for SoundModelDS
% training data normalization)


LLDs.ASP = ASP;


   end