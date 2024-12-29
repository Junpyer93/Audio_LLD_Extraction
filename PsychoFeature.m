
function [PsychoFeatureLLDs] = PsychoFeature(filename, Input, Fs)

%%    (22) Loudness 

% Psychoacoustic feature of sound 

% It can be approximated as the RMS level of the signal in dB.
% Calculated generally by taking a series of frames
% and evaluating the square root of the sum of the squares of the values ​​
% of the windowed sample

L = struct();

loudMtr = loudnessMeter('ChannelWeights',[1 1],'SampleRate',Fs);

%loudnessMeter = Object that computes the loudness, range loudness, and true
%peak of an audio signal 

[L.momentary,L.shortTerm,L.integrated,L.range,L.peak] = loudMtr(Input);

%acousticPerceived loudness of an acoustic signal
L.acoustic = acousticLoudness (Input, Fs);

PsychoFeatureLLDs.Loudness = L;


%% Brightness

% - "Brightness" in analogy with visual brightness
% - Another parameter for perceptually distinguishing sounds
% - Indication of the amount of high-frequency content in
% sound using a measure such as the spectral centroid.
% - Refers to the harmonic content of the tone: Tones with
% stronger higher harmonics tend to be perceived as
% brighter. What is important is the relative strength of the
% harmonics, one should not look for a large number of harmonics but
% the relationship between their amplitude and the fundamental

% The brightness curve shows the evolution of
% brightness throughout the piece of music. High values ​​
% indicate moments in the music where most of the sound energy is high frequency, while low values ​​
% indicate moments where most of the sound energy is low frequency 

% Input:
% - Frame length
% - Cut off Frequency
% - Min RMS

prompt = {'Frame', 'Cutoff', 'MinRMS'};
dlgtitle = 'Brightness';
dims = [1 100];
definput = {'0.02321', '1500', '0.005'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

frame = str2double(answer(1));
cutoff = str2double(answer(2));
minRMS = str2double(answer(3));

outbrightness = mirbrightness(filename, 'Frame',frame, 'CutOff', cutoff, 'MinRMS', minRMS);
PsychoFeatureLLDs.brightness= mirgetdata(outbrightness);

%% Roughness

% - Emphasizes the "sporty" characteristics of the sounds
% produced by car engines.
% - It is a sensation of hearing that is created by the
% relatively rapid changes produced by the modulation frequencies
% of the range between about 15 and 300 Hz. It
% reaches its maximum near modulation frequencies of
% 70 Hz and decreases at higher modulation frequencies.
% - The unit of measurement is the Asper.

% The roughness curve shows the amount of sensory
% dissonance at each successive moment throughout the
% piece of music. This sensory dissonance corresponds
% to the phenomenon of "beating" when several sounds are heard with
% almost the same frequency, but with a few Hz of difference.
% When the roughness is high, the sounds feel harsher and
% contain stranger oscillations.

% S.P: 23.21 ms Windowsize --> 1024 samples, half overlapping

% Input:
%- Frame Length
% - Method
prompt = {'Frame', 'Method'};
dlgtitle = 'Roughness';
dims = [1 100];
definput = {'0.02321', 'Sethares'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

frame = str2double(answer(1));
method = string(answer(2));

outroughness = mirroughness(filename, 'Frame', frame, method);
PsychoFeatureLLDs.roughness = mirgetdata(outroughness);

figure(18)
plot(1:length(PsychoFeatureLLDs.roughness),PsychoFeatureLLDs.roughness,'b'); 
grid on;
title ("Roughness");
xlabel("Frames");
ylabel("Roughness Windowed (Adim)"');


%% Sharpness

% - Measurement of tone color;
% - Sensation value caused by high frequencies of a
% noise;

% - Unit of measurement is "acum". 1 Acum is attributed to a 1 kHz
% narrowband noise of 150 Hz width and
% value of 60 dB. For narrowband noises, sharpness
% increases as the center frequency increases and increases by a
% factor of two for a level increase from 30 to 90 dB.
% - Addition of low-frequency components is done
% to reduce the aggressiveness of the sounds of a particular
% product. This also increases the overall sound sensation.
% However, if the original sound sensation is not
% too high, the reduction in sharpness and therefore,
% aggressiveness may overcompensate for the increase in
% sound sensation with its effects on the overall sound quality.

% 23.21 ms Windowsize ---> 1024 Samples

% Input:
% - Weighting
% - Sound Field
% - Pressure Reference
% - Time Varying

prompt = {'Weighting', 'SoundField', 'PressureReference', 'TimeVarying'};
dlgtitle = 'Sharpness';
dims = [1 100];
definput = {'DIN 45692', 'free', '20e-6', 'false'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

weighting = string(answer(1));
soundfield = string(answer(2));
pressurereference = str2double(answer(3));
timevarying = string(answer(4));
if strcmp(timevarying, 'false')
    timevarying = false;
else
    timevarying = true;
end

PsychoFeatureLLDs.Sharpness = acousticSharpness(Input, Fs, 'Weighting', weighting, 'SoundField', soundfield, ...
                            'PressureReference', pressurereference, 'TimeVarying', timevarying);

%% Fluctuation Strength

% Objective: - Rhythmic estimation based on the calculation of the transformed spectrogram
% modulated by hearing and then estimation of the spectrum
% in each band

% Description: - Similar to Roughness, but reaches its maximum at modulation
% frequencies of about 4 Hz.
% - The input signal in the fluctuation model
% strength is the same as in the roughness model.

% - Unit of measurement "vacil"
% - Crucial role in speech evaluation

% Input:
% - Frame
% - Frame Rate
% - MinRes

prompt = {'Frame', 'FrameRate Window', 'MinRes'};
dlgtitle = 'Fluctuation';
dims = [1 100];
definput = {'0.02321', '10', '0.01'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

frame = str2double(answer(1));
framerate = str2double(answer(2));
minRes = str2double(answer(3));

outfluctuation = mirfluctuation(filename, 'Frame', frame, framerate, 'MinRes', minRes);
PsychoFeatureLLDs.Fluctuation = mirgetdata(outfluctuation);

%% Subjective Duration
% Duration in seconds of each subsequent event.

% It is calculated by detecting the attack and decay phases
% of each event and taking the portion of the curve between the
% start and offset times i.e. over 40% of the
% maximum amplitude between the attack and
% decay times.

%outduration = mirduration(filename);
%PsychoFeatureLLDs.SubDuration = mirgetdata(outduration);

