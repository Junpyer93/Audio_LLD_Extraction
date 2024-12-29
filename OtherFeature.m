function [OtherLLDs] = OtherFeature(auData, fs)

%%    OTHER LOW-LEVEL DESCRIPTORS

% Part of the descriptors not included in the MPEG7 standard  

%% INIZIALIZATION

Fs = fs; %Sample rate (This is the default rate used in the audioread function) 1s = 44100 samples.
Input = auData; %Signal Audio Path
OtherLLDs = struct();

%%    (19)Zero Crossing Rate                (ZCR)      

% Goal: Count the number of times the waveform crosses the zero
% axis.

% Operation: External function used, normalized as per manual
%

OtherLLDs.ZeroCR = 0.5*ZCR(Input)*Fs;

%%    (20) Spectral Rolloff Frequency 

% - Calculate the frequency above which a % of the spectrum amplitude is concentrated.
% - It is a measure of the spectral shape useful for distinguishing % spoken parts from non-spoken parts.

% Lw = 1024, Overlap = Hopsize = 341.

% Input:
% - Threshold
% - Window Length
% - Overlap
% - SpectrumType

prompt = {'Threshold', 'Window Length', 'Overlap', 'SpectrumType'};
dlgtitle = 'Spectral Rolloff Frequency';
dims = [1 100];
definput = {'0.85', '1024', '341', 'magnitude'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

thr = str2double(answer(1));
wlength = str2double(answer(2));
overlap = str2double(answer(3));
SpectrumType = string(answer(4));


% 1. All windows

OtherLLDs.SRP.Windowed = spectralRolloffPoint(Input, Fs, 'Threshold', thr, ...
                                     'Window', hamming(round(wlength)), ...
                                     'OverlapLength', overlap, ...
                                     'SpectrumType', SpectrumType);

% 2. Full sample

% FFT signal
nfft = 2^nextpow2(length(Input));  
spectrum = abs(fft(Input, nfft));  
powerSpectrum = spectrum.^2;       

% Energy Total
totalEnergy = sum(powerSpectrum);

% Calculate the energy threshold for roll-off (e.g. 85%)
thresholdEnergy = thr * totalEnergy;

% Find roll-off frequency
cumulativeEnergy = cumsum(powerSpectrum);
rolloffFreq = find(cumulativeEnergy >= thresholdEnergy, 1, 'first');
rolloffFreq = (rolloffFreq / nfft) * Fs;  

OtherLLDs.SRP.Total = rolloffFreq;                                 

%%    (21) Spectral Flux 

% - Calculate the average change in the amplitude of the signal spectrum
% between adjacent frames.

% - Measure of the local spectral change. It is used to
% separate music from speech.

% High SF values ​​indicate large spectral
% changes between consecutive frames

% SP: - Lw = 0.03; Overlap = 0.01; rectwin; magnitude;

% Input:
% - Overlap
% - Window Length
% - SpectrumType

prompt = {'Overlap', 'Window Length', 'SpectrumType'};
dlgtitle = 'Spectral Flux';
dims = [1 100];
definput = {'341', '1024', 'magnitude'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

wlength = str2double(answer(2));
overlap = str2double(answer(1));
SpectrumType = string(answer(3));

OtherLLDs.flux = spectralFlux(Input,Fs,'OverlapLength', round(overlap),...
                     'Window', hamming(round(wlength)),...
                     'SpectrumType', SpectrumType);


%%    (23-25) MFCC (Mel Frequency Cepstrum Coefficients) 

% Used for speech recognition as opposed to MPEG7 LLDs.

% - Cepstral parameterization coefficients refer to the
% concept of cepstrum .
% - A non-linear scale is used in the
% domain of mel frequency
% Mel: - unit of pitch, a formula is used to convert a
% frequency from Hertz to a Pitch in Mel. (Be careful of the pitch seen above!)
% - The Mel scale is a pitch scale judged by listeners to
% be equal in distance from each other. Below 500 Hz the Mel scale and
% the Hertz scale coincide, above that increasingly larger intervals are
% judged by listeners to produce equal pitch increments.
% - MFCCs are based on the extraction of signal energy
% within frequency bands i.e. a series of triangular filters
% whose central frequencies are spaced according to the Mel scale.
% - The Mel Non-linear scale takes into account human perception mechanisms of
% frequencies, that is, we tend to perceive low frequencies more selectively
% than higher ones.
% - The input signal is divided into overlapping frames of Nw samples, the
% duration of the frame varies from 20 to 40 ms, with an overlap of 50%
% between adjacent frames. In order to minimize the discontinuities of the signal
% at the edges of each frame, a window function defined as
% the Hanning function is used.
% - An FFT is applied to each frame and the absolute value is taken into account
% to obtain the amplitude spectrum. The spectrum is then processed by a
% mel-filter bank. The log-energy spectrum is measured inside each
% band-pass filter, resulting in a reduced representation
% of the spectrum. The cepstral coefficients are obtained through a DCT
% of the log-frequency spectrum.
% - The estimate of the derivative and acceleration of the MFCC is also taken into account to account for the temporal variations of the spectrum, this is done using the delta coefficients.

% Scelte progettuali

% Frame = 30 ms (20ms - 40ms); Overlap = 50%.
% Windowing = Rect (?)
% Mel filters = 24 typ.
% Coefficient cepstral numbers = 12 typ

% Input:
% - Window Length
% - Overlap
% - Number Coefficients
% - Rectification
% - Delta Window Length
% - Log Energy

prompt = {'WindowLength', 'Overlap', 'NumCoeffs', 'Rectification', 'DeltaWindowLength', 'LogEnergy'};
dlgtitle = 'MFCC';
dims = [1 100];
definput = {'0.03', '0.015', '12', 'log', '2', 'Append'};
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);

WindowLength = round(Fs*str2double(answer(1)));
Overlap = round(Fs*str2double(answer(2)));
NumCoeffs = str2double(answer(3));
Rectification = string(answer(4));
DeltaWindowLength = str2double(answer(5));
LogEnergy = string(answer(6));

% channel 1

[coeffsMFCC1,deltaMFCC1,deltaDeltaMFCC1,locMFCC1] = mfcc(Input(:,1), Fs, 'WindowLength', WindowLength, ...
    'OverlapLength', Overlap, 'NumCoeffs', NumCoeffs, ...
    'Rectification', Rectification, 'DeltaWindowLength', DeltaWindowLength, 'LogEnergy', LogEnergy);

MFCC1.coeffsMFCC1 = coeffsMFCC1;
MFCC1.deltaMFCC1 = deltaMFCC1;
MFCC1.deltaDeltaMFCC1 = deltaDeltaMFCC1;
MFCC1.locMFCC1 = locMFCC1;

MFCC.MFCC1 = MFCC1;

% channel 2

if size(Input,2) == 2
   [coeffsMFCC2,deltaMFCC2,deltaDeltaMFCC2,locMFCC2] = mfcc(Input(:,2), Fs, 'WindowLength', WindowLength, ...
    'OverlapLength', Overlap, 'NumCoeffs', NumCoeffs, ...
    'Rectification', Rectification, 'DeltaWindowLength', DeltaWindowLength, 'LogEnergy', LogEnergy);
   MFCC2.coeffsMFCC2 = coeffsMFCC2;
   MFCC2.deltaMFCC2 = deltaMFCC2;
   MFCC2.deltaDeltaMFCC2 = deltaDeltaMFCC2;
   MFCC2.locMFCC2 = locMFCC2;
   MFCC.MFCC2 = MFCC2;
end

OtherLLDs.MFCC = MFCC;



end



                       