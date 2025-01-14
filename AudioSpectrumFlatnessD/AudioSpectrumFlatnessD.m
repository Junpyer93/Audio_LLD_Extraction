function [AudioSpectrumFlatness ,lo_edge, hi_edge, XMLFile ] = AudioSpectrumFlatnessD(audioFile,hopSize,loEdge,hiEdge,writeXML,XMLFile)

% This function describes the spectral flatness measure of the audio signal
% The frequency range is divided into numbands logarithmically spaced bands
% of 1/4 octave width, starting at loEdge up to hiEdge


% Written by Melanie Jackson & Juergen Herre
% Version 2.0 25 July 2001

% Modified by Melanie Jackson 10 December 2001
% Modified by MJ 17 January 2002
% Modified 16/04/2002 by Thibaut Sacreste - add XML generation
% Modified 19/04/2002 by Thibaut Sacreste - changed to be a stand alone function
% Modified 30/04/2002 by Thorsten Kastner - changed function call (loEdge and hiEdge can be set)
%                                         - modified check for loEdge and hiEdge
% Modified 12/06/2002 by Thorsten Kastner - adapted coefficient grouping to variable loEdge and hiEdge
%                                         - values for loEdge and hiEdge will be recalculated if necessary
%--------------------------------------------------------------------
% audioFile is the name of the audio file to process
% 2 types of files can be read: .wav and .au 

% writeXML is a flag for the generation of the XML file
% writeXML=0 -> no generation
% writeXML=1 -> generation

% XMLFile is the name of the XML file to be generated (optional)


%--------------------------------------------------------------------
% Initialisation:

%Questa parte è aggiunta da me

if nargin<4 
    writeXML=0;
end %poichè non abbiamo detto di scrivere in file XML
if nargin<3 
    loEdge = ' ';
    hiEdge = ' ';
end
if nargin<2 
   hopSize='PT07N1000F'; 
end %poichè non abbiamo passato il valore dell'intervallo di tempo tra due time frames successivi
%PT10N1000F E' un formato 


% Read in audio file


[audioData, sr] = audioread(audioFile);

% Descriptors only deal with monaural recordings 
if size(audioData,2)>1 
    % in the wavread function the second dimension contains the number of channels
    audioData = mean(audioData')';
end

% Start calculating descriptors

% Hopsize conversion
% format PT10N1000F

%Questa parte è corretta da me
hop = '';
if (hopSize)
  i = find(hopSize=='N');
  hop = [str2num(hopSize(3:i-1)) str2num(hopSize(i+1:end-1))];
end


standvar = h_mpeg7init(sr,hop);

%Check loEdge and hiEdge

if isempty(loEdge)	% variable defined but no value
  loEdge = 250;
else
  loEdge = 2^(floor(log2(loEdge/1000)*4)/4)*1000 ; % Setting exact value for loEdge, rounding to next lower frequency if necessary
  if (loEdge < 250 ) 
    loEdge = 250;       % No extraction below 250Hz
  end
end 


if isempty(hiEdge)	% Variable defined but no value
  hiEdge = 16000;       % Setting default hiedge
  disp(hiEdge);
else
  if (hiEdge >= sr/2)     % Check if value for hiedge is valid
  hiEdge =  min(hiEdge,sr/2-1); % sr/2-1 : Skipping extraction up to sr/2; Not possible due to 5 percent band overlap
  end
end
hiEdge = 2^(floor(log2(hiEdge/1000)*4)/4)*1000 ; % Setting exact value for hiEdge, rounding to next lower frequency if necessary
if (hiEdge*1.05 >= sr/2) 
    hiEdge = hiEdge/2^0.25 ; 
end %Now it's possible to check if hiEdge is valid
						       %If hiedge plus band overlap greater than sr/2 skip highest band


%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%------------------------------------------------------------------
% STFT with no overlap; window size equal to  hopsize.
% Zero padding of the last few frames will occur, to ensure there is one spectral frame
% for each corresponding power estimate in the power descriptor. 

[fftout,phase] = h_mpeg7getspec(audioData,standvar);

%-----------------------------------
% AudioSpectrumFlatness calculation:

fs = standvar.fs;
N = standvar.FFTsize;

numbands = floor(4*log2(hiEdge/loEdge)); %formula inversa
firstband = round(log2(loEdge/1000)*4);
overlap = 0.05;

disp(numbands);
grpsize = 1;
%ASF = ' ';
%fm = zeros(numbands,20000);
%am = zeros(numbands,20000);

for k = 1:numbands
  % Correzione con overlap
  f_lo = loEdge * (2^((k-1)/4)) * (1-overlap); 
  f_hi = loEdge * (2^((k)/4)) * (1+overlap);
  i_lo = round( f_lo/(fs/N) ) + 1;
  i_hi = round( f_hi/(fs/N) ) + 1;

  % Rounding of upper index according due to coefficient grouping
  if (k+firstband-1 >= 0)                   %Start grouping at 1kHz 
    grpsize = 2^ceil( (k+firstband )/4);
    i_hi = round((i_hi-i_lo+1)/grpsize)*grpsize + i_lo-1 ;
  else
    grpsize = 1;
  end
  tmp = fftout(i_lo:i_hi,:) .^ 2;         % PSD coefficients
  ncoeffs = i_hi - i_lo + 1;
 
  if (k+firstband-1 >= 0)                   % Coefficient grouping
    tmp2 = tmp(1:grpsize:ncoeffs,:);
    for g=2:grpsize
      tmp2 = tmp2 + tmp(g:grpsize:ncoeffs,:) ;
    end
    tmp = tmp2;
  end
  % Actual calculation
  ncoeffs = ncoeffs/grpsize ;
  tmp = tmp + 1e-50;       % avoid underflow for zero signals
  % Ottengo delle matrici dove le righe sono il numero di bande e le
  % colonne sono il numero dei coefficienti
  fm(k,:) = exp( sum(log(tmp))/ncoeffs ); % log processing avoids overflow
  am(k,:) = sum(tmp) / ncoeffs;
  %AudioSpectrumFlatness(k) = (fm(k,:)/am(k,:))'; %Per ogni banda dovrei avere dei valori
  %return
  disp(k);
end
%AudioSpectrumFlatness = ASF';
AudioSpectrumFlatness = (fm./am)';
lo_edge = loEdge;
hi_edge = hiEdge;

%---------------------
%XML generation:

if writeXML
    if ~exist('XMLFile')
      
        XMLFile=h_ASFtoXML(AudioSpectrumFlatness',loEdge,hiEdge,hopSize); 
    else
      XMLFile
        XMLFile=h_ASFtoXML(AudioSpectrumFlatness',loEdge,hiEdge,hopSize,XMLFile);
    end
end
