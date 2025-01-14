function [harmonicRatio,upperLimitOfHarmonicity] = AudioHarmonicityD(auData,totalSampleNum,samplingRate) 

%% 
%% ----------------AudioHarmonicityType-----------------------
%% The function of this subroutine is to describe 
%% the degree of harmonicity of an audio signal 
%% AudioHarmonicityType is a description of the 
%% spread of the log-frequency power spectrum.
%%
%% input:
%% -auData:                          incoming signal
%% -totalSampleNum:                  total Sample Number of signal
%% -samplingRate:                    sampling rate of the signal
%% 
%% output:
%% -harmonicRatio:                   Series of values of the harmonic ratio
%% -upperLimitOfHarmonicity:         Series of values of the UpperLimitOfHarmonicity
%% 

%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)
% === Initialization  ===

hopSize = 1;
sampleNum = fix(samplingRate*0.032);
sampleNum1 = fix(samplingRate*0.04);
frameNum = floor(totalSampleNum/sampleNum1);

% === analysis parameters  ===

windowSize = sampleNum;
window = hanning(windowSize);
fftSize = 2^nextpow2(windowSize);
freqResolution = samplingRate/fftSize;
lowFreqNum = fix(62.5/freqResolution);


upperLimitOfHarmonicity =[];    
   harmonicRatio =[];
%harmonicRatio = zeros(frameNum,fftSize/2-1);
%upperLimitOfHarmonicity = zeros(1, frameNum);
%x =zeros(1,fftSize/2-1);
frameNum
for i = 1:frameNum
   signal = auData(1+(i-1)*sampleNum1:i*sampleNum1);
   disp(i);
   disp(frameNum);
   
   %% === step a: calculate r(k), minj_pos(r(k)), h_k ===
   [minj_pos, H_k] = h_PeriodicSignalDetection(signal,samplingRate);
    if H_k <= 1.e-4
        disp(i);
   	   harmonicRatio(i,:) = zeros(1, sampleNum);
   	   upperLimitOfHarmoncity(i) = 0;
   	return
    else
         signal = auData(1+(i-1)*sampleNum:i*sampleNum);
         [minj_pos, H_k] = h_PeriodicSignalDetection(signal,samplingRate);
         
            if i*(sampleNum+ ceil(minj_pos)) <= totalSampleNum
               signal1 = auData(1+(i-1)*(sampleNum+ceil(minj_pos)):i*(sampleNum+ceil(minj_pos)));
            else
               signal1 = [signal' zeros(1,ceil(minj_pos))]; %zeros(1,ceil(minj_pos)) ho modificato la trasposta
            end
         signal2 = signal1';
      %[harmonicRatio(i,:),upperLimitOfHarmoncity(i)] = harmonicity(signal2,window,samplingRate,minj_pos);
      %% === Memory Allocation ==
	     snCombfiltered = zeros(sampleNum,1);
		
		 %% === Initialization ===
		 harmonicityRatio = [];
		 lowHarmonicity = [];
		 aflim = [];
	     lowPos = [];  
      %% === Calculate HarmonicRatio ===
      %% === step b) calculate DFT of s(i) and s(i)-s(i+j) === 
		 sn = signal2(1:sampleNum);
		 snFFT = h_amplFFT(sn,window);
		 snPower = snFFT.^2;
		 for i=1:sampleNum
             
   		     pos1 = i+minj_pos;
   		     sniplusj = interp1(signal2, pos1); 
   		     snCombfiltered(i) = signal2(i)-sniplusj;
         end
		 snCombFiltFFt = h_amplFFT(snCombfiltered,window);
		 snCombPower = snCombFiltFFt.^2;
         
        

		%% === Calculate power spectra and group the components below 62.5 Hz ===
		%%fftSize = length(snCombFiltFFt);
		%%freqResolution = samplingRate/fftSize;
		%%loFreqNum = fix(62.5/freqResolution);
		 harmonicityRatio =[harmonicityRatio sum(snCombPower(1:lowFreqNum))/sum(snPower(1:lowFreqNum))];
        

		%% === For each frequency f calculate the sum of power beyond that frequency,===
		%% === for both the original and comb-filtered signal, and take their ratio. ===

		for k=lowFreqNum+1:fftSize/2
   		    af = sum(abs(snCombPower(k:fftSize/2)).^2)/sum(abs(snPower(k:fftSize/2)).^2);
   		    aflim = [aflim af];
        end
		harmonicityRatio = [harmonicityRatio aflim];
       

		%% === Starting from fmax(sampleNum)  and moving down in frequency, find the ===
		%% === lowest frequency for which this ratio is smaller than a threshold (0.5). ===
		%% === If that frequency is 0, replace it by 31.25 Hz   ===
        %lowHarmonicity = [];
		for k=length(aflim):-1:1          %#ok<ALIGN>
 			if (aflim(k) < 0.5) %#ok<ALIGN>
   			lowHarmonicity = [lowHarmonicity aflim(k)];       
   			lowPos = [lowPos (k-1)];
            end
        end
        if isempty(lowPos) 
            fmin = 31.25;
        else
        j = lowPos(length(lowPos));
        fmin = j*freqResolution;
        end
        % modificato da me
		%j = lowPos(length(lowPos));
		%if j == 0 && isempty(j)
    	%	fmin = 31.25;
        %else
        %fmin = j*freqResolution;
	    end

	%% === Convert this value to an octave scale based on 1 kHz. ? ===
	freq1KIndex = ceil(1000/freqResolution);
	freqUpper = ceil((samplingRate/2)/freqResolution);
   upperLimitOfHarmonicity =[upperLimitOfHarmonicity fmin];    
   harmonicRatio =[harmonicRatio harmonicityRatio];
   
   %x = strcat(num2str(harmonicityRatio));



   %upperLimitOfHarmonicity(i)= x
   %harmonicRatio(i,:) =y 
%upperLimitOfHarmonicity  

   %	harmonicRatio(i,:) = x;
    %  UpperLimitOfHarmoncity(i) = y;
   end
end

     
     

