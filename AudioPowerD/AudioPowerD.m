function AudioPower_SeriesOfScalar = AudioPowerD(auData,totalSampleNum,samplingRate,scalingRatio,elementNum,weight) 

%% File: AudioPowerType.m
%%
%% ------------- AudioPowerType--------------------
%%
%% The function of this subroutine is to describe the temporally-smoothed 
%% instantaneous power (square of waveform values).
%%
%% 
%%
%% Copyright (c), IME, All Rights Reserved
%%
%% Author: Dong Yan Huang
%% Version: 1.0  Time: 28 October 2000 (N3489)
%% Last Modified: 27 Mar 2001 (N3704)

write_flag = 1; % write_flag = 1 (write DDL in XML file)
                % write_flag = 0 (not write DDL in XML file)
                % write_flag = 2 (write values in XML file)
                
if length(weight)==0
   weight_flag = 0;
else
   weight_flag = 1;
end
hopSize = samplingRate; %1/Fs
sampleNum = 1; %1024
frameNum = floor(totalSampleNum/sampleNum); %(3*44100/1024)

sumElement = sum(elementNum); %3*Fs
AudioPower_SeriesOfScalar = zeros(frameNum, sumElement);
for i = 1:frameNum
   signal = auData(1+(i-1)*sampleNum:i*sampleNum);
   audioPowerData = signal.^2;
   meanValues = h_Mean_SeriesOfScalar(audioPowerData, scalingRatio, elementNum, weight_flag,weight, write_flag);
end


