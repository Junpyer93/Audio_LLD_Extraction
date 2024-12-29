%% ESTRAZIONE DELLE FEATURES DA UN SUONO 

% The goal is the iteration of extracting low-level features from a
% single sound

%% Clear Screen And Workspace

clear; clf; clc;

%% Input Data
                                                                       
prompt = {'Filename:','Start Second(s):', 'End Second(s):', 'Partials(y/n):','Plot(y/n):'};           
dlgtitle = 'Input';                                                                  
dims = [1 50];                                                                       
definput = {'path_audiofile.wav','0','2','y','y'};    
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput,opts);                               

%% Setting Data

filename = string(answer(1));                                              % Folder name
Startsec = cell2mat(answer(2));                                            % Seconds to analyze
Endsec = cell2mat(answer(3));
flagp = string(answer(4));                                                 % Flag for partial window
flagplot = string(answer(5));                                              % Flag for descriptor's plot

infoTrack = audioinfo(filename);                                           % Track info

Fs = infoTrack.SampleRate;                                                 % Sample frequency

%% MPEG7 Extraction and Plotting 

if infoTrack.NumChannels == 2
     LLDs = LLDsfeatureStereo(filename,Fs,Startsec,Endsec,flagp);          % Features extraction
     if strcmp(flagplot, 'y')
      FeaturePlotLowLevelStereo;
     end
else 
     LLDs = LLDsfeatureMono(filename,Fs,sec,flagp);                        % Features extraction
     if strcmp(flagplot, 'y')
      FeaturePlotLowLevelMono;
     end
end

%% Psychoacoustical Extraction

PsychoF = PsychoFeature(answer(1),LLDs.AWF.Input,Fs); 

%% NO-MPEG7 Extraction

OtherLowLevel = OtherFeature(LLDs.AWF.Input,Fs); 



