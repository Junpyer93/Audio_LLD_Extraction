%% Features Extraction From A Dataset

% The objective is the iteration of the extraction of low-level features from an entire dataset

%% Clear Screen and Workspace

clear; 
clf;
clc;

%% Input Data

prompt = {'Folder:', 'Number of Seconds:', 'Partials:', 'Plot:'};
dlgtitle = 'Input';                                                        
dims = [1 35];                                                             
definput = {'path_dataset', '1','y','n'};       
opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle,dims,definput);                          

%% Setting Data 

folder = string(answer(1));
sec = cell2mat(answer(2));                                                  
flagp = string(answer(3));                                                  
flagplot = string(answer(4));

ADS = audioDatastore(folder ,'FileExtensions','.wav');

%% Features Extraction from Dataset

for i=1:1:length(ADS.Files)

disp(i);

filename = string(ADS.Files(i));                   

infoTrack = audioinfo(filename);

Fs = infoTrack.SampleRate;

%% MPEG7 Extraction             

if infoTrack.NumChannels == 2
   LLDs(i) = LLDsfeatureStereo(filename,Fs,sec,flagp);                     
   if strcmp(flagplot, 'y')
      FeaturePlotLowLevelStereo;
   end
else 
   LLDs(i) = LLDsfeatureMono(filename,Fs,sec,flagp);                       
   if strcmp(flagplot, 'y')
      FeaturePlotLowLevelMono;
   end
end

%% PsychoFeature Extraction

LLDs(i).PsychoF = PsychoFeature(ADS.Files(i),LLDs(i).AWF.Input,Fs);

%% NO-MPEG7 Extraction

LLDs(i).OtherLowLevel = OtherFeature(LLDs(i).AWF.Input,Fs); 


end