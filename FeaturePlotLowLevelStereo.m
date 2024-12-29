 %% Feature Plot Low-Level

% Plot Low-level descriptors


%% Audio Waveform Plot (1)

% Convert in dB (dB)
AWF.Input_dB = 20 * log10(abs(LLDs.AWF.Input) + eps); % eps evita log(0)

figure(1);

% Channel 1
subplot(2,2,1); 
plot(LLDs.AWF.t, AWF.Input_dB(:,1), 'b'); 
grid on;
title("Signal Audio dB (Channel 1)");
xlabel("Time (t)");
ylabel("Amplitude (dB)");

subplot(2,2,2); 
plot(LLDs.AWF.t, LLDs.AWF.Input(:,1), 'b');
grid on;
title("Signal Audio (Channel 1)");
xlabel("Time (t)");
ylabel("Amplitude");

% Channel 2
if size(AWF.Input_dB, 2) > 1
    subplot(2,2,3);
    plot(LLDs.AWF.t, AWF.Input_dB(:,2), 'r');
    grid on;
    title("Signal Audio dB (Channel 2)");
    xlabel("Time (t)");
    ylabel("Amplitude (dB)");
else
    subplot(2,2,3); 
    text(0.5, 0.5, 'Canale 2 not exist', 'HorizontalAlignment', 'center');
    axis off;
end

subplot(2,2,4); 
plot(LLDs.AWF.t, LLDs.AWF.Input(:,2), 'b');
grid on;
title("Signal Audio (Channel 2)");
xlabel("Time (t)");
ylabel("Amplitude");



%% Audio Power Plot (2)


% Conversion in dB
audioPowerData1pp_dB = 10 * log10(LLDs.AP.audioPowerData1pp + eps);
audioPowerData2pp_dB = 10 * log10(LLDs.AP.audioPowerData2pp + eps);
audioPowerData1L_dB = 10 * log10(LLDs.AP.audioPowerData1L + eps);
audioPowerData2L_dB = 10 * log10(LLDs.AP.audioPowerData2L + eps);

% Instantaneous power 
figure(2);

subplot(2,2,1);
plot(LLDs.AWF.t, audioPowerData1pp_dB, 'b'); 
grid on;
title("Instantaneous signal power (dB)(Ch1)");
xlabel("Time (t)");
ylabel("Power (dB)");

subplot(2,2,2);
plot(LLDs.AWF.t, LLDs.AP.audioPowerData1pp, 'b'); 
grid on;
title("Instantaneous signal power (Ch1)");
xlabel("Time (t)");
ylabel("Power");

subplot(2,2,3);
plot(LLDs.AWF.t, audioPowerData2pp_dB, 'r');
grid on;
title("Instantaneous signal power (dB)(Ch2)");
xlabel("Time (t)");
ylabel("Power (dB)");

subplot(2,2,4);
plot(LLDs.AWF.t, LLDs.AP.audioPowerData2pp, 'b'); 
grid on;
title("Instantaneous signal power (Ch2)");
xlabel("Time (t)");
ylabel("Power");


% Windowed power
figure(3);

subplot(2,2,1);
plot(1:length(audioPowerData1L_dB), audioPowerData1L_dB, 'b'); 
grid on;
title("Windowed signal power (dB)(Ch1)");
xlabel("Frames");
ylabel("Power (dB)");

subplot(2,2,2);
plot(1:length(LLDs.AP.audioPowerData1L), LLDs.AP.audioPowerData1L, 'b'); 
grid on;
title("Windowed signal power (Ch1)");
xlabel("Frames");
ylabel("Power");

subplot(2,2,3);
plot(1:length(audioPowerData2L_dB), audioPowerData2L_dB, 'r'); 
grid on;
title("Windowed signal power (dB) (Ch1)");
xlabel("Frames");
ylabel("Power (dB)");

subplot(2,2,4);
plot(1:length(LLDs.AP.audioPowerData2L), LLDs.AP.audioPowerData2L, 'b'); 
grid on;
title("Windowed signal power (Ch2)");
xlabel("Frames");
ylabel("Power");

%% Audio Spectrum Envelope Plot (3)

figure(4);
z = LLDs.ASE.AudioSpectrumEnvelope;
[nr,nc] = size(z);
t = 1:nr;
x = 1:nc;
surf(t,x,abs(z'),'EdgeColor','none');   
colorbar;
title ("Audio Spectrum Envelope (Linear Amplitude)");
axis xy; axis tight; colormap(jet); view(0,90);
ylabel("Frequency Bands (Hz)");
xlabel("Time Frames"); 


%% Audio Spectrum Centroid Plot (4)

figure(5);
plot((1:length(LLDs.ASC)),LLDs.ASC,'g'); 
grid on;
title ("Audio Spectrum Centroid");
xlabel ("Time frames");
ylabel ("ASC (Hz)");

%% Audio Spectrum Spread Plot (5)

figure(6);
plot((1:length(LLDs.ASS)),LLDs.ASS,'g'); 
grid on;
title ("Audio Spectrum Spread");
xlabel ("Time frames");
ylabel ("ASS (Hz)");

%% Audio Spectrum Flatness Plot (6)

figure(7);
z = LLDs.ASF;
[nr,nc] = size(z);
t = 1:nr;
x = 1:nc;
surf(t,x,abs(z'),'EdgeColor','none');   
colorbar;
title ("Audio Spectrum Flatness (adim)");
axis xy; axis tight; colormap(jet); view(0,90);
ylabel("Frequency Bands (Hz)");
xlabel("Time Frames");


%% Audio Harmonicity (7)

% Harmonic Ratio

figure(8);
plot((1:length(LLDs.AH.harmonicRatioch1)), LLDs.AH.harmonicRatioch1,'b');
hold on;
plot((1:length(LLDs.AH.harmonicRatioch2)), LLDs.AH.harmonicRatioch2,'r');
grid on;
title ("Harmonic Ratio ");
xlabel("Time frames");
ylabel("Harmonic Ratio (Adim)");
legend('channel1','channel2');

% Upper Limit Of Harmonicity

figure(9);
plot((1:length(LLDs.AH.upperLimitOfHarmonicitych1)), LLDs.AH.upperLimitOfHarmonicitych1,'b');
hold on;
plot((1:length(LLDs.AH.upperLimitOfHarmonicitych2)), LLDs.AH.upperLimitOfHarmonicitych2, 'r');
grid on;
title ("Upper Limit Of Harmonicity");
xlabel("Time frames");
ylabel("ULH (Hz)");
legend('channel1','channel2');

%% Pitch Plot (8)

figure(10);
plot((1:length(LLDs.pitch.ch1)),LLDs.pitch.ch1,'b'); 
hold on;
plot((1:length(LLDs.pitch.ch2)),LLDs.pitch.ch2,'r'); 
grid on;
title ("Audio Fundamental Frequency ");
xlabel ("Time frames");
ylabel ("Pitch(Hz)");
legend('channel1','channel2');

%% Signal Envelope 

figure(11)

subplot(2,1,1);
plot(LLDs.AWF.t,LLDs.LAT.Yu1,"Blue");
hold on;
plot(LLDs.AWF.t,LLDs.LAT.Yl1,"Blue");
title('Signal Envelope Ch1');
xlabel("Time (s)");
ylabel("Sign Env");
hold on;

% Log Attack Time Ch1

lx1 = line([LLDs.LAT.StAtt1 LLDs.LAT.StPAtt1], [0 0]); %% Log Attack Time Ch1
set(lx1, 'Color', 'black', 'LineStyle', '--');
grid on;

% Temporal Centroid Ch1

lxTc1 = line([LLDs.TC.ch1 LLDs.TC.ch1], [-0.3 0.3]); 
set(lxTc1, 'Color', 'black', 'LineStyle', '--');

subplot(2,1,2);
plot(LLDs.AWF.t,LLDs.LAT.Yl2,"Red");
hold on;
plot(LLDs.AWF.t,LLDs.LAT.Yu2,"Red");
title('Signal Envelope Ch2');
hold on;

% Log Attack Time CH2

lx2 = line([LLDs.LAT.StAtt2 LLDs.LAT.StPAtt2], [0 0]); 
set(lx2,'Color','black', 'LineStyle', '--');
grid on;
xlabel("Time (s)");
ylabel("Sign Env");

% Temporal Centroid Ch2

lxTc2 = line([LLDs.TC.ch2 LLDs.TC.ch2], [-0.3 0.3]); %Temporal Centroid ch2
set(lxTc2, 'Color', 'black', 'LineStyle', '--');


%% Harmonic Spectral Centroid Plot (11)

figure(12);
plot((1:length(LLDs.HSC.LHSC)),LLDs.HSC.LHSC);  
grid on;
title ("Harmonic Spectral Centroid");
xlabel ("Time frames");
ylabel ("LHSC (Hz)");

%% Harmonic Spectral Deviation Plot (12)

figure(13);
plot((1:length(LLDs.HSD.LHSD)),LLDs.HSD.LHSD); 
grid on;
title ("Harmonic Spectral Deviation");
xlabel ("Time frames");
ylabel ("LHSD (Hz)");

%% Harmonic Spectral Spread Plot (13)

figure(14);
plot((1:length(LLDs.HSS.LHSS)),LLDs.HSS.LHSS); 
grid on;
title ("Harmonic Spectral Spread");
xlabel ("Time frames");
ylabel ("LHSS (Hz)");

%% Harmonic Spectral Variation Plot (14)

figure(15);

plot((1:length(LLDs.HSV.LHSV)),LLDs.HSV.LHSV); 
grid on;
title ("Harmonic Spectral Variation");
xlabel ("Time frames");
ylabel ("LHSV (Hz)");

%% Plot Spectral Centroid

% Time evolution of Spectral centroid

figure(16);

% channel 1

subplot(2,1,1);

plot((1:length(LLDs.SC.ch1)),LLDs.SC.ch1);
grid on;
title("Spectral Centroid ch1");
xlabel("Time Frames");
ylabel("Spectral Centroid (Hz)");

% channel 2

subplot(2,1,2);

plot((1:length(LLDs.SC.ch2)),LLDs.SC.ch2);
grid on;
title("Spectral Centroid ch2");
xlabel("Time Frames");
ylabel("Spectral Centroid (Hz)");

%% Plot Power Spectrum Sample

% Power DFT Channel 1

Nw1 = length(LLDs.AWF.Input(:,1));                               
Inputsmooth1 = smoothdata(LLDs.AWF.Input(:,1));                  
NFT = 2^nextpow2(Nw1); 

ft = fft(Inputsmooth1,NFT);                                                       
powch1 = abs(ft).^2/NFT;                                                                  
powch1one = powch1(1:NFT/2+1); 
powch1one(2:end-1) = 2*powch1one(2:end-1);

% Power DFT Channel 2

Nw2 = length(LLDs.AWF.Input(:,2));                              
Inputsmooth2 = smoothdata(LLDs.AWF.Input(:,2));                  
NFT = 2^nextpow2(Nw2); 
          

ft = fft(Inputsmooth2,NFT);                                                        
powch2 = abs(ft).^2/NFT;                                                                  
powch2one = powch2(1:NFT/2+1); 
powch2one(2:end-1) = 2*powch2one(2:end-1);

f = Fs*(0:(NFT/2))/NFT;

figure(17)

% Power Spectrum Ch1

subplot(2,1,1);
plot(f,powch1one);
grid on;
title ("Power Spectrum Ch1");
xlabel ("f(Hz)");
ylabel ("Power Spectrum");

% Power Spectrum Ch2

subplot(2,1,2);
plot(f,powch2one);
grid on;
title ("Power Spectrum Ch2");
xlabel ("f(Hz)");
ylabel ("Power Spectrum");

