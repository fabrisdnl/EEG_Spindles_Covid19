clear;
close all;
clc;

addpath();

%% LOAD DATA ------------------------------------------------------------- 1

ctrl_load = load('PROGETTO/CTRL033_nap.mat');
ctrl = ctrl_load.EEG;
icu_load = load('ICU023_nap.mat');
icu = icu_load.EEG;

fs = 250;   % Hz
T = 1/fs;     % s

Ns033 = size(ctrl,1);   % number of samples
Nch033 = size(ctrl,2);  % number of channels
time033 = 0 : T : (Ns033-1)*T; %s

Ns023 = size(icu,1);   % number of samples
Nch023 = size(icu,2);  % number of channels
time023 = 0 : T : (Ns023-1)*T; %s

%% FILTER ---------------------------------------------------------------- 2

slow_band = [9 , 12];
fast_band = [12 , 16];

% SLOW  9 - 12 Hz
fc = 9;     % high pass @ 9Hz
[b1,a1] = butter(4,fc/(fs/2),'high');
EEG_ctrl1 = filtfilt(b1,a1,ctrl);
EEG_icu1 = filtfilt(b1,a1,icu);

fc = 12;    % low-pass  @ 12Hz
[b1,a1] = butter(4,fc/(fs/2),'low');
slow_tracks_033 = filtfilt(b1,a1,EEG_ctrl1);
slow_tracks_023 = filtfilt(b1,a1,EEG_icu1);

% FAST  12 - 16 Hz
fc = 12;    % high pass  @ 12Hz
[b2,a2] = butter(4,fc/(fs/2),'high');08
EEG_ctrl2 = filtfilt(b2,a2,ctrl);
EEG_icu2 = filtfilt(b2,a2,icu);

fc = 16;    % low-pass  @ 16Hz
[b2,a2] = butter(4,fc/(fs/2),'low');
fast_tracks_033 = filtfilt(b2,a2,EEG_ctrl2);
fast_tracks_023 = filtfilt(b2,a2,EEG_icu2);

% Save Data 2D
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/slow_tracks_icu.mat','slow_tracks_033')
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/fast_tracks_icu.mat','fast_tracks_033')
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/slow_tracks_icu.mat','slow_tracks_023')
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/fast_tracks_icu.mat','fast_tracks_023')

%% ------------------------------------------------------------------> PLOT
%CRTL
figure(1)
subplot(311)
plot(time033, ctrl')
title('CTRL033 Raw EEG')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

subplot(312)
plot(time033, slow_tracks_033)
title('CTRL033 - SLOW - Filtered EEG [9-12] Hz')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

subplot(313)
plot(time033, fast_tracks_033)
title('CTRL033 - FAST - Filtered EEG [12-16] Hz')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

%% ICU
subplot(311)
plot(time023, icu')
title('ICU023 Raw EEG')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

subplot(312)
plot(time023, slow_tracks_023)
title('ICU023 - SLOW - Filtered EEG [9-12] Hz')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

subplot(313)
plot(time023, fast_tracks_023)
title('ICU023 - FAST - Filtered EEG [12-16] Hz')
xlabel('time [s]')
ylabel('amplitude [\mu V]')

%% EXTRACT SPINDLEs ------------------------------------------------------ 3

spindles_timing_033 = load('spindles_timing_033.mat');
spindles_timing_023 = load('spindles_timing_023.mat');

instants_033_fast = spindles_timing_033.fast;
instants_033_slow = spindles_timing_033.slow;
instants_023_fast = spindles_timing_023.fast;
instants_023_slow = spindles_timing_023.slow;

duration = 0.5 ; % 500 ms

spindles_samples = duration * fs; % converto la durata dello spindle in numero di campioni ovvero, Ns in 500 ms

spindle_signals_033_slow = zeros(spindles_samples , Nch033 , length(instants_033_slow));
spindle_signals_033_fast = zeros(spindles_samples , Nch033 , length(instants_033_fast));
spindle_signals_023_slow = zeros(spindles_samples , Nch023 , length(instants_023_slow));
spindle_signals_023_fast = zeros(spindles_samples , Nch023 , length(instants_023_fast));

spindle_signals_033_slow = extract_spindles(slow_tracks_033, instants_033_slow, spindles_samples, fs);
spindle_signals_033_fast = extract_spindles(fast_tracks_033, instants_033_fast, spindles_samples, fs);
spindle_signals_023_slow = extract_spindles(slow_tracks_023, instants_023_slow, spindles_samples, fs);
spindle_signals_023_fast = extract_spindles(fast_tracks_023, instants_023_fast, spindles_samples, fs);

% Save Data 3D
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/spindle_033_slow.mat', 'spindle_signals_033_slow');
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/spindle_033_fast.mat', 'spindle_signals_033_fast');
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/spindle_023_slow.mat', 'spindle_signals_023_slow');
% save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/spindle_023_fast.mat', 'spindle_signals_023_fast');

%% AVERAGE --------------------------------------------------------------- 4

average_spindles__signals_033_slow = mean(spindle_signals_033_slow, 2);
average_spindles__signals_033_fast = mean(spindle_signals_033_fast, 2);
average_spindles__signals_023_slow = mean(spindle_signals_023_slow, 2);
average_spindles__signals_023_fast = mean(spindle_signals_023_fast, 2);

%% PWELCH ---------------------------------------------------------------- 5

frequency_range = 0.6:0.1:20;

[pxx_033_slow, f_033_slow] = pwelch(average_spindles__signals_033_slow, size(average_spindles__signals_033_slow, 1), [], frequency_range, fs);
[pxx_033_fast, f_033_fast] = pwelch(average_spindles__signals_033_fast, size(average_spindles__signals_033_fast, 1), [], frequency_range, fs);
[pxx_023_slow, f_023_slow] = pwelch(average_spindles__signals_023_slow, size(average_spindles__signals_023_slow, 1), [], frequency_range, fs);
[pxx_023_fast, f_023_fast] = pwelch(average_spindles__signals_023_fast, size(average_spindles__signals_023_fast, 1), [], frequency_range, fs);

%% PLOT SPECTRA AND DISCARD ---------------------------------------------- 6

% Plot of all

figure(2);

subplot(221);
plot(f_033_slow, pxx_033_slow);
title('Spectrum CTRL - 033 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V^2/Hz)');

subplot(223);
plot(f_033_fast, pxx_033_fast);
title('Spectrum CTRL - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (\mu V^2/Hz)');

subplot(222);
plot(f_023_slow, pxx_023_slow);
title('Spectrum ICU - 023 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V^2/Hz)');

subplot(224);
plot(f_023_fast, pxx_023_fast);
title('Spectrum ICU - 023 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (\mu V^2/Hz)');

%% PLOT UNO ALLA VOLTA      ctrl
i = 10;
j = 23;
figure();
subplot(2, 1, 1);
plot(f_033_slow, sqrt(pxx_033_fast(:,i)*fs)/10);
title('VALID 10 - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]')
hold on 
xline(12);
xline(16);
hold off;
subplot(2, 1, 2);
plot(f_033_fast, sqrt(pxx_033_fast(:,j)*fs)/10);
title('REJECT 23 - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]')
hold on 
xline(12);
xline(16);
hold on;
%% PLOT TUTTI INSIEME       ctrl
a = 41;
for i = 1:a
    figure
plot(f_033_slow,sqrt(pxx_033_slow(:,i)*fs)/10);
title('Spectrum - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]')
hold on 
xline(9);
xline(12);
end
%% PLOT UNO ALLA VOLTA      icu
i = 2;
j = 29;
figure(2);
subplot(2, 1, 1);
plot(f_023_slow, sqrt(pxx_023_slow(:,i)*fs)/10);
title('VALID 2 - 023 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]')
hold on 
xline(9);
xline(12);
hold off;

subplot(2, 1, 2);
plot(f_023_fast, sqrt(pxx_023_slow(:,j)*fs)/10);
title('REJECT 29 - 023 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]');
hold on 
xline(9);
xline(12);
hold off;

%% PLOT TUTTI INSIEME       icu
a = 66;
for i = 1:a
    figure
plot(f_023_slow, sqrt(pxx_023_slow(:,i)*fs)/10);
title('Spectrum - 023 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('amplitude [\mu V]')
hold on 
xline(9);
xline(12);
end

%% Discard

max_peak_033_slow = zeros(size(pxx_033_slow, 2), 1);
idx_max_peak_033_slow = zeros(size(pxx_033_slow, 2), 1);
freq_max_peak_033_slow = zeros(size(pxx_033_slow, 2), 1);
VALID_SPINDLES_SIGNALS_033_SLOW = spindle_signals_033_slow;
VALID_TRIALS_033_SLOW = pxx_033_slow;
discard_033_slow = [];
k = 1;

for i = 1 : size(pxx_033_slow,2)
    max_peak_033_slow(i) = max(pxx_033_slow(:,i));
    idx_max_peak_033_slow(i) = find(pxx_033_slow(:,i) == max_peak_033_slow(i));
    freq_max_peak_033_slow(i) = f_033_slow(idx_max_peak_033_slow(i));
    if freq_max_peak_033_slow(i) < slow_band(1) |  freq_max_peak_033_slow(i) > slow_band(2)
    discard_033_slow(k) = i;
    k = k+1;
    end
end

max_peak_033_fast = zeros(size(pxx_033_fast, 2), 1);
idx_max_peak_033_fast = zeros(size(pxx_033_fast, 2), 1);
freq_max_peak_033_fast = zeros(size(pxx_033_fast, 2), 1);
VALID_SPINDLES_SIGNALS_033_FAST = spindle_signals_033_fast;
VALID_TRIALS_033_FAST = pxx_033_fast;
discard_033_fast = [];
k = 1;

for i = 1 : size(pxx_033_fast,2)
    max_peak_033_fast(i) = max(pxx_033_fast(:,i));
    idx_max_peak_033_fast(i) = find(pxx_033_fast(:,i) == max_peak_033_fast(i));
    freq_max_peak_033_fast(i) = f_033_fast(idx_max_peak_033_fast(i));
    if freq_max_peak_033_fast(i) < fast_band(1) |  freq_max_peak_033_fast(i) > fast_band(2)
        discard_033_fast(k) = i;
        k = k + 1;
    end
end

max_peak_023_slow = zeros(size(pxx_023_slow, 2), 1);
idx_max_peak_023_slow = zeros(size(pxx_023_slow, 2), 1);
freq_max_peak_023_slow = zeros(size(pxx_023_slow, 2), 1);
VALID_SPINDLES_SIGNALS_023_SLOW = spindle_signals_023_slow;
VALID_TRIALS_023_SLOW = pxx_023_slow;
discard_023_slow = [];
k = 1;

for i = 1 : size(pxx_023_slow,2)
    max_peak_023_slow(i) = max(pxx_023_slow(:,i));
    idx_max_peak_023_slow(i) = find(pxx_023_slow(:,i) == max_peak_023_slow(i));
    freq_max_peak_023_slow(i) = f_023_slow(idx_max_peak_023_slow(i));
    if freq_max_peak_023_slow(i) < slow_band(1) |  freq_max_peak_023_slow(i) > slow_band(2)
        discard_023_slow(k) = i;
        k = k + 1;
    end 
end

max_peak_023_fast = zeros(size(pxx_023_fast, 2), 1);
idx_max_peak_023_fast = zeros(size(pxx_023_fast, 2), 1);
freq_max_peak_023_fast = zeros(size(pxx_023_fast, 2), 1);
VALID_SPINDLES_SIGNALS_023_FAST = spindle_signals_023_fast;
VALID_TRIALS_023_FAST = pxx_023_fast;
discard_023_fast = [];
k = 1;

for i = 1 : size(pxx_023_fast,2)
    max_peak_023_fast(i) = max(pxx_023_fast(:,i));
    idx_max_peak_023_fast(i) = find(pxx_023_fast(:,i) == max_peak_023_fast(i));
    freq_max_peak_023_fast(i) = f_023_fast(idx_max_peak_023_fast(i));
    if freq_max_peak_023_fast(i) < fast_band(1) |  freq_max_peak_023_fast(i) > fast_band(2)
        discard_023_fast(k) = i;
        k = k + 1;
    end
end

disp('discard_033_slow'); disp(discard_033_slow)
disp('discard_033_fast'); disp(discard_033_fast)
disp('discard_023_slow'); disp(discard_023_slow)
disp('discard_023_fast'); disp(discard_023_fast)

VALID_SPINDLES_SIGNALS_033_SLOW(: , : , discard_033_slow) = [];
VALID_SPINDLES_SIGNALS_033_FAST(: , : , discard_033_fast) = [];
VALID_SPINDLES_SIGNALS_023_SLOW(: , : , discard_023_slow) = [];
VALID_SPINDLES_SIGNALS_023_FAST(: , : , discard_023_fast) = [];

VALID_TRIALS_033_SLOW(: , discard_033_slow) = [];
VALID_TRIALS_033_FAST(: , discard_033_fast) = [];
VALID_TRIALS_023_SLOW(: , discard_023_slow) = [];
VALID_TRIALS_023_FAST(: , discard_023_fast) = [];

figure(3);
subplot(421)
plot(f_033_slow, VALID_TRIALS_033_SLOW);
title('VALID Spectrum CTRL - 033 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

subplot(423)
plot(f_033_fast, VALID_TRIALS_033_FAST);
title('VALID Spectrum CTRL - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (\mu V_2/Hz)');

subplot(422)
plot(f_033_slow, pxx_033_slow);
title('Spectrum CTRL - 033 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

subplot(424)
plot(f_033_slow, pxx_033_fast);
title('Spectrum CTRL - 033 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

subplot(425)
plot(f_023_slow, VALID_TRIALS_023_SLOW);
title('VALID Spectrum ICU - 023 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

subplot(427)
plot(f_023_fast, VALID_TRIALS_023_FAST);
title('VALID Spectrum ICU - 023 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (\mu V_2/Hz)');

subplot(426)
plot(f_023_slow, pxx_023_slow);
title('Spectrum ICU - 023 Slow Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

subplot(428)
plot(f_023_slow, pxx_023_fast);
title('Spectrum ICU - 023 Fast Spindles');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency \mu V_2/Hz)');

%% AVERAGE AND PLOT ---------------------------------------------------- 7

AVG_VALID_SPINDLES_SIGNALS_033_SLOW = mean(VALID_SPINDLES_SIGNALS_033_SLOW, 3);
AVG_VALID_SPINDLES_SIGNALS_033_FAST = mean(VALID_SPINDLES_SIGNALS_033_FAST, 3);
AVG_VALID_SPINDLES_SIGNALS_023_SLOW = mean(VALID_SPINDLES_SIGNALS_023_SLOW, 3);
AVG_VALID_SPINDLES_SIGNALS_023_FAST = mean(VALID_SPINDLES_SIGNALS_023_FAST, 3);

Ns033_slow = size(AVG_VALID_SPINDLES_SIGNALS_033_SLOW,1);
Ns033_fast = size(AVG_VALID_SPINDLES_SIGNALS_033_FAST,1);
Ns023_slow = size(AVG_VALID_SPINDLES_SIGNALS_023_SLOW,1);
Ns023_fast = size(AVG_VALID_SPINDLES_SIGNALS_023_FAST,1);

time_033_slow = 0 : T : (Ns033_slow-1)*T;
time_033_fast = 0 : T : (Ns033_fast-1)*T;
time_023_slow = 0 : T : (Ns023_slow-1)*T;
time_023_fast = 0 : T : (Ns023_fast-1)*T;

figure(4)
subplot(221)
plot(time_033_slow, AVG_VALID_SPINDLES_SIGNALS_033_SLOW);
title('avg CTRL - 033 SLOW');
xlabel('time [s]');
ylabel('amplitude [\mu V]');
xlim

subplot(223)
plot(time_033_fast, AVG_VALID_SPINDLES_SIGNALS_033_FAST);
title('avg CTRL - 033 FAST');
xlabel('time [s]');
ylabel('amplitude [\mu V]');

subplot(222)
plot(time_023_slow, AVG_VALID_SPINDLES_SIGNALS_023_SLOW);
title('avg ICU - 023 SLOW');
xlabel('time [s]');
ylabel('amplitude [\mu V]');

subplot(224)
plot(time_023_fast, AVG_VALID_SPINDLES_SIGNALS_023_FAST);
title('avg ICU - 023 FAST');
xlabel('time [s]');
ylabel('amplitude [\mu V]');

%% ABSOLUTE , AVERAGE AND TOPPLOT --------------------------------------- 8

abs_spindles_signals_033_slow = abs(AVG_VALID_SPINDLES_SIGNALS_033_SLOW);
abs_spindles_signals_033_fast = abs(AVG_VALID_SPINDLES_SIGNALS_033_FAST);
abs_spindles_signals_023_slow = abs(AVG_VALID_SPINDLES_SIGNALS_023_SLOW);
abs_spindles_signals_023_fast = abs(AVG_VALID_SPINDLES_SIGNALS_023_FAST);

avg_abs_spindles_signals_033_slow = mean(abs_spindles_signals_033_slow, 1);
avg_abs_spindles_signals_033_fast = mean(abs_spindles_signals_033_fast, 1);
avg_abs_spindles_signals_023_slow = mean(abs_spindles_signals_023_slow, 1);
avg_abs_spindles_signals_023_fast = mean(abs_spindles_signals_023_fast, 1);

EEG.chanlocs = readlocs('GSN_204.sfp');

figure(5)
subplot(221)
topoplot(avg_abs_spindles_signals_033_slow,EEG.chanlocs);
colormap('jet');
colorbar
caxis([0 max(avg_abs_spindles_signals_033_slow)]);
title('Topography - CTRL - 033 SLOW');

subplot(223)
topoplot(avg_abs_spindles_signals_033_fast,EEG.chanlocs);
colormap('jet');
colorbar
caxis([0 max(avg_abs_spindles_signals_033_fast)]);
title('Topography - CTRL - 033 FAST');

subplot(222)
topoplot(avg_abs_spindles_signals_023_slow,EEG.chanlocs);
colormap('jet');
colorbar
caxis([0 max(avg_abs_spindles_signals_023_slow)]);
title('Topography - ICU - 023 SLOW');

subplot(224)
topoplot(avg_abs_spindles_signals_023_fast,EEG.chanlocs);
colormap('jet');
colorbar
caxis([0 max(avg_abs_spindles_signals_023_fast)]);
title('Topography - ICU - 023 FAST');

%% BRAINSTORM -----------------------------------------------------------9

% Reduce 3D to 2D
BRAINSTORM_033_SLOW = squeeze(AVG_VALID_SPINDLES_SIGNALS_033_SLOW);
BRAINSTORM_033_FAST = squeeze(AVG_VALID_SPINDLES_SIGNALS_033_FAST);
BRAINSTORM_023_SLOW = squeeze(AVG_VALID_SPINDLES_SIGNALS_023_SLOW);
BRAINSTORM_023_FAST = squeeze(AVG_VALID_SPINDLES_SIGNALS_023_FAST);

% Save Data 2D
 save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/BRAINSTORM_033_SLOW.set', 'BRAINSTORM_033_SLOW');
 save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/BRAINSTORM_033_FAST.mat', 'BRAINSTORM_033_FAST');
 save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/BRAINSTORM_023_SLOW.mat', 'BRAINSTORM_023_SLOW');
 save('/Users/nicve/Desktop/BIO ING/3_Neurorobotics and neurorehabilitation/PROGETTO/BRAINSTORM_023_FAST.mat', 'BRAINSTORM_023_FAST');












