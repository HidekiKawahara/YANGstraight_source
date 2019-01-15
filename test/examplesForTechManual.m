%% Example program for reference manual
% by Hideki Kawahara
% 13/Jan./2019

fs = 44100;
fl = 55;
fh = 1200;
output = designAnalyticWavelet(fs, fl, fh);
figure;
plot(output.wvlt(1).t_axis * fl, abs(output.wvlt(1).w), 'linewidth', 2);
hold all
plot(output.wvlt(1).t_axis * fl, real(output.wvlt(1).w), 'linewidth', 2);
plot(output.wvlt(1).t_axis * fl, imag(output.wvlt(1).w), 'linewidth', 2);
grid on;
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized time (re. 1/fc)');
legend('absolute value', 'real peart', 'imaginary part');
print -depsc defaultShapeExample.eps

%%
figure
fftl = 2^20;
fx = (0:fftl - 1) / fftl * fs;
semilogx(fx / fl, 20 * log10(abs(fft(output.wvlt(1).w, fftl))), 'linewidth', 2);
grid on;
axis([0.01 100 -350 0]);
grid on;
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized frequency (re, fc)');
ylabel('gain (dB)');
print -depsc defaultGainExample.eps

%%

fs = 44100;
fl = 55;
fh = 1200;
channels_oct = 3;
output = designAnalyticWavelet(fs, fl, fh, channels_oct);
fftl = 2^16;
fx = (0:fftl - 1) / fftl * fs;
figure;
for ii = 1:length(output.fc_list)
    semilogx(fx, 20 * log10(abs(fft(output.wvlt(ii).w, fftl))), 'linewidth', 2);
    hold all
end
grid on
axis([20 3000 -150 0]);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('frequency (Hz)');
ylabel('gain (dB)');
print -depsc filterAllocationExample.eps

%%
axis([20 3000 -3 0]);
print -depsc filterAllocationMagExample.eps

%%

fs = 44100;
fl = 55;
fh = 1200;
channels_oct = 3;
mag_list = [1 1.05 1.25 1.5 2 3];
fftl = 2^17;
fx = (0:fftl - 1) / fftl * fs;
gain_figure = figure;
shape_figure = figure;
for ii = 1:length(mag_list)
    mag = mag_list(ii);
    output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag);
    figure(shape_figure);
    plot(output.wvlt(1).t_axis * fl, abs(output.wvlt(1).w), 'linewidth', 2);
    grid on;
    hold all
    drawnow
    figure(gain_figure);
    semilogx(fx / fl, 20 * log10(abs(fft(output.wvlt(1).w, fftl))), 'linewidth', 2);
    hold all
    grid on;
end
figure(shape_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized time (re. 1/fc)');
legend(num2str(mag_list'));
drawnow
print -depsc filteShapeStretchExample.eps

figure(gain_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized frequency (re. fc)');
ylabel('gain (dB)');
legend(num2str(mag_list'));
axis([0.01 100 -350 0])
print -depsc filteGainStretchExample.eps

%%
fs = 44100;
fl = 55;
fh = 1200;
channels_oct = 3;
mag = 1.0;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
fftl = 2^17;
fx = (0:fftl - 1) / fftl * fs;
gain_figure = figure;
shape_figure = figure;
for ii = 1:length(wintype_list)
    output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype_list{ii});
    figure(shape_figure)
    plot(output.wvlt(1).t_axis * fl, abs(output.wvlt(1).w), 'linewidth', 2);
    grid on;
    hold all
    drawnow
    figure(gain_figure)
    semilogx(fx / fl, 20 * log10(abs(fft(output.wvlt(1).w, fftl))), 'linewidth', 2);
    hold all
    grid on;
end
figure(shape_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized time (re. 1/fc)');
legend(wintype_list{:});
drawnow
print -depsc filteShapeEnvExample.eps

figure(gain_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized frequency (re. fc)');
ylabel('gain (dB)');
legend(wintype_list{:}, 'location', 'best');
axis([0.01 100 -180 0])
print -depsc filteGainEnvExample.eps

%%
fs = 44100;
fl = 55;
fh = 1200;
channels_oct = 3;
mag = 1.0;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
duration_result = zeros(length(wintype_list), 1);
for ii = 1:length(wintype_list)
    output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype_list{ii});
    duration_result(ii) = sqrt(sum(abs(output.wvlt(1).w) .^ 2 .* output.wvlt(1).t_axis .^ 2) ...
        / sum(abs(output.wvlt(1).w) .^ 2)) * fl;
end

%%

fs = 44100;
fl = 55;
fh = 1200;
channels_oct = 3;
mag_base = 1.0;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
fftl = 2^17;
fx = (0:fftl - 1) / fftl * fs;
gain_figure = figure;
shape_figure = figure;
for ii = 1:length(wintype_list)
    mag = mag_base / duration_result(ii) * duration_result(1);
    output = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype_list{ii});
    figure(shape_figure)
    plot(output.wvlt(1).t_axis * fl, abs(output.wvlt(1).w), 'linewidth', 2);
    grid on;
    hold all
    drawnow
    figure(gain_figure)
    semilogx(fx / fl, 20 * log10(abs(fft(output.wvlt(1).w, fftl))), 'linewidth', 2);
    hold all
    grid on;
end
figure(shape_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized time (re. 1/fc)');
legend(wintype_list{:});
drawnow
print -depsc filteShapeEnvNExample.eps

figure(gain_figure);
set(gca, 'fontsize', 15, 'linewidth', 2);
xlabel('normalized frequency (re. fc)');
ylabel('gain (dB)');
legend(wintype_list{:}, 'location', 'best');
axis([0.01 100 -180 0])
print -depsc filteGainEnvNExample.eps

%%

fs = 44100;
fo = 100;
fl = 50;
fh = 1200;
duration = 0.071;
tt = (1:round(duration * fs))' / fs;
x = tt * 0;
nto = round(fs / fo);
x(1:nto:end) = 1;
wvltStr = designAnalyticWavelet(fs, fl, fh);
output = waveletAttributesAnalyzer(x, fs, wvltStr);
total_channel = size(output.rawWavelet, 2);
ytickfreq = 100:100:1200;
ytickLoc = interp1(wvltStr.fc_list, 1:total_channel, ytickfreq, 'linear', 'extrap');

figure;imagesc(tt([1 end]),[1 total_channel],angle(output.rawWavelet'));
axis('xy');colormap(hsv)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzPhaseMap.eps

figure; imagesc(tt([1 end]),[1 total_channel],log(max(fl, min(fh, output.inst_freq_map'))));
axis('xy'); colormap(jet)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzInstFreqMap.eps

figure; imagesc(tt([1 end]),[1 total_channel],max(-0.005, min(0.005,output.group_delay_map')));
axis('xy'); colormap(jet)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzGroupDelayMap.eps

dB_map = 10 * log10(output.amp_squared_map');
figure; imagesc(tt([1 end]),[1 total_channel],max(max(dB_map(:))-80, dB_map));
axis('xy'); colormap(jet)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzdBPwrMap.eps

%% visualization with normalization

figure; imagesc(tt([1 end]),[1 total_channel], ...
    log(max(sqrt(0.5), min(sqrt(2), (output.inst_freq_map * diag(1 ./ wvltStr.fc_list))'))));
axis('xy'); colormap(jet)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzInstFreqNrmMap.eps

%
figure; imagesc(tt([1 end]),[1 total_channel], ...
    max(-1, min(1, (output.group_delay_map * diag(wvltStr.fc_list))')));
axis('xy'); colormap(jet)
set(gca, 'fontsize', 13, 'ytick', ytickLoc, 'yticklabel', num2str(ytickfreq'));
xlabel('time (s)')
ylabel('frequency (Hz)');
print -depsc sixt100HzGroupDelayNrmMap.eps

%%
figure;
loglog(wvltStr.fc_list, max(0.1, output.inst_freq_map(1323:10:1764, :)));
hold on; loglog([fl fh], [fl fh], 'k');
grid on;
axis([fl fh fl fh]);
set(gca, 'fontsize', 13);
xlabel('frequency (Hz)');
ylabel('instantaneous frequency (Hz)');
title('sample 1323:10:1764');
print -depsc fixedPointInFrequency.eps

figure;plot(tt, tt+output.group_delay_map(:, 24:42));
hold on;plot(tt([1 end]), tt([1 end]), 'k')
grid on;
axis([tt(1) tt(end) tt(1) tt(end)]);
set(gca, 'fontsize', 13);
xlabel('time (s)')
ylabel('time + group delay (s)')
title('channel 24:42')
print -depsc fixedPointInTime.eps

%%

figure;loglog(wvltStr.fc_list, max(0.1, output.inst_freq_map(1323:10:1764, :) ...
    * diag(1 ./ wvltStr.fc_list)));
grid on;
axis([fl fh 1 / 2 2]);
set(gca, 'fontsize', 13);
xlabel('frequency (Hz)');
ylabel('frequency deviation (ratio)');
title('sample 1323:10:1764');
print -depsc fixedPointInFrequencyN.eps

figure;plot(tt, output.group_delay_map(:, 24:42) * diag(wvltStr.fc_list(24:42)));
grid on;
axis([tt(1) tt(end) -2 2]);
set(gca, 'fontsize', 13);
xlabel('time (s)')
ylabel('normalized group delay (re. 1/fc)')
title('channel 24:42')
print -depsc fixedPointInTimeN.eps

%%

fs = 44100;
fo = 100;
duration = 0.071;
tt = (1:round(duration * fs))' / fs;
x = tt * 0;
nto = round(fs / fo);
x(1:nto:end) = 1;

fl = 50;
fh = 1200;
channels_oct = 12;
mag = 1.25;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
figure('position', [680   218   574   760]);
for ii = 1:length(wintype_list)
    wvltStr = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype_list{ii});
    output = waveletAttributesAnalyzer(x, fs, wvltStr);
    subplot(7,1,ii);
    plot(tt(tt > 0.025 & tt < 0.045), output.inst_freq_map(tt > 0.025 & tt < 0.045, 12:14), ...
        'linewidth', 2);
    set(gca, 'fontsize', 13, 'linewidth', 2);
    tmpdata = output.inst_freq_map(tt > 0.025 & tt < 0.045, 12:14);
    maxy = max(abs(tmpdata(:) - 100));
    text(0.0253, 0.6 * maxy + 100, wintype_list{ii}, 'fontsize', 15);
    axis([0.025 0.045 (maxy * [-1 1] + 100)]);
    grid on;
    drawnow
end
xlabel('time (s)')
print -depsc instFreqWin100Hz.eps

%%

fs = 44100;
fo = 100;
duration = 0.071;
tt = (1:round(duration * fs))' / fs;
x = tt * 0;
nto = round(fs / fo);
x(1:nto:end) = 1;

fl = 50;
fh = 1200;
channels_oct = 12;
mag = 1.25;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
figure('position', [680   218   574   760]);
for ii = 1:length(wintype_list)
    wvltStr = designAnalyticWavelet(fs, fl, fh, channels_oct, mag, wintype_list{ii});
    output = waveletAttributesAnalyzer(x, fs, wvltStr);
    subplot(7,1,ii);
    plot(tt(tt > 0.025 & tt < 0.045), 1000000 * output.group_delay_map(tt > 0.025 & tt < 0.045, 11:13), ...
        'linewidth', 2);
    set(gca, 'fontsize', 13, 'linewidth', 2);
    tmpdata = 1000000 * output.group_delay_map(tt > 0.025 & tt < 0.045, 11:13);
    maxy = max(abs(tmpdata(:)));
    text(0.0253, 0.6 * maxy, wintype_list{ii}, 'fontsize', 15);
    axis([0.025 0.045 maxy * [-1 1]]);
    grid on;
    drawnow
end
xlabel('time (s)')
print -depsc groupDelayWin100Hz.eps

%% example of using sourceAttributesAnalysis

close all
clear variables
%%
dataBaseDir = '/Users/kawahara/Music/CMU_arctic/';
talker_dir = {'cmu_us_bdl_arctic', 'cmu_us_jmk_arctic', 'cmu_us_slt_arctic'};
talkerId = 3;
fileNames = dir([dataBaseDir talker_dir{talkerId} '/orig/*.wav']);
fileId = 3;
[x, fs] = audioread([fileNames(fileId).folder '/' fileNames(fileId).name]);
outStr = sourceAttributesAnalysis(x(:, 1), fs); % speech signal is in ch-1

tsig = (1:length(x)) / fs;
tx = outStr.time_axis_wavelet;
fl = outStr.wvltStrDs.fc_list(1);
fh = outStr.wvltStrDs.fc_list(end);
figure; 
semilogy(outStr.time_axis_wavelet, outStr.fixed_points_freq(:, 1:4), '.', ...
    'linewidth', 2);
hold all;
semilogy(tsig, x(:, 1)/max(abs(x(:, 1))) * 20 + 80, 'k');
grid on;
axis([tx([1 end]) fl fh]);
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('time (s)');
ylabel('frequency (Hz)');
title([talker_dir{talkerId} '  ' fileNames(fileId).name], 'interpreter', 'none');
legend(num2str((1:4)'));
print -depsc foFreqCandExample.eps

figure; 
plot(outStr.time_axis_wavelet, outStr.fixed_points_measure(:, 1:4), '.', ...
    'linewidth', 2);
hold all
plot(tsig, x(:, 1)/max(abs(x(:, 1))) * 4 - 10, 'k');
grid on;
axis([tx([1 end]) -15 50]);
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('time (s)');
ylabel('estimated SNR (dB)');
title([talker_dir{talkerId} '  ' fileNames(fileId).name], 'interpreter', 'none');
legend(num2str((1:4)'));
print -depsc foSNRCandExample.eps

%%

outStr = sourceAttributesAnalysis(x(:, 1), fs); % speech signal is in ch-
fl = outStr.wvltStrDs.fc_list(1);
fh = outStr.wvltStrDs.fc_list(end);
figure;
semilogx(outStr.fixed_points_freq(:, 1:4), outStr.fixed_points_measure(:, 1:4), '.');
grid on;
set(gca, 'fontsize', 14, 'linewidth', 2);
axis([fl fh -15 50]);
xlabel('frequency (Hz)');
ylabel('estimated SNR (dB)');
title([talker_dir{talkerId} '  ' fileNames(fileId).name], 'interpreter', 'none');
print -depsc foCandScatterExample.eps

%%
fl = 100;
fh = 320;
channels_in_octave = 6;
dsOn = 0;
stretching_factor = 1.01;
wintype = 'dpss';
integration_time = 8;
outStr = sourceAttributesAnalysis(x(:, 1), fs, [1 length(x)], fl, fh, ...
            channels_in_octave, dsOn, stretching_factor, wintype, ...
            integration_time);
tsig = (1:length(x)) / fs;
tx = outStr.time_axis_wavelet;
fl = outStr.wvltStrDs.fc_list(1);
fh = outStr.wvltStrDs.fc_list(end);
figure; 
subplot(211)
semilogy(outStr.time_axis_wavelet, outStr.fixed_points_freq(:, 1:4), '.', ...
    'linewidth', 2);
hold all;
semilogy(tsig, x(:, 1)/max(abs(x(:, 1))) * 20 + 80, 'k');
grid on;
axis([tx([1 end]) fl fh]);
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('time (s)');
ylabel('frequency (Hz)');
title([talker_dir{talkerId} '  ' fileNames(fileId).name], 'interpreter', 'none');
legend(num2str((1:4)'), 'location', 'eastoutside');
subplot(212)
plot(outStr.time_axis_wavelet, outStr.fixed_points_measure(:, 1:4), '.', ...
    'linewidth', 2);
hold all
plot(tsig, x(:, 1)/max(abs(x(:, 1))) * 4 - 10, 'k');
grid on;
axis([tx([1 end]) -15 50]);
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('time (s)');
ylabel('estimated SNR (dB)');
title([talker_dir{talkerId} '  ' fileNames(fileId).name], 'interpreter', 'none');
legend(num2str((1:4)'), 'location', 'eastoutside');
print -depsc foFreqSNRASRCandExample.eps
