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
mag_base = 1.25;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
figure('position', [680   218   574   760]);
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(ii);
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
mag_base = 1.25;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
figure('position', [680   218   574   760]);
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(ii);
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
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
stretching_factor = 1.01 / mag_duration_list(7) * mag_duration_list(1);
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

%% Speed test skipping initialization

fs = 44100;
duration = 1;
tt = (1:duration * fs)' / fs;
x = randn(length(tt), 1);
initialoutput = sourceAttributesAnalysis(x, fs);
duration_list = 0.1 * 2 .^ (0:1/2:4);
elapsed_time_mean = zeros(length(duration_list), 2);
elapsed_time_SD = zeros(length(duration_list), 2);
iteration = 20;
for ii = 1:length(duration_list)
    tmp_default = zeros(iteration, 1);
    tmp_skip = zeros(iteration, 1);
    for jj = 1:iteration
        x = randn(round(duration_list(ii) * fs), 1);
        out_default = sourceAttributesAnalysis(x, fs);
        out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
        tmp_default(jj) = out_default.elapsedTime;
        tmp_skip(jj) = out_skip.elapsedTime;
    end
    elapsed_time_mean(ii, :) = [mean(tmp_default), mean(tmp_skip)];
    elapsed_time_SD(ii, :) = [std(tmp_default), std(tmp_skip)];
end
set(gca, 'fontsize', 14, 'linewidth', 2);
title('sixterm');
xlabel('signal duration (s)');
ylabel('elapsed time (s)');
legend(' default', ' skip', 'location', 'southeast')
axis([0.1 2 0.01 0.5]);
print -depsc sixtermEtimeSkip.eps


%%

fs = 44100;
fl = 55;
fh = 1200;
duration = 1;
ch_in_oct = 6;
dsOn = 1;
mag = 1.05;
winType = 'dpss';
tt = (1:duration * fs)' / fs;
x = randn(length(tt), 1);
initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
    ch_in_oct, dsOn, mag, winType);
duration_list = 0.1 * 2 .^ (0:1/2:4);
elapsed_time_mean = zeros(length(duration_list), 2);
elapsed_time_SD = zeros(length(duration_list), 2);
iteration = 20;
for ii = 1:length(duration_list)
    tmp_default = zeros(iteration, 1);
    tmp_skip = zeros(iteration, 1);
    for jj = 1:iteration
        x = randn(round(duration_list(ii) * fs), 1);
        out_default = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, winType);
        out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
        tmp_default(jj) = out_default.elapsedTime;
        tmp_skip(jj) = out_skip.elapsedTime;
    end
    elapsed_time_mean(ii, :) = [mean(tmp_default), mean(tmp_skip)];
    elapsed_time_SD(ii, :) = [std(tmp_default), std(tmp_skip)];
end
set(gca, 'fontsize', 14, 'linewidth', 2);
title('DPSS');
xlabel('signal duration (s)');
ylabel('elapsed time (s)');
legend(' default', ' skip', 'location', 'southeast')
axis([0.1 2 0.01 0.5]);
print -depsc dpssEtimeSkip.eps

%% Speed test for very short segment

fs = 44100;
fl = 55;
fh = 1200;
duration = 0.1;
ch_in_oct = 6;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
iteration = 100;
elapsed_time_median = zeros(length(wintype_list), 2);
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    tmp_default = zeros(iteration, 1);
    tmp_skip = zeros(iteration, 1);
    initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
        ch_in_oct, dsOn, mag, wintype_list{ii});
    for jj = 1:iteration + 1
        x = randn(round(duration * fs), 1);
        out_default = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii});
        out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
        tmp_default(max(1, jj - 1)) = out_default.elapsedTime;
        tmp_skip(max(1, jj - 1)) = out_skip.elapsedTime;
    end
    elapsed_time_median(ii, :) = [median(tmp_default), median(tmp_skip)];
end
figure;barh(elapsed_time_median);grid on;
set(gca, 'fontsize', 14, 'linewidth', 2, 'yticklabel', wintype_list);
xlabel('elapsed time (s)');
print -depsc winAndSpeedInit.eps

%% Speed test for duration and wintype

fs = 44100;
fl = 55;
fh = 1200;
duration = 0.1;
ch_in_oct = 6;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
duration_list = 0.1 * 2 .^ (0:1/2:4);
iteration = 10;
elapsed_time_median = zeros(length(duration_list), 2);
figure;
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    for kk = 1:length(duration_list)
        tmp_default = zeros(iteration, 1);
        tmp_skip = zeros(iteration, 1);
        initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii});
        for jj = 1:iteration
            x = randn(round(duration_list(kk) * fs), 1);
            out_default = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
                ch_in_oct, dsOn, mag, wintype_list{ii});
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tmp_default(jj) = out_default.elapsedTime;
            tmp_skip(jj) = out_skip.elapsedTime;
        end
        elapsed_time_median(kk, :) = [median(tmp_default), median(tmp_skip)];
    end
    loglog(duration_list, elapsed_time_median(:, 2), 'o-', 'linewidth', 2);
    grid on;hold on;
    axis([0.1 2 0.01 0.3]);
    drawnow
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('signal duration (s)');
ylabel('elapsed time (s)');
legend(wintype_list, 'location', 'southeast');
print -depsc winDurationAndSpeedInit.eps

%%

fs = 44100;
fl = 55;
fh = 1200;
duration = 0.1;
ch_in_oct = 6;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
%duration_list = 0.1 * 2 .^ (0:1/2:4);
fl_list = 20 * 2 .^ (0:1 / 2:4);
iteration = 50;
elapsed_time_median = zeros(length(fl_list), 2);
figure;
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    for kk = 1:length(fl_list)
        fl = fl_list(kk);
        tmp_default = zeros(iteration, 1);
        tmp_skip = zeros(iteration, 1);
        initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii});
        for jj = 1:iteration
            x = randn(round(duration * fs), 1);
            out_default = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
                ch_in_oct, dsOn, mag, wintype_list{ii});
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tmp_default(jj) = out_default.elapsedTime;
            tmp_skip(jj) = out_skip.elapsedTime;
        end
        elapsed_time_median(kk, :) = [median(tmp_default), median(tmp_skip)] * 1000;
    end
    loglog(fl_list, elapsed_time_median(:, 2), 'o-', 'linewidth', 2);
    grid on;hold on;
    %axis([0.1 2 0.01 0.2]);
    drawnow
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('freq Low (Hz)');
ylabel('elapsed time (ms)');
legend(wintype_list, 'location', 'northeast');
print -depsc winFlAndSpeedInit.eps

%% speed on density

fs = 44100;
fl = 55;
fh = 1200;
duration = 0.1;
ch_in_oct = 6;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
ch_oct_list = [4 6 12 24 48];
iteration = 50;
elapsed_time_median = zeros(length(ch_oct_list), 2);
figure;
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    for kk = 1:length(ch_oct_list)
        %fl = fl_list(kk);
        ch_in_oct = ch_oct_list(kk);
        tmp_default = zeros(iteration, 1);
        tmp_skip = zeros(iteration, 1);
        initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii});
        for jj = 1:iteration
            x = randn(round(duration * fs), 1);
            out_default = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
                ch_in_oct, dsOn, mag, wintype_list{ii});
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tmp_default(jj) = out_default.elapsedTime;
            tmp_skip(jj) = out_skip.elapsedTime;
        end
        elapsed_time_median(kk, :) = [median(tmp_default), median(tmp_skip)] * 1000;
    end
    loglog(ch_oct_list, elapsed_time_median(:, 2), 'o-', 'linewidth', 2);
    grid on;hold on;
    drawnow
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('channel density (ch/oct)');
ylabel('elapsed time (ms)');
legend(wintype_list, 'location', 'southeast');
print -depsc winChpOctAndSpeedInit.eps

%% SNR caribration check

fs = 44100;
fo = 240;
fl = fo * 2 ^ (-1 / 2);
fh = fo * 2 ^ (1 / 2);
duration = 1;
ch_in_oct = 12;
dsOn = 0;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
integration_time = 15;
iteration = 1;
snr_list = 0:10:80;
snr_median = zeros(length(snr_list), 1);
tt = (1:round(duration * fs))' / fs;
xp = tt * 0;
nto = round(fs / fo);
xp(1:nto:end) = 1;
figure;
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    for kk = 1:length(snr_list)
        snr = snr_list(kk);
        tmp_default = zeros(iteration, 1);
        tmp_skip = zeros(iteration, 1);
        initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii}, integration_time);
        for jj = 1:iteration
            x = xp + randn(round(duration * fs), 1) * std(xp) * 10 ^ (-snr / 20);
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tx = out_skip.time_axis_wavelet;
            tmp_skip(jj) = median(out_skip.fixed_points_measure(tx > 0.2 & tx < duration - 0.2, 1));
        end
        snr_median(kk) = median(tmp_skip);
    end
    plot(snr_list, snr_median, 'o-', 'linewidth', 2);
    grid on;hold on;
    drawnow
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('input SNR (dB)');
ylabel('estimated SNR (dB)');
legend(wintype_list, 'location', 'northwest');
title(['fo:' num2str(fo) ' mag:' num2str(mag_base) ' dsOn:' num2str(dsOn)]);
outFname = ['winSNRcalibration' num2str(dsOn) 'fo' num2str(fo) 'Hz.eps'];
print('-depsc',  outFname);

%%

figure;plot(out_skip.if_smooth_map(tx > 0.2 & tx < duration - 0.2, 7));grid on;
hold all
plot(out_skip.fixed_points_freq(tx > 0.2 & tx < duration - 0.2, 1));grid on;
plot(out_skip.inst_freq_map(tx > 0.2 & tx < duration - 0.2, 7));grid on;

%% Test on SNR and estimated fo

fs = 44100;
fo = 120;
fl = fo * 2 ^ (-1 / 2);
fh = fo * 2 ^ (1 / 2);
duration = 1;
ch_in_oct = 12;
dsOn = 0;
mag_base = 1.05;
wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
integration_time = 15;
iteration = 1;
snr_list = 0:10:80;
fo_err_SD = zeros(length(snr_list), 1);
tt = (1:round(duration * fs))' / fs;
xp = tt * 0;
nto = round(fs / fo);
xp(1:nto:end) = 1;
figure;
for ii = 1:length(wintype_list)
    mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
    x = randn(round(duration * fs), 1);
    for kk = 1:length(snr_list)
        snr = snr_list(kk);
        tmp_default = zeros(iteration, 1);
        tmp_skip = zeros(iteration, 1);
        initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
            ch_in_oct, dsOn, mag, wintype_list{ii}, integration_time);
        for jj = 1:iteration
            x = xp + randn(round(duration * fs), 1) * std(xp) * 10 ^ (-snr / 20);
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tx = out_skip.time_axis_wavelet;
            tmp_skip(jj) = std(out_skip.fixed_points_freq(tx > 0.2 & tx < duration - 0.2, 1)) / fo * 100;
        end
        fo_err_SD(kk) = mean(tmp_skip);
    end
    semilogy(snr_list, fo_err_SD, 'o-', 'linewidth', 2);
    grid on;hold on;
    drawnow
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('input SNR (dB)');
ylabel('SD of relative error (%)');
legend(wintype_list, 'location', 'northeast');
title(['fo:' num2str(fo) ' mag:' num2str(mag_base) ' dsOn:' num2str(dsOn)]);
outFname = ['winSNRFoErrTest' num2str(dsOn) 'fo' num2str(fo) 'Hz.eps'];
print('-depsc',  outFname);

%% Integration time and fo error

fs = 44100;
fo = 240;
fl = fo * 2 ^ (-1 / 2);
fh = fo * 2 ^ (1 / 2);
fo_list = 55 * 2 .^ (0:1/2:3);
duration = 100;
ch_in_oct = 12;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm'};%, 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
integration_time_list = [2 5 10 15 20 30 40];
iteration = 1;
snr = 30;
fo_err_SD = zeros(length(integration_time_list), 1);
fo_raw_SD = zeros(length(fo_list), 1);
tt = (1:round(duration * fs))' / fs;
figure;
for jk = 1:length(fo_list)
    fo = fo_list(jk);
    fl = fo * 2 ^ (-1 / 2);
    fh = fo * 2 ^ (1 / 2);
    xp = tt * 0;
    nto = round(fs / fo);
    xp(1:nto:end) = 1;
    for ii = 1:length(wintype_list)
        mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
        x = randn(round(duration * fs), 1);
        x = xp + randn(round(duration * fs), 1) * std(xp) * 10 ^ (-snr / 20);
        for kk = 1:length(integration_time_list)
            tmp_default = zeros(iteration, 1);
            tmp_skip = zeros(iteration, 1);
            initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
                ch_in_oct, dsOn, mag, wintype_list{ii}, integration_time_list(kk));
            for jj = 1:iteration
                out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
                tx = out_skip.time_axis_wavelet;
                tmp_skip(jj) = std(out_skip.fixed_points_freq(tx > 0.8 & tx < duration - 0.8, 1)) / fo * 100;
                foSeq = out_skip.fixed_points_freq(tx > 0.8 & tx < duration - 0.8, 1);
                foSeq = foSeq - mean(foSeq);
            end
            fo_err_SD(kk) = mean(tmp_skip);
        end
        semilogx(integration_time_list, fo_err_SD, 'o-', 'linewidth', 2);
        grid on;hold on;
        drawnow
    end
    fo_raw_SD(jk) = std(out_skip.inst_freq_map(tx > 0.8 & tx < duration - 0.8, 7)) / fo * 100;
end
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('integration length (ms)');
ylabel('SD of relative error (%)');
legend(num2str(fo_list', '%6.2f'), 'location', 'southwest');
title([' mag:' num2str(mag_base) ' dsOn:' num2str(dsOn) ' ' wintype_list{1}]);
outFname = ['winIntFoErrTest' num2str(dsOn) 'z.eps'];
print('-depsc',  outFname);

%% Integration time and fo traj power spec

fs = 44100;
fo_list = 55 * 2 .^ (0:1/2:3);
duration = 100;
ch_in_oct = 12;
dsOn = 1;
mag_base = 1.05;
wintype_list = {'sixterm'};%, 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
mag_duration_list = [0.4444 0.2831 0.3050 0.3566 0.4020 0.3619 0.3676];
integration_time_list = [2 5 10 15 20 30 40];
iteration = 1;
snr_list = 0:10:80;
snr = 30;
fo_err_SD = zeros(length(integration_time_list), 1);
tt = (1:round(duration * fs))' / fs;
fftl = 32768;
figure;
for jk = 5%1:length(fo_list)
    fo = fo_list(jk);
    fl = fo * 2 ^ (-1 / 2);
    fh = fo * 2 ^ (1 / 2);
    xp = tt * 0;
    nto = round(fs / fo);
    xp(1:nto:end) = 1;
    for ii = 1:length(wintype_list)
        mag = mag_base / mag_duration_list(ii) * mag_duration_list(1);
        x = randn(round(duration * fs), 1);
        x = xp + randn(round(duration * fs), 1) * std(xp) * 10 ^ (-snr / 20);
        for kk = 1:length(integration_time_list)
            tmp_raw = zeros(iteration, 1);
            tmp_skip = zeros(iteration, 1);
            initialoutput = sourceAttributesAnalysis(x, fs, [1 length(x)], fl, fh, ...
                ch_in_oct, dsOn, mag, wintype_list{ii}, integration_time_list(kk));
            out_skip = sourceAttributesAnalysis(x, fs, [1 length(x)], initialoutput);
            tx = out_skip.time_axis_wavelet;
            tmp_skip(jj) = std(out_skip.fixed_points_freq(tx > 0.8 & tx < duration - 0.8, 1)) / fo * 100;
            foSeq = out_skip.fixed_points_freq(tx > 0.8 & tx < duration - 0.8, 1);
            foSeq = (foSeq - mean(foSeq)) / fo;
            fds = fs / out_skip.downSamplinbgRate;
            sgramF = stftSpectrogramStructure(foSeq,fds,2000,100,'nuttallwin');
            fx = sgramF.frequencyAxis;
            semilogx(fx, 10*log10(mean(sgramF.rawSpectrogram,2)), 'linewidth', 2);
            grid on;hold on;
            drawnow
        end
    end
end
fo_raw_seq = out_skip.inst_freq_map(tx > 0.8 & tx < duration - 0.8, 7) / fo;
fo_raw_seq = fo_raw_seq - mean(fo_raw_seq);
sgramF = stftSpectrogramStructure(fo_raw_seq,fds,2000,100,'nuttallwin');
fx = sgramF.frequencyAxis;
semilogx(fx, 10*log10(mean(sgramF.rawSpectrogram,2)), 'k', 'linewidth', 3);
set(gca, 'fontsize', 14, 'linewidth', 2);
xlabel('modulation frequency (Hz)');
ylabel('relative power (dB re. fo)');
legend(num2str(integration_time_list'), 'location', 'southwest');
title([' fo:' num2str(fo, '%5.2f') ' mag:' num2str(mag_base) ' dsOn:' num2str(dsOn) ' ' wintype_list{1} ' SNR:' num2str(snr)]);
axis([1 200 -70 -5]);
outFname = ['winIntFoTrjSpec' num2str(dsOn) 'fo' num2str(fo, '%04.0f') 'Hz.eps'];
print('-depsc',  outFname);

%% CMU ARCTIC
% 09/Jan./2019
% 20/Jan./2019

baseDir = '/Users/kawahara/Music/CMU_arctic/';
talkerDir = {'cmu_us_bdl_arctic','cmu_us_jmk_arctic','cmu_us_slt_arctic'};

wintype_list = {'sixterm', 'hanning', 'hamming', 'blackman', 'nuttall12','kaiser','dpss'};
wtypeId = 1;
mag_effective = 1.05;
channels_oct = 12;
dsOn = 1;
fl = 55;
fh = 1200;
integration_time = 15;

for sourceID = 1:2
    switch sourceID
        case 1
            sourceType = 'Speech';
        case 2
            sourceType = 'EGG';
    end
    figure;
    for tlkrID = 1:length(talkerDir)
        origDir = [baseDir talkerDir{tlkrID} '/orig'];
        dirContents = dir([origDir '/arc*.wav']);
        selected_measure_agg = [];
        for fileId = 1:50%length(dirContents)
            pathname = [origDir '/' dirContents(fileId).name];
            disp(pathname);
            [x, fs] = audioread(pathname);
            output = sourceAttributesAnalysis(x(:, sourceID), fs, [1 length(x)], fl, fh, ...
                channels_oct, dsOn, mag_effective, wintype_list{wtypeId}, integration_time);
            fc_list = output.wvltStrDs.fc_list;
            tmp_measure = output.fixed_points_measure(:, 1:4);
            tmp_measure = tmp_measure(:);
            selected_measure = tmp_measure(~isnan(tmp_measure) & tmp_measure ~= -Inf & tmp_measure ~= Inf);
            selected_measure_agg = [selected_measure_agg; selected_measure(:)];
        end
        sorted_measure = sort(selected_measure_agg);
        sorted_measure = sorted_measure(diff(sorted_measure([1 1:end])) ~= 0);
        prob = (1:length(sorted_measure)) / length(sorted_measure);
        %semilogy(20 * log10(sorted_measure), prob, 'linewidth', 2);
        minM = sorted_measure(1);
        maxM = sorted_measure(end);
        half_bin = 1;
        centerM = minM + half_bin:maxM - half_bin;
        probU = interp1(sorted_measure, prob, centerM + half_bin, 'linear', 'extrap');
        probL = interp1(sorted_measure, prob, centerM - half_bin, 'linear', 'extrap');
        semilogy(centerM, probU-probL);grid on;
        grid on;
        hold all
        drawnow
    end
    set(gca, 'fontsize', 14, 'linewidth', 2)
    xlabel('measure (dB)');
    ylabel('probability in 2 dB bin');
    legend('bdl', 'jmk', 'slt')
    title(['CMU ARCTIC ' sourceType ' tI:' num2str(integration_time) ' ms'], 'interpreter', 'none');
    foutName = ['CMUarctic' sourceType num2str(integration_time) 'msTmp.eps'];
    print('-depsc',  foutName);
end
