%% Test script for basic perfromance of wavelet analysis function
% by Hideki Kawahara
% 16/Dec./2018

fs_list = [4000 6000 8000 11025 12000 16000 22050 32000 44100 48000];
for kk = 1:length(fs_list)
    fs = fs_list(kk);
    low_frequency = 55;
    high_freuency = 1000;
    stretching_factor = 1;
    channels_in_octave = 12;
    fftl = 8192;
    
    wvltStr = designCos6Wavelet(fs, low_frequency, high_freuency, ...
        fftl, stretching_factor, channels_in_octave);
    
    disp(['Maximum length of wavelet:' num2str(length(wvltStr.wvlt(1).w))]);
    %%  Dependency on signal length
    
    duration_list = [0.2 0.4 0.8 1.6 3.2 6.4 12.8 25.6 51.2 102.4];
    el_total = zeros(length(duration_list), 1);
    el_filtering = zeros(length(duration_list), 1);
    el_post = zeros(length(duration_list), 1);
    disp(['fs: ' num2str(fs)]);
    disp('duration total_el filtering_el post_el');
    for ii = 1:length(duration_list)
        tt = (0:1 / fs:duration_list(ii))';
        x = sin(2 * pi * 440 * tt) + 0.2 * sin(2 * pi * 440 * tt) .^ 3 ...
            + 0.25 * sin(2 * pi * 440 * tt) .^ 5 ...
            - 0.3 * sin(2 * pi * 440 * tt) .^ 4 ...
            + 0.4 * sin(2 * pi * 440 * tt) .^ 2;
        output = waveletSourceAnalyzer(x, fs, wvltStr);
        el_total(ii) = output.elapsedTime;
        el_filtering(ii) = output.elapsedTimeForFiltering;
        el_post(ii) = output.elapsedTimeForPostProcess;
        disp(num2str([duration_list(ii) el_total(ii) el_filtering(ii) el_post(ii)]));
    end
    %%
    figure;
    loglog(duration_list, el_total, 'o-');
    grid on; hold all;
    loglog(duration_list, el_filtering, 'o-');
    loglog(duration_list, el_post, 'o-');
    set(gca, 'fontsize', 14);
    legend('total', 'filtering', 'post processing', 'location', 'southeast');
    xlabel('signal duration (s)')
    ylabel('elapsed time (s)')
    title(['fs: ' num2str(fs)]);
    drawnow
end
