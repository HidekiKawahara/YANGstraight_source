function output = sourceInformationAnalysis(x, fs, varargin)
% This function analysis source information using an analytic signal
% with the six-term cosine series envelope
%
% sourceInformationAnalysis(x, fs)
% sourceInformationAnalysis(x, fs, range)
% sourceInformationAnalysis(x, fs, range, outputStruct) % After initialization
% sourceInformationAnalysis(x, fs, range, low_frequency, high_freuency)
% sourceInformationAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave);
% sourceInformationAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn)
% sourceInformationAnalysis(x, fs, range, low_frequency, high_freuency, ...
%     channels_in_octave, dsOn, stretching_factor)
%
% Arguments
%   x       : input signal. One column vector
%   fs      : sampling frequency (Hz)
%   range   : index of the start and end point of the input
%   outputStruct : output of the structure variable generated in the first
%                  run. Using this skips initialization
%   low_frequency : Lower limit of periodicity check (Hz)
%   high_freuency : Higher limit of periodicity check (Hz)
%   channels_in_octave  : Number of filters in each octave
%   dsOn : Switch for internal downsampling. 1:downsample, 0:no downsampling
%   stretching_factor : Envelope stretching factor: 1: isometric in the
%                       time-frequency plane. Larger the longer in time
%
% Output
%   output  : structure variable with the following field
%     estPeriod   : indicator of periodicity 1:periodic 0:random (sd hoc!)
%     wvltStrDs   : structure consisting of used analytic signal bank
%     rawWavelet  : Each column vector consits of each filter output
%     downSamplinbgRate : A mnumber representing downsampling ratio
%     fftlds   : FFT buffer length for debug
%     half_aaf_length : for debug
%     w_aaf           : for debug
%     time_axis_wavelet : time axis for wavelet analysis results
%     signal_time_axis  : time axis for input signal
%     gd_dev_map        : output deviation based on group delay
%     dgd_dev_map       : output deviation based on diffed group delay
%     mix_rev_measure   : mixed periodicity meeasure based on gd and dgd
%     fixed_points_freq : frequency of fo candidates
%     fixed_points_measure : mixed measure of fo candidates
%     fixed_points_amp  : signal amplitude of fo condidates
%     elapsedTime       : total elapsed time (s)

% Designed and implemented by Hideki Kawahara
%
%Copyright 2018 Hideki Kawahara
%
%Licensed under the Apache License, Version 2.0 (the "License");
%you may not use this file except in compliance with the License.
%You may obtain a copy of the License at
%
%    http://www.apache.org/licenses/LICENSE-2.0
%
%Unless required by applicable law or agreed to in writing, software
%distributed under the License is distributed on an "AS IS" BASIS,
%WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%See the License for the specific language governing permissions and
%limitations under the License.

start_tic = tic;
output = struct;
output.narg = nargin;
fftl = 8192;
low_frequency = 55;
high_freuency = 1200;
stretching_factor = 1.0;
channels_in_octave = 24;
isInitialized = 0;
dsOn = 1;
switch nargin
    case {5, 6, 7, 8}
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
        low_frequency = varargin{2};
        high_freuency = varargin{3};
    case 4
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
        outputStruct = varargin{2};
        if ~isstruct(outputStruct)
            help sourceInformationAnalysis
            disp('The last argument should be an initialized wavelet structure');
        else
            isInitialized = 1;
        end
end
switch nargin
    case 2
        ixrange = zeros(2, 1);
        ixrange(1) = 1;
        ixrange(2) = length(x);
    case 3
        ixrange = varargin{1};
        if length(ixrange) ~= 2
            help sourceInformationAnalysis
            disp('Range is the start and the end index of the data');
        end
    case 4
    case 5
    case 6
        channels_in_octave = varargin{4};
    case 7
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
    case 8
        channels_in_octave = varargin{4};
        dsOn = varargin{5};
        stretching_factor = varargin{6};
    otherwise
        help sourceInformationAnalysis
end
ixrange = round(ixrange);
x_trim = x(ixrange(1):ixrange(2));
if ~isInitialized
    downSamplinbgRate = 1;
    if dsOn
        downSamplinbgRate = max(1, floor (fs / (high_freuency * 4)));
    end
    fsd = fs / downSamplinbgRate;
    half_aaf_length = max(1, downSamplinbgRate * 2 - 1);
    w_aaf = hanning(round(2 * half_aaf_length + 1));
    fftlds = fftl / 2;
    wvltStrDs = designCos6Wavelet(fsd, low_frequency, high_freuency, ...
        fftlds, stretching_factor, channels_in_octave);
else
    downSamplinbgRate = outputStruct.downSamplinbgRate;
    wvltStrDs = outputStruct.wvltStrDs;
    fftlds = outputStruct.fftlds;
    half_aaf_length = outputStruct.half_aaf_length;
    w_aaf = outputStruct.w_aaf;
    fsd = fs / downSamplinbgRate;
end
xd = fftfilt(w_aaf/sum(w_aaf), x_trim);
xd = xd(half_aaf_length:downSamplinbgRate:end);

n_channels = length(wvltStrDs.fc_list);
fc_list = wvltStrDs.fc_list;

outputDs = waveletSourceAnalyzer(xd, fsd, wvltStrDs);
%%

n_samples = length(outputDs.rawWavelet(:, 1));
scanning_width = 40; % ms *** MAGIC NUMBER ****
sample_distance = round(scanning_width / 1000 * fsd);
selector = sample_distance:n_samples - sample_distance;
bias = floor(sample_distance / 2);
cum_gd_squared = (outputDs.cum_gd_squared + outputDs.cum_gd_squared(:, [1 1:end - 1])) / 2;
cum_diff_gd_squared = (outputDs.cum_diff_gd_squared + outputDs.cum_diff_gd_squared(:, [1 1:end - 1])) / 2;
gd_dev_map = (cum_gd_squared(selector - bias + sample_distance, :) ...
    - cum_gd_squared(selector - bias, :)) * diag(1 ./ fc_list.^2) / sample_distance;
dgd_dev_map = (cum_diff_gd_squared(selector - bias + sample_distance, :) ...
    - cum_diff_gd_squared(selector - bias, :)) / sample_distance;
mix_rev_measure = exp(-(log(dgd_dev_map) + log(gd_dev_map)) / 2);

if_dev_map = outputDs.inst_freq_map(selector, :) * diag(1 ./ fc_list);
log_if_dev_map = log(max(0.5, if_dev_map));
fixed_points = zeros(n_samples, n_channels);
fixed_points_freq = zeros(n_samples, n_channels);
fixed_points_measure = zeros(n_samples, n_channels);
fixed_points_amp = zeros(n_samples, n_channels);
amplitude = abs(outputDs.rawWavelet);
channel_idx = 1:n_channels;
for ii = selector
    buffer_id = ii - selector(1) + 1;
    tmp_fixp = channel_idx(log_if_dev_map(buffer_id, :) .* log_if_dev_map(buffer_id, [2:end, end]) < 0 & ...
        log_if_dev_map(buffer_id, :) > 0);
    if ~isempty(tmp_fixp)
        fixed_points(ii, 1:length(tmp_fixp)) = tmp_fixp(:)';
        for kk = 1:length(tmp_fixp)
            r = log_if_dev_map(buffer_id, tmp_fixp(kk)) / ...
                (log_if_dev_map(buffer_id, tmp_fixp(kk)) - log_if_dev_map(buffer_id, tmp_fixp(kk) + 1));
            fixed_points_freq(ii, kk) = ...
                fc_list(tmp_fixp(kk)) + (fc_list(tmp_fixp(kk) + 1) - fc_list(tmp_fixp(kk))) * r;
            fixed_points_measure(ii, kk) = (1 - r) * mix_rev_measure(buffer_id, max(1, tmp_fixp(kk) - 1)) ...
                + r * mix_rev_measure(buffer_id, tmp_fixp(kk));
            fixed_points_amp(ii, kk) = (1 - r) * amplitude(buffer_id, max(1, tmp_fixp(kk) - 1)) ...
                + r * amplitude(buffer_id, tmp_fixp(kk));
        end
        [sortedv, sortIdx] = sort(fixed_points_measure(ii, 1:length(tmp_fixp)), 'descend');
        fixed_points_measure(ii, 1:length(tmp_fixp)) = sortedv;
        fixed_points_freq(ii, 1:length(tmp_fixp)) = fixed_points_freq(ii, sortIdx);
        fixed_points_amp(ii, 1:length(tmp_fixp)) = fixed_points_amp(ii, sortIdx);
    end
end
minaa = -1.9731 / 1000;
biasEst = -28.2955;
rev_mes_dB = -10 * log10(fixed_points_measure);
estRandom = rev_mes_dB +  minaa * 1200 * log2(fixed_points_freq) - biasEst + 25;
output.estPeriod = sqrt(1 ./ (1 + 10 .^ (estRandom / 10)));
output.wvltStrDs = wvltStrDs;
output.rawWavelet = outputDs.rawWavelet;
output.downSamplinbgRate = downSamplinbgRate;
output.fftlds = fftlds;
output.half_aaf_length = half_aaf_length;
output.w_aaf = w_aaf;
output.time_axis_wavelet = (1:n_samples) / fsd + ixrange(1) / fs;% + wvltStrDs.wvlt(1).bias / fsd;
output.signal_time_axis = (ixrange(1):ixrange(2)) / fs;
output.gd_dev_map = gd_dev_map;
output.dgd_dev_map = dgd_dev_map;
output.mix_rev_measure = mix_rev_measure;
output.fixed_points_freq = fixed_points_freq;
output.fixed_points_measure = fixed_points_measure;
output.fixed_points_amp = fixed_points_amp;
output.elapsedTime = toc(start_tic);
end