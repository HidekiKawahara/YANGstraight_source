function output = waveletAttributesAnalyzer(x_trim, fs, wvltStr)
% output = waveletAttributesAnalyzer(x_trim, fs, wvltStr)
%
% Source information analysis using an analyitic signal with the envelope
% defined by the six-term cosine seris proposed in ref.[1]
% Note, other analitic signals with different envelope definition can be
% used
% 
% Arguments
%   x_trim  : input signal. One column vector
%   fs      : samplinlg frequency (Hz)
%   wvltStr : structure variable consisting of the analytic signals
%             Fileds used:
%     fc_list  : a vector of carrier frequencies (Hz)
%     wvlt     : a set of structure variables. 
%                Each variable defines an analytic signal. Two fields
%       w      : complex column vector defining an analytic signal
%       bias   : half length of the vector (2 * bias + 1) is the length
%
% Output
%   output : structure variabl with the following fields
%     rawWavelet                : Each column has filtered output
%     inst_freq_map             : Each column has sample-wise inst. freq.
%     inst_freq_map_center      : same i.f. with symmetric calculation
%     group_delay_map           : Each column has sample-wise group delay
%     group_delay_map_center    : same g.d. with symmetric calculation
%     elapsedTimeForFiltering   : Elapsed time for filtering (s)
%     elapsedTimeForPostProcess : Elapsed time for postprocess (s)
%     elapsedTime               : Total elapsed time (s)

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
n_channels = length(wvltStr.fc_list);
n_data = length(x_trim);
rawWavelet = zeros(n_data, n_channels);
n_buffer = n_data;
x_trim = [zeros(wvltStr.wvlt(1).bias, 1);x_trim;zeros(wvltStr.wvlt(1).bias, 1)];

tic
buffer_index = 1:n_data;
for ii = 1:n_channels
    y = fftfilt(wvltStr.wvlt(ii).w, x_trim);
    rawWavelet(:, ii) = y(wvltStr.wvlt(ii).bias + wvltStr.wvlt(1).bias + buffer_index);
end
elapsedTimeForFiltering = toc;
%% 
tic
fc_list = wvltStr.fc_list;
amp_map = abs(rawWavelet);
amp_squared_map = amp_map([2:n_buffer n_buffer], :) .* amp_map;
inst_freq_map = angle(rawWavelet([2:n_buffer n_buffer], :) ./ rawWavelet) * fs / 2 / pi;
inst_freq_map(end, :) = inst_freq_map(end - 1, :);
%inst_freq_map_center = (inst_freq_map .* amp_squared_map + ...
%    inst_freq_map([1 1:end-1], :) .* amp_squared_map([1 1:end-1], :)) ./ ...
%    (amp_squared_map + amp_squared_map([1 1:end-1], :));
fc_list_extend = [fc_list fc_list(end) * fc_list(end) / fc_list(end - 1)];
freq_step = diff(fc_list_extend);
group_delay_map = -angle(rawWavelet(:, [2:end end]) ./ rawWavelet) * diag(1.0 ./ freq_step) / 2 / pi;
group_delay_map(:, end) = group_delay_map(:, end - 1);
%group_delay_map_center = (group_delay_map .* amp_map .* amp_map(:, [2:end end]) ...
%    + group_delay_map(:, [1 1:end-1]) .* amp_map(:, [1 1:end-1]) .* amp_map) ./ ...
%    (amp_map .* amp_map(:, [2:end end]) + amp_map(:, [1 1:end-1]) .* amp_map);
elapsedTimeForPostProcess = toc;
output.rawWavelet = rawWavelet;
output.amp_squared_map = amp_squared_map;
output.inst_freq_map = inst_freq_map;
%output.inst_freq_map_center = inst_freq_map_center;
output.group_delay_map = group_delay_map;
%output.group_delay_map_center = group_delay_map_center;
output.elapsedTimeForFiltering = elapsedTimeForFiltering;
output.elapsedTimeForPostProcess = elapsedTimeForPostProcess;
output.elapsedTime = toc(start_tic);
end