function output = designCos6Wavelet(fs, fl, fh, fftl, mag, channels_oct)
% output = designCos6Wavelet(fs, fl, fh, fftl, mag, channels_oct)
% Generate a set of complex wavelets using the six-term cosine series 
% envelope
%
% Arguments
%   fs    : sampling frequency (Hz)
%   fl    : lower bound of the carrier frequency (Hz)
%   fh    : upper bound of the carrier frequency (Hz)
%   fftl  : FFT buffer length (samples): 2^K, where K: positive integer
%           This parameter is used for showing gain no effects on
%           performance
%   mag   : envelope stretching factor: mag=1 places zero at 0 and 2*fc
%   channels_oct : number of wavelets in one octave
%
% Output
%   output   : structure with the following fields : nch: number of wavelets
%    input_parameters: [1×1 struct] copy of input argument
%             fc_list: [1×nch double] list of carrier frequencies (Hz)
%           centerIdx: Index of the center of each wavelet: IGNORE
%              f_axis: [1×fftl double] frequency axis of the composite filter
%            response: [nresp×1 double] impulse response of the composite filter
%                gain: [fftl×1 double] gain of the composite filter
%                wvlt: [1×nch struct] a set of wavelets
%  note: the composite filter is not used anymore
%
% Copyright 2018 Hideki Kawahara
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

input_parameters.sampling_frequency = fs;
input_parameters.lower_frequency = fl;
input_parameters.higher_frequency = fh;
input_parameters.fft_length = fftl;
input_parameters.stretching_factor = mag;
input_parameters.channels_per_octave = channels_oct;
fx = (0:fftl-1)/fftl*fs;
fc_list = fl * 2.0 .^ (0:1/channels_oct:log2(fh/fl));
fx(fx > fs / 2) = fx(fx > fs / 2) - fs;
ak = [0.2624710164;0.4265335164;0.2250165621;0.0726831633;0.0125124215;0.0007833203];
wvlt = struct;
for ii = 1:length(fc_list)
  fc = fc_list(ii);
  n_to = round(mag * 3 * fs / fc);
  tx_we = (-n_to:n_to)'/ n_to;
  tx = (-n_to:n_to)'/ fs;
  we = cos(tx_we * [0 1 2 3 4 5] * pi) * ak;
  wr = we .* cos(2 * pi * tx * fc);
  wi = we .* sin(2 * pi * tx * fc);
  wvlt(ii).w = (wr+1i*wi)/sum(we);
  wvlt(ii).bias = n_to;
  wvlt(ii).t_axis = tx;
end
wcmp = wvlt(1).w * 0;
centerIdx = wvlt(1).bias + 1;
for ii = 1:length(wvlt)
  tmpIdx = -wvlt(ii).bias:wvlt(ii).bias;
  wcmp(centerIdx + tmpIdx) = wcmp(centerIdx + tmpIdx) + wvlt(ii).w;
end
%--- calibration for linear phase filter
x1khz = exp(1i * 2 * pi * 1000 * (1:fs)' / fs);
y1khz = fftfilt(wcmp, x1khz);
cf = mean(abs(y1khz(3 * centerIdx:end - 3 * centerIdx)));
wcmp = wcmp / cf;
w_gain = 20*log10(abs(fft(wcmp, fftl)));
%---
output.input_parameters = input_parameters;
output.fc_list = fc_list;
output.centerIdx = wvlt(1).bias + 1;
output.f_axis = fx;
output.response = wcmp;
output.gain = w_gain;
output.wvlt = wvlt;
end