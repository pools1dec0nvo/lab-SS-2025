function [X,f] = cftransform(x,i0,n,fs,varargin)
% Evaluate the continuous-time Fourier transform of a sampled
% continuous-time signal.
% x - Vector of input samples
% i0 - Index of sample corresponding to time t = 0
% n — Number of evaluation points
% fs - Sample rate
% varargin - Optional arguments for freqz ('whole')
%
% X - (Sampled) continuous-time Fourier transform
% f - Frequency axis (Hz)

if isempty(i0), i0 = 1; end

[X,f] = freqz(x,1,n,fs,varargin{:});
X = (X/fs).*exp(1i*2*pi*f*(i0-1)/fs);

% Convert frequencies to range [-fs/2, fs/2]
f = fs/(2*pi)*angle(exp(1i*2*pi*f/fs));
end