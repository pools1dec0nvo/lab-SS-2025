function x = seqsin(fs,varargin)
% Generate a sequence of sinusoids with specified frequencies and durations.
% The arguments come in pairs: frequency and duration of each sinusoid.
% A frequency of zero corresponds to a period of silence.

if mod(numel(varargin),2) ~= 0
    error('seqsin: Number of (frequency, duration) arguments must be even')
end


n_in = round(.001*fs);   % Duration, in samples, of the initial taper
n_out = round(.01 *fs);   % Duration, in samples, of the final taper

fade_in = 0.5*(1 + cos(pi/n_in*(-n_in:-1)'));
fade_out = 0.5*(1 + cos(pi/n_out*(-1:-1:-n_out)'));

x = [];
for i = 1:2:numel(varargin)
    freq = varargin{i};
    nsamples = round(varargin{i+1}*fs);

    x = cat(1,x,sin(2*pi*freq/fs*(0:nsamples-1)').* ...
                [fade_in; ones(nsamples-n_in-n_out,1); fade_out]);
end

end