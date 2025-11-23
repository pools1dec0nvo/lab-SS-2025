samplerate = 48e3;
t = (0:1/samplerate:5)';

x = cos(2*pi*100*t);
soundsc(x, samplerate)