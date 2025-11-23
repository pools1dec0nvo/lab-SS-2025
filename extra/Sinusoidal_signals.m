samplerate = 48e3;

%Imagem_1
t= (0:1/samplerate:2)';
x = cos(2*pi*10*t);
plot(t,x) ; grid on;

%Imagem_2
t2= (0:1/samplerate:1)';
x2 = cos(2*pi*1000*t2);
plot(t2,x2) ; grid on;
soundsc(x2,samplerate); % Dá som