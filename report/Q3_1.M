

samplerate = 16e3;
t = (-320:320)'/samplerate;           % time vector

w0 = 2*pi*80;                         % 80 Hz fundamental frequency
T  = 2*pi/w0;                         % period
T1 = T/6;                             % duty cycle parameter (gives 33% duty)


ak = @(k) sin(k*w0*T1)./(k*pi).*exp(-1i*k*w0*T1);
a0 = 2*T1/T;                        

N_vec = [1 2 3 5 10 20 50 100];       

figure; hold on; grid on
for N = N_vec
    k = -N:N;
    coeff = [ak(-N:-1) a0 ak(1:N)];

    % x_N(t) = sum_{k=-N}^N a_k * exp(j*k*w0*t)
    xN = exp(1i*w0*t*k) * coeff.';   

    plot(t, real(xN), 'LineWidth', 1);
end
xlabel('t (s)');
ylabel('$\tilde{x}_N(t)$', 'Interpreter', 'latex');
legend(string(N_vec), 'Location', 'best');
title('Superimposed Fourier approximations $\tilde{x}_N(t)$', 'Interpreter', 'latex');
