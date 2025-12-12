

samplerate = 16e3;
T = 0.4;                              % periodo
w0 = 2*pi/T;                         


t = (-5*T : 1/samplerate : 5*T)';

A_pp = 3;                             % peak-to-peak amplitude
A = A_pp/2;                           % amplitude of centered sawtooth = 1.5
DC = 2.5;                             % DC offset

ak = @(k) -1i*A./(k*pi);
a0 = DC;

N_vec = [1 2 3 5 10 20 50 100];

figure; hold on; grid on;
for N = N_vec
    k = -N:N;
    % vetor coeficientes [a_{-N}, ..., a_{-1}, a_0, a_1, ..., a_N]
    coeff = [ak(-N:-1) a0 ak(1:N)];
    xN = exp(1i*w0*t*k) * coeff.';

    plot(t, real(xN), 'LineWidth', 1);
end

xlabel('t (s)');
ylabel('$\tilde{x}_N(t)$', 'Interpreter', 'latex');
legend(string(N_vec), 'Location', 'best');
title('Fourier para $\tilde{x}_N(t)$ - Sawtooth Wave', 'Interpreter', 'latex');

fprintf('w0 = 2*pi/T = 2*pi/0.4 = %.4f rad/s\n', w0);
fprintf('a0 = %.4f\n', a0);
fprintf('a_{-1} = %.4f + %.4fi\n', real(ak(-1)), imag(ak(-1)));
fprintf('a_3 = %.4f + %.4fi\n', real(ak(3)), imag(ak(3)));

fprintf('\n|a_{-1}| = %.4f, angle(a_{-1}) = %.4f rad = %.2f deg\n', ...
    abs(ak(-1)), angle(ak(-1)), angle(ak(-1))*180/pi);
fprintf('|a_3| = %.4f, angle(a_3) = %.4f rad = %.2f deg\n', ...
    abs(ak(3)), angle(ak(3)), angle(ak(3))*180/pi);
