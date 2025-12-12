

samplerate = 16e3;
t = (-320:320)'/samplerate;           % time vector

w0 = 2*pi*80;                         % freq fundamental
T  = 2*pi/w0;                         % period
T1 = T/6;                             

ak = @(k) sin(k*w0*T1)./(k*pi).*exp(-1i*k*w0*T1);
a0 = 2*T1/T;                         

N_vec = [1 2 3 5 10 20 50 100];       

figure;
num_plots = length(N_vec);

for idx = 1:num_plots
    N = N_vec(idx);
    k = -N:N;

    % vect coeficientes
    coeff = [ak(-N:-1) a0 ak(1:N)];

    xN = exp(1i*w0*t*k) * coeff.';
    xN = real(xN);

    [XN, f] = cftransform(xN, find(t == 0), 1e4, samplerate);

    subplot(ceil(num_plots/2), 2, idx);
    plot(f, abs(XN), 'LineWidth', 1);
    grid on;
    xlabel('f (Hz)');
    ylabel('$|X_N(f)|$', 'Interpreter', 'latex');
    title(['N = ' num2str(N)]);
    xlim([-1000 1000]);  % zoom (opcional)
end
sgtitle('Magnitude spectra of $\tilde{x}_N(t)$ for different N', 'Interpreter', 'latex');

figure; hold on; grid on;
colors = lines(num_plots);
for idx = 1:num_plots
    N = N_vec(idx);
    k = -N:N;

    coeff = [ak(-N:-1) a0 ak(1:N)];
    xN = real(exp(1i*w0*t*k) * coeff.');

    [XN, f] = cftransform(xN, find(t == 0), 1e4, samplerate);

    plot(f, abs(XN), 'Color', colors(idx,:), 'LineWidth', 1);
end
xlabel('f (Hz)');
ylabel('$|X_N(f)|$', 'Interpreter', 'latex');
legend(string(N_vec), 'Location', 'best');
title('Superimposed magnitude spectra of $\tilde{x}_N(t)$', 'Interpreter', 'latex');
xlim([-500 500]);
