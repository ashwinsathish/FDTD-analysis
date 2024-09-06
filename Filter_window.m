clc;
clear;
close all;

% Define parameters
Ts = 0.5e-9;    % pulse width
t0 = 2*Ts;      % pulse delay
N_steps = 1600; % maximum iteration steps
fc = 7.5e9;     % center frequency
bw = 9e9;       % bandwidth 

% Define time vector
t = linspace(0,10*Ts,N_steps);

% Define pulse
fm = bw/2; % modulation frequency
pulse1 = sin(2*pi*fm*t); % sine pulse
pulse2 = exp(-((t-t0)/Ts).^2); % Gaussian pulse
pulse = pulse2.*pulse1;

% Compute Fourier transform
f = linspace(-N_steps/(2*Ts), N_steps/(2*Ts), N_steps);
Y = fftshift(fft(pulse));

% Define filter
f1 = 3e9; % lower cutoff frequency
f2 = 12e9; % upper cutoff frequency
B = 2e9; % transition bandwidth
H = zeros(size(Y));

H(abs(f)>=f1-B & abs(f)<=f1) = 0.5*(1 - cos(pi*(abs(f(abs(f)>=f1-B & abs(f)<=f1))-f1+B)/B));
H(abs(f)>=f1 & abs(f)<=f2) = 1;
H(abs(f)>=f2 & abs(f)<=f2+B) = 0.5*(1 - cos(pi*(abs(f(abs(f)>=f2 & abs(f)<=f2+B))-f2+B)/B));
H(abs(f)>f2+B) = 0;

% Apply filter
Y_filt = H .* Y;

% Plot magnitude of filtered Fourier transform
subplot(3,1,1)
plot(t,pulse)
xlabel('Time (s)')
ylabel('Amplitude')
title('Pulse waveform')

subplot(3,1,2)
plot(f/1e9, abs(Y_filt))
xlabel('Frequency (GHz)')
ylabel('Magnitude')
xlim([1, 15])
ylim([0,0.4])
title('Fourier transform')

subplot(3,1,3)
semilogx(f/1e9, 20*log10(abs(Y_filt)))
xlabel('Frequency (GHz)')
ylabel('Magnitude (dB)')
xlim([1, 15])
ylim([-100,0])
title('Filtered Fourier transform (dB)')
grid on
