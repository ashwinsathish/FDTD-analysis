clc;
close all;
clear;

fs = 1e12; % sampling frequency
t = 0:1/fs:1e-8; % time base

Ts = 4e-11; % pulse width
t0 = 2*Ts; % pulse delay
N_steps = 1e7; % maximum iteration steps
fc = 7.5e9; % center frequency


pulse1 = sin(2*pi*fc*t); 
pulse2 = exp(-((t-t0)/Ts).^2); % Gaussian pulse

% Composite modulated Gaussian pulse
pulse = pulse1.*pulse2;

% Plotting the pulse
figure;
subplot(2,1,1)
plot(t,pulse);
title('Modulated Gaussian Pulse');
xlabel('Time(s)');
ylabel('Amplitude');
xlim([0 0.1e-8])

% FFT of the pulse
X = fftshift(fft(pulse,N_steps));
f = fs*(-N_steps/2:N_steps/2-1)/N_steps; %Frequency Vector

% Finding the 3 dB bandwidth
mag_dB = 20*log10(abs(X));
max_mag = max(mag_dB);
band_3dB = max_mag - 3;

% indices for starting frequency of the 3 dB bandwidth
idx_posfreq = find(f > 0);
idx_start = idx_posfreq(find(mag_dB(idx_posfreq) >= band_3dB, 1, 'first'));
idx_end = find(mag_dB >= band_3dB, 1, 'last');

% starting and ending frequencies of the 3 dB bandwidth
freq_start = f(idx_start);
freq_end = f(idx_end);

min_mag = min(mag_dB(idx_posfreq)); % min magnitude in the positive frequency range

% Normalizing between -1 and 0
mag_normalized = (mag_dB - max_mag);
band_3dB_normalized = (band_3dB - max_mag);

% Ploting
subplot(2,1,2)
plot(f/1e9, mag_normalized);
hold on;

title('Normalized Magnitude of FFT');
xlabel('Frequency (GHz)');
ylabel('Normalized Magnitude');
xlim([0 30]);
ylim([-3 0]);

% Display the starting and ending frequencies of the 3 dB bandwidth
fprintf('Starting freq of 3 dB bandwidth: %.2f GHz\n', freq_start/1e9);
fprintf('Ending freq of 3 dB bandwidth: %.2f GHz\n', freq_end/1e9);
