clc;
clear;
close all;

%----- Medium ----------
length_x = 0.1;
length_y = 0.1;
eps0=8.854e-12;
meu0=4*pi*1e-7;
real_epsr_conc=4; % relative permittivity of concrete = 4 - j0.06
img_epsr_conc=0.06;
real_epsr_rust=10; % relative permittivity of rust = 10 - j0.86
img_epsr_rust=0.86;
eps_conc=eps0*real_epsr_conc;
eps_rust=eps0*real_epsr_rust;
eps_steel=eps0; 
meur=1;
meu=meu0*meur;

%---- Signal -----------
c=3e8;
v=c/sqrt(real_epsr_conc*meur);
freq=12e9;
lamda=v/freq;
lamda_min = lamda/sqrt(real_epsr_conc*meur);

%------ Conductivity values ---------
cond_conc=2*pi*freq*eps0*img_epsr_conc;
cond_rust=2*pi*freq*eps0*img_epsr_rust;

%------ Cell length and time step---------
dx=lamda_min/25;
dy=dx;
dt=1/(v*sqrt(1/(dx^2) + 1/(dy^2)));

%------- Rust conditions ---------
bar_exists = 0;
corr = 0;

%------- Corrosion base values ---------
rs_paper = [12.50,12.18,11.85,11.52,11.18,10.82,10.45]/1000;
rr_paper = [12.50,13.40,14.25,15.05,15.81,16.53,17.23]/1000;
if corr~=0
    for i=1:6
        if corr==5*i
            rs = ceil(rs_paper(i+1)/(2*dx));
            rr = ceil((rr_paper(i+1)-rs_paper(i+1))/(2*dx));    
        end
    end
else
    rs = ceil(rs_paper(1)/(2*dx));
    rr = 0;
end

N_cells_x = round(length_x/dx);
N_cells_y = round(length_y/dy);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);

[X,Y] = meshgrid(1:N_cells_x, 1:N_cells_y);

% Initialize space with concrete
space = ones(N_cells_x, N_cells_y)*eps_conc;

%------- Rebar radius and conditions --------
if(bar_exists==1)
    rebarMask = ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) <= rs^2;
    rustMask = ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) <= (rs+rr)^2 & ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) > rs^2;

    for i = 1:N_cells_x
        for j = 1:N_cells_y
            if rebarMask(i, j)
                space(i, j) = eps_steel;
            elseif rustMask(i, j)
                space(i, j) = eps_rust;
            end
        end
    end
end

% Visualize the space
figure;
imagesc(space);
xlabel('y');
ylabel('x');
if bar_exists == 1
    title(sprintf('Space visualization (%.0f%% corrosion)', corr));
else
    title('Space visualization (no rebar)');
end
axis equal tight;

%------ Multiplying constants --------
const_e1=zeros(N_cells_x,N_cells_y);
const_e2=zeros(N_cells_x,N_cells_y);
const_e3=zeros(N_cells_x,N_cells_y);

%------- Steel, Rust and Concrete region ---------
if(bar_exists==1)
    const_e1(~rebarMask & ~rustMask) = ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
    const_e2(~rebarMask & ~rustMask) = (2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
    const_e3(~rebarMask & ~rustMask) = (2*dt/dy)/((2*eps_conc)+(cond_conc*dt));

    const_e1(rebarMask) = 0;
    const_e2(rebarMask) = 0;
    const_e3(rebarMask) = 0;

    const_e1(rustMask) = ((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
    const_e2(rustMask) = (2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
    const_e3(rustMask) = (2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
else 
    const_e1 = ones(N_cells_x, N_cells_y) * ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
    const_e2 = ones(N_cells_x, N_cells_y) * (2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
    const_e3 = ones(N_cells_x, N_cells_y) * (2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
end

const_hx=dt/(meu*dx);
const_hy=dt/(meu*dy);

const_abc1=(v*dt-dx)/(v*dt+dx);
const_abc2=(2*dx)/(v*dt+dx);
const_abc3=(dx*(v*dt)^2)/(2*(dy^2)*(v*dt+dx));

%------- Initialise E and H arrays ----------
ez=zeros(N_cells_x,N_cells_y);
hx=zeros(N_cells_x,N_cells_y);
hy=zeros(N_cells_x,N_cells_y);

%------- Tangential component around rebar ----------
gz = ones(N_cells_x, N_cells_y);
if bar_exists==1
    for i = 1:N_cells_x
        for j = 1:N_cells_y
            if (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 <= (rs+1)^2 && (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 >= (rs-1)^2
                gz(i,j) = 0;
            end
        end
    end
end

%------ Gaussian pusle ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=2000; % maximum iteration steps

vin = zeros(N_steps,1);
vout = zeros(N_steps,N_cells_y);
vout_avg = zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

ez_past = zeros(N_cells_x,N_cells_y,N_steps);

fc = 7.5e9; % Center frequency
bw = 9e9;

%************** Iteration loop ****************

choice = input("Do you want to continue with the simulation (0-no, 1-yes)? ");

while choice==1
   
for n=1:N_steps         
    time=(n-1)*dt;  
    time_stamp(n,1)=time;
        
%------- Launching the signal ---------
    pulse1 = sin(2*pi*fc*time); % sine pulse
    pulse2 = exp(-((time-t0)/Ts)^2); % Gaussian pulse
    pulse = pulse2.*pulse1;

    ez(1, Mid_cell_y-10:Mid_cell_y+10) = ez(1, Mid_cell_y-10:Mid_cell_y+10) + pulse;
    vin(n) = ez(1, Mid_cell_y);

     %--------- compute H ----------
     i = 1:N_cells_x-1;
     j = 1:N_cells_y-1;
     hx(i,j) = hx(i,j) - const_hy*(ez(i,j+1)-ez(i,j));
     hy(i,j) = hy(i,j) + const_hx*(ez(i+1,j)-ez(i,j));
 

     %--------- compute E ----------
     i = 2:N_cells_x;
     j = 2:N_cells_y;
     ez(i,j) = gz(i,j).*(const_e1(i,j).*ez(i,j) + const_e2(i,j).*(hy(i,j)-hy(i-1,j)) - const_e3(i,j).*(hx(i,j)-hx(i,j-1)));

     %------------ 2nd order ABC ------------- 
     n=n+2;
     
     for j=2:N_cells_y-1
            
          ez(1,j) = - ez_past(2,j,n-2) + const_abc1*(ez(2,j) + ez_past(1,j,n-2)) + const_abc2*(ez_past(2,j,n-1) + ez_past(1,j,n-1)) + const_abc3*(ez_past(2,j+1,n-1) - 2*ez_past(1,j,n-1) + ez_past(2,j-1,n-1) + ez_past(1,j+1,n-1) - 2*ez_past(2,j,n-1) + ez_past(1,j-1,n-1));
          ez(N_cells_x,j) = - ez_past(N_cells_x-1,j,n-2) + const_abc1*(ez(N_cells_x-1,j) + ez_past(N_cells_x,j,n-2)) + const_abc2*(ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j,n-1)) + const_abc3*(ez_past(N_cells_x-1,j+1,n-1) - 2*ez_past(N_cells_x,j,n-1) + ez_past(N_cells_x-1,j-1,n-1) + ez_past(N_cells_x,j+1,n-1) - 2*ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j-1,n-1));
       
     end 

     
%      for i=2:N_cells_x-1
%                         
%           ez(i,1) = - ez_past(i,2,n-2) + const_abc1*(ez(i,2) + ez_past(i,1,n-2)) + const_abc2*(ez_past(i,2,n-1) + ez_past(i,1,n-1)) + const_abc3*(ez_past(i+1,2,n-1) - 2*ez_past(i,1,n-1) + ez_past(i-1,2,n-1) + ez_past(i+1,1,n-1) - 2*ez_past(i,2,n-1) + ez_past(i-1,1,n-1)); % left boundary
%           ez(i,N_cells_y) = - ez_past(i,N_cells_y-1,n-2) + const_abc1*(ez(i,N_cells_y-1) + ez_past(i,N_cells_y,n-2)) + const_abc2*(ez_past(i,N_cells_y-1,n-1) + ez_past(i,N_cells_y,n-1)) + const_abc3*(ez_past(i+1,N_cells_y-1,n-1) - 2*ez_past(i,N_cells_y,n-1) + ez_past(i-1,N_cells_y-1,n-1) + ez_past(i+1,N_cells_y,n-1) - 2*ez_past(i,N_cells_y-1,n-1) + ez_past(i-1,N_cells_y,n-1));
% 
%      end

     ez_past(:,:,n) = ez;
     n=n-2;

     % ------- Recording output -------
     for i=1:N_cells_y
        vout(n,i)=ez(N_cells_x,i);
     end
     vout_avg = mean(vout,2);
%      for j=1:N_cells_y
%         vout(n,i)=ez(N_cells_x,i);
%      end
%      vout_avg_inputside = mean(vout_inputside,2);

     % ------- Plotting --------
     surf(X, Y, ez);
     xlabel('y');
     ylabel('x');
     zlabel('Ez');
     xlim([1 N_cells_x]);
     ylim([1 N_cells_y]);
     zlim([-1 1]);
     if bar_exists == 1
        title(sprintf('Ez(%.0f%% corrosion) at time step %d', corr, n));
     else
        title(sprintf('Ez(no rebar) at time step %d', n));
     end
     pause(0.01);
        
end

break
end

if (choice==1)
%----- FFT parameters ------
fs = 1e12; % sampling frequency
t = 0:1/fs:1e-8; % time base
FFT_steps = 1e7; % maximum iteration steps

% FFT of the pulse
X = fftshift(fft(vout_avg,FFT_steps));
f = fs*(-FFT_steps/2:FFT_steps/2-1)/FFT_steps; %Frequency Vector

% Finding the 3 dB bandwidth
mag_dB = 20*log10(abs(X));

% indices for starting frequency of the 3 dB bandwidth
idx_posfreq = find(f > 0);
[max_mag, max_mag_idx] = max(mag_dB(idx_posfreq));
band_3dB = max_mag - 3;
idx_start = idx_posfreq(find(mag_dB(idx_posfreq) >= band_3dB, 1, 'first'));
idx_end = find(mag_dB >= band_3dB, 1, 'last');
max_freq = f(idx_posfreq(max_mag_idx));

% starting and ending frequencies of the 3 dB bandwidth
freq_start = f(idx_start);
freq_end = f(idx_end);

min_mag = min(mag_dB(idx_posfreq)); % min magnitude in the positive frequency range

% Normalizing between -1 and 0
mag_normalized = (mag_dB - max_mag);
band_3dB_normalized = (band_3dB - max_mag);

% Plotting Magnitude Spectrum
figure;
subplot(2,1,1)
plot(f/1e9, mag_dB);
hold on;

if bar_exists == 0
    title('FFT (no rebar)');
else
    title(sprintf('FFT (%.0f%% corrosion)', corr));
end
xlabel('Frequency (GHz)');
ylabel('Magnitude (dB)');
xlim([0 30]);
%ylim([-3 0]);

% Plotting Phase Spectrum
phase_rad = angle(X);
phase_deg = phase_rad * (180/pi);
max_phase = phase_rad(idx_posfreq(max_mag_idx));

subplot(2,1,2)
plot(f/1e9, phase_deg);
title('Phase of FFT');
xlabel('Frequency (GHz)');
ylabel('Phase (Degrees)');
xlim([0 15]);

% Display the magnitude results
if bar_exists == 0
    fprintf('\nResults for no rebar:\n');
else
    fprintf('\nResults for %.2f%% corrosion\n', corr);
end

% Display the starting and ending frequencies of the 3 dB bandwidth
fprintf('Starting freq of 3 dB bandwidth: %.2f GHz\n', freq_start/1e9);
fprintf('Ending freq of 3 dB bandwidth: %.2f GHz\n', freq_end/1e9);

fprintf('\nFrequency at max magnitude: %.2f GHz\n', max_freq/1e9);
fprintf('Maximum magnitude: %.2f dB\n', max_mag);
fprintf('Phase at max magnitude: %.2f radians\n', max_phase);
end