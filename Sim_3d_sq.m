clc;
clear;
close all;

%----- Medium ----------
length_x = 0.01;
length_y = 0.01;
length_z = 0.01;
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
dx=lamda_min/20;
dy=dx;
dz=dx;
dt=1/(v*sqrt(1/(dx^2) + 1/(dy^2) + 1/(dz^2)));

N_cells_x = round(length_x/dx);
N_cells_y = round(length_y/dy);
N_cells_z = round(length_z/dz);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);
Mid_cell_z = round(N_cells_z/2);

[X,Y,Z] = meshgrid(1:N_cells_x, 1:N_cells_y, 1:N_cells_z);

%------ Multiplying constants --------
const_e1=zeros(N_cells_x,N_cells_y,N_cells_z);
const_ex=zeros(N_cells_x,N_cells_y,N_cells_z);
const_ey=zeros(N_cells_x,N_cells_y,N_cells_z);
const_ez=zeros(N_cells_x,N_cells_y,N_cells_z);

%------- Rust conditions ---------
bar_exists = 1;
corr = 15;

%------- Corrosion base values ---------
% rs_paper = [12.50,12.18,11.85,11.52,11.18,10.82,10.45]/3000;
% rr_paper = [12.50,13.40,14.25,15.05,15.81,16.53,17.23]/3000;
% 
% if corr~=0
%     for i=1:6
%         if corr==5*i
%             rs = ceil(rs_paper(i+1)/(2*dx));
%             rr = ceil((rr_paper(i+1)-rs_paper(i+1))/(2*dx));    
%         end
%     end
% else
%     rs = ceil(rs_paper(1)/(2*dx));
%     rr = 0;
% end

rs = 7; % 9,8,7,6,5
rr = 2; % 0,1,2,3,4


%------- Corrosion simulation---------
const_e1(1:N_cells_x,1:N_cells_y,1:N_cells_z)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_ex(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
const_ey(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
const_ez(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dz)/((2*eps_conc)+(cond_conc*dt));

if bar_exists==1
    
        % Steel region
        const_e1(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=1;
        const_ex(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dx)/eps_steel;
        const_ey(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dy)/eps_steel;
        const_ez(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dz)/eps_steel;

        % Left boundary rust
        const_e1(Mid_cell_x-rs-rr:Mid_cell_x-rs-1, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
        const_ex(Mid_cell_x-rs-rr:Mid_cell_x-rs-1, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
        const_ey(Mid_cell_x-rs-rr:Mid_cell_x-rs-1, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
        const_ez(Mid_cell_x-rs-rr:Mid_cell_x-rs-1, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

        % Right boundary rust
        const_e1(Mid_cell_x+rs+1:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
        const_ex(Mid_cell_x+rs+1:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
        const_ey(Mid_cell_x+rs+1:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
        const_ez(Mid_cell_x+rs+1:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y+rs+rr, :)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

        % Top boundary rust
        const_e1(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y-rs-1, :)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
        const_ex(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y-rs-1, :)=(2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
        const_ey(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y-rs-1, :)=(2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
        const_ez(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y-rs-rr:Mid_cell_y-rs-1, :)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

        % Bottom boundary rust
        const_e1(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y+rs+1:Mid_cell_y+rs+rr, :)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
        const_ex(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y+rs+1:Mid_cell_y+rs+rr, :)=(2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
        const_ey(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y+rs+1:Mid_cell_y+rs+rr, :)=(2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
        const_ez(Mid_cell_x-rs-rr:Mid_cell_x+rs+rr, Mid_cell_y+rs+1:Mid_cell_y+rs+rr, :)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

end

const_abc1z = (v*dt - dx)/(v*dt + dx);
const_abc1y = (v*dt - dy)/(v*dt + dy);
const_abc1x = (v*dt - dx)/(v*dt + dx);

const_hx=dt/(meu*dx);
const_hy=dt/(meu*dy);
const_hz=dt/(meu*dz);

% % Create empty arrays to hold the x, y, z coordinates of points for each material
% x_concrete = []; y_concrete = []; z_concrete = [];
% x_steel = []; y_steel = []; z_steel = [];
% x_rust = []; y_rust = []; z_rust = [];
% 
% % Iterate over all cells in the grid
% for x = 1:N_cells_x
%     for y = 1:N_cells_y
%         for z = 1:N_cells_z
%             % Determine the material in the current cell based on the const_e1 value
%             if const_e1(x, y, z) == ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt))
%                 x_concrete = [x_concrete, x];
%                 y_concrete = [y_concrete, y];
%                 z_concrete = [z_concrete, z];
%             elseif const_e1(x, y, z) == 1
%                 x_steel = [x_steel, x];
%                 y_steel = [y_steel, y];
%                 z_steel = [z_steel, z];
%             elseif const_e1(x, y, z) == ((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt))
%                 x_rust = [x_rust, x];
%                 y_rust = [y_rust, y];
%                 z_rust = [z_rust, z];
%             end
%         end
%     end
% end
% 
% % Plot the points for each material in a different color
% figure;
% scatter3(x_concrete, y_concrete, z_concrete, 'b'); % Concrete in blue
% hold on;
% scatter3(x_steel, y_steel, z_steel, 'k'); % Steel in black
% scatter3(x_rust, y_rust, z_rust, 'r'); % Rust in red
% 
% % Set the plot title and labels
% title('Concrete-Rebar Arrangement');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% 
% % Set the legend
% legend('Concrete', 'Steel', 'Rust');
% 
% % Enable rotation for better view
% rotate3d on;

% Create empty arrays to hold the x, y coordinates of points for each material
x_concrete = []; y_concrete = [];
x_steel = []; y_steel = [];
x_rust = []; y_rust = [];

% Set z to middle slice
z = Mid_cell_z;

% Iterate over all cells in the slice
for x = 1:N_cells_x
    for y = 1:N_cells_y
        % Determine the material in the current cell based on the const_e1 value
        if const_e1(x, y, z) == ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt))
            x_concrete = [x_concrete, x];
            y_concrete = [y_concrete, y];
        elseif const_e1(x, y, z) == 1
            x_steel = [x_steel, x];
            y_steel = [y_steel, y];
        elseif const_e1(x, y, z) == ((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt))
            x_rust = [x_rust, x];
            y_rust = [y_rust, y];
        end
    end
end

% Plot the points for each material in a different color
figure;
scatter(x_concrete, y_concrete, 'b', 'filled'); % Concrete in blue
hold on;
scatter(x_steel, y_steel, 'k', 'filled'); % Steel in black
scatter(x_rust, y_rust, 'r', 'filled'); % Rust in red

% Set the plot title and labels
title('Concrete-Rebar Arrangement (Middle Z Slice)');
xlabel('X');
ylabel('Y');

% Set the legend
legend('Concrete', 'Steel', 'Rust');

% Set the axes equal for a correct aspect ratio
axis equal;


%------- Initialise E and H arrays ----------
ex=zeros(N_cells_x,N_cells_y,N_cells_z);
ey=zeros(N_cells_x,N_cells_y,N_cells_z);
ez=zeros(N_cells_x,N_cells_y,N_cells_z);

hx=zeros(N_cells_x,N_cells_y,N_cells_z);
hy=zeros(N_cells_x,N_cells_y,N_cells_z);
hz=zeros(N_cells_x,N_cells_y,N_cells_z);

%------- Tangential component around rebar ----------
gx = ones(N_cells_x, N_cells_y, N_cells_z);
gy = ones(N_cells_x, N_cells_y, N_cells_z);
gz = ones(N_cells_x, N_cells_y, N_cells_z);

if(bar_exists==1)
    gx(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :) = 0;
    gy(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :) = 0;
    gz(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :) = 0;
end

%------ Gaussian pulse ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=1700; % maximum iteration steps
vin = zeros(N_steps,1);
vout = zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

fc = 7.5e9; % Center frequency
bw = 9e9;

ez_past = ez;
ey_past = ey;
ex_past = ex;

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

    ez(1, Mid_cell_y, Mid_cell_z) = ez(1, Mid_cell_y, Mid_cell_z) + pulse;
    vin(n) = ez(1, Mid_cell_y, Mid_cell_z);

     %--------- compute H ----------
     i = 1:N_cells_x-1;
     j = 1:N_cells_y-1;
     k = 1:N_cells_z-1;
     hx(i,j,k) = hx(i,j,k) + const_hz*(ey(i,j,k+1)-ey(i,j,k)) - const_hy*(ez(i,j+1,k)-ez(i,j,k));
     hy(i,j,k) = hy(i,j,k) + const_hx*(ez(i+1,j,k)-ez(i,j,k)) - const_hz*(ex(i,j,k+1)-ex(i,j,k));
     hz(i,j,k) = hz(i,j,k) + const_hy*(ex(i,j+1,k)-ex(i,j,k)) - const_hx*(ey(i+1,j,k)-ey(i,j,k));

     %--------- compute E ----------
     i = 2:N_cells_x;
     j = 2:N_cells_y;
     k = 2:N_cells_z;
     ex(i,j,k) = gx(i,j,k).*(const_e1(i,j,k).*ex(i,j,k) + const_ey(i,j,k).*(hz(i,j,k)-hz(i,j-1,k)) - const_ez(i,j,k).*(hy(i,j,k)-hy(i,j,k-1)));
     ey(i,j,k) = gy(i,j,k).*(const_e1(i,j,k).*ey(i,j,k) + const_ez(i,j,k).*(hx(i,j,k)-hx(i,j,k-1)) - const_ex(i,j,k).*(hz(i,j,k)-hz(i-1,j,k)));
     ez(i,j,k) = gz(i,j,k).*(const_e1(i,j,k).*ez(i,j,k) + const_ex(i,j,k).*(hy(i,j,k)-hy(i-1,j,k)) - const_ey(i,j,k).*(hx(i,j,k)-hx(i,j-1,k)));

     %-----Mur's 1st order ABC---------
     %X boundaries
     ex(1,:,:) = ex_past(2,:,:) + const_abc1x*(ex_past(2,:,:)-ex_past(1,:,:)); 
     ey(1,:,:) = ey_past(2,:,:) + const_abc1y*(ey_past(2,:,:)-ey_past(1,:,:)); 
     ez(1,:,:) = ez_past(2,:,:) + const_abc1z*(ez_past(2,:,:)-ez_past(1,:,:)); 
 
     ex(N_cells_x,:,:) = ex_past(N_cells_x-1,:,:) + const_abc1x*(ex_past(N_cells_x-1,:,:)-ex_past(N_cells_x,:,:)); 
     ey(N_cells_x,:,:) = ey_past(N_cells_x-1,:,:) + const_abc1y*(ey_past(N_cells_x-1,:,:)-ey_past(N_cells_x,:,:));
     ez(N_cells_x,:,:) = ez_past(N_cells_x-1,:,:) + const_abc1z*(ez_past(N_cells_x-1,:,:)-ez_past(N_cells_x,:,:)); 

     %Y boundaries
     ex(:,1,:) = ex_past(:,2,:) + const_abc1x*(ex_past(:,2,:)-ex_past(:,1,:));
     ey(:,1,:) = ey_past(:,2,:) + const_abc1y*(ey_past(:,2,:)-ey_past(:,1,:));
     ez(:,1,:) = ez_past(:,2,:) + const_abc1z*(ez_past(:,2,:)-ez_past(:,1,:));

     ex(:,N_cells_y,:) = ex_past(:,N_cells_y-1,:) + const_abc1x*(ex_past(:,N_cells_y-1,:)-ex_past(:,N_cells_y,:));
     ey(:,N_cells_y,:) = ey_past(:,N_cells_y-1,:) + const_abc1y*(ey_past(:,N_cells_y-1,:)-ey_past(:,N_cells_y,:));
     ez(:,N_cells_y,:) = ez_past(:,N_cells_y-1,:) + const_abc1z*(ez_past(:,N_cells_y-1,:)-ez_past(:,N_cells_y,:));

     %Z boundaries
     ex(:,:,1) = ex_past(:,:,2) + const_abc1x*(ex_past(:,:,2)-ex_past(:,:,1));
     ey(:,:,1) = ey_past(:,:,2) + const_abc1y*(ey_past(:,:,2)-ey_past(:,:,1));
     ez(:,:,1) = ez_past(:,:,2) + const_abc1z*(ez_past(:,:,2)-ez_past(:,:,1));

     ex(:,:,N_cells_z) = ex_past(:,:,N_cells_z-1) + const_abc1x*(ex_past(:,:,N_cells_z-1)-ex_past(:,:,N_cells_z));
     ey(:,:,N_cells_z) = ey_past(:,:,N_cells_z-1) + const_abc1y*(ey_past(:,:,N_cells_z-1)-ey_past(:,:,N_cells_z));
     ez(:,:,N_cells_z) = ez_past(:,:,N_cells_z-1) + const_abc1z*(ez_past(:,:,N_cells_z-1)-ez_past(:,:,N_cells_z));

     ez_past = ez;
     ey_past = ey;
     ex_past = ex;

     vout(n) = ez(N_cells_x, Mid_cell_y, Mid_cell_z);

     % ------- Plotting --------
     xslice = Mid_cell_x; % Define the x slice plane
     yslice = Mid_cell_y; % Define the y slice plane
     zslice = Mid_cell_z; % Define the z slice plane

     % Plotting the slices of the electric field (Ez)
     figure(2)
     slice(X, Y, Z, ez, xslice, yslice, zslice)
     colorbar
     xlabel('x')
     ylabel('y')
     zlabel('z')
     title(['Electric field (Ez) - Iteration: ', num2str(n)])
     drawnow

end

break
end

if (choice==1)
%----- FFT parameters ------
fs = 1e12; % sampling frequency
t = 0:1/fs:1e-8; % time base
FFT_steps = 1e7; % maximum iteration steps

% FFT of the pulse
X = fftshift(fft(vout,FFT_steps));
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
    title(sprintf('FFT (rs=%.0f, rr=%.0f)', rs, rr));
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
    fprintf('\nResults for rs=%.2f, rr=%.2f\n', rs,rr);
end

% Display the starting and ending frequencies of the 3 dB bandwidth
fprintf('Starting freq of 3 dB bandwidth: %.2f GHz\n', freq_start/1e9);
fprintf('Ending freq of 3 dB bandwidth: %.2f GHz\n', freq_end/1e9);

fprintf('\nFrequency at max magnitude: %.2f GHz\n', max_freq/1e9);
fprintf('Maximum magnitude: %.2f dB\n', max_mag);
fprintf('Phase at max magnitude: %.2f radians\n', max_phase);
end



