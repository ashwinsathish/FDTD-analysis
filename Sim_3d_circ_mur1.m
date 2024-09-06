clc;
clear;
close all;

%----- Medium ----------
length_x = 0.05;
length_y = 0.05;
length_z = 0.05;

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
dx=1e-3;
dy=dx;
dz=dy;
dt=1/(v*sqrt(1/(dx^2) + 1/(dy^2) + 1/(dz^2)));

%------- Rust conditions ---------
bar_exists = 1;
corr = 5;

%------- Corrosion base values ---------
rs_paper = [12.50,12.18,11.85,11.52,11.18,10.82,10.45]/2000;
rr_paper = [12.50,13.40,14.25,15.05,15.81,16.53,17.23]/2000;

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
N_cells_z = round(length_z/dz);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);
Mid_cell_z = round(N_cells_z/2);

[Y, X, Z] = meshgrid(1:N_cells_y, 1:N_cells_x, 1:N_cells_z);

% Initialize space with concrete
space = ones(N_cells_x, N_cells_y, N_cells_z)*eps_conc;

%------- Rebar radius and conditions --------
if(bar_exists==1)
    rebarMask = false(N_cells_x, N_cells_y, N_cells_z);
    rustMask = false(N_cells_x, N_cells_y, N_cells_z);
    for i = 1:N_cells_x
        for j = 1:N_cells_y
            for k = 1:N_cells_z
                distSq = (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2;
                rebarMask(i, j, k) = distSq <= rs^2;
                rustMask(i, j, k) = distSq <= (rs+rr)^2 & distSq > rs^2;
            end
        end
    end

    for i = 1:N_cells_x
        for j = 1:N_cells_y
            for k = 1:N_cells_z
                if rebarMask(i, j)
                    space(i, j, k) = eps_steel;
                elseif rustMask(i, j)
                    space(i, j, k) = eps_rust;
                end
            end
        end
    end
end

% ---- Visualizing the space ----
figure;
z_slice = round(N_cells_z/2);
imagesc(space(:, :, z_slice));
xlabel('y');
ylabel('x');
title(['Steel rebar and rust layer at z = ', num2str(z_slice)]);
axis equal tight;

%------ Multiplying constants --------
const_e1=zeros(N_cells_x,N_cells_y, N_cells_z);
const_ex=zeros(N_cells_x,N_cells_y, N_cells_z);
const_ey=zeros(N_cells_x,N_cells_y, N_cells_z);
const_ez=zeros(N_cells_x,N_cells_y, N_cells_z);

%------- Steel, Rust and Concrete region ---------
if(bar_exists==1)
    const_e1(~rebarMask & ~rustMask) = ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
    const_ex(~rebarMask & ~rustMask) = (2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
    const_ey(~rebarMask & ~rustMask) = (2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
    const_ez(~rebarMask & ~rustMask) = (2*dt/dz)/((2*eps_conc)+(cond_conc*dt));

    const_e1(rebarMask) = 0;
    const_ex(rebarMask) = 0;
    const_ey(rebarMask) = 0;
    const_ez(rebarMask) = 0;

    const_e1(rustMask) = ((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
    const_ex(rustMask) = (2*dt/dx)/((2*eps_rust)+(cond_rust*dt));
    const_ey(rustMask) = (2*dt/dy)/((2*eps_rust)+(cond_rust*dt));
    const_ez(rustMask) = (2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

else 

    const_e1 = ones(N_cells_x, N_cells_y, N_cells_z) * ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
    const_ex = ones(N_cells_x, N_cells_y, N_cells_z) * (2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
    const_ey = ones(N_cells_x, N_cells_y, N_cells_z) * (2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
    const_ez = ones(N_cells_x, N_cells_y, N_cells_z) * (2*dt/dz)/((2*eps_conc)+(cond_conc*dt));

end

const_hx = dt/(meu*dx);
const_hy = dt/(meu*dy);
const_hz = dt/(meu*dz);

const_abc1z = (v*dt - dx)/(v*dt + dx);
const_abc1x = (v*dt - dy)/(v*dt + dy);

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

if bar_exists==1
    for i = 1:N_cells_x
        for j = 1:N_cells_y
            for k = 1:N_cells_z
                if (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 <= (rs+1)^2 && (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 >= (rs-1)^2
                    gx(i,j,k) = 0;
                    gy(i,j,k) = 0;
                    gz(i,j,k) = 0;
                end
            end
        end
    end
end

%------ Gaussian pusle ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=1200; % maximum iteration steps

vin = zeros(N_steps,1);
vout = zeros(N_steps,1);
vout_avg = zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

ez_past = ez;
ex_past = ex;

fc = 7.5e9; % Center frequency
bw = 9e9;

%************** Iteration loop ****************

choice = input("Do you want to continue with the simulation (0-no, 1-yes)? ");

while choice==1

figure;
   
for n=1:N_steps         
    time=(n-1)*dt;  
    time_stamp(n,1)=time;
        
%------- Launching the signal ---------
    pulse1 = sin(2*pi*fc*time); % sine pulse
    pulse2 = exp(-((time-t0)/Ts)^2); % Gaussian pulse
    pulse = pulse2.*pulse1;

    ez(3, Mid_cell_y, Mid_cell_z) = pulse;
    vin(n) = ez(1, Mid_cell_y, Mid_cell_z);

     %--------- compute H ----------
     i = 1:N_cells_x-1;
     j = 1:N_cells_y-1;
     k = 1:N_cells_z-1;

     hx(i,j,k) = hx(i,j,k) + const_hz*(ey(i,j,k+1)-ey(i,j,k)) - const_hy*(ez(i,j+1,k)-ez(i,j,k));
     hy(i,j,k) = hy(i,j,k) + const_hx*(ez(i+1,j,k)-ez(i,j,k)) - const_hz*(ex(i,j,k+1)-ex(i,j,k));
     hz(i,j,k) = hz(i,j,k) + const_hy*(ez(i,j+1,k)-ez(i,j,k)) - const_hx*(ez(i+1,j,k)-ez(i,j,k));

     %--------- compute E ----------
     i = 2:N_cells_x;
     j = 2:N_cells_y;
     k = 2:N_cells_z;

     ex(i,j,k) = gx(i,j,k).*(const_e1(i,j,k).*ex(i,j,k) + const_ey(i,j,k).*(hz(i,j,k)-hz(i,j-1,k)) - const_ez(i,j,k).*(hy(i,j,k)-hy(i,j,k-1)));
     ey(i,j,k) = gy(i,j,k).*(const_e1(i,j,k).*ey(i,j,k) + const_ez(i,j,k).*(hx(i,j,k)-hx(i,j,k-1)) - const_ex(i,j,k).*(hy(i,j,k)-hy(i-1,j,k)));
     ez(i,j,k) = gz(i,j,k).*(const_e1(i,j,k).*ez(i,j,k) + const_ex(i,j,k).*(hy(i,j,k)-hy(i-1,j,k)) - const_ey(i,j,k).*(hy(i,j,k)-hy(i,j-1,k)));

     %----------- 1st order ABC ------------- 
     % X boundaries
%      ez(1,:,:) = ez_past(2,:,:) + const_abc1z*(ez(2,:,:)-ez(1,:,:)); % top boundary in X
%      ex(1,:,:) = ex_past(2,:,:) + const_abc1x*(ex(2,:,:)-ex(1,:,:)); % top boundary in X

%      ez(N_cells_x,:,:) = ez_past(N_cells_x-1,:,:) + const_abc1z*(ez(N_cells_x-1,:,:)-ez(N_cells_x,:,:)); % bottom boundary in X
%      ex(N_cells_x,:,:) = ex_past(N_cells_x-1,:,:) + const_abc1x*(ex(N_cells_x-1,:,:)-ex(N_cells_x,:,:)); % bottom boundary in X

     % Y boundaries
%      ez(:,1,:) = ez_past(:,2,:) + const_abc1z*(ez(:,2,:)-ez(:,1,:)); % top boundary in Y
%      ex(:,1,:) = ex_past(:,2,:) + const_abc1x*(ex(:,2,:)-ex(:,1,:)); % top boundary in Y
% 
%      ez(:,N_cells_y,:) = ez_past(:,N_cells_y-1,:) + const_abc1z*(ez(:,N_cells_y-1,:)-ez(:,N_cells_y,:)); % bottom boundary in Y
%      ex(:,N_cells_y,:) = ex_past(:,N_cells_y-1,:) + const_abc1x*(ex(:,N_cells_y-1,:)-ex(:,N_cells_y,:)); % bottom boundary in Y

     % Z boundaries
%      ez(:,:,1) = ez_past(:,:,2) + const_abc1z*(ez(:,:,2)-ez(:,:,1)); % top boundary in Z
%      ex(:,:,1) = ex_past(:,:,2) + const_abc1x*(ex(:,:,2)-ex(:,:,1)); % top boundary in Z

%      ez(:,:,N_cells_z) = ez_past(:,:,N_cells_z-1) + const_abc1z*(ez(:,:,N_cells_z-1)-ez(:,:,N_cells_z)); % bottom boundary in Z
%      ex(:,:,N_cells_z) = ex_past(:,:,N_cells_z-1) + const_abc1x*(ex(:,:,N_cells_z-1)-ex(:,:,N_cells_z)); % bottom boundary in Z

     ez_past = ez;
     ex_past = ex;

     % ------- Recording output -------
     vout(n)=ez(N_cells_x,Mid_cell_y,Mid_cell_z);

    % ------- Plotting --------
    z_slice = round(N_cells_z/2);
    imagesc(squeeze(ez(:, :, z_slice)));
    xlabel('y');
    ylabel('x');
    title(['Ez at timestep ', num2str(n), ' and z = ', num2str(z_slice)]);
    colorbar;
    clim([-0.5 0.5]);
    axis equal tight;
    pause(0.01);
end

break
end
