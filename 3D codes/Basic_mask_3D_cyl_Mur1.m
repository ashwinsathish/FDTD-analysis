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
eps_steel=1e7; 
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
const_e2=zeros(N_cells_x,N_cells_y,N_cells_z);
const_e3=zeros(N_cells_x,N_cells_y,N_cells_z);
const_e4=zeros(N_cells_x,N_cells_y,N_cells_z);

%------- Rust conditions ---------
bar_exists = 1;
corr = 10;

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

%------- Corrosion simulation---------
const_e1(1:N_cells_x,1:N_cells_y,1:N_cells_z)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_ex(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
const_ey(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
const_ez(1:N_cells_x,1:N_cells_y,1:N_cells_z)=(2*dt/dz)/((2*eps_conc)+(cond_conc*dt));

if bar_exists==1
    
    % Rust region
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

        % Steel region
        const_e1(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=1;
        const_ex(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dx)/eps_steel;
        const_ey(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dy)/eps_steel;
        const_ez(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :)=(dt/dz)/eps_steel;

end

const_hx=dt/(meu*dx);
const_hy=dt/(meu*dy);
const_hz=dt/(meu*dz);

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
% g = ones(N_cells_x, N_cells_y, N_cells_z);
% if(bar_exists==1)
%     g(Mid_cell_x-rs:Mid_cell_x+rs, Mid_cell_y-rs:Mid_cell_y+rs, :) = 0;
% end

%------ Gaussian pulse ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=1600; % maximum iteration steps

vin = zeros(N_steps,1);
vout = zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

fc = 7.5e9; % Center frequency
bw = 9e9;

ez_past = ez;
ey_past = ey;

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

    %------ Update H-field -----------
    i = 1:N_cells_x-1;
    j = 1:N_cells_y-1;
    k = 1:N_cells_z-1;

    hx(i,j,k) = hx(i,j,k) + const_hz*(ey(i,j,k+1)-ey(i,j,k)) - const_hy*(ez(i,j+1,k)-ez(i,j,k));
    hy(i,j,k) = hy(i,j,k) + const_hx*(ez(i+1,j,k)-ez(i,j,k)) - const_hz*(ex(i,j,k+1)-ex(i,j,k));
    hz(i,j,k) = hz(i,j,k) + const_hy*(ez(i,j+1,k)-ez(i,j,k)) - const_hx*(ez(i+1,j,k)-ez(i,j,k));

    % --------- Update E-field ------------
    for i = 2:N_cells_x
        for j = 2:N_cells_y
            for k = 2:N_cells_z
                ex(i,j,k) = (const_e1(i,j,k).*ex(i,j,k) + const_ey(i,j,k).*(hz(i,j,k)-hz(i,j-1,k)) - const_ez(i,j,k).*(hy(i,j,k)-hy(i,j,k-1)));
                ey(i,j,k) = (const_e1(i,j,k).*ey(i,j,k) + const_ez(i,j,k).*(hx(i,j,k)-hx(i,j,k-1)) - const_ex(i,j,k).*(hz(i,j,k)-hz(i-1,j,k)));
                ez(i,j,k) = (const_e1(i,j,k).*ez(i,j,k) + const_ex(i,j,k).*(hy(i,j,k)-hy(i-1,j,k)) - const_ey(i,j,k).*(hx(i,j,k)-hx(i,j-1,k)));
            end
        end
    end

    % ------ Mur's 1st order ABC --------
 
     ez(1,:,:) = ez_past(2,:,:) + const_abc1z*(ez(2,:,:)-ez(1,:,:)); % top boundary in Y
     ey(1,:,:) = ey_past(2,:,:) + const_abc1x*(ey(2,:,:)-ey(1,:,:)); % top boundary in Y
     ez(N_cells_x,:,:) = ez_past(N_cells_x-1,:,:) + const_abc1z*(ez(N_cells_x-1,:,:)-ez(N_cells_x-1,:,:)); % bottom boundary in Y
     ey(N_cells_x,:,:) = ey_past(N_cells_x-1,:,:) + const_abc1x*(ey(N_cells_x-1,:,:)-ey(N_cells_x,:,:)); % bottom boundary in Y

     ez_past = ez;
     ey_past = ey;

     % ------- Recording output -------
     vout(n)=ez(10,Mid_cell_y,Mid_cell_z);

     % ------- Plotting --------
     xslice = Mid_cell_x; % Define the x slice plane
     yslice = Mid_cell_y; % Define the y slice plane
     zslice = Mid_cell_z; % Define the z slice plane

     % Plotting the slices of the electric field (Ez)
     figure(1)
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


