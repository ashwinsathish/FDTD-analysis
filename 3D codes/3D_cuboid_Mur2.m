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
N_cells_z = round(length_z/dz);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);
Mid_cell_z = round(N_cells_z/2);

[X,Y,Z] = meshgrid(1:N_cells_x, 1:N_cells_y, 1:N_cells_z);

% Initialize space with concrete
space = ones(N_cells_x, N_cells_y, N_cells_z)*eps_conc;

%------- Rebar radius and conditions --------
if(bar_exists==1)
    rebarMask = ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) <= rs^2;
    rustMask = ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) <= (rs+rr)^2 & ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) > rs^2;

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

    const_e1 = ones(N_cells_x, N_cells_y) * ((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
    const_ex = ones(N_cells_x, N_cells_y) * (2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
    const_ey = ones(N_cells_x, N_cells_y) * (2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
    const_ez = ones(N_cells_x, N_cells_y) * (2*dt/dz)/((2*eps_conc)+(cond_conc*dt));

end

const_hx = dt/(meu*dx);
const_hy = dt/(meu*dy);
const_hz = dt/(meu*dz);

const_abc1 = (v*dt - dx)/(v*dt + dx);
const_abc2 = (2*dx)/(v*dt + dx);
const_abc3 = (dx*(v*dt)^2)/(2*(dy^2)*(v*dt + dx));
const_abc4 = (dx*(c*dt)^2)/(2*(dx^2)*(c*dt + dx));

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
N_steps=2500; % maximum iteration steps

vin = zeros(N_steps,1);
vout = zeros(N_steps,1);
vout_avg = zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

ez_past = zeros(N_cells_x,N_cells_y,N_cells_z,N_steps);

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

    ez(1, Mid_cell_y-10:Mid_cell_y+10, 1) = ez(1, Mid_cell_y-10:Mid_cell_y+10, 1) + pulse;
    vin(n) = ez(1, Mid_cell_y, 1);

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

     %------------ 2nd order ABC ------------- 
     n=n+2;
     
%      for j = 2:N_cells_y-1
%         for k = 2:N_cells_z-1
%             % x = 1 face
%             ez(1, j, k) = -ez_past(2,j,k,n-2) + const_abc1*(ez(2,j,k) + ez_past(1,j,k,n-2)) + const_abc2*(ez_past(2,j,k,n-1) + ez_past(1,j,k,n-1)) + const_abc3*(ez_past(1,j+1,k,n-1) - 2*ez_past(1,j,k,n-1) + ez_past(1,j-1,k,n-1) + ez_past(2,j+1,k,n-1) - 2*ez_past(2,j,k,n-1) + ez_past(2,j-1,k,n-1)) + const_abc4*(ez_past(1,j,k+3,n-1) - 2*ez_past(1,j,k,n-1) + ez_past(1,j,k-1,n-1) + ez_past(2,j,k+3,n-1) - 2*ez_past(2,j,k+1,n-1) + ez_past(2,j,k-1,n-1));
%        
%             % x = N_cells_x face
%             ez(N_cells_x, j, k) = -ez_past(N_cells_x-1, j, k, n-2) + const_abc1*(ez(N_cells_x-1, j, k) + ez_past(N_cells_x, j, k, n-2)) + const_abc2*(ez_past(N_cells_x-1, j, k, n-1) + ez_past(N_cells_x, j, k, n-1)) + const_abc3*(ez_past(N_cells_x-1, j+1, k, n-1) - 2*ez_past(N_cells_x, j, k, n-1) + ez_past(N_cells_x-1, j-1, k, n-1) + ez_past(N_cells_x, j+1, k, n-1) - 2*ez_past(N_cells_x-1, j, k, n-1) + ez_past(N_cells_x, j-1, k, n-1)) + const_abc4*(ez_past(N_cells_x-1,j,k+3,n-1) - 2*ez_past(N_cells_x-1,j,k,n-1) + ez_past(N_cells_x-1,j,k-1,n-1) + ez_past(N_cells_x,j,k+3,n-1) - 2*ez_past(N_cells_x,j,k+1,n-1) + ez_past(N_cells_x,j,k-1,n-1));
%         end
%      end

     for i = 2:N_cells_x-1
        for k = 2:N_cells_z-1
            % y = 1 face
            ez(i, 1, k) = -ez_past(i, 2, k, n-2) + const_abc1*(ez(i, 2, k) + ez_past(i, 1, k, n-2)) + const_abc2*(ez_past(i, 2, k, n-1) + ez_past(i, 1, k, n-1)) + const_abc3*(ez_past(i+1, 1, k, n-1) - 2*ez_past(i, 1, k, n-1) + ez_past(i-1, 1, k, n-1) + ez_past(i+1, 2, k, n-1) - 2*ez_past(i, 2, k, n-1) + ez_past(i-1, 2, k, n-1)) + const_abc4*(ez_past(i, 1, k+1, n-1) - 2*ez_past(i, 1, k, n-1) + ez_past(i, 1, k-1, n-1) + ez_past(i, 2, k+1, n-1) - 2*ez_past(i, 2, k, n-1) + ez_past(i, 2, k-1, n-1));
            ex(i, 1, k) = -ex_past(i, 2, k, n-2) + const_abc1*(ex(i, 2, k) + ex_past(i, 1, k, n-2)) + const_abc2*(ex_past(i, 2, k, n-1) + ex_past(i, 1, k, n-1)) + const_abc3*(ez_past(i+1, 1, k, n-1) - 2*ez_past(i, 1, k, n-1) + ez_past(i-1, 1, k, n-1) + ez_past(i+1, 2, k, n-1) - 2*ez_past(i, 2, k, n-1) + ez_past(i-1, 2, k, n-1)) + const_abc4*(ez_past(i, 1, k+1, n-1) - 2*ez_past(i, 1, k, n-1) + ez_past(i, 1, k-1, n-1) + ez_past(i, 2, k+1, n-1) - 2*ez_past(i, 2, k, n-1) + ez_past(i, 2, k-1, n-1));

            % y = N_cells_y face
            ez(i, N_cells_y, k) = -ez_past(i, N_cells_y-1, k, n-2) + const_abc1*(ez(i, N_cells_y-1, k) + ez_past(i, N_cells_y, k, n-2)) + const_abc2*(ez_past(i, N_cells_y-1, k, n-1) + ez_past(i, N_cells_y, k, n-1)) + const_abc3*(ez_past(i+1, N_cells_y, k, n-1) - 2*ez_past(i, N_cells_y, k, n-1) + ez_past(i-1, N_cells_y, k, n-1) + ez_past(i+1, N_cells_y-1, k, n-1) - 2*ez_past(i, N_cells_y-1, k, n-1) + ez_past(i-1, N_cells_y-1, k, n-1)) + const_abc4*(ez_past(i, N_cells_y, k+1, n-1) - 2*ez_past(i, N_cells_y, k, n-1) + ez_past(i, N_cells_y, k-1, n-1) + ez_past(i, N_cells_y-1, k+1, n-1) - 2*ez_past(i, N_cells_y-1, k, n-1) + ez_past(i, N_cells_y-1, k-1, n-1));
        end
     end

%     for i = 2:N_cells_x-1
%         for j = 2:N_cells_y-1
%              % z = 1 face
%             ez(i, j, 1) = -ez_past(i, j, 2, n-2) + const_abc1*(ez(i, j, 2) + ez_past(i, j, 1, n-2)) + const_abc2*(ez_past(i, j, 2, n-1) + ez_past(i, j, 1, n-1)) + const_abc3*(ez_past(i+1, j, 1, n-1) - 2*ez_past(i, j, 1, n-1) + ez_past(i-1, j, 1, n-1) + ez_past(i+1, j, 2, n-1) - 2*ez_past(i, j, 2, n-1) + ez_past(i-1, j, 2, n-1)) + const_abc4*(ez_past(i, j+1, 1, n-1) - 2*ez_past(i, j, 1, n-1) + ez_past(i, j-1, 1, n-1) + ez_past(i, j+1, 2, n-1) - 2*ez_past(i, j, 2, n-1) + ez_past(i, j-1, 2, n-1));
% 
%             % z = N_cells_z face
%             ez(i, j, N_cells_z) = -ez_past(i, j, N_cells_z-1, n-2) + const_abc1*(ez(i, j, N_cells_z-1) + ez_past(i, j, N_cells_z, n-2)) + const_abc2*(ez_past(i, j, N_cells_z-1, n-1) + ez_past(i, j, N_cells_z, n-1)) + const_abc3*(ez_past(i+1, j, N_cells_z, n-1) - 2*ez_past(i, j, N_cells_z, n-1) + ez_past(i-1, j, N_cells_z, n-1) + ez_past(i+1, j, N_cells_z-1, n-1) - 2*ez_past(i, j, N_cells_z-1, n-1) + ez_past(i-1, j, N_cells_z-1, n-1)) + const_abc4*(ez_past(i, j+1, N_cells_z, n-1) - 2*ez_past(i, j, N_cells_z, n-1) + ez_past(i, j-1, N_cells_z, n-1) + ez_past(i, j+1, N_cells_z-1, n-1) - 2*ez_past(i, j, N_cells_z-1, n-1) + ez_past(i, j-1, N_cells_z-1, n-1));
%         end
%     end

     ez_past(:,:,:,n) = ez;
     n=n-2;

     % ------- Recording output -------
     vout(n)=ez(N_cells_x,N_cells_y,Mid_cell_z);

    % ------- Plotting --------
    slice(X, Y, Z, ez, [], [], Mid_cell_z); % Display a slice at the mid Z position
    axis([1 N_cells_x 1 N_cells_y 1 N_cells_z]);
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title(['Electric field (ez) at iteration: ', num2str(n)]);
    colorbar;
    drawnow;

end

break
end
