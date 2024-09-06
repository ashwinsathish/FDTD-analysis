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

%----- Circular Rebar Parameters -----
rs = 50; % radius of the steel rebar

N_cells_x = round(length_x/dx);
N_cells_y = round(length_y/dy);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);

[X,Y] = meshgrid(1:N_cells_x, 1:N_cells_y);

% Initialize space with concrete
space = ones(N_cells_x, N_cells_y);

% Create circular steel rebar in the center
for i = 1:N_cells_x
    for j = 1:N_cells_y
        if (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 <= rs^2
            space(i, j) = eps_steel;
        end
    end
end

% Visualize the space
figure;
imagesc(space);
%colormap([0.8, 0.8, 0.8; 0.5, 0.5, 0.5]); % Concrete: light gray, Steel: dark gray
xlabel('y');
ylabel('x');
title('2D rectangular piece of concrete with a steel rebar');
axis equal tight;

%------- Rebar radius and conditions --------
bar_exists=1;
rebarMask = ((X - Mid_cell_x).^2 + (Y - Mid_cell_y).^2) <= rs^2;

%------ Multiplying constants --------
const_e1=zeros(N_cells_x,N_cells_y);
const_e2=zeros(N_cells_x,N_cells_y);
const_e3=zeros(N_cells_x,N_cells_y);

%------- Concrete region ---------
const_e1(~rebarMask)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_e2(~rebarMask)=(2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
const_e3(~rebarMask)=(2*dt/dy)/((2*eps_conc)+(cond_conc*dt));

%------- Steel region ---------
const_e1(rebarMask)=1;
const_e2(rebarMask)=(dt/dx)/eps_steel;
const_e3(rebarMask)=(dt/dy)/eps_steel;

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
for i = 1:N_cells_x
    for j = 1:N_cells_y
        if (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 <= (rs+1)^2 && (i - Mid_cell_x)^2 + (j - Mid_cell_y)^2 >= (rs-1)^2
            gz(i,j) = 0;
        end
    end
end

%------ Gaussian pusle ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=1600; % maximum iteration steps
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

    ez(1, Mid_cell_y) = ez(1, Mid_cell_y) + pulse;
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
     zlim([-0.1 0.1]);
     title('Electric field (Ez)');
     pause(0.01);
        
end

break
end

