%----- Pulse propagation in homogeneous medium-----
% ABC is not applied at the ends, so reflection will be visible

clc;
close all;

%----- Medium ----------
length = 0.1; % length of concrete block = 2m 
eps0=8.854e-12;
meu0=4*pi*1e-7;
real_epsr_conc=4; % relative permittivity of concrete = 4 - j0.06
real_epsr_rust=10; % relative permittivity of rust = 10 - j0.86
img_epsr_conc=0.06;
img_epsr_rust=0.86;
meur=1;
eps_conc=eps0*real_epsr_conc;
eps_rust=eps0*real_epsr_rust;
meu=meu0*meur;

%---- Signal -----
c=3e8;
v=c/sqrt(real_epsr_conc*meur);
freq=12e9;
lamda=v/freq;
lamda_min = lamda/sqrt(real_epsr_conc*meur);

%------ Cell length and time step---------
dz=lamda_min/25;
dt=dz/v;
N_cells=round(length/dz);
Mid_cell=round(N_cells/2); % position to launch the signal

%------ Conductivity values ---------
cond_conc=2*pi*freq*eps0*img_epsr_conc;
cond_rust=2*pi*freq*eps0*img_epsr_rust;

rr=14;

%------ Multiplying constants --------
const_e1=zeros(N_cells,1);
const_e1(1:N_cells)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_e1(Mid_cell-1:Mid_cell+1)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));

const_e2=zeros(N_cells,1);
const_e2(1:N_cells)=(2*dt/dz)/((2*eps_conc)+(cond_conc*dt));
const_e2(Mid_cell-rr:Mid_cell+rr)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

const_h=dt/(meu*dz);
const_abc=(v*dt-dz)/(v*dt+dz);

%------- Initialise E and H arrays ------------
ex=zeros(N_cells,1);
% ex_past=ex;
hy=zeros(N_cells,1);

%------ Sine pulse ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=1600; % maximum iteration steps

fc = 7.5e9; % Center frequency
bw = 9e9;

%------ Recording i/p & o/p --------
vin=zeros(N_steps,1);
vout=zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

ex_past=zeros(N_cells, N_steps);

%************** Iteration loop ****************
   
    for n=1:N_steps         
        time=(n-1)*dt; 
        time_stamp(n,1)=time;
        
%------- Signal launched in the middle cell -----------
        %pulse=sin(2*pi*freq*time);   % sine pulse
         pulse1 = sin(2*pi*fc*time); % sine pulse
         pulse2 = exp(-((time-t0)/Ts)^2); % Gaussian pulse
         pulse = pulse2.*pulse1;


        ex(3)=ex(3)+pulse; % launching pulse
        vin(n,1)=ex(1);

        %--------- compute H ----------
        k=1:N_cells-1;
        hy(k)=hy(k)-const_h*(ex(k+1)-ex(k)); 

        %--------- compute E ----------
        k=2:N_cells; 
        ex(k)=const_e1(k).*ex(k)-const_e2(k).*(hy(k)-hy(k-1)); 

        %------------ 1st order ABC ------------- 
%         ex(1)=ex_past(2)+const_abc*(ex(2)-ex(1)); % left boundary
%         ex(N_cells)=ex_past(N_cells-1)+const_abc*(ex(N_cells-1)-ex(N_cells)); % right boundary
%         ex_past=ex;

        n=n+1;

        ex(1)=ex_past(2,n-1)+const_abc*(ex(2)-ex_past(1,n-1)); % left boundary
        ex(N_cells)=ex_past(N_cells-1,n-1)+const_abc*(ex(N_cells-1)-ex_past(N_cells,n-1)); % right boundary
        ex_past(:,n)=ex; % update ex_past array
        
        n=n-1;


       
        plot(1:N_cells,ex,'-');
        axis([1 N_cells -1 1]);
        pause(0.001);  
        vout(n,1)=ex(N_cells);       
        
    end  % for loop ends

%---------- Plot the results -----------------
close all; 
figure;
plot(1:N_cells,ex);
title('Electric field Ex vs Z');
xlabel('cell number');
ylabel('Ex, V/m');
 
figure;
plot(1:N_cells,hy);
title('Magnetic field Hy vs Z');
xlabel('cell number');
ylabel('Hy, A/m');


