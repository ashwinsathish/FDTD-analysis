%----- Pulse propagation in homogeneous medium-----
% ABC is not applied at the ends, so reflection will be visible

clc;
close all;

%----- Medium ----------
length = 2; % length of concrete block = 2m 
eps0=8.854e-12;
meu0=4*pi*1e-7;
real_epsr_conc=4; % relative permittivity of concrete = 4 - j0.06
img_epsr_conc=0.06;
meur=1;
eps_conc=eps0*real_epsr_conc;
meu=meu0*meur;
eps_steel = 1e7;

%---- Signal -----------
c=3e8;
v=c/sqrt(2*meur);
freq=3e9;
lamda=v/freq;

%------ Cell length and time step --------
dz=lamda/20;
dt=dz/v;
N_cells=round(length/dz);
Mid_cell=round(N_cells/2); % position to launch the signal

%------ Conductivity values ---------
cond_conc=2*pi*freq*eps0*img_epsr_conc;

%------ Multiplying constants --------
const_e1=zeros(N_cells,1);
const_e1(1:N_cells)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_e1(Mid_cell-1:Mid_cell+1)=1;

const_e2=zeros(N_cells,1);
const_e2(1:N_cells)=(2*dt/dz)/((2*eps_conc)+(cond_conc*dt));
const_e2(Mid_cell-1:Mid_cell+1)=dt/(2*eps_steel*dz);

const_h=dt/(meu*dz);
const_abc=(v*dt-dz)/(v*dt+dz); % v = vmax

%------- Initialise E and H arrays ------------
ex=zeros(N_cells,1);
ex_past=ex;
hy=zeros(N_cells,1);

%------ Sine pulse ------------
Ts=10*dt; % pulse width
t0=3*Ts; % pulse delay
N_steps=1200; % maximum iteration steps

%------ Recording i/p & o/p --------
vin=zeros(N_steps,1);
vout=zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

%************** Iteration loop ****************
   
    for n=1:N_steps         
        time=(n-1)*dt; 
        time_stamp(n,1)=time;
        
%------- Signal launched in the middle cell -----------
        %pulse=sin(2*pi*freq*time);   % sine pulse
        pulse = exp(-((time-t0)/Ts)^2)/dz;
        ex(5)=ex(5) + pulse; % launching pulse
        vin(n,1)=ex(1);

        %--------- compute H ----------
        k=1:N_cells-1;
        hy(k)=hy(k)-const_h*(ex(k+1)-ex(k)); 

        %--------- compute E ----------
        k=2:N_cells-1; 
        ex(k)=const_e1(k).*ex(k)-const_e2(k).*(hy(k)-hy(k-1)); 

        %------------ ABC ------------- 
        ex(1)=ex_past(2)+const_abc*(ex(2)-ex(1)); % left boundary
        ex(N_cells)=ex_past(N_cells-1)+const_abc*(ex(N_cells-1)-ex(N_cells)); % right boundary
        ex_past=ex;
        % now ex_past array must be initialised to zero along with ex and hy arrays 
        % const_abc = (v*dt-dz)/(v*dt+dz) must be calculated before the loop
       
        plot(1:N_cells,ex,'-');
        axis([1 N_cells -250 250]);
        pause(0.00001);  
        vout(n,1)=ex(N_cells);       
        
    end  % for loop ends

    % Save in AVI format
    % movie2avi(frame,'C:\Users\Dell\Desktop\gauss_ex.avi');

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
