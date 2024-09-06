%----- Pulse propagation in homogeneous medium-----
% ABC is not applied at the ends, so reflection will be visible

clc;
close all;

%----- Medium ----------
length = 2; 
eps0=8.854e-12;
meu0=4*pi*1e-7;
epsr=2; 
meur=1;
eps=eps0*epsr;
meu=meu0*meur;

%---- Signal -----------
c=3e8;
v=c/sqrt(1*meur);
freq=3e9;
lamda=v/freq;

%------ Cell length and time step---------
dz=lamda/20;
dt=dz/v;
N_cells=round(length/dz);
Mid_cell=round(N_cells/2); % position to launch the signal

%------ Multiplying constants --------
const_e=zeros(N_cells,1);
const_e(1:N_cells)=dt/(eps0*dz);
const_e(Mid_cell:Mid_cell+1)=dt/(eps*dz);
%const_e(Mid_cell-20:Mid_cell+20)=dt/(eps*dz);
%const_e(Mid_cell-50:Mid_cell+50)=dt/(eps*dz);
%const_e(Mid_cell-100:Mid_cell+100)=dt/(eps*dz);
const_h=dt/(meu*dz);
const_abc=(v*dt-dz)/(v*dt+dz);

%------- Initialise E and H arrays ------------
ex=zeros(N_cells,1);
ex_past=ex;
hy=zeros(N_cells,1);

%------ Sine pulse ------------
Ts=10*dt; % pulse width
t0=3*Ts; % pulse delay
N_steps=600; % maximum iteration steps

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
        ex(1)=pulse; % launching pulse
        vin(n,1)=ex(1);

        %--------- compute H ----------
        k=1:N_cells-1;
        hy(k)=hy(k)-const_h*(ex(k+1)-ex(k)); 

        %--------- compute E ----------
        k=2:N_cells-1; 
        ex(k)=ex(k)-const_e(k).*(hy(k)-hy(k-1)); 

        %------------ 1st order ABC ------------- 
        ex(1)=ex_past(2)+const_abc*(ex(2)-ex(1)); % left boundary
        ex(N_cells)=ex_past(N_cells-1)+const_abc*(ex(N_cells-1)-ex(N_cells)); % right boundary
        ex_past=ex;
       
        plot(1:N_cells,ex,'-');
        axis([1 N_cells -300 300]);
        pause(0.01);  
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
