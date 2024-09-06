%----- Pulse propagation in homogeneous medium-----

clc;
close all;

%----- Medium ----------
length = 2; % length of concrete block = 2m 
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

%---- Signal -----------
c=3e8;
v=c/sqrt(real_epsr_conc*meur);
freq=3e9;
lamda=v/freq;
lamda_min=lamda/sqrt(real_epsr_conc*meur);

%------ Cell length and time step --------
dz=lamda_min/25;
dt=dz/v;
N_cells=round(length/dz);
Mid_cell=round(N_cells/2); % position to launch the signal

%------ Conductivity values ---------
cond_conc=2*pi*freq*eps0*img_epsr_conc;
cond_rust=2*pi*freq*eps0*img_epsr_rust;

%------ Multiplying constants --------
rr = 25;
const_e1=zeros(N_cells,1);

const_e1(1:N_cells)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_e1(Mid_cell-rr:Mid_cell+rr)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
% const_e1(Mid_cell-7:Mid_cell-5)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
% const_e1(Mid_cell+5:Mid_cell+7)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));

const_e2=zeros(N_cells,1);

const_e2(1:N_cells)=(2*dt/dz)/((2*eps_conc)+(cond_conc*dt));
const_e2(Mid_cell-rr:Mid_cell+rr)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));
% const_e2(Mid_cell-7:Mid_cell-5)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));
% const_e2(Mid_cell+5:Mid_cell+7)=(2*dt/dz)/((2*eps_rust)+(cond_rust*dt));

const_h=dt/(meu*dz);
const_abc1=(v*dt-dz)/(v*dt+dz);
const_abc2=(2*dz)/(v*dt+dz);

%------ Sine pulse ------------
Ts=0.5e-9; % pulse width
t0=2*Ts; % pulse delay
N_steps=2500; % maximum iteration steps
fc = 7.5e9; % Center frequency

%------- Initialise E and H arrays ------------
ex=zeros(N_cells,1);
hy=zeros(N_cells,1);

%------ Recording i/p & o/p --------
vin=zeros(N_steps,1);
vout=zeros(N_steps,1);
time_stamp=zeros(N_steps,1);

ex_past = zeros(N_cells,N_steps);

%************** Iteration loop ****************
for n=1:N_steps 

        time=(n-1)*dt; 
        time_stamp(n,1)=time;
        
%------- Signal launched in the middle cell -----------

        pulse1 = sin(2*pi*fc*time); % sine pulse
        pulse2 = exp(-((time-t0)/Ts)^2); % Gaussian pulse
        pulse = pulse2.*pulse1;

        ex(3)=ex(3) + pulse; % launching pulse
        vin(n,1)=ex(1);

        %--------- compute H ----------
        k=1:N_cells-1;
        hy(k)=hy(k)-const_h*(ex(k+1)-ex(k)); 

        %--------- compute E ----------
        k=2:N_cells-1; 
        ex(k)=const_e1(k).*ex(k)-const_e2(k).*(hy(k)-hy(k-1)); 

        %------------ 2nd order ABC ------------- 
        n=n+2;

        % left boundary
        ex(1) = - ex_past(2,n-2) + const_abc1*(ex(2) + ex_past(1,n-2)) + const_abc2*(ex_past(2,n-1) + ex_past(1,n-1));
        % right boundary
        ex(N_cells) = - ex_past(N_cells-1,n-2) + const_abc1*(ex(N_cells-1) + ex_past(N_cells,n-2)) + const_abc2*(ex_past(N_cells-1,n-1) + ex_past(N_cells,n-1));
        ex_past(:,n)=ex; % update ex_past array

        n=n-2;
    
        plot(1:N_cells,ex,'-');
        title(sprintf('Ez (rust = %.0f cells) at time step %d', 2*rr, n));
        axis([1 N_cells -1 1]);
        pause(0.001);  
        vout(n,1)=ex(N_cells);       
        
end  % for loop ends

