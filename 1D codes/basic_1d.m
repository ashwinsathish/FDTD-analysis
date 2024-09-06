%------------------ Pulse propagation in homogeneous medium----------------------
%------------------- ABC is not applied at the ends, so reflection will be
%visible

clc;
close all;
%----- Medium ----------
length = 2;
eps0=8.854e-12;
meu0=4*pi*1e-7;
epsr=1;
meur=1;
eps=eps0*epsr;
meu=meu0*meur;
%---- Signal -----------
c=3e8;
v=c/sqrt(epsr*meur);
freq=3e9;
lamda=v/freq;
%------ Cell length and time step---------
dz=lamda/20;
dt=dz/v;
N_cells=length/dz;
Mid_cell=N_cells/2; % position to launch the signal
%------ Multiplying constants --------
const_e=dt/(eps*dz);
const_h=dt/(meu*dz);
%------- Initialise E and H arrays ------------
ex=zeros(N_cells,1);
hy=zeros(N_cells,1);
%------ Gaussian pusle ------------
Ts=10*dt; % pulse width
t0=3*Ts; % pulse delay
N_steps=400; % maximum iteration steps
%************** Iteration loop ****************
   
    for n=1:N_steps;        
        time=(n-1)*dt;              
       
%------- Signal launched in the middle cell -----------
        pulse=sin(2*pi*freq*time);         % Gaussian pulse exp(-((time-t0)/Ts)^2)/dz;              
        ex(Mid_cell)=pulse;    
        %------------------------ compute H -------------------------
        k=1:N_cells-1;
        hy(k)=hy(k)-const_h*(ex(k+1)-ex(k));          
        %------------------------- compute E ------------------------
        k=2:N_cells-1;
        ex(k)=ex(k)-const_e*(hy(k)-hy(k-1));  
       
        plot(1:N_cells,ex,'-');
        axis([1 N_cells -2 2]);
        pause(0.05);  
%         frame(i)=getframe(1);
%         i=i+1;
       
    end  % for loop ends
   



    %---------- Plot the results -----------------
%     close all;
%     plot(1:N_cells,ex);
%     title('Eelectric field Ex vs Z');
%     xlabel('cell number');
%     ylabel('Ex, V/m'); 
%
%    figure;
%    plot(1:N_cells,hy);
%     title('Magnetic field Hy vs Z');
%     xlabel('cell number');
%     ylabel('Hy, A/m');