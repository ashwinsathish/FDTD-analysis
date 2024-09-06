%----- Pulse propagation in homogeneous medium-----
% ABC is not applied at the ends, so reflection will be visible

clc;
close all;

%----- Medium ----------
length_x = 2;
length_y = 2;
eps0=8.854e-12;
meu0=4*pi*1e-7;
real_epsr_conc=1; % relative permittivity of concrete = 4 - j0.06
img_epsr_conc=0;
real_epsr_rust=2; % relative permittivity of rust = 10 - j0.86
img_epsr_rust=0;
eps_conc=eps0*real_epsr_conc;
eps_rust=eps0*real_epsr_rust;
meu=meu0*meur;

%---- Signal -----------
c=3e8;
v=c/sqrt(real_epsr_conc*meur);
freq=3e9;
lamda=v/freq;

%------ Conductivity values ---------
cond_conc=2*pi*freq*eps0*img_epsr_conc;
cond_rust=2*pi*freq*eps0*img_epsr_rust;

%------ Cell length and time step---------
dx=lamda/20;
dy=lamda/20;
dt=1/(v*sqrt(1/(dx^2) + 1/(dy^2)));

N_cells_x = round(length_x/dx);
N_cells_y = round(length_y/dy);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);

[X,Y] = meshgrid(1:N_cells_x, 1:N_cells_y);

%------ Multiplying constants --------
const_e1=zeros(N_cells_x,N_cells_y);
const_e1(1:N_cells_x,1:N_cells_y)=((2*eps_conc)-(cond_conc*dt))/((2*eps_conc)+(cond_conc*dt));
const_e1(50:100,50:100)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));
const_e1(50:100,50:100)=((2*eps_rust)-(cond_rust*dt))/((2*eps_rust)+(cond_rust*dt));

const_ex=zeros(N_cells_x,N_cells_y);
const_ex(1:N_cells_x,1:N_cells_y)=(2*dt/dx)/((2*eps_conc)+(cond_conc*dt));
const_ex(50:100,50:100)=(2*dt/dx)/((2*eps_rust)+(cond_rust*dt));

const_ey=zeros(N_cells_x,N_cells_y);
const_ey(1:N_cells_x,1:N_cells_y)=(2*dt/dy)/((2*eps_conc)+(cond_conc*dt));
const_ey(50:100,50:100)=(2*dt/dy)/((2*eps_rust)+(cond_rust*dt));

const_hx=dt/(meu*dx);
const_hy=dt/(meu*dy);

const_abc1=(v*dt-dx)/(v*dt+dx);
const_abc2=(2*dx)/(v*dt+dx);
const_abc3=(dx*(v*dt)^2)/(2*(dy^2)*(v*dt+dx));

%------- Initialise E and H arrays ------------
ez=zeros(N_cells_x,N_cells_y);
hx=zeros(N_cells_x,N_cells_y);
hy=zeros(N_cells_x,N_cells_y);

%------ Sine pusle ------------
Ts=10*dt; % pulse width
t0=3*Ts; % pulse delay
N_steps=600; % maximum iteration steps

ez_past = zeros(N_cells_x,N_cells_y,N_steps);


%************** Iteration loop ****************
   
    for n=1:N_steps         
        time=(n-1)*dt;               
        
%------- Signal launched in the middle cell -----------
        % pulse = sin(2*pi*freq*time);
        pulse = exp(-((time-t0)/Ts)^2)/dz; % Gaussian pulse
        ez(Mid_cell_x, Mid_cell_y) = ez(Mid_cell_x, Mid_cell_y) + pulse;

        %--------- compute H ----------
        i = 1:N_cells_x-1;
        j = 1:N_cells_y-1;
        hx(i,j) = hx(i,j) - const_hy*(ez(i,j+1)-ez(i,j));
        hy(i,j) = hy(i,j) + const_hx*(ez(i+1,j)-ez(i,j));
 

        %--------- compute E ----------
        i = 2:N_cells_x;
        j = 2:N_cells_y;
        ez(i,j) = const_e1(i,j).*ez(i,j) + const_ex(i,j).*(hy(i,j)-hy(i-1,j)) - const_ey(i,j).*(hx(i,j)-hx(i,j-1));

        %------------ 2nd order ABC ------------- 
        n=n+2;
        % X boundaries
        for j=2:N_cells_y-1

%             ez(1,:) = - ez_past(2,:,n-2) + const_abc1*(ez(2,:) + ez_past(1,:,n-2)) + const_abc2*(ez_past(2,:,n-1) + ez_past(1,:,n-1));
%             ez(N_cells_x,:) = - ez_past(N_cells_x-1,:,n-2) +
%             const_abc1*(ez(N_cells_x-1,:) + ez_past(N_cells_x,:,n-2)) + const_abc2*(ez_past(N_cells_x-1,:,n-1) + ez_past(N_cells_x,:,n-1));
            
            ez(1,j) = - ez_past(2,j,n-2) + const_abc1*(ez(2,j) + ez_past(1,j,n-2)) + const_abc2*(ez_past(2,j,n-1) + ez_past(1,j,n-1)) + const_abc3*(ez_past(2,j+1,n-1) - 2*ez_past(1,j,n-1) + ez_past(2,j-1,n-1) + ez_past(1,j+1,n-1) - 2*ez_past(2,j,n-1) + ez_past(1,j-1,n-1));
            ez(N_cells_x,j) = - ez_past(N_cells_x-1,j,n-2) + const_abc1*(ez(N_cells_x-1,j) + ez_past(N_cells_x,j,n-2)) + const_abc2*(ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j,n-1)) + const_abc3*(ez_past(N_cells_x-1,j+1,n-1) - 2*ez_past(N_cells_x,j,n-1) + ez_past(N_cells_x-1,j-1,n-1) + ez_past(N_cells_x,j+1,n-1) - 2*ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j-1,n-1)); % bottom boundary
       
        end 

        % Y boundaries
        for i=2:N_cells_x-1
            
            % ez(:,1) = - ez_past(:,2,n-2) + const_abc1*(ez(:,2) + ez_past(:,1,n-2)) + const_abc2*(ez_past(:,2,n-1) + ez_past(:,1,n-1)); % left boundary
            % ez(:,N_cells_y) = - ez_past(:,N_cells_y-1,n-2) + const_abc1*(ez(:,N_cells_y-1) + ez_past(:,N_cells_y,n-2)) + const_abc2*(ez_past(:,N_cells_y-1,n-1) + ez_past(:,N_cells_y,n-1)); % right boundary
            
            ez(i,1) = - ez_past(i,2,n-2) + const_abc1*(ez(i,2) + ez_past(i,1,n-2)) + const_abc2*(ez_past(i,2,n-1) + ez_past(i,1,n-1)) + const_abc3*(ez_past(i+1,2,n-1) - 2*ez_past(i,1,n-1) + ez_past(i-1,2,n-1) + ez_past(i+1,1,n-1) - 2*ez_past(i,2,n-1) + ez_past(i-1,1,n-1)); % left boundary
            ez(i,N_cells_y) = - ez_past(i,N_cells_y-1,n-2) + const_abc1*(ez(i,N_cells_y-1) + ez_past(i,N_cells_y,n-2)) + const_abc2*(ez_past(i,N_cells_y-1,n-1) + ez_past(i,N_cells_y,n-1)) + const_abc3*(ez_past(i+1,N_cells_y-1,n-1) - 2*ez_past(i,N_cells_y,n-1) + ez_past(i-1,N_cells_y-1,n-1) + ez_past(i+1,N_cells_y,n-1) - 2*ez_past(i,N_cells_y-1,n-1) + ez_past(i-1,N_cells_y,n-1)); % right boundary

        end

        ez_past(:,:,n) = ez;
        n=n-2;
      

        surf(X, Y, ez);
        xlabel('x');
        ylabel('y');
        zlabel('Ez');
        xlim([1 N_cells_x]);
        ylim([1 N_cells_y]);
        zlim([-50 50]);
        title('Electric field (Ez)');
        pause(0.01);
       
        
    end  % for loop ends
   

    %---------- Plot the results -----------------
%      close all; 
%      figure;
%      plot(1:N_cells,ex);
%      title('Electric field Ex vs Z');
%      xlabel('cell number');
%      ylabel('Ex, V/m');
%  
%      figure;
%      plot(1:N_cells,hy);
%      title('Magnetic field Hy vs Z');
%      xlabel('cell number');
%      ylabel('Hy, A/m');
