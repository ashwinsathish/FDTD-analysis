%----- Pulse propagation in homogeneous medium-----
% ABC is not applied at the ends, so reflection will be visible

clc;
clear;
close all;

%----- Medium ----------
length_x = 2;
length_y = 2;
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
dx=lamda/20;
dy=lamda/20;
dt=1/(v*sqrt(1/(dx^2) + 1/(dy^2)));

N_cells_x = round(length_x/dx);
N_cells_y = round(length_y/dy);

Mid_cell_x = round(N_cells_x/2); % position to launch the signal
Mid_cell_y = round(N_cells_y/2);

[X,Y] = meshgrid(1:N_cells_x, 1:N_cells_y);

%------ Multiplying constants --------
const_ex=dt/(eps*dx);
const_ey=dt/(eps*dy);
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
        pulse = sin(2*pi*freq*time);
        % pulse = exp(-((time-t0)/Ts)^2)/dz; % Gaussian pulse
        ez(Mid_cell_x, Mid_cell_y) = pulse;

        %--------- compute H ----------
        for i = 1:N_cells_x-1
            for j = 1:N_cells_y-1
                hx(i,j) = hx(i,j) - const_hy*(ez(i,j+1)-ez(i,j));
                hy(i,j) = hy(i,j) + const_hx*(ez(i+1,j)-ez(i,j));
            end
        end

        %--------- compute E ----------
        for i = 2:N_cells_x
            for j = 2:N_cells_y-1
                ez(i,j) = ez(i,j) + const_ex*(hy(i,j)-hy(i-1,j)) - const_ey*(hx(i,j)-hx(i,j-1));
            end
        end

        %------------ 2nd order ABC ------------- 
        n=n+2;
%         % X boundaries
%         ez(1,:) = - ez_past(2,:,n-2) + const_abc1*(ez(2,:) + ez_past(1,:,n-2)) + const_abc2*(ez_past(2,:,n-1) + ez_past(1,:,n-1)); % top boundary
%         ez(N_cells_x,:) = - ez_past(N_cells_x-1,:,n-2) + const_abc1*(ez(N_cells_x-1,:) + ez_past(N_cells_x,:,n-2)) + const_abc2*(ez_past(N_cells_x-1,:,n-1) + ez_past(N_cells_x,:,n-1)); % bottom boundary
% 
%         % Y boundaries
%         ez(:,1) = - ez_past(:,2,n-2) + const_abc1*(ez(:,2) + ez_past(:,1,n-2)) + const_abc2*(ez_past(:,2,n-1) + ez_past(:,1,n-1)); % left boundary
%         ez(:,N_cells_y) = - ez_past(:,N_cells_y-1,n-2) + const_abc1*(ez(:,N_cells_y-1) + ez_past(:,N_cells_y,n-2)) + const_abc2*(ez_past(:,N_cells_y-1,n-1) + ez_past(:,N_cells_y,n-1)); % right boundary
             
        for j=2:N_cells_y-1
            
          ez(1,j) = - ez_past(2,j,n-2) + const_abc1*(ez(2,j) + ez_past(1,j,n-2)) + const_abc2*(ez_past(2,j,n-1) + ez_past(1,j,n-1)) + const_abc3*(ez_past(2,j+1,n-1) - 2*ez_past(1,j,n-1) + ez_past(2,j-1,n-1) + ez_past(1,j+1,n-1) - 2*ez_past(2,j,n-1) + ez_past(1,j-1,n-1));
          ez(N_cells_x,j) = - ez_past(N_cells_x-1,j,n-2) + const_abc1*(ez(N_cells_x-1,j) + ez_past(N_cells_x,j,n-2)) + const_abc2*(ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j,n-1)) + const_abc3*(ez_past(N_cells_x-1,j+1,n-1) - 2*ez_past(N_cells_x,j,n-1) + ez_past(N_cells_x-1,j-1,n-1) + ez_past(N_cells_x,j+1,n-1) - 2*ez_past(N_cells_x-1,j,n-1) + ez_past(N_cells_x,j-1,n-1));
       
        end 

     
%         for i=2:N_cells_x-1
%                         
%           ez(i,1) = - ez_past(i,2,n-2) + const_abc1*(ez(i,2) + ez_past(i,1,n-2)) + const_abc2*(ez_past(i,2,n-1) + ez_past(i,1,n-1)) + const_abc3*(ez_past(i+1,2,n-1) - 2*ez_past(i,1,n-1) + ez_past(i-1,2,n-1) + ez_past(i+1,1,n-1) - 2*ez_past(i,2,n-1) + ez_past(i-1,1,n-1)); % left boundary
%           ez(i,N_cells_y) = - ez_past(i,N_cells_y-1,n-2) + const_abc1*(ez(i,N_cells_y-1) + ez_past(i,N_cells_y,n-2)) + const_abc2*(ez_past(i,N_cells_y-1,n-1) + ez_past(i,N_cells_y,n-1)) + const_abc3*(ez_past(i+1,N_cells_y-1,n-1) - 2*ez_past(i,N_cells_y,n-1) + ez_past(i-1,N_cells_y-1,n-1) + ez_past(i+1,N_cells_y,n-1) - 2*ez_past(i,N_cells_y-1,n-1) + ez_past(i-1,N_cells_y,n-1));
% 
%         end

      
        ez_past(:,:,n) = ez;
        n=n-2;
      

        surf(X, Y, ez);
        xlabel('x');
        ylabel('y');
        zlabel('Ez');
        xlim([1 N_cells_x]);
        ylim([1 N_cells_y]);
        zlim([-2 2]);
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
