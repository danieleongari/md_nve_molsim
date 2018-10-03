clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete integration of the Equation of Motion
% for an harmonic oscillator 1d
%
% Kinetic energy K = 0.5 * m * v^2 = 0.5 * p^2 /m where p = mv
% Potential energy, linear spring V = 0.5* k * x^2 = int(Fdx) with F= - k * x
%
% harmonic constant k is conventionally expressed as k = m * omega^2
% because the frequency of the oscillationg system is f=omega/2*pi
%
% Daniele Ongari 22/10/15 (17/10/16)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS SYSTEM
% parameters for spring costant k
k = 1;          % [mass]/[time]^2 = [force]/[length]

% mass of the object
m = 1;          % [mass]

% Starting position and velocity
x0 = 0.5;       % [length]
v0 = 0.3;       % [lenght]/[time]

%% PARAMETERS SIMULATION

dt = 0.01;         % [time]   
t_tot = 20;       % [time]
steps = t_tot/dt; %[-]

chose_integrator=1; %Euler (1), Verlet (2), Leap frog (3), Velocity verlet(4)

%% START SIMULATION
tnow=0;
xnow=x0;
vnow=v0;

for i=1:steps
  %% Integration of the Equation of Motion
      
  % t-dt	t-dt/2		t       t+dt/2      t+dt
  % xold				xnow    		    xnext
  % vold	voldhalf	vnow    vnexthalf   vnext
  % Fold   		  		Fnow                Fnext
  
  	tnext=tnow+dt;

 if chose_integrator==1 % EULER scheme  [WARNING: unavoidable energy drift! You can only reduce it with small dt]
	Fnow=-k*xnow;
	xnext=xnow+vnow*dt+0.5*Fnow/m*dt^2; % +O(dt^3)
	vnext=vnow+Fnow/m*dt;          % +O(dt^2)   
    disp('I will never use Euler again to integrate the EoM!!!')
 end
 
 if chose_integrator==2 % VERLET scheme 
	if i==1
		xold=xnow-vnow*dt;
    end
	Fnow=-k*xnow;
	xnext=2*xnow-xold+Fnow/m*dt^2; % +O(dt^4)
	vnow=(xnext-xold)/(2*dt);      % +O(dt^2)
 end

if chose_integrator==3 % LEAP FROG scheme [WARNING: H=K(v)+V(x) will be wrong!]
	if i==1
		xold=xnow-vnow*dt;
		voldhalf=(xnow-xold)/dt;
    end
	Fnow=-k*xnow;
	vnexthalf=voldhalf+Fnow/m*dt;
	xnext=xnow+vnexthalf*dt;
    
    vnow=vnexthalf; %WRONG! just to show the picture
end

 if chose_integrator==4 % VELOCITY VERLET scheme
    if i==1
        Fnow=-k*xnow; % [mass]*[length]/[time]^2            F=dV/dx
    end
    xnext=xnow+vnow*dt+0.5*Fnow/m*dt^2; % +O(dt^3)
    Fnext=-k*xnext;    
    vnext=vnow + 0.5*(Fnext+Fnow)/m*dt; % +O(dt^3)
 end   
	
    %% Computing Energies and errors
    
    % Energies
    K = 0.5* m * vnow^2;
    V = 0.5* k * xnow^2;
    H = V+K; %Hamiltonian
    %L = V-K; %Lagrangian
    
    % Compute the max error on H
    if i==1
        Hmin=H;
        Hmax=H;
    end
    
    if H>Hmax
        Hmax=H;
    end
    if H<Hmin
        Hmin=H;
    end
        
    %% Plotting graphs
    figure(1)
    subplot(2,1,1); hold on; plot( [tnow tnext] , [xnow xnow], 'b')
    xlabel('time')
    ylabel('space (equil @ x=0)')
    subplot(2,1,2); hold on; plot([tnow tnext],[K K],'g') %Kinetic energy > green 
    subplot(2,1,2); hold on; plot([tnow tnext],[V V],'b') %Potential energy > blue
    subplot(2,1,2); hold on; plot([tnow tnext],[H H],'r') %Hamiltonian > red
    %subplot(2,1,2); hold on; plot([tnow tnext],[L L],'k') %Lagrangian > black
       figure(2)
       title('phasespace')
       xlabel('space')
       ylabel('velocity')
       plot(xnow,vnow,'r-o')
       hold on

    %% Updating 
    tnow=tnext;
    
 	xold=xnow;
    xnow=xnext;  
    
    if chose_integrator==3 %Leap frog
    voldhalf=vnexthalf;
    elseif chose_integrator==2; %Verlet
    else
    vold=vnow;
    vnow=vnext;
    end
  
    if chose_integrator==4 %Velocity Verlet
    Fnow=Fnext;
    end
	
    pause(0.005)
end

%% Final Output

disp(['Simulation Parameters and Results'])
disp ' '
disp(['k = ' num2str(k) ' N/m'])
disp(['m = ' num2str(m) ' kg'])    
disp(['x0 = ' num2str(x0) ' m']) 
disp(['v0 = ' num2str(v0) ' m/s'])
disp ' '

Hmed=(Hmin+Hmax)/2;
Herr=Hmax-Hmin;

disp(['Hmed = ' num2str(Hmed) ' N*m'])
disp(['Herr = ' num2str(Herr) ' N*m'])



