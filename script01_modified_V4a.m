clc, clear all; close all;
nq = 400; % number of grid points % 400 points

[A,x,I] = KS_init_V4(nq);
% A = full(A); %making sparse matrix into full matrix
% I = full(I); %making sparse matrix into full matrix

% Inputs matrix B

% disturbance d (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
Bd  = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% actuator u (Gaussian shape at x_u with sigma_u variance)
x_u = 300; sigma_u = 4; 
Bu  = exp(-((x-x_u)/sigma_u).^2)/sigma_u;


% Outputs matrix C

% measurement y (Gaussian shape at x_y with sigma_y variance)
x_y = 300; sigma_y = 4; 
Cy  = (exp(-((x-x_y)/sigma_y).^2)/sigma_y).' * I;

% output z (Gaussian shape at x_z with sigma_z variance)
x_z = 700; sigma_z = 4; 
Cz  = (exp(-((x-x_z)/sigma_z).^2)/sigma_z).' * I;
   
% Time integration -------------------------------------------------------

% Parameters
dt   = 1;
Tend = 50000;

t  = 0:dt:Tend;
nt = length(t);


% Time-stepper matrices (Crank-Nicholson)
CNi = (speye(nq) - dt/2*A);
CNe = (speye(nq) + dt/2*A);


% % Noise inizialization
d = randn(1,nt);     % unitary variance Gaussian white noise
d = d - mean(d,2);   % enforce zero-mean
d = d / std(d,[],2); % enforce unitary variance

% Variables initialization
v = zeros(nq,nt); % velocity
y = zeros(1 ,nt); % measurement signal
z = zeros(1 ,nt); % output signal

% Initial condition (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
v(:,1) = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% Time-loop
for i = 1:nt-1
    
    % outputs
    y(:,i) = Cy*v(:,i);
    z(:,i) = Cz*v(:,i);
    
    % timestep
    v(:,i+1) = CNi\(CNe*v(:,i) + dt * Bd*d(:,i));
    
%     v(:,i+1) = CNi\(CNe*v(:,i));
    % runtime output
%     runtime_output

end




% % Velocity statistics ----------------------------------------------------
Tsta = 20000;
ista = (Tend-Tsta <= t) & (t <= Tend);

% Root-Mean-Square
v_rms = sqrt( mean(v(:,ista).^2,2) - mean(v(:,ista),2).^2 );

% visualization
% rms_style  = '-b';
% rms_legend = 'P = 100';
% 
% plot_rms1

figure(1001)
semilogy(x,v_rms,'-b','linewidth',2)

set(gca,'YScale','log')
axis([0 800 0.33 0.5*1e2])
xlabel('Streamwise x location'); ylabel('u_{rms}'); grid on
legend('C_{2} = 0.05 (at C_{1} = 0.25 , U = 0.4)')
hold on

%%
nq = 400; % number of grid points % 400 points

[A,x,I] = KS_init_V5(nq);
% A = full(A); %making sparse matrix into full matrix
% I = full(I); %making sparse matrix into full matrix

% Inputs matrix B

% disturbance d (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
Bd  = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% actuator u (Gaussian shape at x_u with sigma_u variance)
x_u = 300; sigma_u = 4; 
Bu  = exp(-((x-x_u)/sigma_u).^2)/sigma_u;


% Outputs matrix C

% measurement y (Gaussian shape at x_y with sigma_y variance)
x_y = 300; sigma_y = 4; 
Cy  = (exp(-((x-x_y)/sigma_y).^2)/sigma_y).' * I;

% output z (Gaussian shape at x_z with sigma_z variance)
x_z = 700; sigma_z = 4; 
Cz  = (exp(-((x-x_z)/sigma_z).^2)/sigma_z).' * I;
   
% Time integration -------------------------------------------------------

% Parameters
dt   = 1;
Tend = 50000;

t  = 0:dt:Tend;
nt = length(t);


% Time-stepper matrices (Crank-Nicholson)
CNi = (speye(nq) - dt/2*A);
CNe = (speye(nq) + dt/2*A);


% % Noise inizialization
d = randn(1,nt);     % unitary variance Gaussian white noise
d = d - mean(d,2);   % enforce zero-mean
d = d / std(d,[],2); % enforce unitary variance

% Variables initialization
v = zeros(nq,nt); % velocity
y = zeros(1 ,nt); % measurement signal
z = zeros(1 ,nt); % output signal

% Initial condition (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
v(:,1) = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% Time-loop
for i = 1:nt-1
    
    % outputs
    y(:,i) = Cy*v(:,i);
    z(:,i) = Cz*v(:,i);
    
    % timestep
    v(:,i+1) = CNi\(CNe*v(:,i) + dt * Bd*d(:,i));
    
%     v(:,i+1) = CNi\(CNe*v(:,i));
    % runtime output
%     runtime_output

end




% % Velocity statistics ----------------------------------------------------
Tsta = 20000;
ista = (Tend-Tsta <= t) & (t <= Tend);

% Root-Mean-Square
v_rms = sqrt( mean(v(:,ista).^2,2) - mean(v(:,ista),2).^2 );

% visualization
% rms_style  = '-b';
% rms_legend = 'P = 100';
% 
% plot_rms1

figure(1001)
semilogy(x,v_rms,'-r','linewidth',2)

set(gca,'YScale','log')
axis([0 800 0.33 0.5*1e2])
xlabel('Streamwise x location'); ylabel('u_{rms}'); grid on
legend('C_{2} = 0.1 (at C_{1} = 0.25 , U = 0.4)')
hold on

%%
nq = 400; % number of grid points % 400 points

[A,x,I] = KS_init_V6(nq);
% A = full(A); %making sparse matrix into full matrix
% I = full(I); %making sparse matrix into full matrix

% Inputs matrix B

% disturbance d (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
Bd  = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% actuator u (Gaussian shape at x_u with sigma_u variance)
x_u = 300; sigma_u = 4; 
Bu  = exp(-((x-x_u)/sigma_u).^2)/sigma_u;


% Outputs matrix C

% measurement y (Gaussian shape at x_y with sigma_y variance)
x_y = 300; sigma_y = 4; 
Cy  = (exp(-((x-x_y)/sigma_y).^2)/sigma_y).' * I;

% output z (Gaussian shape at x_z with sigma_z variance)
x_z = 700; sigma_z = 4; 
Cz  = (exp(-((x-x_z)/sigma_z).^2)/sigma_z).' * I;
   
% Time integration -------------------------------------------------------

% Parameters
dt   = 1;
Tend = 50000;

t  = 0:dt:Tend;
nt = length(t);


% Time-stepper matrices (Crank-Nicholson)
CNi = (speye(nq) - dt/2*A);
CNe = (speye(nq) + dt/2*A);


% % Noise inizialization
d = randn(1,nt);     % unitary variance Gaussian white noise
d = d - mean(d,2);   % enforce zero-mean
d = d / std(d,[],2); % enforce unitary variance

% Variables initialization
v = zeros(nq,nt); % velocity
y = zeros(1 ,nt); % measurement signal
z = zeros(1 ,nt); % output signal

% Initial condition (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
v(:,1) = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% Time-loop
for i = 1:nt-1
    
    % outputs
    y(:,i) = Cy*v(:,i);
    z(:,i) = Cz*v(:,i);
    
    % timestep
    v(:,i+1) = CNi\(CNe*v(:,i) + dt * Bd*d(:,i));
    
%     v(:,i+1) = CNi\(CNe*v(:,i));
    % runtime output
%     runtime_output

end




% % Velocity statistics ----------------------------------------------------
Tsta = 20000;
ista = (Tend-Tsta <= t) & (t <= Tend);

% Root-Mean-Square
v_rms = sqrt( mean(v(:,ista).^2,2) - mean(v(:,ista),2).^2 );

% visualization
% rms_style  = '-b';
% rms_legend = 'P = 100';
% 
% plot_rms1

figure(1001)
semilogy(x,v_rms,'-g','linewidth',2)

set(gca,'YScale','log')
axis([0 800 0.33 0.5*1e2])
xlabel('Streamwise x location'); ylabel('u_{rms}'); grid on
legend('C_{2} = 0.025 (at C_{1} = 0.25 , U = 0.4)')
hold on
