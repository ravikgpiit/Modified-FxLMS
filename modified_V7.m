clc, clear all; 
close all;
fprintf('\nFxLMS algorithm.\n\n')


rand('seed',0);

% Initialization of the plant --------------------------------------------

% LTI system matrix A and x grid
nq = 400; % number of grid points

[A,x,I] = KS_init_V1(nq);


% Inputs matrix B

% disturbance d (Gaussian shape at x_d with sigma_d variance)
x_d = 35; sigma_d = 4; 
Bd  = exp(-((x-x_d)/sigma_d).^2)/sigma_d;

% actuator u (Gaussian shape at x_u with sigma_u variance)
x_u = 400; sigma_u = 4; 
Bu  = exp(-((x-x_u)/sigma_u).^2)/sigma_u;


% Outputs matrix C

% measurement y (Gaussian shape at x_y with sigma_y variance)
x_y = 300; sigma_y = 4; 
Cy  = (exp(-((x-x_y)/sigma_y).^2)/sigma_y).' * I;

% output z (Gaussian shape at x_z with sigma_z variance)
x_z = 700; sigma_z = 4; 
Cz  = (exp(-((x-x_z)/sigma_z).^2)/sigma_z).' * I;




% LQG compensator --------------------------------------------------------

% LQR controller
Wz = 1.0e+0; % z-output weight
Wu = 1.0e+0; % control penalty

% - solution of the control Riccati equation (Eqn. 66)
X = care(A,Bu,Cz'*Wz*Cz,Wu);

% - control gains matrix K (Eqn. 65)
K = - Wu\Bu'*X;


% Kalman filter
Rd = 1.0e+0; % d variance
Rn = 1.0e-2; % n variance

% - solution of the estimation Riccati equation (Eqn. 91)
Re = care(A',Cy',Bd*Rd*Bd',Rn);

% - estimation gains matrix L (Eqn. 90)
L = - Re*Cy'/Rn;




% External description (plant and compensator) ---------------------------
dt   = 1;
Tker = 2000;

tker = 0:dt:Tker;
nker = length(tker);


% Plant
Ad  = expm(A*dt); Bud = dt * Bu;
Czd = Cz;

% - convolution kernel Pzy (Eqn. 48)
Pvu = zeros(nq,nker);

Pvu(:,1) = Bud;
for i = 2:nker
    Pvu(:,i) = Ad*Pvu(:,i-1);
end

Pzu = Czd * Pvu;

% - FIR reduction (Eqn. 50)
toll = 1e-2; % kernel reduction tollerance
ifirzu = (find( abs(Pzu) > max(abs(Pzu))*toll ,1,'first')) : ...
         (find( abs(Pzu) > max(abs(Pzu))*toll ,1,'last'));
   
Pzu_fir = Pzu(ifirzu);


% Compensator
Acd  = expm((A + L*Cy + Bu*K)*dt);

% - convolution kernel Kzy (Eqn. 116)
Kvy = zeros(nq,nker);

Kvy(:,1) = -L * dt;
for i = 2:nker
    Kvy(:,i) = Acd*Kvy(:,i-1);
end

Kuy = K * Kvy;

% - FIR reduction (Eqn. 119)
toll = 1e-2; % kernel reduction tollerance
ifiruy = (find( abs(Kuy) > max(abs(Kuy))*toll ,1,'first')) : ...
         (find( abs(Kuy) > max(abs(Kuy))*toll ,1,'last'));
   
Kuy_fir = Kuy(ifiruy);

% - FIR kernel length
nfiruy = length(Kuy_fir);




% FxLMS ------------------------------------------------------------------
mumax = 1e-6; % max step length




% Actuator displacement --------------------------------------------------
x_u = x_u + 5; 
Bu = exp(-((x-x_u)/sigma_u).^2)/sigma_u;



        
% Time integration -------------------------------------------------------

% Parameters
dt   = 1;
Tend = 50000;

Ton = 4000;

t  = 0:dt:Tend;
nt = length(t);


% Time-stepper matrices (Crank-Nicholson)
CNi = (speye(nq) - dt/2*A);
CNe = (speye(nq) + dt/2*A);


% Noises inizialization
d = randn(1,nt);     % unitary variance Gaussian white noise
d = d - mean(d,2);   % enforce zero-mean
d = d / std(d,[],2); % enforce unitary variance
% d=d+0.5*sin(2*pi*t*50);


n = randn(1,nt);     % unitary variance Gaussian white noise
n = n - mean(n,2);   % enforce zero-mean
n = n./ std(n,[],2) * 1e-1;


% Variables initialization
v = zeros(nq,nt); % velocity
y = zeros(1 ,nt); % measurement signal
z = zeros(1 ,nt); % output signal
u = zeros(1 ,nt); % control signal
u1 = zeros(1 ,nt); % control signal
u2 = zeros(1 ,nt); % control signal

Kuy_lms = zeros(1,nfiruy,nt);   % identified convolution kernel
y_buf   = zeros(1,ifirzu(end)); % measurement signal buffer
 u_buf= zeros(1,ifirzu(end)); % measurement signal buffer
yf_buf  = zeros(1,ifiruy(end)); % filtered measurement signal buffer
uf_buf  = zeros(1,1); % filtered measurement signal buffer
mu      = zeros(1,nt);          % step length
lambda  = zeros(1,nfiruy);      % step direction

Kuy_lms(1,:,1) = Kuy_fir;       % initial guess (LQG)
% Kuy_lms(1,:,1) = zeros(1,500);       % initial guess (LQG)


beta = 0;
% Time-loop
for i = 1:nt-1
    
    % outputs
    y(:,i) = Cy*v(:,i) + n(:,i);
    z(:,i) = Cz*v(:,i);
    
    % control forcing (Eqn. 106)
    if t(i) >= Ton
        u(:,i) = Kuy_lms(:,:,i) * y_buf(ifiruy)';
        u1(:,i)=(beta-1).*u(:,i) ; 
        u2(:,i)=(beta).*u(:,i) ; 
    end
    
    % saturation (uncomment to activate)
    % u(:,i) = min([max([-2,u(:,i)]),2]);
    
    % update measurement buffer
    y_buf  = [y(:,i) y_buf(:,1:end-1)];
    u_buf  = [u2(:,i) u_buf(:,1:end-1)];
    % update filtered measurement buffer
    yf_buf = [Pzu_fir * y_buf(ifirzu)', yf_buf(:,1:end-1)];
    uf_buf = Pzu_fir * u_buf(ifirzu)';
    % timestep
    v(:,i+1) = CNi\(CNe*v(:,i) + dt * Bd*d(:,i) + dt * Bu*u1(:,i));
        
    % FxLMS algorithm
    if t(i) >= Ton
    % - update direction (Eqn. 102)
        lambda = -2 * (uf_buf-z(:,i))* yf_buf(ifiruy);
        
    % - step length (Eqn. 107-108)
        mu(i) = - (uf_buf-z(:,i)) / (lambda * yf_buf(ifiruy)');
        mu(i) = min([mu(i), mumax/var(yf_buf)]);
        
    % - kernel update (Eqn. 101)
        Kuy_lms(1,:,i+1) = Kuy_lms(1,:,i) + mu(i)*lambda;
    else
        Kuy_lms(:,:,i+1) = Kuy_lms(:,:,i);
    end
    
    % runtime output
     runtime_output

end




% Velocity statistics ----------------------------------------------------
Tsta = 20000;
ista = (Tend-Tsta <= t) & (t <= Tend);

% Root-Mean-Square
v_rms = sqrt( mean(v(:,ista).^2,2) - mean(v(:,ista),2).^2 );

% visualization
rms_style  = '-b';
rms_legend = ' \beta = 0';

plot_rms
grid on;
%%
