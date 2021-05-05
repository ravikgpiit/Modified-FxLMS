

%--------------------------------------------------------------------------
clc, clear all; close all;
% close all;
fprintf('\n modified fxlms\n\n')


rand('seed',0);

% Initialization of the plant --------------------------------------------

% LTI system matrix A and x grid
nq = 400; % number of grid points

[A,x,I] = KS_init(nq);


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
% d= 0.5*sin(2*pi*t*50);


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


beta = 0.5;
% Time-loop
for i = 1:nt-1
    
    % outputs
    y(:,i) = Cy*v(:,i) + n(:,i);
    z(:,i) = Cz*v(:,i);
    
    % control forcing (Eqn. 106)
    if t(i) >= Ton
        u(:,i) = Kuy_lms(:,:,i) * y_buf(ifiruy)';
        u1(:,i)=(beta-1).*u(:,i) ; 
        u2(:,i)=beta.*u(:,i) ; 
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
rms_legend = ' \beta = 0.5';

plot_rms

% % kinetic energy of system
%  KE1 = (v.^2)/2;
%  KE = sqrt( mean(KE1(:,ista).^2,2) - mean(KE1(:,ista),2).^2 );
 
%  figure(50)
%  semilogy(v_rms,'r');
%  hold on; 
%  semilogy(KE,'b');
%  
%  figure(51)
%  KE_ucontrol = u.^2;
%  
%  plot(u,'r');hold on; plot(KE_ucontrol);
 
 
  d1 = fft(d);
  tt = 0:1:50000-1;   
  ff_d = (0:length(d1)-1)*50000/length(d1);
  figure(52)
  plot(ff_d,abs(d1/(length(d)/2)))
  xlabel('frequency (d)');
%     axis([0 2000 0 0.05]);

 
  y1 = fft(y);   
  ff_y = (0:length(y1)-1)*50000/length(y1);
  figure(53)
  plot(ff_y,abs(y1/(length(y)/2)))
  xlabel('frequency (y)');
%   axis([0 2000 0 0.5]);
 
  u1 = fft(u);   
  ff_u = (0:length(u1)-1)*50000/length(u1);
  figure(54)
  plot(ff_u,abs(u1/(length(u)/2)))
 xlabel('frequency (u)');
%   axis([0 2000 0 0.5]);

   
  z1 = fft(z);   
  ff_z = (0:length(z1)-1)*50000/length(z1);
  figure(55)
  plot(ff_z,abs(y1/(length(z)/2)))
 xlabel('frequency (z)');
%   axis([0 2000 0 0.5]);
%   Fs = 1000;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L1 = 1500;             % Length of signal
% t1 = (0:L1-1)*T;        % Time vector
% 
%   P2 = abs(Y/L1);
% P1 = P2(1:L1/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs*(0:(L1/2))/L1;

% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')


 
%  figure(56)
%  plot(Tot_KE);
 
 



% %% Figure 24 --------------------------------------------------------------
% figure(24); clf
% 
% plot(tker(ifiruy),Kuy_fir         ,'or',...
%      tker(ifiruy),Kuy_lms(1,:,end),'sb','MarkerSize',5);
% ax = axis; axis([tker(ifiruy([1,end])) ax(3:4)])
% xlabel('i \Delta t'); ylabel('K_{uy}(i)'); grid on
% 
% legend('LQG','FxLMS','Location','SW')




% Kernel time-evolution --------------------------------------------------
% figure; clf
% 
% subplot(5,1,2:3); imagesc(t,tker(ifiruy),squeeze(Kuy_lms)); hold on
%                   plot([Ton Ton], t(ifiruy([1,end])),'--b','Linewidth',2); hold off
%                   axis xy; axis([0 20000 t(ifiruy([1,end]))]);
%                   set(gca,'XTickLabel',[]); ylabel('x')
%                   
%                   caxis([-1,1]*1.5e-2); cax = caxis;
%                  
%                     
% subplot(5,1,1);   caxis(cax);
%                   pos = get(gca,'Position'); pos(4) = pos(4)/3;
%                   set(gca,'Visible','off','Position',pos)
%                   hc = colorbar('NO'); xlabel(hc,'v(x,t)')
%                   set(hc,'Position',pos);
%                   
% subplot(5,1,4:5); semilogy(t,mu,'-k','LineWidth',1); hold on
%                   ax = axis; axis([0 20000 ax(3:4)]);
%                   plot([Ton Ton], ax(3:4),'--b','Linewidth',2); hold on
%                   xlabel('t'), ylabel('\mu(t)'), grid on
% 
   figure
   plot(t,z,'k')
   xlabel('Time'); ylabel('z(t) error signal'); grid on
   
                  