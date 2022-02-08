function [A,x,I] = KS_init_V6(nq)
% clc, clear all;close all;
% nq=400;
%--------------------------------------------------------------------------
%
% [A,x,I,k] = KS_init(nq)
%
%   Linearised Kuramoto-Sivashinsky equation for v(x,t):
%
%       v_t = - V * v_x - 1/R * (P * v_xx + v_xxxx) ;
%
%       b.c.:   v(0,t)     = 0
%               v_x(0,t)   = 0
%               v_x(L,t)   = 0
%               v_xxx(L,t) = 0
%
%   A finite-difference scheme is used to approximate the spatial
%   derivatives
%   The resultant evolution equation for the nodal values is:
%
%       v_t = A * v.
%
% 

%
%--------------------------------------------------------------------------


% KS Parameters for validation with rms data fabbiane
%     R = 1000;
%     P = 35;
%     V = 0.5;    
%     % KS Parameters in general
    R = 0.15; % C1
    P = 0.05; % C2
    V = 0.4;  % U
% Domain length
    L = 800;

    
% Grid
    x  = linspace(0,L,nq+1).';
    dx = x(2)-x(1);
    
% - adding boundary nodes (to be removed)
    x  = [-dx; x; L + [1;2]*dx];
    nx = length(x);

    
% Integral 
    e = ones(nx,1);
    I = spdiags(e*dx,0,nx,nx)/2;

    
% Derivatives
    d  = fd_coeff(4,dx(1));
    nd = (size(d,2)-1)/2;
    
    %inner points
    D0 = speye(nx,nx);
    D1 = spdiags(-e*d(1,:,2),-1-(-nd:nd),nx,nx); %up-wind fd scheme
    D2 = spdiags( e*d(2,:,1),-nd:nd,nx,nx);
    D3 = spdiags( e*d(3,:,1),-nd:nd,nx,nx);
    D4 = spdiags( e*d(4,:,1),-nd:nd,nx,nx);
    %left-boundary
    for j = 1:nd
        D1(j,1:2*nd+1) = d(1,:,nd-j+2);
        D2(j,1:2*nd+1) = d(2,:,nd-j+2);
        D3(j,1:2*nd+1) = d(3,:,nd-j+2);
        D4(j,1:2*nd+1) = d(4,:,nd-j+2);
    end
    D1(nd+1,1:2*nd+1) = d(1,:,1); %because of the up-wind scheme
    %left-boundary
    for j = 1:nd
        D1(nx+1-j,nx+1-(1:2*nd+1)) =-d(1,:,nd-j+2);
        D2(nx+1-j,nx+1-(1:2*nd+1)) = d(2,:,nd-j+2);
        D3(nx+1-j,nx+1-(1:2*nd+1)) =-d(3,:,nd-j+2);
        D4(nx+1-j,nx+1-(1:2*nd+1)) = d(4,:,nd-j+2);
    end
    
   
% Linearized KS operator
    A = - V * D1 - P/R *D2 - 1/R * D4;

    
% Boundary conditions
    bc(1,:) = D0(2,:);  % v    (0,t) = 0
    bc(2,:) = D1(2,:);  % v_x  (0,t) = 0
    bc(3,:) = D1(nx,:); % v_x  (L,t) = 0
    bc(4,:) = D3(nx,:); % v_xxx(L,t) = 0
    
    r = [1:2,nx-1:nx]; % nodes to remove (boundary nodes)
    k = 3:nx-2;        % nodes to keep
    
    BC = speye(nx);
    BC(r,r) = 0;
    BC(r,k) = -bc(:,r)\bc(:,k);
    
    
% Apply boundary conditions
    A = BC*A*BC; A = A(k,k);
    I =    I*BC; I = I(k,k);
                 x = x(k);
    