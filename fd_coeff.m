function d = fd_coeff(n, varargin)
%--------------------------------------------------------------------------

if isempty(varargin)
    dx = 1;
else
    dx = varargin{1};
end

% stencil size
nd = ceil(n/2);
nx = 2*nd + 1;

% initialization
d  = zeros(n,nx,nd+1);

% shifting (j=0 centered scheme, j=1,2... forward)
for j = 0:nd

    % grid points
    x  = ((-ceil(n/2):ceil(n/2))+j)*dx;
    
    F0 = zeros(nx,nx);
    for i = 1:nx
        F0(i,:) = x.^(i-1);
    end
    
    % k-th derivative
    for k = 1:n
        b = zeros(nx,1); b(k+1) = factorial(k);
        d(k,:,j+1) = F0\b;
    end

end