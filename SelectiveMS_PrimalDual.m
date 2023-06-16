function u = SelectiveMS_PrimalDual(im,lambda,mu,eta,gd)
%--------------------------------------------------------------------------
%   This function implements the primal dual scheme for model one in:
%   "On Two Convex Variational Models and Their Iterative Solutions for 
%   Selective Segmentation of Images with Intensity Inhomogeneity"    
%
%   Inputs: 
%       - im: input image
%       - lambda: fitting parameter
%       - mu: L2 regularisation parameter
%       - eta: distance function parameter
%       - gd: geodesic distance function
%
%   Outputs: 
%       - u: minimiser of model one.
% 
%   Code by: Liam Burrows
%   Last updated: 16/06/2023
%--------------------------------------------------------------------------

%Parameters for primal dual scheme
tau = 0.01;
L_squared = 8;
sigma = 1/(tau*L_squared);
theta = 0.9; % theta = [0,1]



[n,m]= size(im);
px = zeros(n,m); py = px;
u = im; uhat = im;

uker = 0*im;
uker(1,1) = 4;uker(1,2)=-1;uker(2,1)=-1;
uker(end,1)=-1;uker(1,end)=-1;  

% lhs of fft
uker = 1+tau*lambda+(mu*tau)*fft2(uker);

for iters=1:500
    
    px_arg = px + sigma*Dx(uhat);
    py_arg = py + sigma*Dy(uhat);
    
    norm_p = sqrt(px_arg.^2 + py_arg.^2);
    
    px = px_arg./(max(norm_p,1));
    py = py_arg./(max(norm_p,1));
    
    
    u0 = u;
    u_arg = u - tau*(Dxt(px)+Dyt(py));
    %u = (u_arg + tau*lambda*im)./(1+tau*lambda);
    rhs = lambda.*tau.*im - eta*tau*gd + u_arg;
    u = real(ifft2(fft2(rhs)./uker));
    
    
    uhat = u + theta*(u-u0);
    
    err=norm(u-u0,'fro')/norm(u,'fro');
    
    if mod(iters,10)==0
        disp(['iterations: ' num2str(iters) '!  ' 'error is:   ' num2str(err)]);
    end
    
    % check the stopping criterion
    if err<10^(-4)
        break;
    end
    
end



end



function d = Dx(u)
%--------------------------------------------------------------------------
%   This program is computing the gradient operator according to x-axis 
%   using backward difference.
%    
%   Usage: d = Dx(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: gradient u according to x-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------


[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,2:cols) = u(:,2:cols)-u(:,1:cols-1);

% use periodic boundary
d(:,1) = u(:,1)-u(:,cols);

% %use reflective boundary
% d(:,1) = 0;
% d(:,1) = d(:,2);
end

function d = Dxt(u)
%--------------------------------------------------------------------------
%   This program is computing the transpose of the gradient operator 
%   according to x-axis using backward difference.
%    
%   Usage: d = Dxt(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: transpose gradient u according to x-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(:,1:cols-1) = u(:,1:cols-1)-u(:,2:cols);

% use periodic boundary
d(:,cols) = u(:,cols)-u(:,1);

% %use reflective boundary
% d(:,cols)=0;
% d(:,cols)=d(:,cols-1);
end

function d = Dy(u)
%--------------------------------------------------------------------------
%   This program is computing the gradient operator according to y-axis 
%   using backward difference.
%    
%   Usage: d = Dy(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: gradient u according to y-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(2:rows,:) = u(2:rows,:)-u(1:rows-1,:);

% use periodic boundary
d(1,:) = u(1,:)-u(rows,:);

% %use reflective boundary
% d(1,:) = 0;
% d(1,:) = d(2,:);
end

function d = Dyt(u)
%--------------------------------------------------------------------------
%   This program is computing the transpose of the gradient operator 
%   according to y-axis using backward difference.
%    
%   Usage: d = Dyt(u);
%
%   Inputs: 
%       - u: 2d data
%
%   Outputs: 
%       - d: transpose gradient u according to y-axis
% 
%   Code by: Xiaohao Cai
%   Last updated: 10/11/2012 
%--------------------------------------------------------------------------

[rows,cols] = size(u); 
d = zeros(rows,cols);
d(1:rows-1,:) = u(1:rows-1,:)-u(2:rows,:);

% use periodic boundary
d(rows,:) = u(rows,:)-u(1,:);

% %use reflective boundary
% d(rows,:)=0;
% d(rows,:)=d(rows-1,:);
end