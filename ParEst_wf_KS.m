%%%% Integral Method Kuramoto-Sivashinsky Parameter Estimation %%%%

%% PSEUDO-CODE
%{

INPUTS 
------
filename : string for path to -v7.3 .mat file containing data
N_d : number of integration domains
D = [fr, Dt] : size of integration domain

OUTPUTS
-------
par : estimated parameters in physical units
ksi : estimated parameters in non-dimensinoal units
res : residual of estimation (utility unclear)

Initialize
    -size of local domain
    -any tracking variables
    -library and target

Fill terms
    -advection term
    -laplacian term
    -rayleigh term
    -time-derivative term

Regression
    -invert library onto target
    -calculate norm of b-\Theta\xi

Reminders
    -disp to-do list

weight_full(k)
    -inputs: k = [kx,ky,kt], order of derivative(s)
    -output: 3D array corresponding to local sub-domain

weight_poly(x,m,k)
    -inputs: absiscae, polynomial order, derivative order
    -output: additive windowing polynomial

%}
%%
function [ksi,res,Q,q0] = ParEst_wf_KS(filename,N_d,D,if_track,if_symreg)

%% INITIALIZE
% Define matfile
traj = matfile(filename);

% Get size of velocity fields
[Lx,Lt] = size(traj,'uu');

% Grid densities
dt = traj.dt;
dx = traj.dx; 

% % Size of local domain
clearvars -global
global var
var.Dx = D(1);
var.Dt = D(2);

% Create subdomain
var.x = linspace(-1,1,var.Dx+1); 
var.t = linspace(-1,1,var.Dt+1);

% Define variable conversions
S_x = 2/(dx*var.Dx);
S_t = 2/(dt*var.Dt);

% Time sampling scheme
P = zeros(2,N_d);
P(1,:) = randi([1,Lx-var.Dx],N_d,1);
P(2,:) = randi([1,Lt-var.Dt],N_d,1);

% Initialize Target and Library
q0 = zeros(N_d*N_h(1)*N_h(2),1);
Q = zeros(length(q0),8);

%% FILL TERMS

n_lib = 0;
n_track = 10;
            
% Pre-make derivatives of windowing functions
dA00 = weight_full([0,0]);
dA01 = weight_full([0,1]);
dA10 = weight_full([1,0]);
dA20 = weight_full([2,0]);
dA40 = weight_full([4,0]);
dA30 = weight_full([3,0]);

for np = 1:N_d  

    n_lib = n_lib + 1;

    if if_track && n_lib == n_track
        if n_lib < 100
            disp(['Library Row # : ',num2str(n_lib)])
            n_track = n_track + 10;
        elseif n_lib < 1000
            disp(['Library Row # : ',num2str(n_lib)])
            n_track = n_track + 100;
        else
            disp(['Library Row # : ',num2str(n_lib)])
            n_track = n_track + 1000;
        end
    end

    % Indices for integration domain
    rx = P(1,np):(P(1,np)+var.Dx);
    rt = P(2,np):(P(2,np)+var.Dt);

    % Velocity fields on integration domain
    U = traj.uu(rx,rt); 

    % Target
    B = U.*dA01*S_t; 
    q0(n_lib,1) = trapz(var.x,trapz(var.t,B,2),1);

    % Advection Term 
    th1 = -(1/2)*U.^2.*dA10*S_x; 
    Q(n_lib,1) = trapz(var.x,trapz(var.t,th1,2),1);

    % Laplacian Term
    th2 = U.*dA20*S_x^2;
    Q(n_lib,2) = trapz(var.x,trapz(var.t,th2,2),1);

    % Biharmonic Term
    th3 = U.*dA40*S_x^4; 
    Q(n_lib,3) = trapz(var.x,trapz(var.t,th3,2),1);

    % Linear Term
    th4 = U.*dA00;
    Q(n_lib,4) = trapz(var.x,trapz(var.t,th4,2),1);

    % First Order Derivative
    th5 = U.*dA10*S_x;
    Q(n_lib,5) = trapz(var.x,trapz(var.t,th5,2),1);

    % Third Order Derivative
    th6 = U.*dA30*S_x^3;
    Q(n_lib,6) = trapz(var.x,trapz(var.t,th6,2),1);

    % Quadratic Term
    th7 = U.^2.*dA00;
    Q(n_lib,7) = trapz(var.x,trapz(var.t,th7,2),1);

    % Cubic Term
    th8 = U.^3.*dA00;
    Q(n_lib,8) = trapz(var.x,trapz(var.t,th8,2),1);

end
            


%% REGRESSION
% Parameters
if if_symreg
    ksi = SINDy(Q, q0); % sparsify library
else
    ksi = Q(:,1:3) \ q0;
end

% How good was the estimation
res = norm(q0 - Q*ksi);

%% REMINDERS

disp(' ')
disp(' ')

end
%% --------------------------------------------------------------- SINDy()
function Xi = SINDy (Theta, dXdt)
%{
Compute sparse regression on dX = Theta * Xi
Regression technique used: sequential least squares

Modified procedure from:
S. H. Rudy, S. L. Brunton, J. L. Proctor, J. N. Kutz, Data-driven 
discovery of partial differential equations. Sci. Adv. 3, e1602614 (2017)
%}

Xi = Theta \ dXdt;

gamma = 0.05;
lambda = gamma*mean(abs(dXdt)); % threshold to determine term as "small"
for i = 1:20

  product = zeros(size(Xi)); 
  [~,w] = size(Theta);
  for p_ind = 1:w
    product(p_ind) = Xi(p_ind)*mean(abs(Theta(:,p_ind)));
  end

  smallinds = (abs(product) < lambda);
  Xi(smallinds) = 0;                        % set negligible terms to 0
  for ind = 1:size(dXdt,2)   
    biginds = ~smallinds(:,ind);
    Xi(biginds,ind) = Theta(:,biginds) \ dXdt(:,ind);
  end
end
    
end
%% Assorted functions
function W = weight_full(k)
global var
%{
Assemble the 1D weight functions into the full weight

k = [kx,ky,kt]: order of derivative(s)
%}

wx = weight_poly(var.x,4,k(1));
wt = weight_poly(var.t,1,k(2));

[wT,wX] = meshgrid(wt,wx);

W = wX.*wT;

end
function p = weight_poly(x,m,k)
%{
Polynomial piece of weighting function used to satisfy BC

A = d^k/dx^k[ (x^2 - 1)^m ]

x: independent variable
m: power of base function
k: order of derivative
%}

a = zeros(m*2 + 1,1); % initial coefficent vector
for l = 0:m
    a(2*l+1) = (-1)^(m-l)*nchoosek(m,l); % set polynomial coefficients
end 

c = zeros(2*m+1,1); % final coefficient vector
for n = 0:(2*m - k)
    c(n+1) = a(n+1+k)*factorial(n+k)/factorial(n);
end

p = 0;
for n = 0:(2*m-k)
    p = p + c(n+1)*x.^n; % final windowing function
end


end