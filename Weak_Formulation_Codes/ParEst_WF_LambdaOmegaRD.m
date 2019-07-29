%%%% Lambda-Omega RD Weak Formulation System Identification %%%%%
function [ksi,res,Q,q0] = ParEst_WF_LambdaOmegaRD(filename,N_d,D,if_track)
%{
INPUTS:
-------
filename                -- string containing path to matfile w/ data
N_d                     -- number of integration domains to saple
D = [Dx,Dt]             -- size of integration domain in space/time
if_track                -- bool for showing progress in command window

OUTPUTS:
--------
ksi     -- non dimensional estimated parameters
res     -- residual = mean(abs(Q*ksi-q0))
Q       -- Assembled library
q0      -- time derivative term
%}

%% INITIALIZE
% Define matfile
traj = matfile(filename);

% Get size of velocity fields
[~,Lx,Lt] = size(traj,'U_t');

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

% Initialize integration domain location
P = zeros(3,1);

% Initialize Target and Library
q0 = zeros(N_d,2);
Q = zeros(length(q0),10,2);

%% FILL TERMS
n_lib = 0;
n_track = 10;
            
% Pre-make derivatives of weight functions
dA000 = weight_full([0,0,0]);
dA001 = weight_full([0,0,1]);
dA200 = weight_full([2,0,0]);
dA020 = weight_full([0,2,0]);
        
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

    % Choose random point
    P(1) = randi([1,Lx-var.Dx]);
    P(2) = randi([1,Lx-var.Dx]);
    P(3) = randi([1,Lt-var.Dt]);

    % Indices for integration domain
    rx = P(1):(P(1) + var.Dx);
    ry = P(2):(P(2) + var.Dx);
    rt = P(3):(P(3) + var.Dt);

    % U in integration domain
    U = traj.U_t(ry,rx,rt); 
    V = traj.V_t(ry,rx,rt);

    %%%%% U %%%%%
    % du/dt
    B = -U.*dA001*S_t; 
    q0(n_lib,1) = trapz(var.x,trapz(var.x,trapz(var.t,B,3),2),1);

    % \nabla^2 u
    th1 = U.*(dA020 + dA200)*S_x^2;
    Q(n_lib,1,1) = trapz(var.x,trapz(var.x,trapz(var.t,th1,3),2),1);

    % u
    th2 = U.*dA000;
    Q(n_lib,2,1) = trapz(var.x,trapz(var.x,trapz(var.t,th2,3),2),1);

    % u^2
    th3 = U.^2.*dA000;
    Q(n_lib,3,1) = trapz(var.x,trapz(var.x,trapz(var.t,th3,3),2),1);

    % u^3
    th4 = U.^3.*dA000;
    Q(n_lib,4,1) = trapz(var.x,trapz(var.x,trapz(var.t,th4,3),2),1);

    % v
    th5 = V.*dA000;
    Q(n_lib,5,1) = trapz(var.x,trapz(var.x,trapz(var.t,th5,3),2),1);

    % v^2
    th6 = V.^2.*dA000;
    Q(n_lib,6,1) = trapz(var.x,trapz(var.x,trapz(var.t,th6,3),2),1);

    % v^3
    th7 = V.^3.*dA000;
    Q(n_lib,7,1) = trapz(var.x,trapz(var.x,trapz(var.t,th7,3),2),1);

    % uv
    th8 = U.*V.*dA000;
    Q(n_lib,8,1) = trapz(var.x,trapz(var.x,trapz(var.t,th8,3),2),1);

    % u^2v
    th9 = U.^2.*V.*dA000;
    Q(n_lib,9,1) = trapz(var.x,trapz(var.x,trapz(var.t,th9,3),2),1);

    % uv^2
    th10 = U.*V.^2.*dA000;
    Q(n_lib,10,1) = trapz(var.x,trapz(var.x,trapz(var.t,th10,3),2),1);

    %%%%% V %%%%%
    % dv/dt
    B = -V.*dA001*S_t;
    q0(n_lib,2) = trapz(var.x,trapz(var.x,trapz(var.t,B,3),2),1);

    % \nabla^2 v
    th1 = V.*(dA020 + dA200)*S_x^2;
    Q(n_lib,1,2) = trapz(var.x,trapz(var.x,trapz(var.t,th1,3),2),1);

    % u
    th2 = U.*dA000;
    Q(n_lib,2,2) = trapz(var.x,trapz(var.x,trapz(var.t,th2,3),2),1);

    % u^2
    th3 = U.^2.*dA000;
    Q(n_lib,3,2) = trapz(var.x,trapz(var.x,trapz(var.t,th3,3),2),1);

    % u^3
    th4 = U.^3.*dA000;
    Q(n_lib,4,2) = trapz(var.x,trapz(var.x,trapz(var.t,th4,3),2),1);

    % v
    th5 = V.*dA000;
    Q(n_lib,5,2) = trapz(var.x,trapz(var.x,trapz(var.t,th5,3),2),1);

    % v^2
    th6 = V.^2.*dA000;
    Q(n_lib,6,2) = trapz(var.x,trapz(var.x,trapz(var.t,th6,3),2),1);

    % v^3
    th7 = V.^3.*dA000;
    Q(n_lib,7,2) = trapz(var.x,trapz(var.x,trapz(var.t,th7,3),2),1);

    % uv
    th8 = U.*V.*dA000;
    Q(n_lib,8,2) = trapz(var.x,trapz(var.x,trapz(var.t,th8,3),2),1);

    % u^2v
    th9 = U.^2.*V.*dA000;
    Q(n_lib,9,2) = trapz(var.x,trapz(var.x,trapz(var.t,th9,3),2),1);

    % uv^2
    th10 = U.*V.^2.*dA000;
    Q(n_lib,10,2) = trapz(var.x,trapz(var.x,trapz(var.t,th10,3),2),1);
        
end % np
            
%% REGRESSION
% Parameters
ksi_u = SINDy(Q(:,:,1), q0(:,1));
ksi_v = SINDy(Q(:,:,2), q0(:,2));

fields = ["Laplacian";"u";"u2";"u3";"v";"v2";"v3";"uv";"u2v";"uv2"];

ksi.u = cell2struct(num2cell(ksi_u),fields);
ksi.v = cell2struct(num2cell(ksi_v),fields);

% How good was the estimation
res_u = mean(abs(q0(:,1) - Q(:,:,1)*ksi_u));
res_v = mean(abs(q0(:,2) - Q(:,:,2)*ksi_v));

res = [res_u, res_v];

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
lambda = gamma*mean(abs(dXdt)); 

for i = 1:5

  product = zeros(size(Xi)); 
  [~,w] = size(Theta);
  for p_ind = 1:w
    product(p_ind) = mean(abs(Xi(p_ind)*Theta(:,p_ind)));
  end

  smallinds = product < lambda;
  Xi(smallinds) = 0;                % set negligible terms to 0
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

wx = weight_poly(var.x,2,k(1));
wy = weight_poly(var.x,2,k(2));
wt = weight_poly(var.t,1,k(3));

[wX,wY,wT] = meshgrid(wx,wy,wt);

W = wX.*wY.*wT;

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
