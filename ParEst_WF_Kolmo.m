%%%% Integral Method Q2D-NS Parameter Estimation %%%%
%% ----------------------------------------------------------------------
%                           PREAMBLE
% -----------------------------------------------------------------------
%{
INPUTS 
------
filename        : string for path to -v7.3 .mat file containing data
    - U_t,V_t   : non-dim velocity fields
    - dx        : non-dim spatial grid spacing
    - dt_mks    : dimensional temporal grid spacing (seconds)
    - Ts_mks    : time-scale in seconds
N_d             : number of integration domains
N_h = [Q R S] 
    - S         : max # harmonics in time (>=1)
    - Q,R       : max # harmonics in x,y resp. (>=0)
D = [fr, Dt] 
    - fr        : spatial-area of int domain as fraction of full area
    - Dt        : time-length of integration domain in grid points
track           : boolean to track library filling in cmd window

OUTPUTS
-------
ksi : estimated parameters in non-dimensinoal units
res : residual of estimation (utility unclear)
Q   : matrix containing integral evaluations of gov eq terms
q0  : vector containing integral evaluation of time-derivative term
P   : 3xN_d matrix containing integration domain locations 

SECTIONS
--------
Initialize
    -size of local domain
    -how many harmonics to sample
    -any tracking variables
    -library and target

Construct Library
    -advection term
    -laplacian term
    -rayleigh term
    -time-derivative term

Regression
    -invert library onto target or use SINDy
    -calculate norm of b-\Theta\xi

Reminders
    -disp to-do list

FUNCTIONS
---------
Construct 3D Weight Function
    - use General Liebniz Rule to find the k^th oder product rule of base
        weight and harmonics
    - use meshgrid to combine 3 1D composite fields into final 3D weight
        function

Construct 1D Base Weight Function
    - find base poly coefficients based on exponent m
    - convert poly coefficients based on derivative order k
    - compute full polynomial

Construct 1D Weight harmonics
    - Calculate derivatives of sin(pi*q*x) in 1D

%}
%%
function [ksi,res,Q,q0,P] = ParEst_WF_Kolmo(filename,N_d,N_h,D,track)
%% ----------------------------------------------------------------------
%                           INITIALIZE
% -----------------------------------------------------------------------
% Define matfile
traj = matfile(filename);

% Load time-scale
Ts_mks = traj.Ts_mks;

% Get size of velocity fields
[Ly,Lx,Lt] = size(traj,'U_t');

% Grid densities
dt_mks = traj.dt_mks;
dx = traj.dx; 

% Size of local domain
fr = D(1);
clearvars -global
global var
var.Dx = round(fr*Lx); % sizes of the integration domain
var.Dy = round(fr*Lx);
var.Dt = D(2);

% Sampling scheme
P = zeros(3,N_d); % Starting corner of integration domain
P(1,:) = randi([1,Lx-var.Dx],N_d,1);
P(2,:) = randi([1,Ly-var.Dy],N_d,1);
P(3,:) = randi([1,Lt-var.Dt],N_d,1);

% Define derivative conversions
S_x = 2/(dx*var.Dx);
S_y = 2/(dx*var.Dy);
S_t = (2*Ts_mks)/(dt_mks*var.Dt);

% Create subdomain
var.x = linspace(-1,1,var.Dx); 
var.y = linspace(-1,1,var.Dy);
var.t = linspace(-1,1,var.Dt);

% Initialize Target and Library
q0 = zeros(N_d*N_h(1)*N_h(2)*N_h(3),1);
q1 = q0; q2 = q0; q3 = q0;

%% ----------------------------------------------------------------------
%                           CONSTRUCT LIBRARY
% -----------------------------------------------------------------------
n_lib = 0;
n_track = 10;

for np = 1:N_d  

    % Indices for integration domain
    rx = P(1,np):(P(1,np)+var.Dx-1);
    ry = P(2,np):(P(2,np)+var.Dy-1);
    rt = P(3,np):(P(3,np)+var.Dt-1);

    % Velocity fields on integration domain
    U = traj.U_t(ry,rx,rt);
    V = traj.V_t(ry,rx,rt);

    for q = 0:N_h(1)
    for r = 0:N_h(2)
    for s = 1:N_h(3)
                
            n_lib = n_lib + 1;

            if track && n_lib == n_track
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
            
            % Make wave numbers global for use in weight_full()
            var.q = q;
            var.r = r;
            var.s = s;

            % Make derivatives of windowing functions
            dA100 = weight_full([1,0,0]); dA010 = weight_full([0,1,0]);           
            dA101 = weight_full([1,0,1]); dA011 = weight_full([0,1,1]);            
            dA110 = weight_full([1,1,0]);
            dA020 = weight_full([0,2,0]); dA200 = weight_full([2,0,0]);                         
            dA210 = weight_full([2,1,0]); dA120 = weight_full([1,2,0]);              
            dA030 = weight_full([0,3,0]); dA300 = weight_full([3,0,0]);                                   
            % --> the 3 digits correspond to derivatives, NOT to harmonics

            % Time-derivative
            b = (V.*dA101*S_x - U.*dA011*S_y)*S_t;
            q0(n_lib,1) = trapz(var.x,trapz(var.y,trapz(var.t,b,3)));

            % Advection Term (incompressible)
            th1 = U.*V.*(dA020*S_y^2 - dA200*S_x^2) + ...
                  (U.^2 - V.^2).*dA110*S_x*S_y;
            q1(n_lib,1) = trapz(var.x,trapz(var.y,trapz(var.t,th1,3)));   

            % Laplacian Term
            th2 = U.*(dA210*S_x^2*S_y + dA030*S_y^3) - ...
                V.*(dA300*S_x^3 + dA120*S_x*S_y^2);
            q2(n_lib,1) = trapz(var.x,trapz(var.y,trapz(var.t,th2,3)));

            % Rayleigh Term
            th3 = V.*dA100*S_x - U.*dA010*S_y;
            q3(n_lib,1) = trapz(var.x,trapz(var.y,trapz(var.t,th3,3)));           

    end % s
    end % r
    end % q
end % np

% Fill Library
Q = [q1 q2 q3];

%% ----------------------------------------------------------------------
%                               REGRESSION
% -----------------------------------------------------------------------
% Compute parameters
ksi = Q \ q0;

% Compute residual
res = mean(abs(q0 - Q*ksi));

%% ----------------------------------------------------------------------
%                               REMINDERS
% -----------------------------------------------------------------------
disp(' ')
disp(' ')

end

%% -----------------------------------------------------------------------
%                       FULL 3D WINDOWING FUNCTION
% ------------------------------------------------------------------------
function W = weight_full(k)
global var
%{
Combine base windows and harmonics in 1D, then assemble into 3D
k = [kx,ky,kt]: order of derivative(s)
%}

% General Liebniz Rule 
wx = 0;
for n = 0:k(1)
    p = weight_poly(var.x,3,k(1)-n); 
    g = weight_harm(var.x,n,var.q);
    wx = wx + nchoosek(k(1),n)*p.*g;  
end
wy = 0;
for n = 0:k(2)
    p = weight_poly(var.y,3,k(2)-n);
    g = weight_harm(var.y,n,var.r);
    wy = wy + nchoosek(k(2),n)*p.*g;
end
wt = 0;
for n = 0:k(3)
    p = weight_poly(var.t,0,k(3)-n);
    g = weight_harm(var.t,n,var.s);
    wt = wt + nchoosek(k(3),n)*p.*g;
end

[wX,wY,wT] = meshgrid(wx,wy,wt);        % Make 1D components 3D

W = wX.*wY.*wT;                         % combine all components

end
%% -----------------------------------------------------------------------
%                    1D BASE WEIGHT FUNCTION
% ------------------------------------------------------------------------
function p = weight_poly(x,m,k)
%{
Polynomial piece of weighting function used to satisfy BC

A = d^k/dx^k[ (x^2 - 1)^m ]

x: independent variable
m: power of base function
k: order of derivative
%}
a = zeros(m*2 + 1,1);                       
for l = 0:m
    a(2*l+1) = (-1)^(m-l)*nchoosek(m,l);    % set polynomial coefficients
end 
c = zeros(2*m+1,1);                         
for n = 0:(2*m - k)                         % take k^th derivative
    c(n+1) = a(n+1+k)*factorial(n+k)/factorial(n);
end
p = 0;
for n = 0:(2*m-k)
    p = p + c(n+1)*x.^n;                    % final windowing function
end
end

%% -----------------------------------------------------------------------
%                     WEIGHT FUNCTION HARMONICS
% ------------------------------------------------------------------------
function g = weight_harm(x,k,q)
% g_j = sin(q*pi*x) 
if mod(k,2)                                 % odd k
    n = (1/2)*(k - 1);
    g = (-1)^n*(pi*q)^k*cos(pi*q*x);
else                                        % even k
    n = k/2; 
    g = (-1)^n*(pi*q)^k*sin(pi*q*x);
end
if q == 0 && k == 0
    g = ones(size(x));
elseif q == 0 && k > 0
    g = zeros(size(x));
end
end
