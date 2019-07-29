% kursiv.m - solution of Kuramoto-Sivashinsky equation by ETDRK4 scheme
%
%   u_t = -u*u_x - u_xx - u_xxxx, periodic BCs on [0,32*pi]
%   computation is based on v = fft(u), so linear term is diagonal
%   compare p27.m in Trefethen, "Spectral Methods in MATLAB", SIAM 2000
%   AK Kassam and LN Trefethen, July 2002

function [uu,tt,R,w] = KS_integrate(u0,tmax,par,if_plot,if_save,out_str)
% Spatial grid and initial condition:
Nx = 128;                         % no. of grid points in space
Lx = 22;                         % physical length of domain
x = Lx*(1:Nx)'/Nx;
dx = x(2) - x(1);
if isempty(u0)
u0 = cos(2*pi*x/Lx).*(1+sin(2*pi*x/Lx));       % initial condition
end
u = u0;
v = fft(u);

% Precompute various ETDRK4 scalar quantities:
h = 2/10;                           % time-step
k = 2*pi*[0:Nx/2-1 0 -Nx/2+1:-1]'/Lx;    % wave numbers
L = par(1)*k.^2 - par(2)*k.^4;      % Fourier multipliers
E = exp(h*L); E2 = exp(h*L/2);
M = 16;                             % no. of points for complex means
r = exp(1i*pi*((1:M)-.5)/M);        % roots of unity
LR = h*L(:,ones(M,1)) + r(ones(Nx,1),:);
Q  = h*real(mean(           (exp(LR/2)-1)./LR              ,2));
f1 = h*real(mean(  (-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3  ,2));
f2 = h*real(mean(    (2+LR+exp(LR).*(-2+LR))./LR.^3    ,2));
f3 = h*real(mean(  (-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3  ,2));

% Main time-stepping loop:
uu = u; tt = 0;
% tmax = 100;                         % time-length
Nt = 4*tmax;                         % no. of grid points in time
nmax = round(tmax/h); nplt = floor((tmax/Nt)/h);
g = -0.5i*k;
for n = 1:nmax
    t = n*h;
    Nv = g.*fft(real(ifft(v)).^2);
    a = E2.*v + Q.*Nv;
    Na = g.*fft(real(ifft(a)).^2);
    b = E2.*v + Q.*Na;
    Nb = g.*fft(real(ifft(b)).^2);
    c = E2.*a + Q.*(2*Nb-Nv);
    Nc = g.*fft(real(ifft(c)).^2);
    v = E.*v + Nv.*f1 + 2*(Na+Nb).*f2 + Nc.*f3;
    u = real(ifft(v));
    
    if mod(n,nplt)==0 && out_str == "full"        
        uu = [uu,u];
        tt = [tt,t];
    end
end

if out_str == "final"
    uu = u;
end



% Plot results:
if(if_plot && out_str == "full")
%     surf(tt,x,uu), shading interp, lighting phong, axis tight
%     view([-90 90]), colormap(autumn); set(gca,'zlim',[-5 50])
%     light('color',[1 1 0],'position',[-1,2,2])
%     material([0.30 0.60 0.60 40.00 1.00]);
    figure
    subplot(2,1,1)
    imagesc(flipud(uu))
    xlabel('t'),ylabel('x')
    subplot(2,1,2)
    dt = tt(2) - tt(1);
    Tau = 200;
    [R,w] = recurrence(uu,dt,Tau);
    imagesc(1:tmax,1:Tau,flipud(R))
    xlabel('t'),ylabel('\tau')
    
end

% Save results
if(if_save && out_str == "full")
    dt = tt(2)-tt(1);
    [Nx,Nt] = size(uu);
    save_name = ['KS_',num2str(par(1)),'_',num2str(par(2)),'.mat'];
    save(save_name,'tt','uu','par','Nx','Nt','dx','dt','-v7.3')
    
end

end
function [R,w] = recurrence(U_t,dt,Tau)
dtau = 1/dt;
[Lx,Lt] = size(U_t);
lt = length(1:dtau:(dtau*Tau));
R = zeros(lt,Lt);
w = zeros(lt,Lt);
% s1 = figure;
for t = 1:Lt
    n_tau = 0;
    for tau = 1:dtau:min(dtau*Tau,t-1)
        n_tau = n_tau + 1;
        V = fftshift(fft(U_t(:,t-tau)));
        n = [0 -Lx/2+1:Lx/2-1]';
        for l = 0:20:Lx-1
        th = (l/Lx)*(2*pi);        
        D = exp(-1i*n*th);
        U_tau = real(ifft(ifftshift(D.*V)));
        R_int = norm(U_t(:,t) - U_tau);
           
        if l == 0
            R(n_tau,t) = R_int;
        else
            R(n_tau,t) = min(R_int,R(n_tau,t));
            if R(n_tau,t) == R_int
                w(n_tau,t) = l;
            end
        end
        end
    end
end
end
