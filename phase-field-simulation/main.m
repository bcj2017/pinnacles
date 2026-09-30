%--------------------------------------------------------------------------
% Phase field model for simulating ice melting in warm water with the 
% Boussinesq approximation.
% --------------------------------------------------------------------------

close all
addpath('functions')
target = 'frames'; %output directory
if exist(target,'dir'); rmdir(target,'s'); end
mkdir(target)

% Pinnacle sizes: 5, 10, 30, 50, 100, 170, 300, 500

global hz rc zc Nr Nz hp hm hc nr %#ok<*GVMIS> 
% Things I often change
m = .05; %interface mobility
pinnacleSize = 100; % Quantities we run on: 5 10 30 50 100 170 300 500
% Perturbations: blunt, sharp, unperturbed
perturbationVal = 0; % 0 = unperturbed, -1 = blunt, 1 = sharp
perturbationHeightCm = .25;

T_inf = 20; %C, far-field temperature of 20 degrees Celsius
dt = 1e-8; %time step %CHANGE: originally 5e-7
tplt = 20*100*dt; %plot time %CHANGE: originally 100*dt; should be same evenly spaced time steps now
tsave = 20*500*dt; %save time %CHANGE: originally 500*dt; should be same evenly spaced time steps now

% Dependent on pinnacleSize: Nz (# grid points), lr (dimensionless refined region),
% H (height in dimensionless domain)
% lr chosen so that about 1mm of leeway between ice and coarse region
H = .8; H0 = 8 * 1e-2;
if pinnacleSize == 5
    Nz = 1900; lr = 0.1; H = 0.8; % 50
elseif pinnacleSize == 10
    Nz = 1900; lr = 0.11; H = 0.8; % 50
elseif pinnacleSize == 30
    Nz = 1800; lr = 0.15; H = 0.8; % 55
elseif pinnacleSize == 50
    Nz = 1800; lr = 0.17; H = 0.8; % 55
elseif pinnacleSize == 100
    Nz = 1600; lr = 0.2; H = 0.8; % 66 refined micron spacing
elseif pinnacleSize == 170
    % Nz = 1500; lr = 0.2; 
    Nz = 1500; lr = .23;
elseif pinnacleSize == 300
    % Nz = 1500; lr = 0.2; H = 0.6; % 89
    Nz = 1500; lr = .26;
elseif pinnacleSize == 500
    % Nz = 1500; lr = 0.206; H = 0.55; % 100
    Nz = 1500; lr = .3;
end
% H0 = 8 * 1e-2; % initial pinnacle height in cm
L = H0 / H;

%% Constants contributing to dimensionless constants (Bobae) from 02/27/2024
g = 9.8; %m/s^2, gravity
beta = 7.68 * 1e-6; %1/C^2, thermal expansion coefficient
T_0 = 0; %C, melting temperature
Tstar = 4;
nu = 1.57*1e-6; %m^2/s at 4C, kinematic viscosity
kappa_T = 1.3*1e-7; %m^2/s at 4C, thermal diffusivity
c_p = 4.205 * 1e3; %at 4C, J / kg K, specific heat capacity
latent_heat = 3.34*1e5; %J / kg, latent heat of fusion, engineering toolbox


%% System parameters
Lr = 2; %radial box size
Lz = 1; %vertical box size
Ra = g * beta * (T_inf - T_0)^2 * L^3 / (nu * kappa_T);
Pr = nu / kappa_T;
St = c_p * (T_inf - T_0) / latent_heat;
disp(['Rayleigh number: ',num2str(Ra)]);
disp(['Prandtl number: ',num2str(Pr)]);
disp(['Stefan number: ',num2str(St)]);

t0 = 0; %initial time
tf = 1e-1; %maximum time
Nt = floor((tf-t0)/dt); %number of time steps
nplt = floor(tplt/dt); %plot interval 
nsave = floor(tsave/dt); %save interval
nt = 0; nf = 0; ns = 0; nb = 0; %counters

Vf = 0.1; %final volume

%% Discretization

% Spatial discretization
Nr = floor(Lr/Lz)*Nz; % " " in r
nr = floor((lr/Lr)*Nr); %number of grid points in refined region %I THINK we want to change here... and then edit everywhere else
r1 = linspace(0,Lr,Nr+1)'; h1 = r1(2)-r1(1); %r grid
r2 = linspace(r1(nr),Lr,floor((Lr-r1(nr))/(20*h1))+1)'; h2 = r2(2) - r2(1); %coarse grid %CHANGED
r = [r1(1:nr);r2(2:end)]; rc = avg(r); %radial grid
z = linspace(0,Lz,Nz+1)'; zc = avg(z); %vertical grid
Nr = length(rc); %updated number of radial cells

% Grid spacing
h = diff(rc); %radial spacing
hp = [h;h(end)]; %extended right
hm = [h(1);h]; %extended left
hc = avg(hp); %centered
hz = z(2)-z(1); %vertical spacing

[rr,zz] = meshgrid(rc,zc); rr = rr'; zz = zz';
rrr = [-flipud(rr);rr]; zzz = [zz;zz]; %mirror grid
      
% Phase-field parameters (calibrated for stability)
del = hz/2; %interface thickness
a = (1/St)/del; %coupling coefficient
eta = dt; %velocity forcing factor
Cdp = a*St/del; %dp coefficient in phase eq, 1/delta^2
Cdg = 0.25/del^2; %dg coefficient in phase eq, 1/4delta^2
Cphi = (1/St)/(2*dt); %dphidt coefficient in temp eq
Tw = (T_0-Tstar)/(T_inf-T_0); %wall temperature. CHANGE 2026-01-19. 

%% Initialize fields and operators
% Time stepping operators and LU decomposition
[Lu,Lv,Ltm,Ltp,Lphim,Lphip,Lumh,Luph,Lvmh,Lvph,Lp,Lq] = calcops(Nr,Nz,dt,Pr,m);
LU = lustruct(Lumh); LV = lustruct(Lvmh); LP = lustruct(Lp);
LT = lustruct(Ltm); LPhi = lustruct(Lphim);

[u,v,T,phi] = initialize(rr,zz,Nr,Nz,Tstar,T_inf,H,L,pinnacleSize,perturbationVal,perturbationHeightCm); %initial conditions
fum1 = zeros(Nr-1,Nz); fvm1 = zeros(Nr,Nz-1);
p = zeros(Nr,Nz); phim1 = phi;

% Refined region of simulation (ideally slightly larger than pinnacle
  % x-domain)
  refinedRegioncm = lr*L*1e2;
  disp(['Refined region: ',num2str(refinedRegioncm),' cm.']);
  % Micron spacing of refined region
  refinedSpacing = (rc(2)-rc(1))*L*1e2*1e4;
  disp(['r Refined spacing: ',num2str(refinedSpacing),' microns = ',num2str(refinedSpacing*1e-4),' cm.'])
  % Micron spacing of unrefined region
  unrefinedSpacing = (rc(end)-rc(end-1))*L*1e2*1e4;
  disp(['r Unrefined spacing: ',num2str(unrefinedSpacing),' microns = ',num2str(unrefinedSpacing*1e-4),' cm.'])
  % Micron spacing in z direction
  zSpacing = (zc(2)-zc(1))*L*1e2*1e4;
  disp(['z spacing: ',num2str(zSpacing),' microns = ',num2str(zSpacing*1e-4),' cm.'])
  % Total number of data points
  disp(['Total number of data points = ',num2str(Nr),'*',num2str(Nz),' = ',num2str(Nr*Nz)])

V0 = 2*pi*sum(sum((1 - phi).*rr.*hp*hz)); %initial volume
Vsave = [0,V0]; %volume storage

Xsave = {}; %boundary storage

% Velocity boundary conditions (no-slip)
uN = zeros(Nr-1,1); vN = zeros(Nr,1);
uS = zeros(Nr-1,1); vS = zeros(Nr,1);
uW = zeros(1,Nz);   vW = zeros(1,Nz-1);
uE = zeros(1,Nz);   vE = zeros(1,Nz-1);

%% Begin main loop
while nt <= Nt
  
  tic
  t = dt*nt;

  % Phase field predictor
  fphi = m*(Cdp*(T-Tw).*dp(phi) - Cdg*dg(phi));%nonlinear part %2026 CHANGE: Tw is now < 0
  rhs = mult(Lphip,phi) + dt*fphi; %rhs function
  phip1 = lusolve(LPhi,rhs); %predictor step
  phim = 0.5*(phi + phip1); %mid-point phase-field
  pen = (1/eta)*(1 - phim).^2; %penalization factor

  % Temperature predictor
  [Tr,~] = grad([T(1,:);T;T(end,:)]);
  [~,Tz] = grad([T(:,1),T,T(:,end)]);
  uTr = avg([uW;u;uE].*Tr);
  vTz = avg(([vS v vN].*Tz)')';
  fT = -(uTr + vTz + Cphi*(phip1 - phim1).*dp(phi)); %nonlinear part
  rhs = mult(Ltp,T) + dt*fT; %rhs function
  Tp1 = lusolve(LT,rhs); %predictor step
  Tm = 0.5*(T + Tp1); %mid-point temperature

  Fb = Pr*Ra*(avg(Tm')').^2; %buoyancy force % CHANGE 2025-10-06: quadratic
           
  % Half-step velocity  
  [pr,pz] = grad(p);
      %x component: 5/8 stuff due to Adams-Bashforth implementation
      %probably
  [ur,~] = grad([uW;u;uE]);
  [~,uz] = grad([2*uS-u(:,1) u 2*uN-u(:,end)]);
  uur = u.*avg(ur); % u \dot \nabla u
  vuz = avg((avg([vS v vN]).*uz)')'; % same for z component
  fu = -(uur + vuz) + Pr*(-pr + mult(Lu,u)) - avg(pen).*u; %nonlinear part
  uph = u + (dt/8)*(5*fu - fum1); fum1 = fu; %half step
      %z component:
  [vr,~] = grad([v(1,:);v;2*vE-v(end,:)]);
  [~,vz] = grad([vS v vN]);
  uvr = avg(avg([uW;u;uE]')'.*vr);
  vvz = v.*avg(vz')';
  fv = -(uvr + vvz) + Pr*(-pz + mult(Lv,v)) - avg(pen')'.*v + Fb; %nonlinear part
  vph = v + (dt/8)*(5*fv - fvm1); fvm1 = fv; %half step

  % Full step velocity    (including pressure update through pressure correction section)
  % x component. uph = u^* half step velocity
  [ur,~] = grad([uW;uph;uE]);
  [~,uz] = grad([2*uS-uph(:,1) uph 2*uN-uph(:,end)]);
  uur = uph.*avg(ur);
  vuz = avg((avg([vS vph vN]).*uz)')';
  fu = -(uur + vuz) - Pr*pr - avg(pen).*uph; %nonlinear part
  rhs = mult(Luph,u) + dt*fu; %rhs function
  up1 = lusolve(LU,rhs); %full step
  % z component
  [vr,~] = grad([vph(1,:);vph;2*vE-vph(end,:)]);
  [~,vz] = grad([vS vph vN]); 
  uvr = avg(avg([uW;uph;uE]')'.*vr);
  vvz = vph.*avg(vz')';
  fv = -(uvr + vvz) - Pr*pz - avg(pen')'.*vph + Fb; %nonlinear part
  rhs = mult(Lvph,v) + dt*fv; %rhs function
  vp1 = lusolve(LV,rhs); %full step

  % Projection step
  ue = [uW;up1;uE];
  ve = [vS vp1 vN];
  ur = grad(ue) + avg(ue)./rr;
  ur(1,:) = 1/(2*hc(1))*(3*ue(2,:) - ue(1,:));
  [~,vz] = grad(ve);

  rhs = -(ur + vz)/(Pr*dt); %divergence
  q = lusolve(LP,rhs); %pressure corrector
  [qr,qz] = grad(q); %pressure gradient
  
  up1 = up1 - (Pr*dt)*qr; %project u
  vp1 = vp1 - (Pr*dt)*qz; %project v
  p = p + (q - 0.5*dt*mult(Lq,q)); %update pressure

  % Phase field corrector
  fphi = m*(Cdp*(Tm-Tw).*dp(phim) - Cdg*dg(phim)); %nonlinear part %2026 CHANGE: theta_0 now relevant, added -Tw term
  rhs = mult(Lphip,phi) + dt*fphi; %rhs function
  phip1 = lusolve(LPhi,rhs); %full step

  % Temperature corrector
  [Tr,~] = grad([Tm(1,:);Tm;Tm(end,:)]);
  [~,Tz] = grad([Tm(:,1),Tm,Tm(:,end)]);
  uTr = avg([uW;0.5*(u+up1);uE].*Tr);
  vTz = avg(([vS 0.5*(v+vp1) vN].*Tz)')';
  fT = -(uTr + vTz + Cphi*(phip1 - phim1).*dp(phim)); %nonlinear part
  rhs = mult(Ltp,T) + dt*fT; %rhs function
  Tp1 = lusolve(LT,rhs); %full step

  % Prepare for next time step
  u = up1; v = vp1; T = Tp1; phim1 = phi; phi = phip1;
  V = 2*pi*sum(sum((1 - phi).*rr.*hp*hz)); %volume

  el = toc;
  clc; fprintf(' t = %1.4e\n V/V0 = %1.4f\n loop time = %1.4fs\n',t,V/V0,el);
  if isnan(T(1)) || isinf(T(1)); fprintf('\n Unstable...\n'); break; end
  if V/V0 < Vf; break; end

  % Save
    if mod(nt,nsave) == 0
        filename = sprintf('%s/s-%d',target,ns); ns = ns+1;
        save(filename,'u','v','T','phi');

        save('output.mat'); % save output file to reinitialize in new simulations
        save('timeData.mat','dt','tplt','tsave','Nz','m','kappa_T','a','nplt','nsave','rc','zc','r1','nr','hz','T_0','T_inf','del','St');
    end

  % Plot
  if mod(nt,nplt) == 0 % NO PNG SAVE
    pphi = [flipud(phi);phi]; %mirror field
    X = contour(rrr,zzz,pphi,1,'LineWidth',2,'Color',0.5*[1 1 1]); %boundary
    Xsave{nb+1} = X(:,2:end); nb = nb+1; %#ok<SAGROW>
    if nt > dt % saves after the first one
        save('alldata.mat', 'Xsave'); % less frequent, but useful for reinitializing to know the exact stopping point!
    end
  end



 nt = nt + 1;
    
end
