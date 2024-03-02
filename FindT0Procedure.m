%% Find T0 - BC for kinetic issue
% This procedure is OK for thin plasma limit!!!

%% Set constants
% Physical constants
c = 3 * 10^10;              % Light speed in [cm/s]
eVtoErg = 1.6e-12;          % Convertion coeff from [eV] to [Erg]
m = 938.27 * 10^6;          % Mass of proton [eV]

% Spacial grid
a = 50;                     % Plasma radius in [cm]
b = 250;                    % Expander radius in [cm]
r = linspace(0,a,35);       % Grid inside the plasma [cm]
R = linspace(a,b,35);       % Grin outside the plasma [cm]

% Plasma parameters
mp = m;                     % Ions mass [eV]
np = ones(size(r)) * 3e11;  % Ions density [cm^{-3}]
Tp = ones(size(r)) * 100;   % Ions temperature [eV]

% Gas parameters
mg = 2 * m;                 % Gas mass [eV]

% Diff cross section of elastic scatering
diffCross = (1 / (4 * pi)) * 3.6e-15;          % Gas-Ions elastic cross section [cm^2]             

% Wall conditions
Twall = 0.026;
nwall = 1e14;
Pwall = nwall * (Twall * eVtoErg);

% Velocity grid parameters
Nv = 300;

%% Compute T0
% Rough search on the graph. This part does't optimize well
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Trange = [1.1 1.15 1.2 1.25];

% Need to find zero of function F with some accuracy!!!
F = ComputeT0(1, Trange, Twall, Pwall, kappa, np, Tp, mp, mg, diffCross, r, a, b, Nv);

% The T0 which was fuond!!!
T0 = 1.125;

%% Compute n0
% After we find T0. We can get n0
% Compute DF in "cold" and "hot" region for given T0
[DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, 1, T0, mg, diffCross, r, a, Nv, false);
nCold = Ncold(np, Tp, mp, 1, T0, mg, diffCross, r, a, Nv);
[DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, Nv, false);

positiveIdxStart = fix(Nv/2) + 1;

% Compute Prr
PrrHot = 0;
for i = 1:Nv
    PrrHot = PrrHot + sum(DFhot(:,positiveIdxStart:Nv,i,end) .* VrHot(:,positiveIdxStart:Nv).^2, "all");
end
PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;

PrrCold = 0;
for i = 1:Nv
    PrrCold = PrrCold + sum(DFCold(:,positiveIdxStart:Nv,i,end) .* VrCold(:,positiveIdxStart:Nv).^2, "all");
end
PrrCold = PrrCold * mg * eVtoErg / c^2 * gridStepCold^3;

Prr = PrrCold + PrrHot;

% The n0 which was fuond!!!
n0 = ComputeN0(Prr, Pwall, 1);

%% Reference answer to check
% with this params:
%
% Pwall = 4.14
% diffCross = 1e-16
% np = 1e11
% Tp = 100
% mp - proton mass
% mg = 2mp
% a = 50
% b = 250
%
% Answer:
%
% T0 = 0.466        !!!
% n0 = 4.1058e+12   !!!

% H colisions
% T0 = 1.86
% n0 = 7.587e+11

% H2 colisions
% T0 = 1.125
% n0 = 1.0976e+12
