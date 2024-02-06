%% Set constants
% Physical consts
c = 3 * 10^10;                                      % Speed of light [cm/s]
eVtoErg = 1.6e-12;                                  % Convertion coef from [eV] to [Erg]
m = 938.27 * 10^6;                                  % Mass of proton [eV]

% Spacial grid
a = 50;                     % Plasma radius in [cm]
b = 250;                    % Expander radius in [cm]
r = linspace(0,a,35);       % Grid inside the plasma [cm]
R = linspace(a,b,35);       % Grin outside the plasma [cm]

% Plasma parameters
mp = m;                     % Ions mass [eV]
np = 3e13;                  % Ions density [cm^{-3}]
Tp = 100;                   % Ions temperature [eV]

% Gas parameters
mg = 2 * m;                 % Gas mass [eV]

% Diff cross section of elastic scatering
diffCross = 1e-16;          % Gas-Ions elastic cross section [cm^2]             

% Wall conditions
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Twall = 0.026;
Pwall = 4.14;

% Velocity grid size
Nv = 300;

%% Section 1 Title
Trange = [0.6 0.61 0.62 0.63 0.64 0.65 0.66 0.67];

F = ComputeT0DensePlasma(1, Trange, Twall, Pwall, kappa, np, Tp, mp, mg, diffCross, a, b, Nv);

% The T0 which was fuond!!!
T0 = 0.637;

%% Section 2 Title
% Description of second code block

beta23 = (2^(1/3) + 2^(-2/3));
beta = beta23^(3/2);
localCoef = 5^(2/3);
nEff = 0.5 * 1 / (5 * beta) * (3/4) * ( sqrt(pi) * erf(localCoef * beta23) - ...
    2 * localCoef * beta23 * exp(-localCoef * beta23) );
[DFhot, gridStepHot, VrHot, VsqrHot] = DFDensePlasma(np, Tp, mp, nEff, T0, mg, diffCross, a, Nv, false);

PrrHot = 0;
for i = 1:Nv
    PrrHot = PrrHot + sum(DFhot(:,fix(Nv/2):Nv,i,end) .* VrHot(:,fix(Nv/2):Nv).^2, "all");
end
PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;

Prr = PrrHot;

n0 = ComputeN0(Prr, Pwall, 1);

% n0 = 1.4479e+12!!
