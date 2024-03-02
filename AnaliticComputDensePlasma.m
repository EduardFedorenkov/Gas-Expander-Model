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
% n0 = 1.2786e+12;            % Gas dencity [cm^{-3}]
% T0 = 2.36;                 % Gas temperature [eV]

n0 = 1.7802e+12;            % H2
T0 = 1.42;                  % H2

% Diff cross section of elastic scatering
diffCross = (1 / (4 * pi)) * 3.6e-15;          % Gas-Ions elastic cross section [cm^2]     

% Wall conditions
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Twall = 0.026;
nwall = 1e14;
Pwall = nwall * (Twall * eVtoErg);

% Velocity grid size
Nv = 300;

%% Section 1 Title
[da, IntOnePrr, IntTwoPrr] = AnaliticComputDensePlas(mg, mp, T0, np, Tp, diffCross, a, Nv);

beta23 = (2^(1/3) + 2^(-2/3));
beta = beta23^(3/2);
localCoef = 5^(2/3);
nEff = 0.5 * n0 / (5 * beta) * (3/4) * ( sqrt(pi) * erf(localCoef * beta23) - ...
    2 * localCoef * beta23 * exp(-localCoef * beta23) );

% IntOnePrr
% IntTwoPrr * (da / (2*a))

%% Section 2 Title

T0_expect = Twall * (1 + a*log(b/a)*sqrt(mg*mp)/(mp+mg)*sqrt(Tp/Twall) * nwall* diffCross * 2.6 * 9/5 * sqrt(2))^(2/3);

n_gap = T0_expect / Twall;
