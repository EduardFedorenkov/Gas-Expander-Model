%% Set constants
% Physical consts
c = 3 * 10^10;                                      % Speed of light [cm/s]
eVtoErg = 1.6e-12;                                  % Convertion coef from [eV] to [Erg]
m = 938.27 * 10^6;                                  % Mass of proton [eV]

% Spacial grid
a = 50;                     % Plasma radius in [cm]
b = 250;                    % Expander radius in [cm]
r = linspace(0,a,100);       % Grid inside the plasma [cm]
R = linspace(a,b,100);       % Grin outside the plasma [cm]

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

%% Compute DF in "cold" and "hot" region

nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv);

beta23 = (2^(1/3) + 2^(-2/3));
beta = beta23^(3/2);
localCoef = 5^(2/3);
nEff = 0.5 * n0 / (5 * beta) * (3/4) * ( sqrt(pi) * erf(localCoef * beta23) - ...
    2 * localCoef * beta23 * exp(-localCoef * beta23) );
[DFhot, gridStepHot, VrHot, VsqrHot] = DFDensePlasma(np, Tp, mp, nEff, T0, mg, diffCross, a, Nv, true);

%% Compute n_hot and n_cold
% Description of second code block

figure(3);
x = r/b;
plot(x,nCold,'LineWidth',2,'Color','Blue');
hold on

%% Compute n ouside the plasma

positiveIdxStart = fix(Nv/2) + 1;

PrrHot = 0;
qrHot = 0;
for i = 1:Nv
    PrrHot = PrrHot + sum(DFhot(:,positiveIdxStart:Nv,i,end) .* VrHot(:,positiveIdxStart:Nv).^2, "all");
    qrHot = qrHot + sum(DFhot(:,positiveIdxStart:Nv,i,end) .* VrHot(:,positiveIdxStart:Nv) .* VsqrHot(:,positiveIdxStart:Nv,i), "all");
end
PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;
qrHot = qrHot * mg * eVtoErg / c^2 * 0.5 * gridStepHot^3;

Prr = PrrHot;
qrr = qrHot;

NOut = Prr ./((Twall * eVtoErg)^(3/2) - qrr * a / kappa .* log(R./b)).^(2/3);

xout = R/b;
figure(3);
plot(xout,NOut,'LineWidth',2,'Color','Blue');
title('H^+, H_2: n_p = 3 * 10^{13} см^{-3}, T_p = 100 эВ, n_0 = 1.8 * 10^{12} см^{-3}, T_0 = 1.42 эВ');
xlabel('r/b') ;
ylabel('n [cm^{-3}]');


%% Compute electrons flux

VT0 = sqrt(2 * T0 / mg) * c;
j_in = n0 * VT0 / (sqrt(pi) * 2);
j_out = 0;

for i = 1:Nv
    j_out = j_out + sum(DFhot(:,positiveIdxStart:Nv,i,end) .* VrHot(:,positiveIdxStart:Nv), "all");
end
j_out = j_out * gridStepHot^3;

dif = (j_in - j_out)/j_in;

% H2 - 50.3% ионизуется
% H2 - n0 = 1.7802e+12;
% H2 - T0 = 1.42;

% H - 50.4% ионизуется
