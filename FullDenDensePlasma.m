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
np = 1e14;                  % Ions density [cm^{-3}]
Tp = 100;                   % Ions temperature [eV]

% Gas parameters
mg = 2 * m;                 % Gas mass [eV]
n0 = 1.4479e+12;            % Gas dencity [cm^{-3}]
T0 = 0.637;                 % Gas temperature [eV]

% Diff cross section of elastic scatering
diffCross = 1e-16;          % Gas-Ions elastic cross section [cm^2]             

% Wall conditions
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Twall = 0.026;
Pwall = 4.14;

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
semilogy(x,nCold,'LineWidth',2);
hold on
title('Gas dencity inside the plasma column');
xlabel('r/a') ;
ylabel('n [cm^{-3}]');
lgd = legend;
lgd.FontSize = 18;

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
semilogy(xout,NOut,'LineWidth',2);
legend('n_{out}');
lgd1 = legend;
lgd1.FontSize = 18;


%% Compute electrons flux

VT0 = sqrt(2 * T0 / mg) * c;
j_in = n0 * VT0 / (sqrt(pi) * 2);
j_out = 0;

for i = 1:Nv
    j_out = j_out + sum(DFhot(:,positiveIdxStart:Nv,i,end) .* VrHot(:,positiveIdxStart:Nv), "all");
end
j_out = j_out * gridStepHot^3;

dif = j_in - j_out;
