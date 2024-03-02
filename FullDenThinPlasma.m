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
np = ones(size(r)) * 3e11;  % Ions density [cm^{-3}]
Tp = ones(size(r)) * 100;   % Ions temperature [eV]

% Gas parameters
mg = 2 * m;                 % Gas mass [eV]
n0 = 1.0976e+12;            % Gas dencity [cm^{-3}]
T0 = 1.125;                 % Gas temperature [eV]

% T0 = 1.86;
% n0 = 7.587e+11;

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

[DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv, false);
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv);
[DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, Nv, false);

%% Compute n_hot and n_cold
% Description of second code block

nHot = zeros(1,length(r));
for i = 1:length(r)
    nHot(i) = sum(DFhot(:,:,:,i), "all");
end
nHot = nHot * (gridStepHot)^3;
figure(3);
x = r/b;
semilogy(x,nCold,'LineWidth',2);
hold on
semilogy(x,nHot,'LineWidth',2);
semilogy(x,nCold + nHot,'LineWidth',2,'Color','Yellow');

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

PrrCold = 0;
qrCold = 0;
for i = 1:Nv
    PrrCold = PrrCold + sum(DFCold(:,positiveIdxStart:Nv,i,end) .* VrCold(:,positiveIdxStart:Nv).^2, "all");
    qrCold = qrCold + sum(DFCold(:,positiveIdxStart:Nv,i,end) .* VrCold(:,positiveIdxStart:Nv) .* VsqrCold(:,positiveIdxStart:Nv,i), "all");
end
PrrCold = PrrCold * mg * eVtoErg / c^2 * gridStepCold^3;
qrCold = qrCold * mg * eVtoErg / c^2 * 0.5 * gridStepCold^3;

Prr = PrrCold + PrrHot;
qrr = qrHot + qrCold;

NOut = Prr ./((Twall * eVtoErg)^(3/2) - qrr * a / kappa .* log(R./b)).^(2/3);

xout = R/b;
figure(3);
semilogy(xout, NOut,'LineWidth',2,'Color','Yellow');
legend('концентрация "холодного" H_2','концентрация "горячего" H_2','полная концентрация H_2');
title('H^+, H_2: n_p = 3 * 10^{11} см^{-3}, T_p = 100 эВ, n_0 = 1.1 * 10^{12} см^{-3}, T_0 = 1.125 эВ');
xlabel('r/b') ;
ylabel('n [cm^{-3}]');
lgd = legend;
lgd.FontSize = 14;