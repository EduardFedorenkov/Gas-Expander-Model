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
np = ones(size(r)) * 1e11;  % Ions density [cm^{-3}]
Tp = ones(size(r)) * 100;   % Ions temperature [eV]

% Gas parameters
mg = 2 * m;                 % Gas mass [eV]
n0 = 4.1058e+12;            % Gas dencity [cm^{-3}]
T0 = 0.466;                 % Gas temperature [eV]

% Diff cross section of elastic scatering
diffCross = 1e-16;          % Gas-Ions elastic cross section [cm^2]             

% Wall conditions
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Twall = 0.026;
Pwall = 4.14;

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
semilogy(x,nCold + nHot,'LineWidth',2);
legend('n_{cold}','n_{hot}','n_{full}');
title('Gas dencity inside the plasma column');
xlabel('r/a') ;
ylabel('n [cm^{-3}]');
lgd = legend;
lgd.FontSize = 18;

%% Compute n ouside the plasma

PrrHot = 0;
qrHot = 0;
for i = 1:Nv
    PrrHot = PrrHot + sum(DFhot(:,fix(Nv/2):Nv,i,end) .* VrHot(:,fix(Nv/2):Nv).^2, "all");
    qrHot = qrHot + sum(DFhot(:,fix(Nv/2):Nv,i,end) .* VrHot(:,fix(Nv/2):Nv) .* VsqrHot(:,fix(Nv/2):Nv,i), "all");
end
PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;
qrHot = qrHot * mg * eVtoErg / c^2 * 0.5 * gridStepHot^3;

PrrCold = 0;
qrCold = 0;
for i = 1:Nv
    PrrCold = PrrCold + sum(DFCold(:,fix(Nv/2):Nv,i,end) .* VrCold(:,fix(Nv/2):Nv).^2, "all");
    qrCold = qrCold + sum(DFCold(:,fix(Nv/2):Nv,i,end) .* VrCold(:,fix(Nv/2):Nv) .* VsqrCold(:,fix(Nv/2):Nv,i), "all");
end
PrrCold = PrrCold * mg * eVtoErg / c^2 * gridStepCold^3;
qrCold = qrCold * mg * eVtoErg / c^2 * 0.5 * gridStepCold^3;

Prr = PrrCold + PrrHot;
qrr = qrHot + qrCold;

NOut = Prr ./((Twall * eVtoErg)^(3/2) - qrr * a / kappa .* log(R./b)).^(2/3);

xout = R/b;
figure(3);
semilogy(xout,NOut,'LineWidth',2);
legend('n_{out}');
lgd1 = legend;
lgd1.FontSize = 18;