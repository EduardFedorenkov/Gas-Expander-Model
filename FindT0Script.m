%Константы задачи
c = 3 * 10^10;
mp = 938.27 * 10^6;
mg = 2 * mp;
a = 50;
b = 250;
r = linspace(0,a,35);
R = linspace(a,b,35);
n0 = 6e11;
T0 = 0.466;
Twall = 0.026;
Pwall = 4.14;
%np = (1 - (r / a).^2) * 1e14;
np = ones(size(r)) * 1e11;
Tp = ones(size(r)) * 100;
diffCross = 1e-16;
Nv = 300;

eVtoErg = 1.6e-12;

[DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, 1, T0, mg, diffCross, r, a, Nv);
nCold = Ncold(np, Tp, mp, 1, T0, mg, diffCross, r, a, Nv);
[DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, Nv);

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

n00 = ComputeN0(Prr, Pwall, 1);

% kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
% Trange = [0.45 0.46 0.465 0.47 0.475 0.48 0.485 0.5];
% 
% F = ComputeT0(1, Trange, Twall, Pwall, kappa, np, Tp, mp, mg, diffCross, r, a, b, Nv);

% При Pwall = 4.14:
% T0 = 0.466        !!!
% n0 = 4.1058e+12   !!!