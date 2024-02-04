%Константы задачи
c = 3 * 10^10;
eVtoErg = 1.6e-12;
mp = 938.27 * 10^6;
mg = 2 * mp;
a = 50;
b = 250;
r = linspace(0,a,35);
R = linspace(a,b,35);
n0 = 5e+12;
T0 = 0.5;
Twall = 0.026;
np = ones(size(r)) * 1e11;
Tp = ones(size(r)) * 100;
diffCross = 1e-16;
kappa = 5/9/diffCross/sqrt(mg * eVtoErg / c^2);
Nv = 300;

[DFCold, gridStepCold, VrCold, VsqrCold] = PlotDFCold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv, false);
nCold = Ncold(np, Tp, mp, n0, T0, mg, diffCross, r, a, Nv);
[DFhot, gridStepHot, VrHot, VsqrHot] = PlotDFHot(np, Tp, mp, T0, nCold, mg, diffCross, r, a, Nv, false);
nHot = zeros(1,length(r));
for i = 1:length(r)
    nHot(i) = sum(DFhot(:,:,:,i), "all");
end
nHot = nHot * (gridStepHot)^3;
figure(26);
x = r/a;
plot(x,nCold,'LineWidth',2);
hold on
plot(x,nHot,'LineWidth',2);
plot(x,nCold + nHot,'LineWidth',2);
legend('n_{cold}','n_{hot}','n_{full}');
lgd = legend;
lgd.FontSize = 18;
% PrrHot = 0;
% qrHot = 0;
% for i = 1:Nv
%     PrrHot = PrrHot + sum(DFhot(:,fix(Nv/2):Nv,i,end) .* VrHot(:,fix(Nv/2):Nv).^2, "all");
%     qrHot = qrHot + sum(DFhot(:,fix(Nv/2):Nv,i,end) .* VrHot(:,fix(Nv/2):Nv) .* VsqrHot(:,fix(Nv/2):Nv,i), "all");
% end
% PrrHot = PrrHot * mg * eVtoErg / c^2 * gridStepHot^3;
% qrHot = qrHot * mg * eVtoErg / c^2 * 0.5 * gridStepHot^3;
% 
% PrrCold = 0;
% qrCold = 0;
% for i = 1:Nv
%     PrrCold = PrrCold + sum(DFCold(:,fix(Nv/2):Nv,i,end) .* VrCold(:,fix(Nv/2):Nv).^2, "all");
%     qrCold = qrCold + sum(DFCold(:,fix(Nv/2):Nv,i,end) .* VrCold(:,fix(Nv/2):Nv) .* VsqrCold(:,fix(Nv/2):Nv,i), "all");
% end
% PrrCold = PrrCold * mg * eVtoErg / c^2 * gridStepCold^3;
% qrCold = qrCold * mg * eVtoErg / c^2 * 0.5 * gridStepCold^3;
% 
% Prr = PrrCold + PrrHot;
% qrr = qrHot + qrCold;
% 
% NOut = Prr ./((Twall * eVtoErg)^(3/2) - qrr * a / kappa .* log(R./b)).^(2/3);
% 
% xout = R/b;
% figure(27);
% plot(xout,NOut,'LineWidth',2);
% legend('n_{out}');
% lgd1 = legend;
% lgd1.FontSize = 18;