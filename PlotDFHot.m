% Function PlotDFHot make plot of hot part of gas distibution function inside the plasma
%   Input parameters:
% np :
%       type        : 1D double precision array
%       units       : cm^(-3)
%       description : plasma ion's density
% Tp :
%       type        : 1D double precision array
%       units       : eV
%       description : plasma ion's temperature
% mp :
%       type        : double precision
%       units       : eV
%       description : plasma ion's mass
% T0 :
%       type        : double precision
%       units       : eV
%       description : The half of the Maxwell DF temperature penetrate
%                     inside the plasma column
% Ncold :
%       type        : 1D double precision array
%       units       : cm^(-3)
%       desciption  : Density of cold part of DF inside plasma
% mg :
%       type        : double precision
%       units       : eV
%       description : gas mass
% diffCross :
%       type        : double precision
%       units       : cm^2
%       description : gas-plasma ellastic collisions differential
%                     cross section. Asumming that diff cross section is
%                     constant!!!
% r  :
%       type        : 1D double precision array
%       units       : cm
%       description : space grid inside the plasma column. r(end) = a.
% a  : 
%       type        : double precision
%       units       : cm
%       description : Radius of plasma column
% Nv :  type        : integer
%       units       : unitless
%       description : The number of points in the velocity grid
%   Output parameters:
% DF:
%       type        : 4D double precision array
%       units       : s^3 * cm^(-6)
%       description : DF of the gas hot part. Size of [Nv x Nv x Nv x Nr],
%                     where Nr - number of points in space grid
% All input 1D arrays must be the same size!!!
function DF = PlotDFHot(np, Tp, mp, T0, Ncold, mg, diffCross, r, a, Nv)
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp);
    VT0 = sqrt(2 * T0 / mg);
    Q = diffCross ./ (2 * alpha^3 * sqrt(pi) * VTp.^2) .* np;
    VelGridBoundary = eps * VT0 / max(VTp);                       % V < Veps domen is no need for this solution

    UGrid = linspace(-eps,eps,Nv);
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;

    P = zeros(Nv, Nv, Nv,length(r));
    DF0 = zeros(Nv, Nv, Nv);
    DF = zeros(Nv, Nv, Nv, length(r)-1);
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary;
    end

    intNcold0 = trapz(r,Ncold);
    intNcold = zeros(Nv, Nv, length(r));
    x = (r / a);
    for k = 2:length(r)
        for j = 1:Nv
            for i = 1:Nv
                sinPhiSqr = (r(k) / a)^2 * UGrid(i)^2 / (UGrid(j)^2 + UGrid(i)^2);
                if UGrid(j) <= 0
                    integrationDomen = (x.^2 > max([sinPhiSqr, x(k)^2])) ./ sqrt(1 - sinPhiSqr ./ x.^2);
                    intNcold(i, j, k) = trapz(x, integrationDomen .* Ncold);
                else
                    integrationDomen = (x.^2 > sinPhiSqr) ./ sqrt(1 - sinPhiSqr ./ x.^2);
                    intNcold(i, j, k) = trapz(x, integrationDomen .* Ncold);
                    if x(k)^2 > sinPhiSqr
                        integrationDomen2 = (x <= x(k)) .* integrationDomen;
                        intNcold(i, j, k) = intNcold(i, j, k) + ...
                            trapz(x, integrationDomen2 .* Ncold);
                    end
                end
            end
        end
    end

    for j = 1:length(r)
        for i = 1:Nv
            P(:,:,i,j) = (2 * alpha) * Q(j) ./ sqrt(UOrtSqr + Uz(i)^2) .* exp(- (UOrtSqr + Uz(i)^2) ./ (4 * alpha^2));
        end
        P(:,:,:,j) = P(:,:,:,j) .* boundaryFactor;
    end

    for i = 1:Nv
        DF0(:,:,i) = P(:,:,i,1) ./ sqrt(UOrtSqr) * intNcold0;
    end
    DF0 = DF0 * a;

    DF(:,:,:,1) = DF0;
    for k = 2:length(r)
        for j = 1:Nv
            DF(:,:,j,k) = P(:,:,j,k) ./ sqrt(UOrtSqr(:,:)) .* intNcold(:,:,k);
        end
    end
    DF = DF * a;
    
    figure(1);
    pcolor(UGrid, UGrid, DF0(:,:,fix(Nv/2)));
    shading flat;
    shading interp;
    title('F(r = 0)'); 
    xlabel('v_r/v_{Tp}'); 
    ylabel('v_\phi/v_{Tp}');

    for i = 2:length(r)
        figure(i);
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),i));
        shading flat;
        shading interp;
        title('F(r = a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
    end

end

