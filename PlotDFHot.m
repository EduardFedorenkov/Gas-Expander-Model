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
function [DF, gridStep] = PlotDFHot(np, Tp, mp, T0, Ncold, mg, diffCross, r, a, Nv)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = diffCross ./ (2 * alpha^3 * sqrt(pi) * VTp.^2) .* np;
    VelGridBoundary = eps * VT0 / max(VTp);                       % V < Veps domen is no need for this solution

    UGrid = linspace(-1.1*eps,1.1*eps,Nv);
    gridStep = (UGrid(2) - UGrid(1)) * max(VTp);
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;

    P = zeros(Nv, Nv, Nv,length(r));
    DF = zeros(Nv, Nv, Nv, length(r));
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary;
    end

    intNcold = zeros(Nv, Nv, length(r));
    for k = 1:length(r)
        for j = 1:Nv
            for i = 1:Nv
                sinPhiSqr = (r(k) / a)^2 * UGrid(i)^2 / (UGrid(j)^2 + UGrid(i)^2);
                if UGrid(j) > 0
                    intNcold(i,j,k) = NcoldRtoL(Ncold, r, a, r(k), sinPhiSqr, 1, length(r));
                else
                    intNcold(i,j,k) = NcoldRtoL(Ncold, r, a, r(k), sinPhiSqr, 0, length(r));
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

    for k = 1:length(r)
        for j = 1:Nv
            DF(:,:,j,k) = P(:,:,j,k) ./ sqrt(UOrtSqr(:,:)) .* intNcold(:,:,k);
        end
    end

    for i = 1:length(r)
        figure(i);
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),i));
        shading flat;
        shading interp;
        title('F(r = a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
    end

end

