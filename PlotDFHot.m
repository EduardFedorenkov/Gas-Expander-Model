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
function [DF, gridStep, Vr, Vsqr] = PlotDFHot(np, Tp, mp, T0, Ncold, mg, diffCross, r, a, Nv, isPlot)
    c = 3 * 10^10;
    eps = 3 / sqrt(2);
    alpha = mp / (mg + mp);
    VTp = sqrt(2 * Tp / mp) * c;
    VT0 = sqrt(2 * T0 / mg) * c;
    Q = diffCross ./ (2 * alpha^3 * sqrt(pi) * VTp.^2) .* np;
    VelGridBoundary = eps * VT0 / max(VTp);                       % V < Veps domen is no need for this solution

    UGrid = linspace(-eps * 1.5,eps * 1.5,Nv);
    gridStep = (UGrid(2) - UGrid(1)) * max(VTp);
    y0 = ones(1,Nv);
    Uz = UGrid;
    Ur = y0'*UGrid;
    Vr = Ur * max(VTp);
    Uphi = UGrid'*y0;
    UOrtSqr = Ur.^2 + Uphi.^2;
    Vsqr = zeros(Nv, Nv, Nv);
    for i = 1:Nv
        Vsqr(:,:,i) = UOrtSqr(:,:) + UGrid(i)^2;
    end
    Vsqr = Vsqr * max(VTp)^2;

    P = zeros(Nv, Nv, Nv,length(r));
    DF = zeros(Nv, Nv, Nv, length(r));
    boundaryFactor = zeros(Nv,Nv,Nv);
    for i = 1:Nv
        boundaryFactor(:,:,i) = (UOrtSqr + Uz(i)^2) >= VelGridBoundary^2;
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
    DF = DF / max(VTp);

    plotBoard = 3.325 * VT0 / max(VTp);

    if (isPlot)
       tiledlayout(3,3);

        ax1 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),1));
        shading flat;
        shading interp;
        title(ax1,'f_g(r = 0)'); 
        xlabel(ax1,'v_r/v_{Tp}'); 
        ylabel(ax1,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax2 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),2));
        shading flat;
        shading interp;
        title(ax2,'f_g(r = 0.125a)'); 
        xlabel(ax2,'v_r/v_{Tp}'); 
        ylabel(ax2,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax3 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),3));
        shading flat;
        shading interp;
        title(ax3,'f_g(r = 0.25a)'); 
        xlabel(ax3,'v_r/v_{Tp}'); 
        ylabel(ax3,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax4 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),4));
        shading flat;
        shading interp;
        title(ax4,'f_g(r = 0.375a)'); 
        xlabel(ax4,'v_r/v_{Tp}'); 
        ylabel(ax4,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax5 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),5));
        shading flat;
        shading interp;
        title(ax5,'f_g(r = 0.5a)'); 
        xlabel(ax5,'v_r/v_{Tp}'); 
        ylabel(ax5,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax6 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),6));
        shading flat;
        shading interp;
        title(ax6,'f_g(r = 0.625a)'); 
        xlabel(ax6,'v_r/v_{Tp}'); 
        ylabel(ax6,'v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax7 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),7));
        shading flat;
        shading interp;
        title(ax7,'f_g(r = 0.75a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax8 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),8));
        shading flat;
        shading interp;
        title(ax8,'f_g(r = 0.875a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log');
    
        ax9 = nexttile;
        pcolor(UGrid, UGrid, DF(:,:,fix(Nv/2),9));
        shading flat;
        shading interp;
        title(ax9,'f_g(r = a)'); 
        xlabel('v_r/v_{Tp}'); 
        ylabel('v_\phi/v_{Tp}');
        colormap('jet');
        hold on
        rectangle('Position',[-1,-1,2,2]*plotBoard,'Curvature',1,...
            'LineWidth',2,'EdgeColor','w','FaceColor','None')
        hold off
        set(gca,'colorscale','log'); 
    end
end

