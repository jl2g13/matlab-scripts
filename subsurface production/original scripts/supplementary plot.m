% Supplementary plot: Cross correlation surface (3m) to depth integrated
% PP (0-100m) using MONTHLY PP files

% expect XCorr to be higher than for 5day files

yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(14); % top 60 m

    
% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,12); % 12 months
moddpp = NaN(lontotal,lattotal,depth,12);

pos = 1;
% reading in monthly files
for m = 1:12
    if m < 10
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dD.nc', yr, yr, m);
    else
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dD.nc', yr, yr, m);
    end


    fnamesD = dir(dnomD);
    
    fprintf('Reading month %d\n', m);    
    

    
    % (1) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD.name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    
    pos = pos + 1;

end

modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC/m3/d

surfPP = modppC(:,:,1,:);
intPP = nansum(modppC,3) / depth ; % mgC/m3/d


pos = 1;

 for j = 1:size(modppC,2) % for every 1/4 lattitude
        for i = 1:size(modppC,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = intPP(i,j,:);
            b = surfPP(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAn(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading lat = %d of 1021\n', j)
 end
 
 
 
 
 
% Plotting

supplot = reshape(corAn,1442,1021); 
 
% General plotting preamble
    % Setting up M_Map model
    ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
    xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
    yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

    % removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
    land = xx(1,1);
    xx(xx == land) = NaN; xx(xx == 0) = NaN;
    land = yy(1,1);
    yy(yy == land) = NaN; yy(yy == 0) = NaN;
    
    
figure(1);
odvpal(50);
% Issue - plotting nans in water depth < 100 m (depth = 14) and in
% central basin (perhaps PP = 0 at <= 100 m depth -> NaN so corAn -> NaN)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-100 -100],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,supplot); shading flat; hold on; % XCORR 
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis([0 1]); cb = colorbar; ylabel(cb,'unitless');
title('Annual monthly XCORR','FontWeight','bold')

print -depsc -painters xcorr_monthly_PP_0.8-1.eps
