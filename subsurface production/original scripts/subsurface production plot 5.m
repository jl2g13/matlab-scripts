% subsurface production
% plot 5: spatial plots of
%             i) Pregrowth (Apr) DIN integrated over 0-20 m
%             ii) Ice minimum (Sept) PAR at 20 m


% Started jl2g13 (02/09/14)


% i) DIN plot
% ===========

% (1) Reading in monthly model variables - CHl, PP, PAR, SST    
    
% Only considering 2005
yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(4); % top 20m
months = double(12); % 1 year

    
% Preallocating arrays for model output to be read into
modN = NaN(lontotal,lattotal,depth,months);
mld = NaN(lontotal,lattotal,months);

pos = 1;
% reading in monthly files
for m = 1:12
    if m < 10
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dP.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dT.nc', yr, yr, m);
        dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dI.nc', yr, yr, m);        
    else
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dP.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dT.nc', yr, yr, m);
        dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dI.nc', yr, yr, m);
    end

    fnamesP = dir(dnomP);
    fnamesT = dir(dnomT);
    fnamesI = dir(dnomI);
    
    fprintf('Reading month %d\n', m);    
    
    % (1) Dissolved Inorganic Nitrogen   
    
    fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP.name);
    fprintf('Reading %s\n', fnameP);

    % (a) Depth-resolved

    % Reading in nitrogen model field - DIN
    t1 = ncread(fnameP,'DIN',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land) - DIN might actually be zero
    modN(:,:,:,pos) = t1;
    
    
     % (2) Mixed layer depth
           
    fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT.name);
    fprintf('Reading %s\n', fnameT);
 
     % Reading in mixed layer depth model field
    t1 = ncread(fnameT,'somxl010',[1 latstart 1], [lontotal lattotal 1]); % mixing layer depth -> Tfile,somixhgt (not so meaningful for monthly averages)
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    mld(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
    
    % (3) For Sea Ice Cover
        
    % Reading in sea ice cover model field
    t1 = ncread(fnameT,'soicecov',[1 latstart 1], [lontotal lattotal 1]); % sea ice cover is a 2D model field, has no depth dimension
    t1 = squeeze(t1);
    % t1 (t1 == 0) = NaN; % removing zero values (land)

    icecov(:,:,pos) = t1;
    
     % (6) For Sea Ice Thickness

    % Extracting model sea ice thickness
    fnameI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesI.name);
    fprintf('Reading %s\n', fnameI);

    % Reading in sea ice thickness model field
    t1 = ncread(fnameI,'iicethic',[1 latstart 1], [lontotal lattotal 1]); % sea ice thickness is a 2D model field, has no depth dimension
    t1 = squeeze(t1);
    % t1 (t1 == 0) = NaN; % removing zero values (land)

    icethk(:,:,pos) = t1;

    
    pos = pos + 1;

end

icethk (icethk == 9.969209968386869e+36) = NaN;

clearvars t1 % clear up memory


% Plotting 'pregrowth' DIN over 

% Setting up M_Map model
ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

% removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
land = xx(1,1);
xx(xx == land) = NaN; xx(xx == 0) = NaN;
land = yy(1,1);
yy(yy == land) = NaN; yy(yy == 0) = NaN;



% check when  mld deepest (before Apr growth) -> Mar (note mld deeper in Apr)
figure;
odvpal(50);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,mld(:,:,12)); shading flat; hold on; % month = march
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 100]); c = colorbar; ylabel(c,'m');
title('MLD','FontWeight','bold')


% plotting DIN int 0-20 for march
dep_nemo % function that reads in nemo's depth grid
intN = squeeze(nansum(modN,3) * (nemo_dep(4)  / size(modN,3))); % units -> mmolN/m2. nemo_dep(4) is depth of column (~20m)





% option 1: 0.1 and 3 m ice contours for Jan (tried xm contour for Jan -
% Feb - Mar but ice margin doesn't move much over this period -> moves
% DURING growing season
figure;
odvpal(50);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,intN(:,:,2)); shading flat; hold on; % month = march
m_contour(xx,yy,icethk(:,:,1),[2.9],'edgecolor','k','LineWidth',1.5); hold on % Jan 3m icethk contour
m_contour(xx,yy,icethk(:,:,1),[0.2],'edgecolor','b','LineWidth',1.5); hold on % Feb 3m icethk contour 
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 200]); c = colorbar; ylabel(c,'mmolN/m2');
title('DIN int 0-20m Mar','FontWeight','bold');


print -depsc -painters DIN_20-0m_Mar_ice_overlay.eps


% check whether spatial distribution of DIN in march maps best onto Mar ice
% cover or thickness

figure;
odvpal(50);
subplot(2,1,1);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,icecov(:,:,1)); shading flat; hold on; % month = march
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 1]); c = colorbar; ylabel(c,'mmolN/m2');
title('Ice cover Oct','FontWeight','bold');

subplot(2,1,2);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,icethk(:,:,1)); shading flat; hold on; % month = march
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 5]); c = colorbar; ylabel(c,'mmolN/m2');
title('Ice thk Dec','FontWeight','bold');


clearvars -except xx yy icecov







% ii) PAR plot
% ============

% (1) Read in model PAR and Chl
% -----

pos = 1;

lattotal = double(1021);
lontotal = double(1442);
depth = double(4); % top 20 m
yr = 2005;
months = 12; % 1 year


% Preallocating arrays for model output to be read into
totalshortwave = zeros(lontotal,lattotal,months);
modchn = zeros(lontotal,lattotal,depth,months);
modchd = zeros(lontotal,lattotal,depth,months);

    
    % Getting list of all files for this year
for m = 1:12
    if m < 10
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dD.nc', yr, yr, m);
    else
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dD.nc', yr, yr, m);
    end
        


    fnamesP = dir(dnomP);
    fnamesD= dir(dnomD);
    
    fprintf('Reading month %d\n', m); 
    
    % Initiating extraction loop
            
        % (1) For Chlorophyll
        
        fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP.name);
        fprintf('Reading %s\n', fnameP);
        
        % (i) For non-diatom chlorophyll
        t1 = ncread(fnameP,'CHN',[1 1 1 1], [lontotal lattotal depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN;
        
        modchn(:,:,:,pos) = t1;
        
        % (ii) For diatom chlorophyll
        t1 = ncread(fnameP,'CHD',[1 1 1 1], [lontotal lattotal depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN;
        
        modchd(:,:,:,pos) = t1;
        
        
        
      % (2) For Photosynthetically Active Radiation
        
        fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD.name);
        fprintf('Reading %s\n', fnameD);
        
        % Reading in surface PAR model field
        t1 = ncread(fnameD,'MED_QSR',[1 1 1], [lontotal lattotal 1]);
        t1 (t1 == 0) = NaN;
        % don't squeeze because want depth dim for multiplication below
        
        totalshortwave(:,:,pos) = t1;
        
        pos = pos + 1;
        
end

% Freeing up some memory
modchl = modchn + modchd;
clearvars modchn modchd t1;


% (2) Attenuate model par with depth (Levy 2001)
% -----

% Constants for light attenuation - Levy2001 Table 3
rpig = 0.7;
xkr0 = 0.225;
xkrp = 0.037;
xlr = 0.629;
xlg = 0.674;
xkg0 = 0.0232;
xkgp = 0.074;
surface_dim = 1; % par a surface-only field

dep_nemo; % for model grid depth

zpar = totalshortwave * (0.43 / 2); % Levy2001 eq A23. Divide by two because two-band model.
zpar = reshape(zpar,lontotal,lattotal,surface_dim,months); % adding depth dim back in (=1, making explicit)

zparr = zpar; % par for first depth level = par of model output
zparg = zpar;

modchl (modchl == 0) = realmin; % realmin: smallest +ve floating-point number (Fortran: TINY). because can't have log(zpig=0)
zpig = modchl ./ rpig; % attenuates due to pigments not biomass - Levy2001 A22


clearvars totalshortwave modchl zpar % freeing up some more memory


posz = 2;
for z = 2:size(zpig,3)
    
    thick = nemo_dep(z)-nemo_dep(z-1); % depth box thickness
    
    zkr  = xkr0 + xkrp .* exp(xlr .* log(zpig(:,:,(z-1),:))); % red wavelengths attenuation coefficient - Levy2001 eq A20.
    zkg  = xkg0 + xkgp .* exp(xlg .* log(zpig(:,:,(z-1),:))); % green wavelengths attenuation coefficient - Levy2001 eq A21
    
    zparr(:,:,posz,:) = zparr(:,:,(z-1),:) .* exp( -zkr * thick);  % attenuation of red light with depth - Levy2001 eq A24
    zparg(:,:,posz,:) = zparg(:,:,(z-1),:) .* exp( -zkg * thick);  % attenuation of green light with depth - Levy2001 eq A25

    par_J = zparr + zparg; % Levy2001 eq A26
    fprintf('Read level %d\n', z)

    posz = posz + 1;
    
    % you need zparr for each depth in the loop because it's used to
    % calculate next depth - need to store it somehow -> this is why you do
    % line 96 and make zpar size(z dim) = 64
end
clearvars zparg zparr zkr zkg zpig

h = 6.63e-34; c = 3e8; length = 550e-9; % 'length' a crude approximation (depends on light spectrum with depth)
par_E = par_J / ((h*c/length) * 6.022D17); % J/m2/s -> umolphotons/m2/s
% par_J = par_E * ((h*c/length) * 6.022D17)


% NB PAR_J(:,:,1,:) IS NOT RIGHT, I.E. SURFACE PAR FIELD NOT RIGHT - ERROR
% IN CODE

% plotting Sept par_J at 20 m depth
t1 = par_J(:,:,4,9); 

% plotting jul-sep av par_J at 20 m 

% t1 = par_J(:,:,4,7:9) / 3; % 3 = size par_J,4 summed over - to retain units 

% plotting increase (difference) in light from Mar -> Jul
% t1 = par_J(:,:,4,7) - par_J(:,:,4,3);

figure;
odvpal(50);
m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,t1); shading flat; hold on;
m_contour(xx,yy,icecov(:,:,9),[0.95],'edgecolor','k','LineWidth',1.5); hold on % sep 0.9 icecov contour
m_contour(xx,yy,icecov(:,:,9),[0.8],'edgecolor','r','LineWidth',1.5); hold on % sep 0.5 icecov contour
 m_contour(xx,yy,icecov(:,:,9),[0.92],'edgecolor','y','LineWidth',1.5); hold on % sep 0.5 icecov contour
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 5]); c = colorbar; ylabel(c,'W/m2');
title('PAR at 20m Sep av','FontWeight','bold')

print -depsc -painters par_20m_Sep_icethk_overlay.eps

