% Final figure 2

% To do: set white space to be  removeed from figure automatically


% Created by jl2g13 (04/09/13) from oringinal scripts 'subsurface
% production plot 1' and subsurface production plot 5'

% Plot has 3 panels as follows
% a. PP 20-60m prop (as figure 1c)
% b. DIN 0-20m Mar with Jan-Mar ice 0.2 and 3 Mar ice thickness contours 
% c. PAR 20m Sep with 0.8, 0.92 and 0.95 Sep ice cover contours



% Reading in DIN and ice contour fields
% -------------------------------------

yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(10); % top 60m
months = double(12); % 1 year

    
% Preallocating arrays for model output to be read into
modN = NaN(lontotal,lattotal,depth,months);
moddpp = NaN(lontotal,lattotal,depth,months);
modndpp = NaN(lontotal,lattotal,depth,months);
icecov = NaN(lontotal,lattotal,months);
icethk = NaN(lontotal,lattotal,months);

pos = 1;
% reading in monthly files
for m = 1:12
    if m < 10
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dD.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dT.nc', yr, yr, m);
        dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dI.nc', yr, yr, m);        
    else
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dD.nc', yr, yr, m);
        dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dT.nc', yr, yr, m);
        dnomI = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dI.nc', yr, yr, m);
    end

    fnamesP = dir(dnomP);
    fnamesD = dir(dnomD);
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
    
    
    
    % (2) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD.name);
    fprintf('Reading %s\n', fnameD);

    % (a) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (b) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land)

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
    
    
    % (3) For Sea Ice Cover
    
    fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT.name);
    fprintf('Reading %s\n', fnameT);
    
    % Reading in sea ice cover model field
    t1 = ncread(fnameT,'soicecov',[1 latstart 1], [lontotal lattotal 1]); % sea ice cover is a 2D model field, has no depth dimension
    t1 = squeeze(t1);
    % t1 (t1 == 0) = NaN; % removing zero values (land)

    icecov(:,:,pos) = t1;
    
    
    
    % (4) For Sea Ice Thickness

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
modndpp (modndpp == 9.969209968386869e+36) = NaN; % what on earth is this flag? - netcdf fill value, probably occurs when read in problem
moddpp (moddpp == 9.969209968386869e+36) = NaN;

modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC

clearvars modndpp moddpp




% Creating DIN field: integrated 0-20 for march
% ---------------------------------------------

modN = modN(:,:,1:4,:); % only want top 20m (but read in 60 because needed for PP prop field)
dep_nemo % function that reads in nemo's depth grid
intN = squeeze(nansum(modN,3) * (nemo_dep(4)  / size(modN,3))); % units -> mmolN/m2. nemo_dep(4) is depth of column (~20m)

clearvars t1 modN



% creating PP prop 20-60m field
% -----------------------------

modppC_an = (nansum(modppC,4)) .* (365/months); % Annual PP (integrating over all months in year). days in yr/no files. mgC/m3/yr
   
% Water depth filter (60m)
t1 = modppC(:,:,:,1) * 0; % only want filter at top 60 m -> not over whole grid. nb. n=10, is z~60m. 
t2 = sum(isfinite(t1),3); % chl = finite in water, = NaN below seafloor
kbathy = t2; % model water depth in terms of k model levels
kbathy (kbathy == 0) = NaN; % from here just times whatever model field you want to filter by kbathy
kbathy (kbathy >= 1) = 1; % when sums over finite values each x-y water column gets value 1-10 if > 60m deep (would be 10 if didn't set == 0 -> NaN at read in above)

clearvars t1 t2
% clearvars modppC

pp_surf = squeeze(nansum(modppC_an(:,:,1:4),3) * (20/4) / 1000); % surface PP = integrate top 4 boxs ~ top 20m. depth-interval/noboxes. gC/m2/yr
pp_sub = squeeze(nansum(modppC_an(:,:,5:10),3) * (40/6) / 1000); % deep PP = integrate boxes 5-10 ~ 20-60 m depth. gC/m2/yr

prop_t = pp_sub ./ (pp_surf + pp_sub);
prop = prop_t .* kbathy;

clearvars pp_surf pp_sub kbathy modppC_an prop_t








% Reading in light field and creating PAR(z) field
% ------------------------------------------------


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

% convert par field to einsteins
% h = 6.63e-34; c = 3e8; length = 550e-9; % 'length' a crude approximation (depends on light spectrum with depth)
% par_E = par_J / ((h*c/length) * 6.022D17); % J/m2/s -> umolphotons/m2/s
% par_J = par_E * ((h*c/length) * 6.022D17)


% NB PAR_J(:,:,1,:) IS NOT RIGHT, I.E. SURFACE PAR FIELD NOT RIGHT - ERROR
% IN CODE

% plotting Sept par_J at 20 m depth
t1 = par_J(:,:,4,9); 












% ======
% PLOT 2
% ======

% M_map plotting preamble
ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

% removing spurious spikes in N Atl when plotting (beware true zeros, quick fix)
land = xx(1,1);
xx(xx == land) = NaN; xx(xx == 0) = NaN;
land = yy(1,1);
yy(yy == land) = NaN; yy(yy == 0) = NaN;

fig = figure;
odvpal(50);

% PP 20-60m prop
% positionVector1 = [0.1, 0.8, 0.6, 0.3]; % first two [left bottom width height] (values 0-1) . left-bottom positions bottom left corner of figure, width-height to set size
% subplot('Position',positionVector1);
subplot(3,1,1)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-500 -500],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,prop); shading flat; hold on; % boxes 5--9 ~ 20-60 m depth 
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis([0 0.6]); cb = colorbar; ylabel(cb,'dimensionless');
%     title('Annual proportion 20-60 m','FontWeight','bold')
title('a.','position',[-0.5 0.6],'FontWeight','bold')


% DIN Mar
% positionVector1 = [0.1, 0.5, 0.6, 0.3];
% subplot('Position',positionVector1);
subplot(3,1,2)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,intN(:,:,2)); shading flat; hold on; % month = march
m_contour(xx,yy,icethk(:,:,1),[2.9],'edgecolor','k','LineWidth',0.5); hold on % Jan 3m icethk contour
m_contour(xx,yy,icethk(:,:,1),[0.2],'edgecolor','b','LineStyle','--','LineWidth',0.5); hold on % Feb 3m icethk contour 
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 200]); c = colorbar; ylabel(c,'mmolN/m2');
title('b.','position',[-0.5 0.6],'FontWeight','bold')

% PAR Sep
% positionVector1 = [0.1, 0.1, 0.6, 0.3];
% subplot('Position',positionVector1);
subplot(3,1,3)

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,t1); shading flat; hold on;
m_contour(xx,yy,icecov(:,:,9),[0.95],'edgecolor','k','LineWidth',0.5); hold on % sep 0.9 icecov contour
m_contour(xx,yy,icecov(:,:,9),[0.92],'edgecolor','b','LineStyle','--','LineWidth',0.5); hold on % sep 0.5 icecov contour
m_contour(xx,yy,icecov(:,:,9),[0.8],'edgecolor','g','LineStyle','-.','LineWidth',0.5); hold on % sep 0.5 icecov contour
m_coast('patch',[.7 .7 .7],'edgecolor','k');
caxis([0 5]); c = colorbar; ylabel(c,'W/m2');
title('c.','position',[-0.5 0.6],'FontWeight','bold')

tightfig; % remove white space around edges of figure

set(gcf,'PaperPositionMode','auto') % have to manually drag matlab figure window to desired fig size before executing these lines
print -depsc -painters -r0 fig2_test.eps


% % 2) check in ghostscript and print
% lpr -Pspartacus filename % print to spartacus
% lpq -Pspartacus % just to check it's printing
