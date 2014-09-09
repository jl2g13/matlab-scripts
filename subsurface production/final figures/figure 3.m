% Final figure 3

% Created by jl2g13 (04/09/13) from original scripts 'subsurface production
% plot 1'

% Plot has 1 panel as follows
% a. PP Xcorr (5day) for surface (3m) vs depth integrated production (all)






% (1) Reading in 5-day model variable - PP   (need 5-day files for XCORR)

% Only considering 2005
yr = 2005;

lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(1); % 14 = read in top 100 m - need for XCORR (could read in depth-int PP var and depth = 1)
pos = 1;

   
dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr); % get a list of D (Diagnostic) files for this year
fnamesD = dir (dnomD);
[num, x] = size(fnamesD);


% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,num);
moddpp = NaN(lontotal,lattotal,depth,num);

    
    % Initiating extraction loop
    
for f = 1:1:num % looping from file 1 to num (73) - for each year


    % (1) Primary Production

    fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
    fprintf('Reading %s\n', fnameD);

    % (i) For non-diatom production

    t1 = ncread(fnameD,'PRN3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRN is depth-integrated so a 2D model field
    t1 = squeeze(t1);
     t1 (t1 == 0) = NaN; % removing zero values (land). some actual zeros
%     in central basin

    modndpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

    % (ii) For diatom production
    t1 = ncread(fnameD,'PRD3',[1 latstart depth 1], [lontotal lattotal depth 1]); % PRD is depth-integrated so a 2D model field
    t1 = squeeze(t1);
     t1 (t1 == 0) = NaN; % removing zero values (land). some actual zeros
%     in central basin

    moddpp(:,:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)
    
   
    % check values reasonable and most PP in 60 m: read in depth-integrated PP
    
    t1 = ncread(fnameD,'PRD',[1 latstart 1], [lontotal lattotal 1]);
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; 
    
    prd(:,:,pos) = t1;
    
    t1 = ncread(fnameD,'PRN',[1 latstart 1], [lontotal lattotal 1]);
    t1 = squeeze(t1);
    t1 (t1 == 0) = NaN; 
    
    prn(:,:,pos) = t1;
    
    pos = pos + 1;
    
end

modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC

modsurfppC = squeeze((modndpp + moddpp) * 12.011 * 6.625); % in mgC
modppC_int = (prn + prd) * 12.011 * 6.625;


clearvars modndpp moddpp prn prd







% ========
% PLOT 2    
% ========

pos = 1;

 for j = 1:size(modppC,2) % for every 1/4 lattitude
        for i = 1:size(modppC,1) % for every 1/4 longditude
            
            % Annual cross correlation
            a = modppC_int(i,j,:);
            b = modsurfppC(i,j,:);
            [c,lags] = xcorr(a,b,0,'coeff'); % zero lag coeff, lags = 0 
            corAn(:,pos) = c;
            
            pos = pos + 1;
            
        end
        fprintf('Reading lat = %d of 1021\n', j)
 end

plotXcorr = reshape(corAn,1442,1021);    
   
% M_map plotting preamble
ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);
    

figure;
odvpal(50);

m_proj('stereographic','lat',90,'long',30,'radius',30)
m_elev('contour',[-100 -100],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
m_pcolor(xx,yy,plotXcorr); shading flat; hold on; % XCORR 
m_coast('patch',[.7 .7 .7],'edgecolor','k'); 
caxis([0 1]); cb = colorbar; ylabel(cb,'Primary production cross-correlation (dimensionless)');
%     title('Annual XCORR (surface vs depth integrated','FontWeight','bold')

set(gcf,'Position',[0 0 4 3]); % set figure size


print -depsc -painters fig3_test.eps
