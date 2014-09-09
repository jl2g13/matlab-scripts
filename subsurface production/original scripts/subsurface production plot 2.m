% Subsurface production
%    Plot 2: bar plot of 0-20 and 20-60 PP and DIN by month


    
% (1) Reading in monthly model variables - CHl, PP, PAR, SST    
    
% Only considering 2005
yr = 2005; 
    
lattotal = double(1021);
latstart = double(1);
lontotal = double(1442);
depth = double(14); % top 100 m

    
% Preallocating arrays for model output to be read into
modndpp = NaN(lontotal,lattotal,depth,12); % 12 months
moddpp = NaN(lontotal,lattotal,depth,12);

pos = 1;
% reading in monthly files
for m = 1:12
    if m < 10
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dD.nc', yr, yr, m);
%         dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm0%dT.nc', yr, yr, m);
    else
        dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dP.nc', yr, yr, m);
        dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dD.nc', yr, yr, m);
%         dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%dm%dT.nc', yr, yr, m);
    end


    fnamesD = dir(dnomD);
    fnamesP = dir(dnomP);
    
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
    
    
    
    % (2) Dissolved Inorganic Nitrogen   
    
    fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP.name);
    fprintf('Reading %s\n', fnameP);

    % (a) Depth-resolved

    % Reading in nitrogen model field - DIN
    t1 = ncread(fnameP,'DIN',[1 latstart 1 1], [lontotal lattotal depth 1]);
    t1 (t1 == 0) = NaN; % removing zero values (land) - DIN might actually be zero
    modN(:,:,:,pos) = t1;

    
    pos = pos + 1;

end

modndpp (modndpp == 9.969209968386869e+36) = NaN; % what on earth is this flag? - netcdf fill value, probably occurs when read in problem
moddpp (moddpp == 9.969209968386869e+36) = NaN;
modN (modN == 9.969209968386869e+36) = NaN; 


modppC = (modndpp + moddpp) * 12.011 * 6.625; % in mgC/m3/d

clearvars modndpp moddpp t1 % clear up memory



% =======
% PLOT 2
% =======


% integrate over depth
modppC_20 = squeeze(nansum(modppC(:,:,1:4,:),3) ./ 4); % integrating over 0-20 m depth
modN_20 = nansum(modN(:,:,1:4,:),3) ./ 4;

modppC_60 = nansum(modppC(:,:,5:9,:),3) ./ 5; % integrating over 20-60 m depth
modN_60 = nansum(modN(:,:,5:9,:),3) ./ 5;




% % % Check 0-20 against 20-100 m
% % ------------
%  modppC_60 = squeeze(nansum(modppC(:,:,5:14,:),3) ./ 10); % integrating over 20-60 m depth
%  modN_60 = nansum(modN(:,:,5:14,:),3) ./ 10;
%  
%  clearvars modN modppC





% intergrate over area (x,y)
modppC_20 = squeeze(nansum(nansum(modppC_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mgC/m3/d)
modppC_60 = squeeze(nansum(nansum(modppC_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442); 

modN_20 = squeeze(nansum(nansum(modN_20(:,831:end,:),2),1)) ./ ((1021-830) * 1442); % integrating over lat and lon of Arctic (lat > 831 ~ 66° N = Arctic circ ('~' because grid tripolar)) (mmolN/m3)
modN_60 = squeeze(nansum(nansum(modN_60(:,831:end,:),2),1)) ./ ((1021-830) * 1442); 




% plotting bar graph


% check plot

% bar(1:12, [modppC_20 modppC_60], 0.5, 'stack');
% legend('PP 0-20','PP 20-60')
% title('PP')
% 
% bar(1:12, [modN_20 modN_60], 0.5, 'stack');
% legend('DIN 0-20','DIN 20-60')
% title('DIN')

% actual plot


NumStacksPerGroup = 2; % N and PP
NumGroupsPerAxis = 12; % months
NumStackElements = 2; % depth windows


% labels to use on tick marks for groups 
% groupLabels = { 'J'; 'F'; 'M'; 'A'; 'M'; 'J'; 'J'; 'A'; 'S'; 'O'; 'N'; 'D' }; % Months
groupLabels = { 1:12 }; % Months
stackData = reshape([modN_60 modppC_60 modN_20 modppC_20], 12,2,2); % dim ordering needs to be (NumGroupsPerAxis,NumStacksPerGroup,NumStackElements) = [v1_s1 v2_s1 v1_s2 v2_s2 v1_s3 v2_s3 etc]

plotBarStackGroups(stackData, groupLabels); 
set(gca,'FontSize',9) 
set(gcf,'Position',[100 100 720 650]) 
grid off 
set(gca,'Layer','top') % put grid lines on top of stacks

% title('Surface and subsurface primary production','FontWeight','bold'); % Add title and axis labels
xlabel('Month','fontsize',10);
ylabel({'Primary production (mgC/m3/d)','DIN (mmolN/m3)'},'fontsize',10); % cell array to split title over two lines

% set ylabel horizontal - need to get 'ylabel' handle
% set(get(gca,'YLabel'),'Rotation',0)


% legend('DIN 60','Din 20','PP 60', 'PP 20')


print -depsc -painters fig2.eps

% Calculate proportion of PP that happens 0-20 and 20-60 by month
zint_20 = 20; % integrate over depth interval
zint_60 = 40;
% zint_60 = 80; % check 0-20 against 20-100
PP_proport = (modppC_20 * zint_20) ./ (modppC_60 * zint_60) % mgC/m2/d


% % for 0-20 vs 0-100 m
% 
% PP_proport =
% 
%    4.066977780976219
%    3.578788773569942
%    3.271714876624915
%    3.124595195433111
%    2.866966265969533
%    2.120026501247386
%    1.656189944270219
%    1.551430546971300
%    2.325358628413311
%    4.649134342280561
%    5.515067811102148
%    4.168191495086438


% May - Sept
grow = (2.867 + 2.120 + 1.656 + 1.551 + 2.325) / 5;












% =======
% PLOT 5
% =======

% Primary production by sector

% splitting production into sectors of Arrigo2008b
    % sectors used (by eye from Arrigo2008b Fig 1, using Google Earth)
    % Greenland 315-15E
    % Barents 15-55E
    % Kara 55-105E
    % Laptev 105-150E
    % East Siberian Sea 150-180E
    % Chuchki 180-200E
    % Beaufort 200-260E
    % Baffin 260-315E
    
    % Tripolar grid - regrid model output then cut it by sector






