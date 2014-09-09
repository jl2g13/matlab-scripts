% testing script for under-ice PP

% Reads in variables, uses these to look at under ice by creating a mask
% then calculating monthly and annual total and under ice production, both
% globally and for north of the Arctic circle only. Plots bar of under,
% total and fraction-under production; also plots colorplot of spatial
% distributions of under ice production, both monthly and annual. (jl2g13, 15/02/14)

% EDIT (jl2g13, 26/02/14): padded arrays with NaN at padtime (line 74-75) rather
% than padding with -5 then changing -5 = NaN below the loop. (line 92-96).

for y = 1:1:1
    yr = 2004 + y;
    
    lattotal = double(1021);
    latstart = double(1);
    lontotal = double(1442);
    pos = 1;

    % combined
    for m = 1:12
        if m < 10
            dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05D.nc', yr, yr, m);
            dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d0%d*d05T.nc', yr, yr, m);
        else
            dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05D.nc', yr, yr, m);
            dnomT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/ORCA025-N201_%d%d*d05T.nc', yr, yr, m);
        end

        fnamesD = dir(dnomD);
        fnamesT = dir(dnomT);
        [num,x] = size(fnamesD);
        fprintf('Reading month %d\n', m);

            for f = 1:num % ie for each of the files within that given month - read in modpp

                % (1) Depth-integrated primary production

                % extract single file    
                fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
                fprintf('Reading %s\n', fnameD);


                % (i) For non-diatom production
                t1 = ncread(fnameD,'PRN',[1 latstart 1], [lontotal lattotal 1]); % PRN is depth-integrated so a 2D model field
                t1 (t1 == 0) = NaN; % removing zero values (land)
                modndpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)


                % (ii) For diatom production
                t1 = ncread(fnameD,'PRD',[1 latstart 1], [lontotal lattotal 1]); % PRD is depth-integrated so a 2D model field
                t1 (t1 == 0) = NaN; % removing zero values (land)
                moddpp(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)



                % (2) Sea Ice cover

                % extract single file    
                fnameT = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesT(f).name);
                fprintf('Reading %s\n', fnameT);


                % Reading in sea ice cover model field
                t1 = ncread(fnameT,'soicecov',[1 latstart 1], [lontotal lattotal 1]); % sea ice cover is a 2D model field, has no depth dimension
                modice(:,:,pos) = t1; % lattotal, lon, no of 5-day means (73/yr)

                pos = pos + 1;

            end

            % Padding - making arrays for each individual month the same size
            % so that they can be sequentially appended to form a single 4D
            % array. Pad value = NaN (ice may be zero)
            modndpp = padarray(modndpp, [0 0 (7-num)], NaN, 'post');
            moddpp = padarray(moddpp, [0 0 (7-num)], NaN, 'post');
            modice = padarray(modice, [0 0 (7-num)], NaN, 'post');

            % assign current month to array that will store all months for a given var
            modpp(:,:,:,m) = modndpp + moddpp; % t1 is non-diatom pp, t2 is diatom pp
            modicecov(:,:,:,m) = modice;


            % modpp = padarray(modpp,[0 0 (7-num)],'post');      % to be able to append each month to the same array requires that each appended array is the same size as the matrix it's being appended to - 7 is the most files that appear in a given month (December)
            clearvars modndpp moddpp modice % clear all variables read in within f-loop (otherwise when come to assign modice on second,third etc loops it can't assign if the num of files in the previous month isn't the same as the month it's overwriting with
            pos = 1; % must reset for the next month (loop): modndpp and modpp then overwrite preceeding months matricies (but only if num(monthN) >= num(month(N-1)) - this is why all model vars read in within f loop must be cleared

    end

%    EDIT: 
%     % setting padding to NaN
%     modpp (modpp == -10) = NaN; % checked modpp for (:,:,7,1) which should be fully padded (Jan only has 6 files) and all elements = -10 not -5
%     modicecov (modicecov == -10) = NaN;

    % converting pp units from to mgC (from mmol/N)
    modppC = modpp * 12.011 * 6.625;


    
    

    % Looking at under ice production

    % setting up ice mask and calculating under-ice primary production
    icemask = modicecov > 0.10;    % http://www.mathworks.co.uk/help/matlab/ref/logicaloperatorselementwise.html (see example 2)
    undericePP = modppC .* icemask; % returns PP values which are under ice cover >= 10% (of area of grid square)
    undericePP (undericePP == 0) = NaN;

    
    
    % Finding monthly productions

    % Maintaining time units whilst summing to get month bins (keep at mgC/m2/day)
    lenM = [31 28 31 30 31 30 31 31 30 31 30 31];
    filesM = [6 5 7 6 6 6 6 6 6 6 6 7];
    acirc = 831; % 66°N = acirc in model grid coordinates
    
    pos = 1;
        for i = 1:12
            % Global and Arctic under-ice production for each month

            t1 = undericePP(:,acirc:end,:,i); 
            uaM(:,:,pos) = nansum(t1,3) * lenM(i) / filesM(i);

            t2 = undericePP(:,:,:,i); % to plot in M_Map I want the whole world, not just Arctic
            uM(:,:,pos) = nansum(t2,3) * lenM(i) / filesM(i);

            % Global and Arctic total production for each month

            t3 = modppC(:,acirc:end,:,i);
            taM(:,:,pos) = nansum(t3,3) * lenM(i) / filesM(i);

            t4 = modppC(:,:,:,i);
            tM(:,:,pos) = nansum(t4,3) * lenM(i) / filesM(i);

            
            
            % Pan-Arctic production values - integrating over x and y

            t5 = taM(:,:,i); %(i) for total Arctic 
            totalM(pos) = sum(t5(:)) * lenM(i) * 1.4056e13 / ((lattotal - acirc + 1) * lontotal) ; % changing space units /m2 to pan-Arctic (Arctic area in m2 / Arctic model-grid elements)

            t6 = uaM(:,:,i);    %(ii) just under ice production
            underM(pos) = sum(t6(:)) * lenM(i) * 1.4056e13 / ((lattotal - acirc + 1) * lontotal);

            
            
            % Ice cover overlay for under ice spatial distribution plots
            t7 = modicecov(:,:,:,i);
            iceM(:,:,pos) = (nansum(t7,3)) / lenM(i);
       
            fprintf('Read month %d\n', i);

            pos = pos + 1;
        end


    % Annual Arctic and Global production (not x-y integrated) for 
    % (i) under ice production
    uaAn = nansum(uaM,3) * 365 / 12; % changing time units /day to /yr
    uAn = nansum(uM,3) * 365 / 12;

    % (ii) total production
    taAn = nansum(taM,3) * 365 / 12;
    tAn = nansum(tM,3) * 365 / 12;



    % Pan-Arctic production (x-y integrated) for total and under ice
    totalAn = sum(totalM(:)) * 365 / 12 / (sum(lenM)/12); % units /day to /year
    underAn = sum(underM(:)) * 365 / 12 / (sum(lenM)/12);

    under = [underM underAn] / 1e15 ; % mgc to TgC
    total = [totalM totalAn] / 1e15; % mgC to TgC

    frac = [(underM ./ totalM) (underAn ./ totalAn)];
    
    
    % Creating ice overlays for m_map plots
    
    pos  = 1;
    for i = 1:12

       pos = pos + 1; 
    end
   
    

    % Plots

    % (1) Bar plots to show fraction of production that is under ice

        figure(y);
        subplot(1,3,1)
        bar(under);
        title('Under-ice PP','FontWeight','bold'); ylabel('TgC'); xlabel('monthly (1-12) and annual production');

        subplot(1,3,2)
        bar(total);
        title('Total PP','FontWeight','bold'); ylabel('TgC'); xlabel('monthly (1-12) and annual production');

        subplot(1,3,3)
        figure;
        bar(frac)
        title('Monthly proportions of primary production under ice >10%, 2005','FontWeight','bold');
        ylabel('Primary production fraction under ice'); xlabel('month (1-12) and annual average');
        
        % print -depsc -painters bar_fraction.eps

   % (2) M-Map plots to show spatial distributions of under ice production

        % Preamble - setting plotting coordinates
        lattotal = double(1021);
        lontotal = double(1442);

        ft='/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09P.nc';
        xx=ncread(ft,'nav_lon', [1 1], [lontotal lattotal]);
        yy=ncread(ft,'nav_lat', [1 1], [lontotal lattotal]);

        land = xx(1,1);
        xx(xx == land) = NaN; xx(xx == 0) = NaN;
        land = yy(1,1);
        yy(yy == land) = NaN; yy(yy == 0) = NaN;

        
        % Monthly under-ice productions
        
        pUnder = uM / 1000; % converting from mgC/m2/time to gC/m2/time
        pUnder(:,:,13) = uAn; % appending annual under ice to monthly under ice for plotting 
        pUnder (pUnder == 0) = NaN;
        
        pIce = iceM;
        pIce(:,:,13) = iceM(:,:,9); % appending september ice overlay to end of monthly overlays (so plot annual under-ice pp with a sept 10% ice cover contour overlay)
            
        
        mnom = char('May','Jun','Jul','Aug','Sep');

        
        figure(1); odvpal(10);
        suptitle('Monthly and annual primary production under-ice (>10%), 2005') 
        for q = 5:1:9      % plots figure 1-12 as months(i), figure 13 as annual
            subplot(2,3,(q-4))
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,pUnder(:,:,q)); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,pIce(:,:,q),[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis([0 10]); cb = colorbar; ylabel(cb,'(gC/m2/month)') % for figures 1-12 time is 'month' for 13 it's 'year'
            tsrt = sprintf('%s', mnom((q-4),:));
            title(tsrt,'FontWeight','bold')  
        end     
        
            % Appending annual subplot
            subplot(2,3,6)
            m_proj('stereographic','lat',90,'long',30,'radius',30)
            m_elev('contour',[-500 -500],'edgecolor','r'); hold on
            m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
            m_pcolor(xx,yy,pUnder(:,:,13)); shading flat; hold on;
            [cs,h]=m_contour(xx,yy,pIce(:,:,13),[0.1],'edgecolor','k','LineWidth',1); hold on % 10%  ice cover contour 
            m_coast('patch',[.7 .7 .7],'edgecolor','k');
            caxis auto; cb = colorbar; ylabel(cb,'(gC/m2/yr)') % for figures 1-12 time is 'month' for 13 it's 'year'
            title ('Annual','Fontweight','bold')   
            
            % print -depsc -painters spatial_dist_landscape.eps



end