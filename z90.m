% kd490

% Global kd490 (4km, 8day) for summer 2013. z90 plotted for Arctic for June
% 2013 average.                       

% ft='/noc/msm/scratch/ariane/jl2g13/kd490/erdVHk4908day_c3d2_0590_6755.nc'; % June 06/06/13 to 06/07/13
% ft='/noc/msm/scratch/ariane/jl2g13/kd490/erdVHk4908day_d990_db97_0206.nc'; % July 06/07/13 to 06/08/13
% ft='/noc/msm/scratch/ariane/jl2g13/kd490/erdVHk4908day_8e14_f716_6238.nc';  % Aug 06/08/13 to 06/09/13 
 ft='/noc/msm/scratch/ariane/jl2g13/kd490/erdVHk4908day_9255_c9b1_6f9d.nc'; % Sep
% ft='/noc/msm/scratch/ariane/jl2g13/kd490/erdVHk4908day_0612_c9f0_c8d4.nc'; % Oct


kd490 = ncread(ft,'k490', [1 1 1], [8640 4320 4]); % change '5' to match len time dim // lon x lat x time (4km, 8 day = 360° x 65-90°N x 360 days)
kd490 (kd490 == -3.28e4) = NaN; % missing data fill value
z90 = 1 / kd490;

figure; % 6-10 Month global 
subplot(2,1,1); pcolor(kd490(:,:,1)); shading flat; cb = colorbar; ylabel(cb, '(/m)'); title('kd490');
subplot(2,1,2); pcolor(z90(:,:,1)); shading flat; cb = colorbar; ylabel(cb, '(m)'); title('z90');

% % Month average global
%     % retaining land mask
%     kd490av = nansum(kd490,3) / 5;
%     kd490av (kd490av == 0) = NaN;
% 
%     z90av = nansum(z90,3) / 5;
%     z90av (z90av == 0) = NaN;
% figure;
% subplot(2,1,1); pcolor(kd490av); shading flat; colorbar; title('kd490');
% subplot(2,1,2); pcolor(z90av); shading flat; colorbar; title('z90');


% Plotting Arctic

      t1 = permute(z90(:,:,1), [2 1]);
    
%     % Monthly average
%     t1 = permute(nansum(z90,3), [2 1]); % need to divide by number of files
%     t1 (t1 == 0) = NaN;

    % creating grid for M_Map
        x1 = (-180 + (1/48)):(1/24):(180 - (1/48)); % data at 1/24° ~= 4 km
        y1 = (-90 + (1/48)):(1/24):(90 - (1/48));
        [x2,y2] = meshgrid(x1, y1);


    % No data in the Arctic for the entire data file (one month - 5 8-day averages) 
        figure;
        m_proj('stereographic','lat',90,'long',30,'radius',30)
        m_elev('contour',[-500 -500],'edgecolor','r'); hold on
        m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-');
        m_pcolor(x2,y2,t1); shading flat; hold on;
        m_coast('patch',[.7 .7 .7],'edgecolor','k');
        caxis([0 60]); cb = colorbar; ylabel(cb,'Under ice PP (mgC/m2/month)')
        title ('z90,early August 2013','Fontweight','bold')   

    
    
    
    
figure;
    subplot(4,1,1)
    pcolor(z90(:,:,20)); shading flat; colorbar;
    
    subplot(4,1,2)
    pcolor(z90(:,:,25)); shading flat; colorbar;

    subplot(4,1,3)
    pcolor(z90(:,:,30)); shading flat; colorbar;
    
    subplot(4,1,4)
    pcolor(z90(:,:,35)); shading flat; colorbar;