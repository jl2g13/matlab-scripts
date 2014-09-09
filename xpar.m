% Model light field

% The script reads in the model suface light field and then attenuates it
% with depth using the model surface chlorophyll field (jl2g13 13/03/14).
% Based on fortran script trcopt_medusa.F90. Model uses 

% (1) Read in model PAR and Chl
% -----

pos = 1;

lattotal = double(1021);
latstart = double(825);
lontotal = double(1442);
depth = 64;

for y = 1:1
    yr = 2004 + y;
    
    % Getting list of all files for this year
    dnomP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05P.nc', yr);
    fnamesP = dir(dnomP);
    
    dnomD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/*d05D.nc', yr);
    fnamesD = dir(dnomD);
    
    [num, x] = size(fnamesP);
    
    % Preallocating arrays for model output to be read into
    totalshortwave = zeros(lontotal,(lattotal-latstart),num);
    modchn = zeros(lontotal,(lattotal-latstart),depth,num);
    modchd = zeros(lontotal,(lattotal-latstart),depth,num);

    
    % Initiating extraction loop
    
    for f = 1:1:num
        
        % (1) For Chlorophyll
        
        fnameP = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesP(f).name);
        fprintf('Reading %s\n', fnameP);
        
        % (i) For non-diatom chlorophyll
        t1 = ncread(fnameP,'CHN',[1 latstart 1 1], [lontotal (lattotal-latstart) depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN;
        
        modchn(:,:,:,pos) = t1;
        
        % (ii) For diatom chlorophyll
        t1 = ncread(fnameP,'CHD',[1 latstart 1 1], [lontotal (lattotal-latstart) depth 1]);
        t1 = squeeze(t1);
        t1 (t1 == 0) = NaN;
        
        modchd(:,:,:,pos) = t1;
        
        
        
      % (2) For Photosynthetically Active Radiation
        
        fnameD = sprintf('/noc/altix/scratch/omfman/ORCA025-MEDUSA-N201/%d/%s', yr, fnamesD(f).name);
        fprintf('Reading %s\n', fnameD);
        
        % Reading in surface PAR model field
        t1 = ncread(fnameD,'MED_QSR',[1 latstart 1], [lontotal (lattotal-latstart) 1]);
        t1 (t1 == 0) = NaN;
        % don't squeeze because want depth dim for multiplication below
        
        totalshortwave(:,:,pos) = t1;
        
        pos = pos + 1;
        
    end
end

% Freeing up some memory
modchl = modchn + modchd;
clearvars modchn modchd t1


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

dep_nemo; % for model grid depth

zpar = totalshortwave * (0.43 / 2); % Levy2001 eq A23. Divide by two because two-band model.
zpar = reshape(zpar,1442,196,1,73); % adding depth dim back in (=1, making explicit)

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
par_J = par_E * ((h*c/length) * 6.022D17)

% parJ = parE * ((h*c/length) * 6.022D17);

% Example plots

% (1) Par at surface, 10, 20, 46 and 100 m in summer
figure(1);
subplot(2,2,1); pcolor(par_J(:,:,3,35)); shading flat; cb = colorbar; title('par_J at 10 m')
subplot(2,2,2); pcolor(par_J(:,:,4,35)); shading flat; cb = colorbar; title('par_J at 20 m')
subplot(2,2,3); pcolor(par_J(:,:,8,35)); shading flat; cb = colorbar; title('par_J at 46 m')
subplot(2,2,4); pcolor(par_J(:,:,14,35)); shading flat; cb = colorbar; title('par_J at 100 m')

figure(2);
subplot(2,2,1); pcolor(par_E(:,:,3,35)); shading flat; cb = colorbar; title('par_E at 10 m')
subplot(2,2,2); pcolor(par_E(:,:,4,35)); shading flat; cb = colorbar; title('par_E at 20 m')
subplot(2,2,3); pcolor(par_E(:,:,8,35)); shading flat; cb = colorbar; title('par_E at 46 m')
subplot(2,2,4); pcolor(par_E(:,:,14,35)); shading flat; cb = colorbar; title('par_E at 100 m')

% individual figures:
    figure(3);
    pcolor(par_J(:,:,1,35)); shading flat; cb = colorbar; title('summer PAR directly below surface'); ylabel(cb, 'par_J (W/m2)')
    pcolor(par_J(:,:,3,35)); shading flat; cb = colorbar; title('par_J at 10 m')
    pcolor(par_J(:,:,8,35)); shading flat; cb = colorbar; title('par_J at 46 m'); caxis([0 0.2])
    pcolor(par_J(:,:,14,35)); shading flat; cb = colorbar; title('par_J at 100 m'); caxis([0 0.01])
    
    

% (2) Par(z) for an arbirary location (170W,77N, Chukchi Sea shelf break)

% NB at the moment I think the y axis is plotted using box edges (z=65) and
% profile is in the centre of boxes (z=64).

figure(2);
profile_J = squeeze(par_J(391,(910-latstart),:,35));
plot(profile(1:10),nemo_dep(1:10)); xlabel('par_J (W/m2)'); ylabel('Depth (m)'); title ('par_J at 77N,177W','FontWeight','bold'); set(gca,'YDir','reverse');

% profile_E = squeeze(par_E(391,(910-latstart),:,35));
% plot(profile_E(1:10),nemo_dep(1:10)); xlabel('par_E (umolphotons/m2/s)'); ylabel('Depth (m)'); title ('par_E at 77N,177W','FontWeight','bold'); set(gca,'YDir','reverse');