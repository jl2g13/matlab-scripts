%{
Vertical Nitrate Fluxes in the Arctic

Started by jl2g13 on 6/8/2014 - basic outline plan of script
Continued by jl2g13 on 8/9/2014 - started work on gridding nitrate data


Assume no vertical eddy fluxes due to baroclinic instability of boundary
current, instead eddies from surface fronts (subduction of isopycnals;
Timmermans2008).

Useful initial plots
- number of nitrate data per 1x1° grid cell - grid nitrate then plot number
of data as p_colour for each grid cell
- depth of nitracline for each month, depth of mixed layer for each month
(cf Tagliabue 1c)


========================================================================
READING IN OBSERVATIONAL FIELDS
========================================================================

===================
NITRATE
===================
%}

fid = fopen('/noc/msm/scratch/ariane/jl2g13/ArcticNuts/Codispoti_Arctic_Nutrients_Submission_11-11-2010.csv');
nCol = 17; % number of columns in csv file, counted in textwrangler
titles = textscan(fid,'%s',nCol,'delimiter',','); % all strings
dataFormats = ['%s %d %s %s ' repmat('%f ',1,13)];
data = textscan(fid,dataFormats,'delimiter',',','CollectOutput',1); % collecting data into cell array


date = datevec(data{3}(:,2));
lat = data{4}(:,5);
lon = data{4}(:,6);
NO3 = data{4}(:,11); % 12th column in 3rd cell of 'data' cell array (bigN = nitrate + nitrite, following Codispotti2013).
depth = data{4}(:,7);
stations = data{2}; % stations field 1 longer than all others - not sure why
year = data{4}(:,1);
% SHOULD USE COL bigN = 12 (not NO3), JUST USING NO3 AT BUILD TIME BECAUSE
% CONTAINS MORE DATA SO EASIER TO TEST SCRIPT BUILD


%{
Interpolating observations to prescribed grid positions
---------------

Interpolating the position data in the observation file (lat,lon,depth) to
values that the final grid will have: 1x1° and 15 depths (for depths see
z_grid vector initialisation).

%}

latRounded = round(lat); % lat to nearest degree
lonRounded = round(lon); % lon to nearest degree

zRoundTargets = [0 10 20 30 40 50 60 70 80 90 100 125 150 175 200]; % depths to interpolate observations onto (selected by looking at obs depth sampling regime)
depthRounded = interp1(zRoundTargets,zRoundTargets,depth,'nearest'); % 'nearest' means measured depth assigned to nearest z_roundTarget
% plot(depth, depthRounded); % deviation from y=x indicates how large an error this interpolation introduces

nitrate = [latRounded lonRounded depthRounded NO3]; % round pos to nearest degree
nitrate (nitrate == -999) = NaN; % changing missing data flag to NaN



%{
for each lattitude and longditude pair (i.e. for each grid square),
search through the entire nitrate database; for each data point that
matches the given lat-lon pair of that iteration, assign it to nitrate
grid.Grid squares with no data will remain NaN
%}

% t1 = (-180 + (1/2)):1:(180 - (1/2)); % 1° grid
% t2 = (-90 + (1/2)):1:(90 - (1/2));
% [lonnew, latnew] = meshgrid(t1, t2);
t1 = (-180):1:(180); % 1° grid
t2 = (-90):1:(90);
[lonnew, latnew] = meshgrid(t1, t2);
% NOTE: REMOVED HALVES - NEED TO BE CAREFUL OF COLLECTION INTO BINS BELOW

zGridDim = size(zRoundTargets,2); % extract depth grid resolution
nitrateGrid = NaN(181,361,zGridDim,500); % creating empty nitrate grid 1x1° for plotting. 22 depths, 10 m spacing (nitrate obs = 0-209)
                                           % although don't know size of (nitrateGrid,4) loop below much faster if preallocate 
                                           

for x = 1:size(lonnew,2) % lonnew = 180 x 360
    for y = 1:size(latnew,1) % latnew = 180 x 360
        pos = 1;
        gridloc = [lonnew(1,x) , latnew(y,1)]; % location of grid point for this iteration
        
        % for each x-y position, cycle through the entire nitrate dataset
        % looking for matches in position
            
        for i = 1:size(nitrate,1)
            if isnan(nitrate(i,4)) == 1, continue; end % if nitrate measurement is missing for this location no need to assign it (reassigns a NaN)
            
            dataPoint = nitrate(i,:); 
            
            Nlat = dataPoint(1); % nitrate data lat
            Nlon = dataPoint(2); % nitrate data lon
            
            if (Nlon == gridloc(1)) && (Nlat == gridloc(2))
              z  = dataPoint(3);
              zPos = find(ismember(zRoundTargets,z)); % ismember: searches through the depth grid 'z_roundTarges' to find elements (only one) that match(es) the depth 'z' of the data point (for this iteration). Find then gives the index (position) of this element in the depth-grid vector
              
              if (isnan(nitrateGrid(y,x,zPos,pos)) == 0); % if array position is not NaN then it's not empty (data in it), so create a new array to store this measurement in
                    pos = pos + 1; % pos only gets iterated if a nitrate measurement already exists at that position (x,y,z).
              end
              
              %{
              note that at this stage (8/9/2014) pos is most likely always
              'yr' so averaging along this dimension would give fluxes
              averaged over the entire dataset (1957-present). (Else filter
              the dataset immediately after it's read in, taking just
              1998-2005; this is better if enough data.)
              %}
              
              nitrateGrid(y,x,zPos,pos) = dataPoint(4); % dataPoint(4) = NO3
              
%               fprintf('Nitrate data %d here\n', i);
              
              
%             else
%                    fprintf('Nitrate data %d from elsewhere\n', i);
            end
        end
    end
    fprintf('Assigned all lattitudes for lon %d\n', x);
end
%{
test loop:
a = nitrate(:,4);
b = ~isnan(a);     
c = find(b,3,'first'); % first 3 elements that aren't NaN = 375,376,377
nitrate = nitrate(375:377,:); % produce a 3x4 test matrix
y = 171; x = 186; % correspond to lat,lon of the nitate data in the test matrix (gridloc = (5,80))

if (a ~= NaN)
fprintf('yes\n')
end

%}

save('/vertical_nitrate_fluxes/nitrate_data_read_in')

% calculate and plot the nitrate data coverage: number of measurements per 1x1° grid square
dataLocations = ~isnan(nitrateGrid);
dataCoverage = sum(sum(dataLocations,4),3);

dataCoverage (dataCoverage == 0) = NaN;

figure(1);
m_proj('stereo','lat',90,'rad',30,'lon',-60) % lon=-60 sets orientation of projection (60°W is bottom centre point)
m_elev('contour',[],'edgecolor','r'); hold on
m_grid('xtick',6,'xticklabels',[],'tickdir','out','ytick',[70 80],'yticklabels',[],'linest','-','box','on','linestyle','--');
m_pcolor(lonnew,latnew,dataCoverage); shading flat; hold on;
% m_coast('patch',[.7 .7 .7],'edgecolor','k'); % fill land. patch specifiers = standard matlab ones
m_coast('line','linewidth',2,'color','k'); % line coastline. line specifiers = standard matlab ones
% m_gshhs_i('line','linewidth',1,'color','k') % high resolution coastline
set(gcf,'Color','w') % sets figure background to white
caxis([0 10]); cb = colorbar; ylabel(cb,'nitrate data coverage')




             




%{
===================
10M WIND SPEED
===================
NCEP reanalysis: http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.surface.html
(2.5° winds)
    ECMWF do a 1/4° resolution reanalysis since 2005.

QuickSCAT measured wind speeds ftp://ftp.ifremer.fr./ifremer/cersat/products/gridded/MWF/L3/QuikSCAT/


Produces monthly averages of 10 m u_speed and v_speed wind interpolated
onto 1x1° grid, from the period 1981-2010.
%}

% Reading in wind speed variables

% U-wind (lonxlatxt = 144x73x798 - 2.5° 1948-2014.5?)
fnameU = '/noc/msm/scratch/ariane/jl2g13/10mwind/U-wind/monthly/uwnd.mon.mean.nc';
uWind = ncread(fnameU,'uwnd');
latWind = ncread(fnameU,'lat');
lonWind = ncread(fnameU,'lon');
tWind = ncread(fnameU,'time');
% V-wind 
fnameV = '/noc/msm/scratch/ariane/jl2g13/10mwind/V-wind/monthly/vwnd.mon.mean.nc';
vWind = ncread(fnameV,'vwnd');



%{
Producing monthly averages at 1x1° resolution

    Assuming that t_wind = 798 is 1948-present (June2014), taking simple monthly means using interval 1981-2010 - should do something better than this eventually
%}

Jan1981 = (1981-1948) * 12;
Dec2010 = (2010-1948) * 12 + 11; % + 11 because want Dec (not Jan)
avU = uWind(:,:,Jan1981:Dec2010);
avV = vWind(:,:,Jan1981:Dec2010);

% ordering data by month and year
timeInterval = 30;

for yr = 1:timeInterval
    for m = 1:12
    uwndMon(:,:,m,yr) = avU(:,:,m*yr);
    vwndMon(:,:,m,yr) = avV(:,:,m*yr);
    end
end

% producing monthly means over 30 yrs
uwndAvs = squeeze(sum(uwndMon,4) ./ timeInterval);
vwndAvs = squeeze(sum(vwndMon,4) ./ timeInterval);


% create 1x1° grid
t1 = (-180 + (1/12)):1:(180 - (1/12)); % 1° grid
t2 = (-90 + (1/12)):1:(90 - (1/12));
[lonnew, latnew] = meshgrid(t1, t2);

% interpolate 2.5x2.5° reanalysis to 1x1°
for m = 1:12
    vwndInt(:,:,m) = interp2(lonWind,latWind,vwndAvs(:,:,m)',lonnew,latnew); % linear by default; note transpose of v_wind
    uwndInt(:,:,m) = interp2(lonWind,latWind,uwndAvs(:,:,m)',lonnew,latnew);
end


%{
===================
T-S
===================
PHC MLD climatology: http://psc.apl.washington.edu/nonwp_projects/PHC/Climatology.html

Data file spirals under the globe, starting at (-89.5S,0.5E) runs eastward (.5E to
359.5E). Repeated for every depth, starting at surface. See http://psc.apl.washington.edu/nonwp_projects/PHC/Data3.html for details

data depth levels: 0 10 20 30 50 75 100 125 150 200 250 300 400 500 600 700
800 900 1000 1100 1200 1300 1400 1500
%}

% Reading in monthly Temp
% ---------
Temp = zeros(360,180,24,12); % preallocating temperature array

for m = 1:12
    if m < 10, fname  = sprintf('/noc/msm/scratch/ariane/jl2g13/PHC_MLD/TempMonths/Temp0%d_p3.obj', m);
    else, fname = sprintf('/noc/msm/scratch/ariane/jl2g13/PHC_MLD/TempMonths/Temp%d_p3.obj',m);
    end
    
    fid  = fopen (fname);
        
    data = textscan(fid,'%f');
    t1 = cell2mat(data); % recast cell as matrix
    t1 = reshape(t1,360,180,24); % data file formated as 1x1° grid with 24 depth levels ('spirals' so lon before lat in reshape)
    t1 (t1 == -99) = NaN; % land mask
    
    Temp(:,:,:,m) = t1; % x,y,z,t(mon)
      
end
            
       
% Reading in monthly S
% ---------
Sal = zeros(360,180,24,12); % preallocating temperature array

for m = 1:12
    if m < 10, fname  = sprintf('/noc/msm/scratch/ariane/jl2g13/PHC_MLD/SaltMonths/Salt0%d_p3.obj', m);
    else, fname = sprintf('/noc/msm/scratch/ariane/jl2g13/PHC_MLD/SaltMonths/Salt%d_p3.obj',m);
    end
    
    fid  = fopen (fname);
        
    data = textscan(fid,'%f');
    t1 = cell2mat(data); % recast cell as matrix
    t1 = reshape(t1,360,180,24); % data file formated as 1x1° grid with 24 depth levels ('spirals' so lon before lat in reshape)
    t1 (t1 == -99) = NaN; % land mask
    
    Sal(:,:,:,m) = t1; % x,y,z,t(mon)
      
end

% figure; pcolor(permute(Sal(:,:,1,11),[2 1])); shading flat; caxis([20 40]); colorbar % SSS in November


% linear approximation equation of state for seawater - alternative is TEOS-10 (formerly IEOS-80), has ~48 constants

rho_ref = 1035; % ocean reference density (kg/m3) - SHOULD THIS BE 1035?
beta = 8e-4; % haline expansion coefficient 
alpha = 2e-4; % thermal expansion coefficient
T0 = -2; % ref temp, degrees C - NEED TO CHECK THIS IS APPROPRIATE
S0 = 30; % ref salinity, dimensionless (conductivity ratio) - need to check this is appropriate (OA phys notes = 34.5 for global ocean, 30 better for Arctic?)

rho = rho_ref * (1 - alpha(T-T0) + beta(S-S0)); % NEXT calculate density at each x,y and each depth 




% ========================================================================
% CALCULATING VERTICAL NITRATE FLUX
% ========================================================================


% ===================
% PART 1 - DIFFUSIVE FLUX
% ===================

% (1) Find depth of nitracline

% (2) Find mixed layer depth

% (3) let k = 1e-5;



% ======> bin nitrate data by month and regrid onto 1x1° grid

% need N_mld - nitrate concentration at the base of the mixed layer (for
% F_d and F_Ek)

k = 1e-5; % diffusivity (Shaw2014)
dN = N_nitracline - N_mld; % nitrate change
dz = z_nitracline - z_mld; % depth change

F_dif = k * (dN ./ dz) * grid_area; % nitrate diffusive flux. m2/s * (mmol/m3 / m) ->  (mmol/s)


% ===================
% PART 2 - EKMAN FLUX
% ===================

Cd = 1.3e-3; % the drag coefficient. Kara2007 Fig 1
p_air = 1.2; % density of air (kg/m3).
f = 1.432e-4; % coriolis parameter at 80°N - should calculate for each meridional 1° grid square

tau_x = Cd * p_air * (u_speed^2); % zonal wind stress
tau_y = Cd * p_air * (v_speed^2); % meridional wind stress


rho_ref = 1035; % ocean reference density (kg/m3). -- variable already initialised above
dx = 111e3; % zonal length scale - 1° grid so 111 km (m)
dy = 111e3; % meridional length scale - 1° grid (m)

w_ek = ((tau_y ./ (rho_ref * f)) ./ dx) - ((tau_x ./ (rho_ref * f)) ./ dy); % depth-integrated upward Ekman velocity (m2/s)


grid_area = dx * dy;

F_Ek = w_ek * N_mld * grid_area; % nitrate ekman flux. m2/s * mmol/m3 * m2 -> mmolm/s - extra 'm' because depth integrated flux


% ===================
% PART 3 - WINTER ENTRAINMENT PRECHARGING
% ===================


% ======> regrid monthly T-S data to 1x1° then calcualed mld from it by calcualting
% density profiles for each grid square then use density difference (start
% with linear equation of state) - could use salinity difference directly

% Need maximum mixed layer depth for the year for each grid cell
% Need average nitrate concentration over corresponding depth range

N_ent = ((MLD_max * N_wMLD) - N_dis) * grid_area; % mixed layer * it's average nitrate concentration - disentrainment during MLD shoaling. m * mmol/m3 * m2 -> mmol (no s because rate irrelevant - treated as single instantaneous pulse)


