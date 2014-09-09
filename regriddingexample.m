clear; startup;

% ======================================================================
% SETUP GRID
% ======================================================================
% select file to work from
% inom = char('MSA003/1996/ORCA025-M003_1996m01D.nc'); 
inom = char('/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m09D.nc');

% read in longitude and latitude
% lon025 = nc_varget(inom, 'nav_lon'); % obviously I'm using an older version of Matlab
% lat025 = nc_varget(inom, 'nav_lat');
lon025 = ncread(inom, 'nav_lon'); % obviously I'm using an older version of Matlab
lat025 = ncread(inom, 'nav_lat');

% load in PAR
% t1 = nc_varget(inom, 'MED_QSR');
t1 = ncread(inom, 'MED_QSR');
land = t1(1,1); % set land mask
t1(t1 == land) = NaN; 
par = t1 * 0.43;

% create new grid
t1 = (-180 + (1/12)):(1/6):(180 - (1/12)); % 1/6° grid
t2 = (-90 + (1/12)):(1/6):(90 - (1/12));
[lonnew, latnew] = meshgrid(t1, t2);

% ======================================================================
% LOAD VARIABLE
% ======================================================================
% load up PAR 
% par = zeros([1021 1442 12]);
% for m = 1:1:12
%   if m < 10, inom = sprintf('/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m0%dD.nc', m);
%   else, inom = sprintf('/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m%dD.nc', m); end
%   t1 = nc_varget(inom, 'MED_QSR');
%   land = t1(1,1);
%   t1(t1 == land) = NaN;
%   par(:,:,m) = t1 * 0.43;
% end

for m = 1:1:12
  if m < 10, inom = sprintf('/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m0%dD.nc', m);
    else inom = sprintf('/noc/altix2/scratch/omfman/ORCA025-N201/means/1990/ORCA025-N201_1990m%dD.nc', m);
  end
  t1 = ncread(inom, 'MED_QSR');
  land = t1(1,1);
  t1(t1 == land) = NaN;
  par(:,:,m) = t1 * 0.43;
end




% some NaN values in this data are actually zero - find these and replace
% them with zeros
t1 = isfinite(par);
t2 = sum(t1, 3);
mask = t2 * 0;
mask(t2 == 0) = NaN;

for m = 1:1:12
  t1 = par(:,:,m);
  t1(isnan(t1)) = 0;
  par(:,:,m) = t1 + mask;
end



% ======================================================================
% REGRID VARIABLE
% ======================================================================
% convert PAR to this new grid
%
% the conversion is always the same because the grid is the same, so one
% can save time by first calculating how to do the conversion once, and
% then applying the result of this to all the fields; this saves the cost
% of repeatedly recalculating how to do the conversion; obviously this only
% works if the grid is always the same; note: a large negative value is
% used over land so that it can be easily identified after regridding - if
% this isn't done, tinterp will produce PAR values over land (not
% necessarily a bad thing - but disconcerting!)
%
% part 1 - determine the conversion calculation
xx = lon025; yy = lat025;
xx3 = [(xx - 360) xx (xx + 360)];
yy3 = [yy yy yy];

t0 = par(:,:,1);
t0(isnan(t0)) = -1e3; t1 = [t0 t0 t0];
t2 = isfinite(t1);

get_loc = t2;
px = xx3(t2); py = yy3(t2); p = double([px py]);
t = delaunay(p); % delaunayn may only work this way in Matlab v2009a
% delaunay: returns the vertices as numerical output. (for points in 4D)
% DelaunayTri returns a triangulation object (for points in 2-3D)

% i = PointLocation (t, p); 


% part 2 - perform this on all PAR fields
parnew = zeros([1080 2160 12]);
for m = 1:1:12
  t0 = par(:,:,m); t0(isnan(t0)) = -1e6;
  t1 = [t0 t0 t0]; vals = t1(get_loc);
  % t = PointLocation (t, p); % alternative to tsearch/tinterp
  t2 = tinterp(p, t, vals, lonnew, latnew); % tinterp is not a standard Matlab command (already copied Andrew's function, may not work in v2012)
  % see tinterp fix for later versions on matlab page for: http://www.mathworks.co.uk/matlabcentral/fileexchange/13183-tinterp-an-alternative-to-griddata
  t2(t2 < 0) = NaN;
  parnew(:,:,m) = t2;
  fprintf('- done month %d\n', m);
end
