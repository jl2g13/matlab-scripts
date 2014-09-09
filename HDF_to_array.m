% =========================
% EXTRACT CHL DATA FROM SINGLE HDF FILE
% =========================
%
% Read in Chl from a HDF file into a Matlab array
%
% HDF file source: http://orca.science.oregonstate.edu/1080.by.2160.8day.hdf.chl.seawifs.php
%   hdf files are gridded data in an Equidistant Cylindrical projection.
%   SeaWiFS Chl is calculated using OC4: http://oceancolor.gsfc.nasa.gov/DOCS/MSL12/master_prodlist.html/#chlor_a

% [Chl, PAR, SST data for seawifs and modis(modis and gsm) accessible at: http://www.science.oregonstate.edu/ocean.productivity/inputData.php]

%==========
% HDF PAR files are here: http://orca.science.oregonstate.edu/1080.by.2160.8day.hdf.par.seawifs.php
% HDF VGPM_PP files are here: http://orca.science.oregonstate.edu/1080.by.2160.8day.hdf.par.seawifs.php
%
% An introduction to the HDF file type is here: http://www.swa.com/meteorology/hdf/tutorial/Intro.htm
%
% In-Situ (ARCS-PP) data is here: http://psc.apl.washington.edu/cgi-bin/PPobs/PPobs.cgi?page=Data&id=0.859687454373674
%==========
%
% Original instructions followed: http://www.atmos.umd.edu/~gcm/usefuldocs/hdf_netcdf/HDF_matlab.html
% The generic code is:
% 
% fileID = hdfsd('start','filename', 'read')
% [numdatasets, numdescripters, status] = hdfsd('fileinfo', fileID)
% DescripterN = hdfsd('readattr', 'fileID', n)
% DatasetID = hdfsd('select', fileID, DataIndex)
% [name, numdim, dimvector, datatype, numdescriptors, status] = hdfsd('getinfo',DatasetID)
% DatasetVar = hdfsd('readdata', datasetID, startvector, stridevector, endvector)


fileID = hdfsd('start','chl.2000214.hdf', 'read'); % Open the HDF file and assign it a file id (computer stores a numerical file ID which can be used interchangeably with the ID assigned here)
    % File formats: yyyyddd of the file name is the first day of the 8-day average window in
    % year day format (ddd)

[numdata, numdescr] = hdfsd('fileinfo', fileID); % Extract information about the file description (returns number of datasets and number of attributes; numdescr = no of attributes)


%   EOS Chl HDF files have 5 (n.b. index from zero) attributes by trial and error using * these are:
%         attribute 0 = Chl
%         attribute 1 = Start date (Julian time)
%         attribute 2 = End date (Julian time)
%         attribute 3 = Start date (Calender)
%   attribute 4 = End date (Calender)
% 
%  * descriptor0 = hdfsd('readattr', CloudID, 0) etc.


chlid = hdfsd('select', fileID, 0); % Assign Chl dataset ID

[name,numdim,dimvector,type,numdescr] = hdfsd('getinfo', chlid); % Extract dataset information

% Specify how much of chl data to read in
startvector = [0 0] % Start reading HDF file at the beginning of Chl
endvector = dimvector % Finish reading HDF file at end of Chl
stride = [] % specifies number of data steps to read at (n=1 here, i.e. read all)


Chl = hdfsd('readdata', chlid, startvector, stride, endvector);


Chl = double(Chl);% Converting Chl date to double precision (needed to do any operations on
% it in Matlab


% Missing data flag is -9.9990

ii=find(Chl < 0); Chl(ii)=NaN; % removing negative Chl values (missing data etc)


% Surface plot of Chl as a f(latt and long) - surf(x,y,z) is surf(Chl,lat,lon)
    % Latt = 1080, Lon = 2160 (so both at 1/6 degree)
    % Chl in units mg/m3

surf(Chl);
xlabel('lattitude')
ylabel('longditude')
zlabel('Chl concentration (mg/m3)')


% =========================
% EXTRACT PAR DATA FROM SINGLE HDF FILE
% =========================
% Data source: http://orca.science.oregonstate.edu/1080.by.2160.8day.hdf.par.seawifs.php
% PAR is in E/m2/d.
% Code is same as thast implemented above for Chl


fileID = hdfsd('start','par.2007105.hdf', 'read');

[numdata, numdescr] = hdfsd('fileinfo', fileID); % 1 datafile with 4 attributes

parid = hdfsd('select', fileID, 0); % PAR (like chl) is stored in the first attribute (it has 2dim, lat and lon = 1080, 2160)

[name,numdim,dimvector,type,numdescr] = hdfsd('getinfo', parid);

startvector = [0 0]
endvector = dimvector
stride = []

PAR = hdfsd('readdata', parid, startvector, stride, endvector);

PAR = double(PAR);

ii=find(PAR < 0); PAR(ii)=NaN; % Missing data flag is -9.990 (as for Chl)

figure;
surf(PAR)
xlabel('lattitude')
ylabel('longditude')
zlabel('PAR(E/m2/d)')


% =========================
% EXTRACT CHL OR PAR DATA FROM SINGLE HDF FILE
% =========================
% Moving towards a single script

% function [array_out] = hdf_to_array('filename')
%
% 
% fileID = hdfsd('start','filename', 'read');
% 
% [numdata, numdescr] = hdfsd('fileinfo', fileID); % 1 datafile with 4 attributes
% 
% varid = hdfsd('select', fileID, 0); % variable (PAR or CHL) is stored in the first attribute (it has 2dim, lat and lon = 1080, 2160)
% 
% [name,numdim,dimvector,type,numdescr] = hdfsd('getinfo', varid);
% 
% startvector = [0 0]
% endvector = dimvector
% stride = []
% 
% var = hdfsd('readdata', varid, startvector, stride, endvector);
% 
% var = double(var);
% 
% ii=find(var < 0); var(ii)=NaN; % Missing data flag is -9.990 (as for Chl)
% 
% figure;
% surf(var)
% xlabel('lattitude')
% ylabel('longditude')
% zlabel('var (PAR or Chl(E/m2/d or mg/m3))')
