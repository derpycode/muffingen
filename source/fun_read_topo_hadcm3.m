function [topo,mask]  = fun_read_topo_hadcm3(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN TOPO ********************************************** %
% *********************************************************************** %
%
% open netCDF file and load variables
% NOTE: add legacy / traceabile variable reading
%
% str input KEY:
% str(1).nc == par_nc_axes_name
% str(2).nc == par_nc_topo_name
% str(3).nc == par_nc_atmos_name
% str(4).nc == par_nc_ocean_name
% str(5).nc == par_nc_coupl_name
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(2).nc '.nc'],'nowrite');
% load TOPOGRAPHY
if strcmp(str(1).exp,'xbomg'),
    varid  = netcdf.inqVarID(ncid,'bathsmooth');
else
    varid  = netcdf.inqVarID(ncid,'depthdepth');
end
topo(:,:) = netcdf.getVar(ncid,varid);
topo = double(topo);
if strcmp(str(1).exp,'xbomg'),
    topo = topo';
else
    % NOTE: format needed is: [LAT,LON] (rows x columns) orientation
    %       => flip up-down & transpose to give a readable ASCII array
    topo = flipud(topo');
end
% set bathymetry = negative topography
topo = -topo;
% filter out land (mark as invalid)
topo(find(topo==99999))=NaN;
% load OCEAN MASK
if strcmp(str(1).exp,'xbomg'),
    varid  = netcdf.inqVarID(ncid,'depthmask');
else
    varid  = netcdf.inqVarID(ncid,'lsm');
end
mask(:,:) = netcdf.getVar(ncid,varid);
mask = double(mask);
if strcmp(str(1).exp,'xbomg'),
    mask = mask';
    % convert depth mask to mask mask
    mask(find(mask>0)) = 1;
else
    % NOTE: format needed is: [LAT,LON] (rows x columns) orientation
    %       => flip up-down & transpose to give a readable ASCII array
    mask = flipud(mask');
    % convert land to ocean mask
    mask = 1.0 - mask;
end
% close netCDF file
netcdf.close(ncid);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
