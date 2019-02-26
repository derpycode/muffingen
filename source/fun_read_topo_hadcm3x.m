function [topo] = fun_read_topo_hadcm3x(str)
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
% str(4).nc == par_nc_mask_name
% str(5).nc == par_nc_coupl_name
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(2).nc '.nc'],'nowrite');
% load TOPOGRAPHY
varid  = netcdf.inqVarID(ncid,'depthdepth');
topo(:,:) = netcdf.getVar(ncid,varid);
topo = double(topo);
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
%       => flip up-down & transpose to give a readable ASCII array
topo = flipud(topo');
% set bathymetry = negative topography
topo = -topo;
% filter out land (mark as invalid)
topo(find(topo==99999))=NaN;
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
