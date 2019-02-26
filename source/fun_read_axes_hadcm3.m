function [axisloncm,axislonce,axislatcm,axislatce] = fun_read_axes_hadcm3(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN AXES INFO ***************************************** %
% *********************************************************************** %
%
% str input KEY:
% str(1).nc == par_nc_axes_name
% str(2).nc == par_nc_topo_name
% str(3).nc == par_nc_atmos_name
% str(4).nc == par_nc_mask_name
% str(5).nc == par_nc_coupl_name
%
% NOTE: lon grid starts at 0
%
% *** OCEAN GRID ******************************************************** %
%
% open netCDF file
% NOTE: read axes data from pdclann file
%       (complete grid info not present in topo file)
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(2).nc '.nc'],'nowrite');
% latitude
% read latitude mid-point axis -- c
varid     = netcdf.inqVarID(ncid,'latitude');
axislatcm = netcdf.getVar(ncid,varid);
%%%axislatcm = flipud(axislatcm);
% create latitude edge axis - c
axislatce = [axislatcm-180/144/2; 90.0];
% longitude
% read longitude mid-point axis -- c
varid  = netcdf.inqVarID(ncid,'longitude');
axisloncm = netcdf.getVar(ncid,varid);
% create longitude edge axis - c
axislonce = [axisloncm-360/288/2; 360.0-360/288/2];
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

