function [axisloncm,axislonce,axislatcm,axislatce,axislonam,axislonae,axislatam,axislatae] = fun_read_axes_rockee(str)
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
% str(6).nc == par_nc_ocean_name
%
% NOTE: lon grid starts at 0
%
% *** OCEAN GRID ******************************************************** %
%
% open netCDF file
% NOTE: read axes data from ocean file
%       (both centre and edge vectors available)
ncid = netcdf.open([str(1).path '/' str(1).nc '.nc'],'nowrite');
% latitude
% read latitude mid-point axis -- c
varid     = netcdf.inqVarID(ncid,'lato');
axislatcm = netcdf.getVar(ncid,varid);
% read latitude edge axis - c
varid     = netcdf.inqVarID(ncid,'lato2');
axislatce = netcdf.getVar(ncid,varid);
% add northern-most edge
axislatce = [axislatce; 90.0];
% longitude
% read longitude mid-point axis -- c
varid  = netcdf.inqVarID(ncid,'lono');
axisloncm = netcdf.getVar(ncid,varid);
% read longitude edge axis - c
varid     = netcdf.inqVarID(ncid,'lono2');
axislonce = netcdf.getVar(ncid,varid);
% add western-most edge
axislonce = [-180.0; axislonce];
% close netCDF file
netcdf.close(ncid);
%
% *** ATMOSPHERE GRID *************************************************** %
%
axislonam = axisloncm;
axislonae = axislonce;
axislatam = axislatcm;
axislatae = axislatce;
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %