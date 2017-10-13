function [axisloncm,axislonce,axislatcm,axislatce,axislonpm,axislonpe,axislatpm,axislatpe,axislonam,axislonae,axislatam,axislatae]  = fun_read_axes_foam(str)
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
% str(4).nc == par_nc_ocean_name
% str(5).nc == par_nc_coupl_name
%
%  0.0     2.81    5.62
%
%               P
%
%   *****************     90.0
%   *       |
%   c-------c---p---      89.3
%   *       |
%   *       |
%   *       |
%   c-------c---p---      87.9
%   *       |
%   *       |
%   *       |
%
%   *       |               |       *
%           |               |       *
%   *       |               |       *
%   c-------c---p--       --c---p---*   P
%   *       |               |       *
%           |               |       *
%   *       |               |       *
%   c-------c------       --c-------*
%   *       |               |       *
%   *******************   **c************
%
%                        357.2   360.0
%
% *** OCEAN c GRID ****************************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(1).nc '.nc'],'nowrite');
% latitude (c grid)
varid     = netcdf.inqVarID(ncid,'MYLAT');
axislatcm = netcdf.getVar(ncid,varid);
% create latitude edge axis - c
% NOTE: ensure limits are +/- 90
axislatce = axislatcm - 180.0/length(axislatcm)/2.0;
axislatce = [axislatce; axislatce(end)+180.0/length(axislatcm)];
axislatce(1)   = -90.0;
axislatce(end) = 90.0;
% longitude (c grid)
varid     = netcdf.inqVarID(ncid,'MYLON');
axisloncm = netcdf.getVar(ncid,varid);
% create longitude edge axis - c
axislonce = axisloncm - 360.0/length(axisloncm)/2.0;
axislonce = [axislonce; axislonce(end)+360.0/length(axisloncm)];
% close netCDF file
netcdf.close(ncid);
%
% *** OCEAN p [coupled] GRID ******************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(5).nc '.nc'],'nowrite');
% latitude (p grid)
% NOTE: the lat grid is offset outside of -90 to 90
varid     = netcdf.inqVarID(ncid,'lat');
axislatpm = netcdf.getVar(ncid,varid);
% create latitude edge axis - p
% NOTE: ensure limits are +/- 90
axislatpe = axislatpm - 180.0/length(axislatpm)/2.0;
axislatpe = [axislatpe; axislatpe(end)+180.0/length(axislatpm)];
axislatpe(1)   = -90.0;
axislatpe(end) = 90.0;
% longitude (p grid)
% NOTE: the lon grid is offset outside of 0 to 360
varid     = netcdf.inqVarID(ncid,'lon');
axislonpm = netcdf.getVar(ncid,varid);
% create longitude edge axis - p
% NOTE: allow the lon grid offset
axislonpe = axislonpm - 360.0/length(axislonpm)/2.0;
axislonpe = [axislonpe; axislonpe(end)+360.0/length(axislonpm)];
% close netCDF file
netcdf.close(ncid);
%
% *** ATMOSPHERE GRID *************************************************** %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(3).nc '.nc'],'nowrite');
% latitude (c grid)
varid     = netcdf.inqVarID(ncid,'lat');
axislatam = netcdf.getVar(ncid,varid);
% create latitude edge axis - c
% NOTE: ensure limits are +/- 90
axislatae = axislatam - 180.0/length(axislatam)/2.0;
axislatae = [axislatae; axislatae(end) + 180.0/length(axislatam)];
axislatae(1)   = -90.0;
axislatae(end) = 90.0;
% longitude (c grid)
varid     = netcdf.inqVarID(ncid,'lon');
axislonam = netcdf.getVar(ncid,varid);
% create longitude edge axis - c
axislonae = axislonam - 360.0/length(axislonam)/2.0;
axislonae = [axislonae; axislonae(end) + 360.0/length(axislonam)];
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
