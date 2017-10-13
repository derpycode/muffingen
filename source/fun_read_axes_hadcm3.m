function [axisloncm,axislonce,axislatcm,axislatce,axislonpm,axislonpe,axislatpm,axislatpe,axislonam,axislonae,axislatam,axislatae]  = fun_read_axes_hadcm3(str)
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
% NOTE: the UM uses:
%       c-grid for tracers, fluxes, etc.
%       p-grid for wind products
%
%  0.0     3.75    7.5
%
%   c*******c*******c                   90.0
%   *   |       |
%   *---t-------t----                   88.75
%   *   |       |
%   c   |   c   |   c                   87.5
%   *   |       |
%   *---t-------t----                   86.25
%   *   |       |
%   c   |   c   |   c                   85.0
%   *   |       |               |   *
%
%   *   |       |               |   *
%   c   |   c   |   c       c   |   *
%   *   |       |               |   *
%   *---p-------p------   ------p---*
%   *   |       |               |   *
%   c   |   c   |   c       c   |   *
%   *   |       |               |   *
%   *---p-------p------   ------p---*
%   *   |       |               |   *
%   c*******c*******c**   **c********
%
% NOTE: in the UM:
%       latitude    == c grid ( +90.0   to  -90.0  ) -- n=73
%       latitude_1  == p grid ( +88.75  to  -88.75 ) -- n=72
%       longitude   == c grid (  +0.0   to +356.25 ) -- n=96
%       longitude_1 == p grid (  +1.875 to +358.125) -- n=96
% NOTE: lon grid starts at 0
%
% *** OCEAN GRID ******************************************************** %
%
% open netCDF file
% NOTE: read axes data from pdclann file
%       (complete grid info not present in topo file)
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(1).nc '.nc'],'nowrite');
% latitude
% read latitude mid-point axis -- c
varid     = netcdf.inqVarID(ncid,'latitude');
axislatcm = netcdf.getVar(ncid,varid);
axislatcm = flipud(axislatcm);
% read latitude mid-point axis -- p
varid  = netcdf.inqVarID(ncid,'latitude_1');
axislatpm = netcdf.getVar(ncid,varid);
axislatpm = flipud(axislatpm);
% create latitude edge axis - c
axislatce = [-90.0; axislatpm; 90.0];
% create latitude edge axis - p
axislatpe = [axislatcm];
% longitude
% read longitude mid-point axis -- c
varid  = netcdf.inqVarID(ncid,'longitude');
axisloncm = netcdf.getVar(ncid,varid);
% read longitude mid-point axis -- p
varid  = netcdf.inqVarID(ncid,'longitude_1');
axislonpm = netcdf.getVar(ncid,varid);
% create longitude edge axis - c
axislonce = [-1.875; axislonpm];
% create longitude edge axis - p
axislonpe = [axisloncm; 360.0];
% close netCDF file
netcdf.close(ncid);
%
% *** ATMOSPHERE GRID *************************************************** %
%
% set a grid == c grid
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

