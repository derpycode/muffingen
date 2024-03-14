function [grid_orog] = fun_read_orog_hadcm3x(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN OROGRAPHY ***************************************** %
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
% str(6).nc == par_nc_ocean_name
% str(7).nc == par_nc_biome_name
% str(8).nc == par_nc_orog_name
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(8).nc '.nc'],'nowrite');
% load TOPOGRAPHY
varid  = netcdf.inqVarID(ncid,'ht');
grid_orog(:,:) = netcdf.getVar(ncid,varid);
grid_orog = double(grid_orog);
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
%       => flip up-down & transpose to give a readable ASCII array
% NOTE: for some reason this field does not need flipud-ing ... ??
grid_orog = grid_orog';
% filter out land (mark as invalid)
% NOTE: assume that zero is ocean, land-sea corrections made in muffingen.m
grid_orog(find(grid_orog==0.0))=NaN;
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
