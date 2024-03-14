function [grid_orog] = fun_read_orog_foam(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN OROGRAPHY ***************************************** %
% *********************************************************************** %
%
% open netCDF file and load variables
%
% str input KEY:
% str(1).nc == par_nc_axes_name ('topo')
% str(2).nc == par_nc_topo_name ('topo')
% str(3).nc == par_nc_atmos_name ('atmos')
% str(4).nc == par_nc_mask_name ('topo')
% str(5).nc == par_nc_coupl_name
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(2).nc '.nc'],'nowrite');
% load TOPOGRAPHY
varid  = netcdf.inqVarID(ncid,'TOPO');
grid_orog(:,:) = netcdf.getVar(ncid,varid);
grid_orog = double(grid_orog);
% flip array around diagonal to give (lon,lat) array orientation
grid_orog = flipud(grid_orog');
% filter out ocean (mark as invalid)
% NOTE: assume neg topo is ocean, land-sea corrections made in muffingen.m
grid_orog(find(grid_orog<0.0))=NaN;
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
