function [topo,mask]  = fun_read_topo_foam(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN TOPO ********************************************** %
% *********************************************************************** %
%
% open netCDF file and load variables
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
% load topo from netCDF
varid  = netcdf.inqVarID(ncid,'TOPO');
topo(:,:) = netcdf.getVar(ncid,varid);
topo = double(topo);
% flip array around diagonal to give (lon,lat) array orientation
topo = flipud(topo');
% set bathymetry = negative topography
topo = -topo;
% filter out land (mark as invalid)
topo(find(topo<=0.0))   = NaN;
% derive mask
mask = topo;
mask(find(isnan(mask))) = 0;
mask(find(mask>0.0))    = 1;
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
