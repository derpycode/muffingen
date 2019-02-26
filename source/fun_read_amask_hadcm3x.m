function [mask]  = fun_read_amask_hadcm3x(str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN ATM MASK ****************************************** %
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
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(4).nc '.nc'],'nowrite');
% load OCEAN MASK
varid  = netcdf.inqVarID(ncid,'lsm');
mask(:,:) = netcdf.getVar(ncid,varid);
mask = double(mask);
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
%       => transpose to give a readable ASCII array
% NOTE: atmosphere mask is inverted as compared to ocean mask ... !!
%       (so no flip up-down)
mask = mask';
% convert land to ocean mask
mask = 1.0 - mask;
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
