function [grid_albd]  = fun_read_albd_rockee(gilone,gilate,str)
%
%%

% *********************************************************************** %
% *** READ AND RETURN ALBEDO ******************************************** %
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
% set up monthly netCDF stings
str_month = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'];
% set up albedo array
grid_albd = zeros(length(gilone),length(gilate));
% LOOP
for t=1:12
    % open netCDF file
    ncid = netcdf.open([str(1).path '/' str_month(t,:) str(1).exp '.' str(3).nc '.nc'],'nowrite');
    % planetary albedo -- units of %
    varid  = netcdf.inqVarID(ncid,'plan_alb');
    albd(:,:) = netcdf.getVar(ncid,varid);
    % close netCDF file
    netcdf.close(ncid);
    % create annual averages (and correct for % albedo units)
    grid_albd  = grid_albd  + 0.01*albd(:,:)/12.0;
end
% return albedo
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
grid_albd = flipud(grid_albd');
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
