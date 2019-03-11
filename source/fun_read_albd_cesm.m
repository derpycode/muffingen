function [grid_albd]  = fun_read_albd_cesm(str)
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
% str(4).nc == par_nc_ocean_name
% str(5).nc == par_nc_coupl_name
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(3).nc '.nc'],'nowrite');
% incoming (TOA) solar radiation
% NOTE: SOLIN(time=12, lat=96, lon=144)
varid  = netcdf.inqVarID(ncid,'SOLIN');
solin(:,:,:) = netcdf.getVar(ncid,varid);
% outgoing solar radiation (TOA) 
% NOTE: FSUTOA(time=12, lat=96, lon=144)
varid  = netcdf.inqVarID(ncid,'FSUTOA');
fsutoa(:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
% create annual averages
grid_solin  = 0.0*solin(:,:,1);
grid_fsutoa = 0.0*fsutoa(:,:,1);
for t=1:12
    grid_solin  = grid_solin   + solin(:,:,t)/12.0;
    grid_fsutoa = grid_fsutoa  + fsutoa(:,:,t)/12.0;
end
% remove zeros ...
grid_solin(find(grid_solin<=0.0))   = NaN;
grid_fsutoa(find(grid_fsutoa<=0.0)) = NaN;
% calculate albedo
% NOTE: reflected outgoing = grid_fsutoa
grid_albd = grid_fsutoa./grid_solin;
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
