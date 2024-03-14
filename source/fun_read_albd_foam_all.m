function [grid_albd_planet,grid_albd_surface,grid_albd_cloud]  = fun_read_albd_foam_all(str)
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
% NOTE: SOLIN(time=12, lat=40, lon=48)
varid = netcdf.inqVarID(ncid,'SOLIN');
solin(:,:,:) = netcdf.getVar(ncid,varid);
% net (TOA) solar radiation at top (incoming minus outgoing)
% NOTE: FSNT(time=12, lat=40, lon=48)
varid = netcdf.inqVarID(ncid,'FSNT');
fsnt(:,:,:) = netcdf.getVar(ncid,varid);
% downward solar flux at surface (clear sky)
% NOTE: FSDSC(time=12, lat=40, lon=48)
varid = netcdf.inqVarID(ncid,'FSDSC');
fsdsc(:,:,:) = netcdf.getVar(ncid,varid);
% downward solar flux at surface (cloudy sky)
% NOTE: FSDSC(time=12, lat=40, lon=48)
varid = netcdf.inqVarID(ncid,'FSDS');
fsds(:,:,:) = netcdf.getVar(ncid,varid);
% net solar flux at surface 
% NOTE: FSDSC(time=12, lat=40, lon=48)
varid = netcdf.inqVarID(ncid,'FSNS');
fsns(:,:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
% create annual averages
grid_solin = 0.0*solin(:,:,1);
grid_fsnt  = 0.0*fsnt(:,:,1);
grid_fsdsc = 0.0*fsdsc(:,:,1);
grid_fsds  = 0.0*fsds(:,:,1);
grid_fsns  = 0.0*fsns(:,:,1);
for t=1:12
    grid_solin = grid_solin+ solin(:,:,t)/12.0;
    grid_fsnt  = grid_fsnt + fsnt(:,:,t)/12.0;
    grid_fsdsc  = grid_fsdsc + fsdsc(:,:,t)/12.0;
    grid_fsds  = grid_fsds + fsds(:,:,t)/12.0;
    grid_fsns  = grid_fsns + fsns(:,:,t)/12.0;
end
% remove zeros ...
grid_solin(find(grid_solin<=0.0)) = NaN;
grid_fsnt(find(grid_fsnt<=0.0))   = NaN;
grid_fsdsc(find(grid_fsdsc<=0.0)) = NaN;
grid_fsds(find(grid_fsds<=0.0))   = NaN;
grid_fsns(find(grid_fsns<=0.0))   = NaN;
%
% calculate (planetary) albedo
% NOTE: reflected outgoing = (grid_solin-grid_fsnt)
grid_albd_planet = (grid_solin-grid_fsnt)./grid_solin;
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
grid_albd_planet = flipud(grid_albd_planet');
%
% calculate (atmospheric) albedo
% NOTE: (FSDSC-FSDS)/SOLIN
% NOTE: difference btween 'clear-sky' and 'cloud sky' downward shortwave
grid_albd_cloud = (grid_fsdsc-grid_fsds)./grid_solin;
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
grid_albd_cloud = flipud(grid_albd_cloud');
%
% calculate (surface) albedo
% NOTE: (FSDS-FSNS)/FSDS
grid_albd_surface = (grid_fsds-grid_fsns)./grid_fsds;
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
grid_albd_surface = flipud(grid_albd_surface');
%
%

%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
