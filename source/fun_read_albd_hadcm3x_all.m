function [gi_albd_planet gi_albd_surface gi_albd_cloud gi_albd_cloud_sw]  = fun_read_albd_hadcm3x_all(str)
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
str_path = [str(1).path '/' str(1).exp];
if ~(exist([str_path '/' str(3).nc '.nc'],'file') == 2),
    if (exist([str_path '/' str(3).nc '.zip'],'file') == 2),
        unzip([str_path '/' str(3).nc '.zip'],str_path);
        fprintf('       - Extracted zipped netCDF file.\n')
    else
        disp(['       * ERROR: file ' [str_path '/' str(3).nc '.nc'] ' cannot be found (nor a zipped version)']);
        disp(['--------------------------------------------------------']);
        disp([' ']);
        diary off;
        return;
    end
end
% open netCDF file
ncid = netcdf.open([str_path '/' str(3).nc '.nc'],'nowrite');
% load and process data
% NOTE: albedop(time=18, lat=73, lon=96)
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
%       => flip up-down & transpose to give a readable ASCII array
% NOTE: for some reason this field does not need flipud-ing ... ??
% NOTE: to derive alpha(c) from alpha(p) and alpha(s), assume:
%       alpha(p) = alpha(c) + (1 - alpha(c))*alpha(s)
%       => alpha(c)/(1 - alpha(c)) = alpha(p)/alpha(s)
%       => (1 - alpha(c))/alpha(c) = alpha(s)/alpha(p)
%       => 1/alpha(c) - 1 = alpha(s)/alpha(p)
%       => alpha(c) = 1 / (1 + alpha(s)/alpha(p))
%
% ALBEDO -- PLANETARY
varid  = netcdf.inqVarID(ncid,'albedop');
albd(:,:,:) = netcdf.getVar(ncid,varid);
gi_albd_planet = double(albd(:,:,13))';
% ALBEDO -- SURFACE
varid  = netcdf.inqVarID(ncid,'albedos');
albd(:,:,:) = netcdf.getVar(ncid,varid);
gi_albd_surface = double(albd(:,:,13))';
% ALBEDO -- CLOUD (diagnosed from above)
gi_albd_cloud = 1 ./ (1.0 + gi_albd_surface./gi_albd_planet);
%
% NOTE: to derive cloud albedo from SW fluxes:
%       alpha(c) = 1.0 - SWdown(sur)/SWdown(TOA)
% SW (net solar) -- TOA down
varid  = netcdf.inqVarID(ncid,'netsolar');
sw(:,:,:) = netcdf.getVar(ncid,varid);
gi_sw_toa = double(sw(:,:,13))';
% SW -- surface down
varid  = netcdf.inqVarID(ncid,'albedos');
sw(:,:,:) = netcdf.getVar(ncid,varid);
gi_sw_surface = double(sw(:,:,13))';
% ALBEDO -- CLOUD (diagnosed from above)
gi_albd_cloud_sw = 1.0 - gi_sw_surface./gi_sw_toa;
%
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
