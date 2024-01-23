function [gi_albd_planet,gi_albd_surface,gi_albd_cloud]  = fun_read_albd_hadcm3x_all(str)
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
%
% ALBEDO -- PLANETARY
% NOTE: albedop = upSol_mm_s3_TOA./downSol_mm_TOA;
varid  = netcdf.inqVarID(ncid,'albedop');
albd(:,:,:) = netcdf.getVar(ncid,varid);
gi_albd_planet = double(albd(:,:,13))';
% ALBEDO -- SURFACE
% NOTE: albedos = (downSol_Seaice_mm_s3_srf-solar_mm_s3_srf)./downSol_Seaice_mm_s3_srf;
varid  = netcdf.inqVarID(ncid,'albedos');
albd(:,:,:) = netcdf.getVar(ncid,varid);
gi_albd_surface = double(albd(:,:,13))';
% close netCDF file
netcdf.close(ncid);
%
% ALBEDO -- CLOUD 
% open netCDF files
ncid = netcdf.open([str_path '/' str(1).nc '.nc'],'nowrite');
% load and process data
% NOTE: VAR(lat=73, lon=96) ==> annual means
% NOTE: no hadCM3 output -- needs manual calculation
varid  = netcdf.inqVarID(ncid,'downSol_Seaice_mm_s3_srf');
downSol_srf(:,:) = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'upSol_mm_s3_TOA');
upSol_TOA(:,:) = netcdf.getVar(ncid,varid);
varid  = netcdf.inqVarID(ncid,'downSol_mm_TOA');
downSol_TOA(:,:) = netcdf.getVar(ncid,varid);
% close netCDF file
netcdf.close(ncid);
% transpose
downSol_srf = downSol_srf';
upSol_TOA = upSol_TOA';
downSol_TOA = downSol_TOA';
%
% calculate upSol_srf (reflected from surface)
upSol_srf = gi_albd_surface.*downSol_srf;
% calculate diff between TOA outgoing solar and surface-reflected 
% outgoing solar to obtain cloud-reflected fraction of outgoing solar energy
upSol_Cloud = upSol_TOA - upSol_srf;
% calculate cloud albedo and 
gi_albd_cloud = (upSol_Cloud./downSol_TOA);
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
