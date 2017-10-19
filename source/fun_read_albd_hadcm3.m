function [grid_albd]  = fun_read_albd_hadcm3(str)
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
% ALBEDO
% NOTE: albedop(time=18, lat=73, lon=96)
varid  = netcdf.inqVarID(ncid,'albedop');
albd(:,:,:) = netcdf.getVar(ncid,varid);
grid_albd = double(albd(:,:,13));
% NOTE: format needed is: [LAT,LON] (rows x columns) orientation
%       => flip up-down & transpose to give a readable ASCII array
% NOTE: for some reason this field does not need flipud-ing ... ??
grid_albd = grid_albd';
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
