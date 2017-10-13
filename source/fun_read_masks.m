function [maskt,maskv]  = fun_read_masks(par_path,par_expid,par_gcm)
%
%%

% *********************************************************************** %
% *** READ AND RETURN MASKS ********************************************* %
% *********************************************************************** %
%
% BLAH
switch par_gcm
    case 'hadcm3l'
% open netCDF file
        ncid = netcdf.open([par_path '/' par_expid '/' par_expid '.qrparm.omask.nc'],'nowrite');
        % load "Mask for Ocean Temperature" (73x96)
        varid  = netcdf.inqVarID(ncid,'maskt');
        mask_tmp = netcdf.getVar(ncid,varid);
        maskt = double(mask_tmp);
        maskt(:,:,2:end) = [];
        maskt(find(maskt(:,:) == -99999.0)) = NaN;
        maskt = flipdim(maskt',1);
        % load "Mask for Ocean Velocity" (72x96)
        varid  = netcdf.inqVarID(ncid,'masku');
        mask_tmp = netcdf.getVar(ncid,varid);
        maskv = double(mask_tmp);
        maskv(:,:,2:end) = [];
        maskv(find(maskv(:,:) == -99999.0)) = NaN;
        maskv = flipdim(maskv',1);
    case 'foam'
% open netCDF file
        ncid = netcdf.open([par_path '/' par_expid '.topo.nc'],'nowrite');
        %%%
end
% close netCDF file
netcdf.close(ncid);
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
