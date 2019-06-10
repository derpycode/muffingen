function [] = make_grid_winds_rockee(gilonm,gilone,gilatm,gilate,gimask,golonm,golone,golatm,golate,gomask,str,optplots);
% make_grid_winds_rockee
%
%   *********************************************************
%   *** RE_GRID WIND PRODUCTS                             ***
%   *********************************************************
%
%
%   wstress is organised :   level 1   tau_x at u point
%                            level 2   tau_x at v point
%                            level 3   tau_y at u point
%                            level 4   tau_y at v point
%
%   wvelocity is organised : level 1   x velocity at grid point
%                            level 2   y velocity at grid point
%
%   wspeed is organised :    level 1   speed at grid point
%
%   str input KEY:
% str(1).nc == par_nc_axes_name
% str(2).nc == par_nc_topo_name
% str(3).nc == par_nc_atmos_name
% str(4).nc == par_nc_mask_name
% str(5).nc == par_nc_coupl_name
% str(6).nc == par_nc_ocean_name
%
%   *********************************************************
%
%   ***********************************************************************
%%

% *********************************************************************** %
% *** INITIALIZE ******************************************************** %
% *********************************************************************** %
%
% determine output grid size (remember: [rows columns])
[jmax, imax] = size(gomask);
% set up monthly netCDF stings
str_month = ['JAN'; 'FEB'; 'MAR'; 'APR'; 'MAY'; 'JUN'; 'JUL'; 'AUG'; 'SEP'; 'OCT'; 'NOV'; 'DEC'];
% set up averaged input arrays
giwtauu     = zeros(length(gilonm),length(gilatm));
giwtauv     = zeros(length(gilonm),length(gilatm));
giwvelu     = zeros(length(gilonm),length(gilatm));
giwvelv     = zeros(length(gilonm),length(gilatm));
giwspd_wsma = zeros(length(gilonm),length(gilatm));
%
% *********************************************************************** %

% *********************************************************************** %
% *** load input data *************************************************** %
% *********************************************************************** %
%
% NOTE: wind stress in units of: g/m*s^2
%       vs. normally:            N/m^2
%       N m-2 == kg m s-2 m-2
%             == kg m-1 s-2
%       ==> divide g/m*s^2 units by 1000.0 (to give kg/ms^2)
%       (1 N	=	1 kg	?	1 m/s2)
%
% *** load and process data -- wind stress ****************************** %
%
% LOOP
for t=1:12
    % open netCDF file
    ncid = netcdf.open([str(1).path '/' str_month(t,:) str(1).exp '.' str(3).nc '.nc'],'nowrite');
    % U COMPON OF MOMENTUM SRF DRAG
    varid  = netcdf.inqVarID(ncid,'tauus');
    loctauu(:,:) = netcdf.getVar(ncid,varid);
    % V COMPON OF MOMENTUM SRF DRAG
    varid  = netcdf.inqVarID(ncid,'tauvs');
    loctauv(:,:) = netcdf.getVar(ncid,varid);
    % close netCDF file
    netcdf.close(ncid);
    % create annual averages
    giwtauu  = giwtauu  + loctauu(:,:)/1000.0/12.0;
    giwtauv  = giwtauv  + loctauv(:,:)/1000.0/12.0;
end
% re-orientate
giwtauu = flipud(giwtauu');
giwtauv = flipud(giwtauv');
%
% *** load data -- wind velocity **************************************** %
%
% LOOP
for t=1:12
    % open netCDF file
    ncid = netcdf.open([str(1).path '/' str_month(t,:) str(1).exp '.' str(3).nc '.nc'],'nowrite');
    % U COMPONENT OF SURFACE AIR WIND
    varid  = netcdf.inqVarID(ncid,'usurf');
    locwvelu(:,:) = netcdf.getVar(ncid,varid);
    % V COMPONENT OF SURFACE AIR WIND
    varid  = netcdf.inqVarID(ncid,'vsurf');
    locwvelv(:,:) = netcdf.getVar(ncid,varid);
    % close netCDF file
    netcdf.close(ncid);
    % create annual averages
    giwvelu  = giwvelu  + locwvelu(:,:)/12.0;
    giwvelv  = giwvelv  + locwvelv(:,:)/12.0;
end
% re-orientate
giwvelu = flipud(giwvelu');
giwvelv = flipud(giwvelv');
%
% *** process data -- wind speed **************************************** %
%
% LOOP
for t=1:12
    % open netCDF file
    ncid = netcdf.open([str(1).path '/' str_month(t,:) str(1).exp '.' str(3).nc '.nc'],'nowrite');
    % SURFACE WIND SPEED
    varid  = netcdf.inqVarID(ncid,'wsurf');
    locwspd(:,:) = netcdf.getVar(ncid,varid);
    % close netCDF file
    netcdf.close(ncid);
    % create annual averages
    giwspd_wsma  = giwspd_wsma  + locwspd(:,:)/12.0;
end
% re-orientate
giwspd_wsma = flipud(giwspd_wsma');
%
% *** set wind stress mask ********************************************** %
%
% copy wind stress mask from ocean mask
gimaskp = gimask;
% % re-orientate
% gimaskp = flipud(gimaskp');
% derive mask values (ocean == 1, land == NaN)
gimaskp(find(gimaskp~=1.0))  = NaN;
%
% *********************************************************************** %

% *********************************************************************** %
% *** RE-GRID *********************************************************** %
% *********************************************************************** %
%
% *** set up GOLDSTEIN grid ********************************************* %
%
% Need to determine wind stress at u and v points of the Arakawa C
% grid in GENIE ...
%
%   ----v----
%   |       |
%   |   c   u
%   |       |
%   ---------
% 
% remember: [rows columns] == [j i]
% grid boudaries are as follows:
% GENIE c-grid: (golate,golone)   == jmax+1 x imax+1
% GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% GENIE v-grid: (golatve,golonve) == jmax+1 x imax
%
% create GENIE u and v edge axes
% latitude
golatue = golate;
golatve = [golatm 90.0];
% longitude
golonue = [golonm golonm(end)+360.0/imax];
golonve = golone;
% create GENIE NaN mask
gm = gomask;
gm(find(gm == 0)) = NaN;
%
% *** Put on a GOLDSTEIN grid ******************************************* %
%
% NOTE: don't forget to flip the re-gridded orientation back around again
%
% plot raw wind stress
if (optplots), plot_2dgridded(flipud(giwtauu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_u.IN'],['wind stress in -- u']); end
if (optplots), plot_2dgridded(flipud(giwtauv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_v.IN'],['wind stress in -- v']); end
% 
% -> wind stress @ u point
% NOTE: GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN u-grid\n');
[gowtauuu,gofwtauuu] = make_regrid_2d(gilone,gilate,(gimaskp.*giwtauu)',golonue,golatue,false);
gowtauuu(find(isnan(gowtauuu))) = 0.0;
gowtauuu = gowtauuu'; 
gowtauuu = gomask.*gowtauuu;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATu.out'],['wind stress out -- x @ u']); end
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN u-grid\n');
[gowtauvu,gofwtauvu] = make_regrid_2d(gilone,gilate,(gimaskp.*giwtauv)',golonue,golatue,false);
gowtauvu(find(isnan(gowtauvu))) = 0.0;
gowtauvu = gowtauvu';
gowtauvu = gomask.*gowtauvu;
if (optplots), plot_2dgridded(flipud(gm.*gowtauvu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATu.out'],['wind stress out -- y @ u']); end
% 
% -> wind stress @ v point
% NOTE: GENIE v-grid: (golatve,golonve) == jmax+1 x imax
% apply p-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN v-grid\n');
[gowtauuv,gofwtauuv] = make_regrid_2d(gilone,gilate,(gimaskp.*giwtauu)',golonve,golatve,false);
gowtauuv(find(isnan(gowtauuv))) = 0.0;
gowtauuv = gowtauuv'; 
gowtauuv = gomask.*gowtauuv;
if (optplots), plot_2dgridded(flipud(gm.*gowtauuv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATv.out'],['wind stress out -- x @ v']); end
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN v-grid\n');
[gowtauvv,gofwtauvv] = make_regrid_2d(gilone,gilate,(gimaskp.*giwtauv)',golonve,golatve,false);
gowtauvv(find(isnan(gowtauvv))) = 0.0;
gowtauvv = gowtauvv';
gowtauvv = gomask.*gowtauvv;
if (optplots), plot_2dgridded(flipud(gm.*gowtauvv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATv.out'],['wind stress out -- y @ v']); end
% 
% -> wind velocity
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: DO NOT apply a mask (as the output is used globally by the EMBM)
% NOTE: input grid is the FOAM atmosphere grid
% replace NaNs with zeros
% u
fprintf('       - Regridding wind velocity (x) to GOLDSTEIN c-grid\n');
if (optplots), plot_2dgridded(flipud(giwvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.IN'],['wind velocity in -- x']); end
[gowvelu,gofwvelu] = make_regrid_2d(gilone,gilate,giwvelu',golone,golate,false);
gowvelu(find(isnan(gowvelu))) = 0.0;
gowvelu = gowvelu'; 
gofwvelu = gofwvelu';
if (optplots), plot_2dgridded(flipud(gowvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.OUT'],['wind velocity out -- x']); end
% v
fprintf('       - Regridding wind velocity (y) to GOLDSTEIN c-grid\n');
if (optplots), plot_2dgridded(flipud(giwvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.IN'],['wind velocity in -- y']); end
[gowvelv,gofwvelv] = make_regrid_2d(gilone,gilate,giwvelv',golone,golate,false);
gowvelv(find(isnan(gowvelv))) = 0.0;
gowvelv = gowvelv';
gofwvelv = gofwvelv';
if (optplots), plot_2dgridded(flipud(gowvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.OUT'],['wind velocity out -- y']); end
%
% -> wind speed
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: create only no-mask version
% NOTE: only one grid is used 
%       (no seperate atmopsheric grid or gilonae, gilatae axis vectors)
% replace NaNs with zeros
% apply GENIE mask to output wind speed
fprintf('       - Regridding wind speed to GOLDSTEIN c-grid\n');
% wspeed_wsma
plot_2dgridded(flipud(giwspd_wsma),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_wsma.IN'],['wind speed in']);
[gowspdall,gofwspdall] = make_regrid_2d(gilone,gilate,giwspd_wsma',golone,golate,false);
gowspdall = gowspdall';
gofwspdall = gofwspdall';
if (optplots), plot_2dgridded(flipud(gowspdall),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd_wsma.OUTALL'],['wind speed out -- no mask']); end
wspeed_wsma = gomask.*gowspdall;
wspeed_wsma(find(isnan(wspeed_wsma))) = 0.0;
%
% *** Copy to output arrays ********************************************* %
%
% -> wind stress
wstress(:,:,1) = flipud(gowtauuu); %g_taux_u
wstress(:,:,2) = flipud(gowtauuv); %g_taux_v
wstress(:,:,3) = flipud(gowtauvu); %g_tauy_u
wstress(:,:,4) = flipud(gowtauvv); %g_tauy_v
% -> wind velocity
wvelocity(:,:,1) = flipud(gowvelu);
wvelocity(:,:,2) = flipud(gowvelv);
%
% *********************************************************************** %

% *********************************************************************** %
% *** SAVE DATA ********************************************************* %
% *********************************************************************** %
%
% Save regridded data to file
% Taux at u point (g_taux_u == gowtauuu)
outname = [str(2).dir '/' str(2).exp '.taux_u.dat'];
c = wstress(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (u point) data to %s\n',outname);
% Taux at v point (g_taux_v == gowtauuv)
outname = [str(2).dir '/' str(2).exp '.taux_v.dat'];
c = wstress(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (v point) data to %s\n',outname);
% Tauy at u point (g_tauy_u == gowtauvu)
outname = [str(2).dir '/' str(2).exp '.tauy_u.dat'];
c = wstress(:,:,3); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (u point) data to %s\n',outname);
% Tauy at v point (g_tauy_v == gowtauvv)
outname = [str(2).dir '/' str(2).exp '.tauy_v.dat'];
c = wstress(:,:,4); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (v point) data to %s\n',outname);
% U wind velocity 
outname = [str(2).dir '/' str(2).exp '.wvelx.dat'];
c = wvelocity(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written u wind speed data to %s\n',outname);
% V wind speed
outname = [str(2).dir '/' str(2).exp '.wvely.dat'];
c = wvelocity(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written v wind speed data to %s\n',outname);
% Save 2-D ASCII wind speed scalar for BIOGEM
outname = [str(2).dir '/' str(2).exp '.windspeed_wsma.dat'];
a = wspeed_wsma;
save(outname,'a','-ascii');
fprintf('       - Written BIOGEM windspeed data to %s\n',outname);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
% *********************************************************************** %
