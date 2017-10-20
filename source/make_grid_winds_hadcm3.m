function [wstress,wvelocity,wspeed] = make_grid_winds_hadcm3(giloncm,gilonce,gilatcm,gilatce,gilonpm,gilonpe,gilatpm,gilatpe,gimask,golonm,golone,golatm,golate,gomask,str);
% make_grid_winds_hadcm3
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
%   str(1).nc == par_nc_axes_name
%   str(2).nc == par_nc_topo_name
%   str(3).nc == par_nc_atmos_name
%   str(4).nc == par_nc_ocean_name
%   str(5).nc == par_nc_coupl_name
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
%
% *********************************************************************** %

% *********************************************************************** %
% *** SET UP GRID ******************************************************* %
% *********************************************************************** %
%
% *** OPEN netCDF DATA FILE ********************************************* %
%
% open netCDF file
ncid = netcdf.open([str(1).path '/' str(1).exp '/' str(5).nc '.nc'],'nowrite');
% read netCDf information
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
%
% *** load data ********************************************************* %
%
% NOTE: flip array around diagonal to give (lon,lat) array orientation
%       & ensure double
%
%  0.0     3.75    7.5
%
%   c*******c*******c   90.0
%   *   |       |
%   *---p-------p----   88.75
%   *   |       |
%   c   |   c   |   c   87.5
%   *   |       |
%   *---p-------p----   86.25
%   *   |       |
%   c   |   c   |   c   85.0
%
% c grid == concentrations, vertical fluxes etc
%        -> edges are: gilone,gilate
% t grid == all wind products
%        -> edges are: gilonm,gilatm
%
% load TAUX
varid  = netcdf.inqVarID(ncid,'taux_mm_hyb');
loctau(:,:) = netcdf.getVar(ncid,varid);
giwtauu = double(loctau');
% load TAUY
varid  = netcdf.inqVarID(ncid,'tauy_mm_hyb');
loctau(:,:) = netcdf.getVar(ncid,varid);
giwtauv = double(loctau');
% load U
varid  = netcdf.inqVarID(ncid,'u_mm_10m');
locwvel(:,:) = netcdf.getVar(ncid,varid);
giwvelu = double(locwvel');
% load V
varid  = netcdf.inqVarID(ncid,'v_mm_10m');
locwvel(:,:) = netcdf.getVar(ncid,varid);
giwvelv = double(locwvel');
%         % load wind speed scalar
%         varid  = netcdf.inqVarID(ncid,'wind_mm_10m');
%         old_wspd(:,:) = netcdf.getVar(ncid,varid);
%         new_wspd = double(old_wspd);
%         wspd = maskv.*flipdim(new_wspd',1);
% calculate wind speed
giwspd2 = (giwvelu.^2 + giwvelv.^2).^0.5;
%
% *** close netCDF file ************************************************* %
%
netcdf.close(ncid);
%
% *** create land-sea masks for both c and t grids ********************** %
%
% c grid == concentrations, vertical fluxes etc
%        -> edges are: gilonce,gilace
% p grid == wind products
%        -> edges are: gilonpe,gilape
%
fprintf('       - Creating wind product input mask\n');
% create ocean (only) mask -- c grid
gimaskc = gimask;
gimaskc(find(gimaskc == 0)) = NaN;
% create ocean (only) mask -- p grid
% NOTE: edges are: gilonm,gilatm
% *assumption* :: ALL p-grid points surrounding an invalid c-grid point
%                 are also invalid
% NOTE: need search both polar boundaries seperately
% NOTE: also, search only Western boundary
%
n_lonc = length(giloncm);
n_latc = length(gilatcm);
n_lonp = length(gilonpm);
n_latp = length(gilatpm);
%
gimaskp = zeros(n_latp,n_lonp);
gimaskp = gimaskp + 1;
%
%   c   |   c   |   c       c   |   *
%   *   |       |               |   *
%   *---p-------p------   ------p---*
%   *   |       |               |   *
%   c   |   c   |   c       c   |   *
%   *   |       |               |   *
%   *---p-------p------   ------p---*
%   *   |       |               |   *
%   c*******c*******c**   **c********
%
% first -- Western boundary
lonc = 1;
latc = 1;
if isnan(gimaskc(latc,lonc)),
    gimaskp(latc,n_lonc) = NaN;
    gimaskp(latc,lonc) = NaN;
end
for latc=2:n_latc-1,
    if isnan(gimaskc(latc,lonc)),
        gimaskp(latc-1,n_lonc) = NaN;
        gimaskp(latc-1,lonc) = NaN;
        gimaskp(latc,n_lonc) = NaN;
        gimaskp(latc,lonc) = NaN;
    end
end
latc = n_latc;
if isnan(gimaskc(latc,lonc)),
    gimaskp(latc-1,n_lonc) = NaN;
    gimaskp(latc-1,lonc) = NaN;
end
% then -- all remaining longitude points
for lonc=2:n_lonc,
    latc = 1;
    if isnan(gimaskc(latc,lonc)),
        gimaskp(latc,lonc-1) = NaN;
        gimaskp(latc,lonc) = NaN;
    end
    for latc=2:n_latc-1,
        if isnan(gimaskc(latc,lonc)),
            gimaskp(latc-1,lonc-1) = NaN;
            gimaskp(latc-1,lonc) = NaN;
            gimaskp(latc,lonc-1) = NaN;
            gimaskp(latc,lonc) = NaN;
        end
    end
    latc = n_latc;
    if isnan(gimaskc(latc,lonc)),
        gimaskp(latc-1,lonc-1) = NaN;
        gimaskp(latc-1,lonc) = NaN;
    end
end
%
plot_2dgridded(flipud(gimaskc),999.0,'',[[str(2).dir '/' str(2).exp] '.maskc.IN'],['c-grid mask']);
plot_2dgridded(flipud(gimaskp),999.0,'',[[str(2).dir '/' str(2).exp] '.maskp.IN'],['p-grid mask']);
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
% grid doudaries are as follows:
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
plot_2dgridded(flipud(giwtauu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_u.IN'],['wind stress in -- u']);
plot_2dgridded(flipud(giwtauv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_v.IN'],['wind stress in -- v']);
% 
% -> wind stress @ u point
% NOTE: GENIE u-grid: (golatue,golonue) == jmax   x imax+1
% apply P-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN u-grid\n');
[gowtauuu,gofwtauuu] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauu)',golonue,golatue,false);
gowtauuu(find(isnan(gowtauuu))) = 0.0;
gowtauuu = gowtauuu'; 
gowtauuu = gomask.*gowtauuu;
plot_2dgridded(flipud(gm.*gowtauuu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATu.out'],['wind stress out -- x @ u']);
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN u-grid\n');
[gowtauvu,gofwtauvu] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauv)',golonue,golatue,false);
gowtauvu(find(isnan(gowtauvu))) = 0.0;
gowtauvu = gowtauvu';
gowtauvu = gomask.*gowtauvu;
plot_2dgridded(flipud(gm.*gowtauvu),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATu.out'],['wind stress out -- y @ u'])
% 
% -> wind stress @ v point
% NOTE: GENIE v-grid: (golatve,golonve) == jmax+1 x imax
% apply P-grid mask to input wind stress grid
% replace NaNs with zeros
% apply GENIE mask to output wind stress
% u
fprintf('       - Regridding wind stress (x) to GOLDSTEIN v-grid\n');
[gowtauuv,gofwtauuv] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauu)',golonve,golatve,false);
gowtauuv(find(isnan(gowtauuv))) = 0.0;
gowtauuv = gowtauuv'; 
gowtauuv = gomask.*gowtauuv;
plot_2dgridded(flipud(gm.*gowtauuv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_xATv.out'],['wind stress out -- x @ v']);
% v
fprintf('       - Regridding wind stress (y) to GOLDSTEIN v-grid\n');
[gowtauvv,gofwtauvv] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwtauv)',golonve,golatve,false);
gowtauvv(find(isnan(gowtauvv))) = 0.0;
gowtauvv = gowtauvv';
gowtauvv = gomask.*gowtauvv;
plot_2dgridded(flipud(gm.*gowtauvv),999.0,'',[[str(2).dir '/' str(2).exp] '.wtau_yATv.out'],['wind stress out -- y @ v'])
% 
% -> wind velocity
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% NOTE: DO NOT apply P-grid mask to input wind velocity grid
%       (as the output is used globally by the EMBM)
% replace NaNs with zeros
% u
fprintf('       - Regridding wind velocity (x) to GOLDSTEIN c-grid\n');
plot_2dgridded(flipud(giwvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.IN'],['wind velocity in -- x']);
[gowvelu,gofwvelu] = make_regrid_2d(gilonpe,gilatpe,giwvelu',golone,golate,false);
gowvelu(find(isnan(gowvelu))) = 0.0;
gowvelu = gowvelu'; 
gofwvelu = gofwvelu';
plot_2dgridded(flipud(gowvelu),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_x.OUT'],['wind velocity out -- x']);
% v
fprintf('       - Regridding wind velocity (y) to GOLDSTEIN c-grid\n');
plot_2dgridded(flipud(giwvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.IN'],['wind velocity in -- y']);
[gowvelv,gofwvelv] = make_regrid_2d(gilonpe,gilatpe,giwvelv',golone,golate,false);
gowvelv(find(isnan(gowvelv))) = 0.0;
gowvelv = gowvelv';
gofwvelv = gofwvelv';
plot_2dgridded(flipud(gowvelv),999.0,'',[[str(2).dir '/' str(2).exp] '.wvel_y.OUT'],['wind velocity out -- y']);
%
% -> wind speed
% NOTE: GENIE c-grid: (golate,golone) == jmax+1 x imax+1
% create both output grids with and without the (p-grid) mask applied to
% the input wind speed grid
% replace NaNs with zeros
% apply GENIE mask to output wind speed
fprintf('       - Regridding wind speed to GOLDSTEIN c-grid\n');
plot_2dgridded(flipud(giwspd2),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd.IN'],['wind speed in']);
[gowspd2all,gofwspd2all] = make_regrid_2d(gilonpe,gilatpe,giwspd2',golone,golate,false);
gowspd2all = gowspd2all';
gofwspd2all = gofwspd2all';
[gowspd2,gofwspd2] = make_regrid_2d(gilonpe,gilatpe,(gimaskp.*giwspd2)',golone,golate,false);
gowspd2 = gowspd2';
gofwspd2 = gofwspd2';
plot_2dgridded(flipud(gowspd2all),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd.OUTALL'],['wind speed out -- without mask']);
wspeed = gomask.*gowspd2;
wspeed(find(isnan(wspeed))) = 0.0;
plot_2dgridded(flipud(gm.*wspeed),999.0,'',[[str(2).dir '/' str(2).exp] '.wspd.OUT'],['wind speed out -- with mask']);
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
outname = [str(2).dir '/' str(2).exp '_taux_u.dat'];
c = wstress(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (u point) data to %s\n',outname);
% Taux at v point (g_taux_v == gowtauuv)
outname = [str(2).dir '/' str(2).exp '_taux_v.dat'];
c = wstress(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau u (v point) data to %s\n',outname);
% Tauy at u point (g_tauy_u == gowtauvu)
outname = [str(2).dir '/' str(2).exp '_tauy_u.dat'];
c = wstress(:,:,3); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (u point) data to %s\n',outname);
% Tauy at v point (g_tauy_v == gowtauvv)
outname = [str(2).dir '/' str(2).exp '_tauy_v.dat'];
c = wstress(:,:,4); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written tau v (v point) data to %s\n',outname);
%
% X wind velocity 
outname = [str(2).dir '/' str(2).exp '_wvelx.dat'];
c = wvelocity(:,:,1); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written u wind speed data to %s\n',outname);
% Y wind speed
outname = [str(2).dir '/' str(2).exp '_wvely.dat'];
c = wvelocity(:,:,2); b = permute(c, [2 1]); a = reshape(b, [imax*jmax 1]);
save (outname, 'a', '-ascii');
fprintf('       - Written v wind speed data to %s\n',outname);
% Save 2-D ASCII wind speed scalar for BIOGEM
outname = [str(2).dir '/' str(2).exp '_windspeed.dat'];
a = wspeed;
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
