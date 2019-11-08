function [] = muffingen(POPT)
% muffingen
%
%   ***********************************************************************
%   *** muffingen [MAKE MUFFIN CONFIG] ************************************
%   ***********************************************************************
%
%   muffingen(POPT)
%   takes 1 argument:
%   POPT [STRING] (e.g., 'muffingen_settings_BLANK')
%   --> the string for the configuration parameter file
%   --> if an empty (i.e., '') value is passed to this parameter
%       then the default parameter set is used (muffingen_settings_BLANK)
%
%   NOTE: grid/matrix orientations are :: [LAT,LON] (rows x columns)
%         (i.e., the orientation as you would view the raw data)
%
%   NOTE: for plotting grids (plot_2dgridded)
%         arrays need to be flipped up-down
%         (to un-do the visually-correct orientation used as standard)
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/01/04: renamed and generalized yoolarate_hadcm3
%   14/02/09: changed to new re-gridding routine
%             also updated plotting
%             set opengl rendering to 'neverselect' to fix plotting bug
%   15/04/30: moved experiment name into parameter file
%   16/03/10: time to sort sh*t out ...
%             added new parameter file parameters
%             (sub-gridding)
%   16/08/23: no really ... this time ...
%             (and complete re-build and simplification of the function)
%   16/09/09: a little mroe getting-it-together ...
%   17/02/??: lots of work ...
%   17/02/23: removed most input parameters (now in option file)
%             & renamed function
%   17/02/24: added paths
%             added user-selectable output dir for output
%   17/02/26: edits to results path
%             alt option for 'GCM' ('k1')
%             improved code formatting and messaging
%             updated island ordering (largest now occurring first as #1)
%             & creation of .psiles grid
%   17/02/27: fixed paths calculation
%   17/02/28: added visual output of paths
%   17/03/05: worked on winds ...
%   17/03/10: more winds ...
%   17/03/13: added SEDGEM component
%             removed unused first passed parameter
%   17/03/18: added user modification of mask
%             added used modification of k1 (topography)
%             *** VERSION 0.1 *********************************************
%   17/03/20: adjusted a few parameter option names
%             fixed user-editing of mask and k1
%             added user-editing of borders
%   17/03/22: developed text file (k1 or mask) input
%   17/03/23: fixed paths output
%             revised & added EXAMPLES (but need cleaning up)
%   17/03/24: further path revision
%   17/03/24: path generation fixing ... finished(!)
%   17/03/26: added logging
%             added saving of config file parameter lines
%             fixed bug in path cleanup
%             added checking for tripple junctions in paths
%             *** VERSION 0.2 *********************************************
%   17/03/31: added FOAM option code
%             split apart functions for HadCM3 vs. FOAM
%             renamed GCM option hadcm3l -> hadcm3
%   17/04/02: completed FOAM GCM option (wind product re-gridding)
%             *** VERSION 0.3 *********************************************
%   17/04/18: sorted out blank option
%             added (non GCM) albedo profile generation
%             revised messaging
%             *** VERSION 0.31 ********************************************
%   17/04/19: revised config file output (added SEDGEM and ROKGEM lines)
%             added calculation of solar constant and salinity (from age)
%             *** VERSION 0.32 ********************************************
%   17/06/05: fixed missing albedo parameter output ...
%             *** VERSION 0.33 ********************************************
%   17/06/07: adjusted mask filtering option behavior
%             *** VERSION 0.34 ********************************************
%   17/07/21: added orbital parameters
%             added options for prescribing specific netCDF filenames
%   17/07/22: revised file path and name passing and usage
%             *** VERSION 0.35 ********************************************
%   17/07/23: added BIOGEM parameters to base config output
%             *** VERSION 0.36 ********************************************
%   17/07/27: added simple ('roof') runoff generator
%             *** VERSION 0.37 ********************************************
%   17/07/28: small modifications to grid filtering [fun_grid_topo_filter]
%             removal fo flagging of border topology issues in border #1
%             (becasue it is not used ot generate a path from)
%             [find_grid_borders_check]
%             *** VERSION 0.4 *********************************************
%   17/07/30: added plot to random runoff generation
%             *** VERSION 0.41 ********************************************
%   17/07/31: fixed bug in roofing runoff scheme
%             *** VERSION 0.42 ********************************************
%   17/08/04: added air-sea gas exchange parameter re-scaling (HadCM3 FOAM)
%             *** VERSION 0.43 ********************************************
%   17/08/15: edited output format for 2D albedo
%   17/08/16: added minimum ocean k value parameter
%             (e.g. for flat-bottom ocean)
%             *** VERSION 0.5 *********************************************
%   17/08/22: added paleo Ca/Mg
%             *** VERSION 0.51 ********************************************
%   17/08/24: added max age (100 Ma) for paleo Ca/Mg
%             *** VERSION 0.52 ********************************************
%   17/08/25: fixed bug in written-out file input paths
%             *** VERSION 0.53 ********************************************
%   17/10/13: added data path
%             *** VERSION 0.54 ********************************************
%   17/10/16: changed '\' -> '/'
%             ('/' valid under windows MATLAB (as well as normal Win '\')
%             fixed parsing of which return (to avoid occurrence of '//')
%             *** VERSION 0.55 ********************************************
%             NOW MacOSX FRIENDLY!!!
%   17/10/19: added zip-file extraction for HadCM3L 'sed' netCDF file
%             *** VERSION 0.56 ********************************************
%   17/10/20: added par_wor_name (world name) length check
%             revised wind velocity component naming (now x and y)
%             (consistent with tau components)
%             *** VERSION 0.57 ********************************************
%   17/10/22: added option to control assignment of zonal wind stress
%             added calibrated air-sea gas exchange scaling parameter
%             for zonal wind stress
%             *** VERSION 0.58 ********************************************
%   17/11/29: fixed erroreous default longitude of perihelion
%             *** VERSION 0.59 ********************************************
%   17/12/27: minor comment edits, plus corrected default parameter set
%   18/02/13: changed log file name
%             moved order of <RE-GRID TOPO> (later in sequence)
%             *** VERSION 0.60 ********************************************
%   18/04/19: fixes for pole-to-pole continents
%             *** VERSION 0.61 ********************************************
%   18/04/19: added check and warning when trying to turn land cell to ocn,
%             but when there is no ocean depth value available
%             *** VERSION 0.62 ********************************************
%   18/09/19: added minimum wind stress value [make_grid_winds_zonal.m]
%             *** VERSION 0.63 ********************************************
%   18/10/17: minor messaging changes
%             *** VERSION 0.64 ********************************************
%   19/01/29: added a parameter to enable n ocean layers to be selected,
%             but with maximum depth scaled such that surface ocean layer
%             is the same (this is the parameter -- par_sur_D)
%             *** VERSION 0.65 ********************************************
%   19/02/21: fixed missing pass of par_sur_D in 2nd make_genie_grid call
%             *** VERSION 0.66 ********************************************
%   19/02/25: made distinction between HadCM3 and HadCM3L
%             (including moving sed grid generation earlier)
%             added option not to create and show plots [MUCH FASTER!!]
%             rationalized stage numbering
%             *** VERSION 0.67 ********************************************
%   19/02/28: added option for zonal rather than GCM winds for GCM configs
%             fixed issues with hadcm3/hadcm3l distinction
%             name more consistent wind product file naming convention
%             (replacing '_' with '.')
%             *** VERSION 0.68 ********************************************
%   19/03/03: changed detault of surface layer depth scaling [par_sur_D]
%             to zero (to disable non conventional grid creation)
%             NOTE: scale depth for a 5000 m 16-level ocean, is 80.840706 m
%             *** VERSION 0.69 ********************************************
%   19/03/11: added option for CESM input
%             also tidied up the HadCM3/HadCM3l options a little
%             *** VERSION 0.70 ********************************************
%   19/03/11: fixed FOAM nc name bug
%             *** VERSION 0.71 ********************************************
%   19/03/12: further back-compatability
%             added new option for specifying how wind speed is
%             averaged/calculated
%             *** VERSION 0.72 ********************************************
%   19/03/13: added new mask option for modifying FINAL mask
%             *** VERSION 0.73 ********************************************
%   19/03/17: finalized CESM air-sea gas transfer coefficient settings
%             *** VERSION 0.74 ********************************************
%   19/05/15: fix for incorrect (sign, NaN) SEDGEM topo
%             and shitty SEDGEM mask output
%             *** VERSION 0.75 ********************************************
%   19/06/10: added (crude) functionality for ROCKE-3D GCM input
%             NOTE: no air-sea gas exchange calibration yet
%                   (or any testing ...)
%             *** VERSION 0.76 ********************************************
%   19/06/19: generalized out grid in make_regrid_2d
%             *** VERSION 0.77 ********************************************
%   19/07/04: initialization bug-fix
%             *** VERSION 0.78 ********************************************
%   19/07/08: added an alt 'k1' format -- 'k2'
%             (as per k1, but without the 'borders' / data buffer)
%             *** VERSION 0.79 ********************************************
%   19/08/14: minor netCDF name bug-fixes
%             *** VERSION 0.80 ********************************************
%   19/10/23: change to how missing and empty parameter options are handled
%             for GCM netCDF file names
%             *** VERSION 0.81 ********************************************
%   19/11/08: added mask editing modifications
%             *** VERSION 0.82 ********************************************
%   ***********************************************************************
%%

% *********************************************************************** %
% *** INITIALIZE MUFFINGEN ********************************************** %
% *********************************************************************** %
%
% *** initialize muffingen ********************************************** %
%
disp([' ']);
disp(['>>> INITIALIZING ...']);
% set function name
str_function = 'muffingen';
% set version!
par_muffingen_ver = 0.82;
% set date
str_date = [datestr(date,11), datestr(date,5), datestr(date,7)];
% close existing plot windows
% end any active journaling
% NOTE: don't clear variable space here ...
close all;
ans = get(0,'Diary');
if strcmp(ans,'on'), diary off; end
% load options file
if isempty(POPT), POPT='muffingen_settings_BLANK'; end
eval(POPT);
%
% *** backwards compatability ******************************************* %
% 
% zonal wind-stress generaton parameter
if ~exist('par_tauopt','var'), par_tauopt = 0; end
% surface layer reference thickness (m)
if ~exist('par_sur_D','var'),  par_sur_D  = 0.0; end
% minimum k level
if ~exist('par_min_k','var'),  par_min_k  = 1; end
% create plots?
if ~exist('opt_plots','var'),  opt_plots  = true; end
% debug?
if ~exist('opt_debug','var'),  opt_debug  = false; end
% make zonal winds (rather than re-grid GCM)?
if ~exist('opt_makezonalwind','var'), opt_makezonalwind = false; end
% set dummy mask nc name
% NOTE: backwards compatability provided by
%       subsequent default variable name setting
%if ~exist('par_nc_mask_name','var'), par_nc_mask_name = ''; end
% mask mask!
if ~exist('par_mask_mask_name','var'), par_mask_mask_name = ''; end
%
% *** check / filter options ******************************************** %
%
% age parameter
if ~exist('par_age','var'),
    par_age = 0.0;
    par_age_emty = true;
else
    par_age_emty = false;
end
% process GCM name string
if strcmp(par_gcm,'hadcm3l'), par_gcm = 'hadcm3l'; end
if strcmp(par_gcm,'HadCM3L'), par_gcm = 'hadcm3l'; end
if strcmp(par_gcm,'HADCM3L'), par_gcm = 'hadcm3l'; end
if strcmp(par_gcm,'hadcm3'),  par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HadCM3'),  par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HADCM3'),  par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'um'),      par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'UM'),      par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'FOAM'),    par_gcm = 'foam';   end
if strcmp(par_gcm,'CESM'),    par_gcm = 'cesm';   end
if strcmp(par_gcm,'ROCKEE'),  par_gcm = 'rockee'; end
if strcmp(par_gcm,'K1'),      par_gcm = 'k1';     end
if strcmp(par_gcm,'.k1'),     par_gcm = 'k1';     end
if strcmp(par_gcm,'.K1'),     par_gcm = 'k1';     end
if strcmp(par_gcm,'K2'),      par_gcm = 'k2';   end
if strcmp(par_gcm,'MASK'),    par_gcm = 'mask';   end
if strcmp(par_gcm,'dat'),     par_gcm = 'mask';   end
if strcmp(par_gcm,'.dat'),    par_gcm = 'mask';   end
if strcmp(par_gcm,''),        par_gcm = 'blank';  end
if strcmp(par_gcm,'BLANK'),   par_gcm = 'blank';  end
if strcmp(par_gcm,'none'),    par_gcm = 'blank';  end
if strcmp(par_gcm,'NONE'),    par_gcm = 'blank';  end
% adjust options accroding to input (GCM) type
switch par_gcm
    case {'hadcm3','hadcm3l','foam','cesm','rockee'}
    case {'k1','mask','k2'}
    otherwise
        opt_makeall=false;
        opt_user=true;
end
% deal with meta selections
if opt_makeall
    opt_makemask=true;    %
    opt_maketopo=true;    %
    opt_makeocean=true;   %
    opt_makerunoff=true;  %
    opt_makewind=true;    %
    opt_makealbedo=true;  %
    opt_makeseds=true;    %
    opt_filtermask=true;  %
    opt_filtertopo=true;  %
end
if opt_makeseds
    opt_maketopo=true;
end
% rename variables
% NOTE: for some earlier code consistency ... now redundant ...
imax = par_max_i;
jmax = par_max_j;
kmax = par_max_k;
str_nameout = par_wor_name;
% initialize optional netCDF filenames
if ~exist('par_nc_topo_name','var'),  par_nc_topo_name  = ''; end
if ~exist('par_nc_mask_name','var'),  par_nc_mask_name  = ''; end
if ~exist('par_nc_axes_name','var'),  par_nc_axes_name  = ''; end
if ~exist('par_nc_atmos_name','var'), par_nc_atmos_name = ''; end
if ~exist('par_nc_ocean_name','var'), par_nc_ocean_name = ''; end
if ~exist('par_nc_coupl_name','var'), par_nc_coupl_name = ''; end
% -> set default variable names
switch par_gcm
    case {'hadcm3','hadcm3l'}
        if isempty(par_nc_topo_name),  par_nc_topo_name  = [par_expid '.qrparm.omask']; end
        if isempty(par_nc_mask_name),  par_nc_mask_name  = [par_expid '.qrparm.mask']; end
        if isempty(par_nc_axes_name),  par_nc_axes_name  = [par_expid 'a.pdclann']; end
        if isempty(par_nc_atmos_name), par_nc_atmos_name = [par_expid '_sed']; end
        if isempty(par_nc_ocean_name), par_nc_ocean_name = [par_expid '.qrparm.mask']; end
        if isempty(par_nc_coupl_name), par_nc_coupl_name = [par_expid 'a.pdclann']; end
    case ('foam')
        if isempty(par_nc_topo_name),  par_nc_topo_name  = 'topo'; end
        if isempty(par_nc_mask_name),  par_nc_mask_name  = par_nc_topo_name; end
        if isempty(par_nc_axes_name),  par_nc_axes_name  = par_nc_topo_name; end
        if isempty(par_nc_atmos_name), par_nc_atmos_name = 'atmos'; end
        if isempty(par_nc_ocean_name), par_nc_ocean_name = ''; end
        if isempty(par_nc_coupl_name), par_nc_coupl_name = ''; end
    case {'cesm'}
        if isempty(par_nc_topo_name),  par_nc_topo_name  = 'climo'; end
        if isempty(par_nc_mask_name),  par_nc_mask_name  = par_nc_topo_name; end
        if isempty(par_nc_axes_name),  par_nc_axes_name  = par_nc_topo_name; end
        if isempty(par_nc_atmos_name), par_nc_atmos_name = par_nc_topo_name; end
        if isempty(par_nc_ocean_name), par_nc_ocean_name = par_nc_topo_name; end
        if isempty(par_nc_coupl_name), par_nc_coupl_name = par_nc_topo_name; end
    case ('rockee')
        if isempty(par_nc_topo_name),  par_nc_topo_name  = ''; end
        if isempty(par_nc_mask_name),  par_nc_mask_name  = par_nc_topo_name; end
        if isempty(par_nc_atmos_name), par_nc_atmos_name = ''; end
        if isempty(par_nc_ocean_name), par_nc_ocean_name = ''; end
        if isempty(par_nc_axes_name),  par_nc_axes_name  = par_nc_ocean_name; end
        if isempty(par_nc_coupl_name), par_nc_coupl_name = ''; end
end
% set default annual wind averaging to not based on monthly winds
switch par_gcm
    case {'hadcm3'}
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = 'uvaa'; end
    case {'hadcm3l'}
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = 'uvaa'; end
    case ('foam')
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = 'uvaa'; end
    case {'cesm'}
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = 'wsma'; end
    case {'rockee'}
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = 'wsaa'; end
    otherwise
        if ~exist('par_wspeed_avstr','var'), par_wspeed_avstr = ''; end
end
%
% *** initialize I/O **************************************************** %
%
% add library path to muffingen functions
% NOTE: find where muffin lives ...
%       remove its name (+ '.m' extension) from the returned path ...
%       add relative source path to it
tmp_path = which(str_function);
tmp_path = tmp_path(1:end-length(str_function)-3);
addpath([tmp_path '/' par_dpath_source]);
addpath([tmp_path '/' 'DATA']);
% set full file name string
if ~isempty(par_gcm),
    str_file = [str_function '.' par_gcm '.' str_nameout];
else
    str_file = [str_function '.' str_nameout];
end
% set/create output directory
if opt_outputdir,
    str_dirout = uigetdir(par_pathout);
else
    if isempty(par_pathout)
        str_dirout = [pwd '/' par_wor_name];
    else
        if ~(exist(par_pathout,'dir') == 7), mkdir(par_pathout); end
        str_dirout = [pwd '/' par_pathout '/' par_wor_name];
    end
    if ~(exist(str_dirout,'dir') == 7), mkdir(str_dirout); end
end
%
% *** create strings object ********************************************* %
%
str = struct('gcm', {}, 'exp', {}, 'path', {}, 'dir', {}, 'nc', {});
str = setfield(str, {1}, 'gcm', par_gcm);
str = setfield(str, {1}, 'exp', par_expid);
str = setfield(str, {2}, 'exp', str_nameout);
str = setfield(str, {1}, 'path', par_pathin);
str = setfield(str, {2}, 'path', par_pathout);
str = setfield(str, {1}, 'dir', '');
str = setfield(str, {2}, 'dir', str_dirout);
str = setfield(str, {2}, 'nc', par_nc_topo_name);
str = setfield(str, {4}, 'nc', par_nc_mask_name);
str = setfield(str, {1}, 'nc', par_nc_axes_name);
str = setfield(str, {3}, 'nc', par_nc_atmos_name);
str = setfield(str, {6}, 'nc', par_nc_ocean_name);
str = setfield(str, {5}, 'nc', par_nc_coupl_name);
str = setfield(str, {1}, 'wspd', par_wspeed_avstr);
str = setfield(str, {1}, 'mask', par_mask_mask_name);
%
% *** initialize reporting ********************************************** %
%
% initialize muffingen step number
n_step = 0;
% create copy of .m file options
copyfile([POPT '.m'],[str(2).dir '/' POPT '_' str_date '.m'])
% start logging
% NOTE: delete any existing (current date) log file in the directory first
str_log = [str(2).dir '/' 'log.' str_function '.' str_date '.txt'];
if (exist(str_log,'file') == 2), delete(str_log); end
diary(str_log)
% display warm and comforting welcoming message
disp([' ']);
disp(['------------------------------------------------------------']);
disp(['   Hello! Welcome to muffingen v',num2str(par_muffingen_ver)]);
disp(['   We are going to make a GREAT model configuration!']);
disp(['------------------------------------------------------------']);
disp([' ']);
%
% *********************************************************************** %

% *********************************************************************** %
% *** GET SH*T DONE ***************************************************** %
% *********************************************************************** %
%
% %%% SUMMARY OF STEP SEQUENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: not all of these steps are carried out, depending on the options
% (initialize)
% CONFIRM OPTIONS
% CREATING GENIE GRID
% READ AXES INFORMATION
% LOAD TOPO & MASK DATA
% RE-GRID MASK
% FILTER MASK
% ADJUST MASK -- USER!
% RE-GRID TOPO
% RE-GRID VERTICALLY
% ADJUST TOPO -- AUTOMATIC BATHYMETRY FILTERING
% ADJUST TOPO -- USER!
% CALCULATE RUNOFF & COMPLETE k1 FILE
% IDENTIFY ISLANDS
% UPDATE ISLANDS AND ISLAND PATHS
% GENERATE ISLAND PATHS
% GENERATE PSI ISLANDS
% GENERATE SEDIMENT GRID
% RE-GRID WIND SPEED/STRESS DATA
% LOAD ALBEDO DATA
% RE-GRID & PROCESS ALBEDO
% GENERATE CONFIG FILE PARAMETER LINES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *** CONFIRM OPTIONS *************************************************** %
%
n_step = n_step+1;
disp(['>   ' num2str(n_step) '. CHECKING PRIMARY OPTIONS ...']);
% check world name
if (length(par_wor_name) ~= 8)
    disp(['       * ERROR: World name (par_wor_name) must be 8 characters long.']);
    disp(['--------------------------------------------------------']);
    disp([' ']);
    diary off;
    return;
end
% check GCM options
switch str(1).gcm
    case {'hadcm3','hadcm3l','foam','cesm','rockee'}
        disp(['       * GCM == ' str(1).gcm ' (OK)']);
    case {'k1','mask','k2'}
        disp(['       * GENIE grid will be loaded directly from k1, mask, or k2 text file: ' str(1).exp]);
    case {'blank'}
        disp(['       * A blank template grid will be generated: ' str(1).exp]);
    otherwise
        disp(['       * ERROR: Unknown GCM name or input format. Maybe it is LOSCAR? Good luck to you! (Bye bye)']);
        disp(['--------------------------------------------------------']);
        disp([' ']);
        diary off;
        return;
end
%
% *** SET UP OUTPUT GRID ************************************************ %
%
n_step = n_step+1;
disp(['>   ' num2str(n_step) '. CREATING GENIE GRID ...']);
% create GENIE grid
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de,par_max_D] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea,par_sur_D);
disp(['       - GENIE grid generated.']);
%
% *** LOAD GRID (AXES) DATA ********************************************* %
%
n_step = n_step+1;
disp(['>   ' num2str(n_step) '. READING AXES INFORMATION ...']);
%
switch str(1).gcm
    case {'hadcm3','hadcm3l','foam','cesm','rockee'}
        % read axes
        if strcmp(str(1).gcm,'hadcm3')
            % NOTE: axes need to be re-generated later (for winds etc.)
            [gi_loncm,gi_lonce,gi_latcm,gi_latce] = fun_read_axes_hadcm3(str);
        elseif strcmp(str(1).gcm,'hadcm3l')
            [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_hadcm3x(str);
        elseif strcmp(str(1).gcm,'foam')
            [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_foam(str);
        elseif strcmp(str(1).gcm,'cesm')
             [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_cesm(str);
        elseif strcmp(str(1).gcm,'rockee')
             [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_rockee(str);
        else
            disp(['       * ERROR: Unknown error.']);
            disp([' ']);
            diary off;
            return;
        end
        disp(['       - Axis info read.']);
    otherwise
        % DO NOTHING
        disp(['         (Nothing to load.)']);
end
%
% *** LOAD TOPO & MASK DATA ********************************************* %
%
% NOTE: mask is defined with '1' for ocean
%
n_step = n_step+1;
disp(['>   ' num2str(n_step) '. READING MASK & TOPO GRIDS ...']);
%
switch str(1).gcm
    case {'hadcm3','hadcm3l','foam','cesm','rockee'}
        % read topo
        if (strcmp(str(1).gcm,'hadcm3') || strcmp(str(1).gcm,'hadcm3l'))
            [gi_mask] = fun_read_omask_hadcm3x(str);
            [gi_topo] = fun_read_topo_hadcm3x(str);
        elseif strcmp(str(1).gcm,'foam')
            [gi_topo,gi_mask] = fun_read_topomask_foam(str);
        elseif strcmp(str(1).gcm,'cesm')
            [gi_topo,gi_mask] = fun_read_topomask_cesm(str);
        elseif strcmp(str(1).gcm,'rockee')
            [gi_topo,gi_mask] = fun_read_topomask_rockee(str);
        end
        disp(['       - Mask & topo info read.']);
        % plot input mask & topo
        if (opt_plots), plot_2dgridded(flipud(gi_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_in'],['mask in']); end
        if (opt_plots), plot_2dgridded(flipud(gi_topo),6000.0,'',[[str_dirout '/' str_nameout] '.topo_in'],['topo in']); end
    case {'k1','mask','k2'}
        % load topo directly
        [go_k1,go_mask,imax,jmax] = fun_read_k1(str);
        disp(['       - k1 read.']);
        % re-create GENIE grid with derived (/updated?) grid dimensions
        % NOTE: the value of kmax is taken from the config file
        %      (while imax and jmax are deduced from the file)
        [go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea,par_sur_D);
        disp(['       - GENIE grid re-generated.']);
    otherwise
        % DO NOTHING
        disp(['         (Nothing to load.)']);
end
%
% *** RE-GRID MASK ****************************************************** %
%
n_step = n_step+1;
disp(['>   ' num2str(n_step) '. RE-GRIDING MASK ...']);
%
switch str(1).gcm
    case {'hadcm3','hadcm3l','foam','cesm','rockee'}
        % initial re-gridding of mask
        % NOTE: need to transpose around [gi_mask] to have correct input format
        %       to make_regrid_2d
        %       similarly, output needs to be transposed back again
        % NOTE: pass edges of c-grid
        [go_mask,go_tmp] = make_regrid_2d(gi_lonce,gi_latce,gi_mask',go_lone,go_late,opt_debug);
        disp(['       - Mask re-gridded.']);
        go_mask = go_mask';
        go_fmask = go_mask;
        % create mask (<>= par_A_frac_threshold fractional area thresold)
        % NOTE: 1.0 == 100% ocean
        go_mask(find(go_mask>=par_A_frac_threshold)) = 1.0;
        go_mask(find(go_mask<par_A_frac_threshold))  = 0.0;
        % calculate respective fractional areas
        [si_farea,si_farearef] = fun_grid_calc_ftotarea(gi_mask,gi_lonce,gi_latce);
        [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
        disp(['       * Original land area fraction    = ', num2str(1.0-si_farea)]);
        disp(['       * Re-gridded land area fraction  = ', num2str(1.0-so_farea)]);
    case {'k1','mask','k2'}
        disp(['         (Nothing to do ... k1/mask/k2 file already loaded.)']);
        go_fmask = zeros(jmax,imax) + 1;
    otherwise
        go_mask = zeros(jmax,imax) + 1;
        disp(['       - Blank mask created (nothing to re-grid).']);
        go_fmask = go_mask;
end
% plot & save initial mask re-grid
if (opt_plots), plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.RAW'],['mask out -- RAW re-gridded']); end
%
% *** FILTER MASK ******************************************************* %
%
n_step = n_step+1;
% filter mask if requested
% NOTE: when loading in a default 'k1' file, best to skip this step
% set VERSION 0 (raw)
grid_ver = 0;
str_ver = num2str(grid_ver);
%
if opt_makemask && (opt_filtermask || (par_min_oceann > 0)),
    %
    disp(['>   ' num2str(n_step) '. FILTERING MASK ...']);
    %
    if opt_filtermask,
        %
        % BASIC MASK FILTERING
        %
        % filter out single cell embayments iteratively
        % initialize
        go_mask_fills_cnt = 1;
        % LOOP >>>
        while (go_mask_fills_cnt > 0)
            % increment VERSION
            grid_ver = grid_ver + 1;
            str_ver = num2str(grid_ver);
            % fill in 4-surrounded cells
            go_mask_fills = find_grid_4cells(go_mask);
            % fill in 3-surrounded cells (excepting open ocean bordering ones)
            go_mask_fills = go_mask_fills + find_grid_3cells(go_mask);
            % update mask
            go_mask = go_mask - go_mask_fills;
            % count number of cell changes
            go_mask_fills_cnt = sum(sum(go_mask_fills));
        end
        % <<< LOOP
        fprintf('       - Single cell embayments filtered out.\n')
        % plot mask
        if (opt_plots), plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]); end
        %
    end
    %
    if opt_filtermask,
        %
        % ADJUST MASK -- CLEAN UP POLAR CONNECTIONS
        %
        % increment VERSION
        grid_ver = grid_ver + 1;
        str_ver = num2str(grid_ver);
        % open up narrow polar connections
        go_mask_fills = find_grid_poleopen(go_mask);
        % update mask
        go_mask = go_mask + go_mask_fills;
        %
        fprintf('       - Polar connections cleaned up.\n')
        % plot mask
        if (opt_plots), plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]); end
        %
    end
    %
    if (par_min_oceann > 0),
        %
        % SMALL WATER BODY FILTERING MASK FILTERING
        %
        [go_oceans,n_oceans,i_oceans] = find_grid_oceans(go_mask);
        % plot oceans!
        if (opt_plots), plot_2dgridded(flipud(go_oceans),999,'',[[str_dirout '/' str_nameout] '.ocean_out.INIT'],['oceans out -- INITIAL']); end
        % increment VERSION
        grid_ver = grid_ver + 1;
        str_ver = num2str(grid_ver);
        % clean up small water bodies
        [go_mask,go_oceans,n_oceans] = find_grid_oceans_update(go_mask,go_oceans,n_oceans,par_min_oceann);
        fprintf('       - Small water bodies cleaned up.\n')
        % plot mask
        if (opt_plots), plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]); end
        %
    end
    %
    % calculate new fractional area
    [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
    disp(['       * Revised land area fraction = ', num2str(1.0-so_farea)]);
    %
end
%
% *** ADJUST MASK -- USER! ********************************************** %
%
n_step = n_step+1;
%
if opt_user
    %
    disp(['>   ' num2str(n_step) '. USER EDITING OF MASK ...']);
    %
    [go_mask]  = fun_grid_edit_mask(go_mask,go_fmask);
    % increment VERSION
    grid_ver = grid_ver + 1;
    str_ver = num2str(grid_ver);
    % plot mask
    if (opt_plots), plot_2dgridded(flipud(go_mask),2,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]); end
    % calculate new fractional area
    [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
    disp(['       * Revised land area fraction = ', num2str(1.0-so_farea)]);
    %
    fprintf('       - User-editing complete.\n')
    %
end
%
% *** CREATE FINAL MASK ************************************************* %
%
n_step = n_step+1;
%
disp(['>   ' num2str(n_step) '. CREATE FINAL MASK ...']);
%
% load and apply mask mask to ocean mask!
if ~isempty(str(1).mask)
    loc_mask = fun_read_mask(str);    
    go_mask(find(loc_mask == 1))  = 1;
    go_mask(find(loc_mask == -1)) = 0;
end
%
% create GENIE NaN mask
go_masknan = go_mask;
go_masknan(find(go_masknan == 0)) = NaN;
%
% plot final mask
if (opt_plots), plot_2dgridded(flipud(go_mask),99999.0,'',[[str_dirout '/' str_nameout] '.mask_out.FINAL'],['mask out -- FINAL version']); end
%
% save mask
fprint_2DM(go_mask(:,:),[],[[str_dirout '/' str_nameout] '.mask_out.FINAL.dat'],'%4.1f','%4.1f',true,false);
fprintf('       - .mask_out.FINAL.dat saved\n')
%
% calculate new fractional area
[so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
disp(['       * Final land area fraction   = ', num2str(1.0-so_farea)]);
%
% *** RE-GRID TOPO ****************************************************** %
%
n_step = n_step+1;
%
if opt_maketopo
    %
    disp(['>   ' num2str(n_step) '. RE-GRIDING TOPOGRAPHY ...']);
    %
    switch str(1).gcm
        case {'hadcm3','hadcm3l','foam','cesm','rockee'}
            % initial re-gridding of topo
            % NOTE: need to transpose around [gi_topo] to have correct input format
            %       to make_regrid_2d
            %       similarly, output needs to be transposed back again
            % NOTE: pass edges of c-grid
            [go_topo,go_ftopo] = make_regrid_2d(gi_lonce,gi_latce,gi_topo',go_lone,go_late,opt_debug);
            disp(['       - Topography re-gridded.']);
            go_topo = go_topo';
            go_ftopo = go_ftopo';
        case {'k1','k2'}
            % convert k1 to depth
            [go_topo] = fun_conv_k1(go_de,go_k1);
            disp(['         (Nothing to re-grid -- convert k1 file data.)']);
        otherwise
            go_topo = par_max_D*go_mask;
            disp(['         (Nothing to re-grid -- set uniform ocean depth.)']);
    end
    % plot & save initial topo re-grid
    if ( ~strcmp(str(1).gcm,'k1') || ~strcmp(str(1).gcm,'k2') )
        if (opt_plots), plot_2dgridded(flipud(go_topo),99999.0,'',[[str_dirout '/' str_nameout] '.topo_out.RAW'],['topo out -- RAW']); end
    end
    %
end
%
% *** RE-GRID VERTICALLY ************************************************ %
%
n_step = n_step+1;
%
if opt_maketopo,
    %
    disp(['>   ' num2str(n_step) '. RE-GRIDING OCEAN BATHYMETRY ...']);
    %
    switch str(1).gcm
        case {'k1','k2'}
            disp(['         (Nothing to re-grid as k1 file already loaded.)']);
        otherwise
            % convert depth into k levels (and create k1 grid)
            [go_k1] = find_grid_k(par_min_Dk,go_dm,go_de,go_mask,go_topo);
            fprintf('       - Bathymetry re-gridding complete.\n')
    end
    % filter min k value
    go_k1(find(go_k1 < par_min_k)) = par_min_k;
    % plot initial k1 re-grid
    if (opt_plots), plot_2dgridded(flipud(go_k1),89.0,'',[[str_dirout '/' str_nameout] '.k1_out.RAW'],['k1 out -- RAW re-gridded']); end
    %
end
%
% *** ADJUST TOPO -- AUTOMATIC BATHYMETRY FILTERING ********************* %
%
n_step = n_step+1;
%
% carry out basic automatic topo filtering
if opt_maketopo && opt_filtertopo,
    %
    disp(['>  ' num2str(n_step) '. FILTERING BATHYMETRY ...']);
    %
    [go_k1] = fun_grid_topo_filter(go_k1);
    fprintf('       - Topography filtered.\n')
    % plot adjusted k1 re-grid
    if (opt_plots), plot_2dgridded(flipud(go_k1),89.0,'',[[str_dirout '/' str_nameout] '.k1_out.FILTERED'],['k1 out -- auto filtered']); end
end
%
% *** ADJUST TOPO -- USER! ********************************************** %
%
n_step = n_step+1;
%
if (opt_maketopo && opt_user)
    % filter for ocean but not depth info (k > kmax)
    % search for ocean in the mask that has a 'land' value ...
    % set to the shallowest k1 level and a nominal middepth of that level
    loc_dry = intersect(find(go_mask == 1),find(go_k1 > par_max_k));
    if ~isempty(loc_dry)
        go_k1(loc_dry) = par_max_k;
        go_topo(loc_dry) = -go_dm(par_max_k);
    end
    %
    disp(['>  ' num2str(n_step) '. USER EDITING OF TOPO ...']);
    % user-editing! what can go wrong?
    [go_k1] = fun_grid_edit_k1(go_k1,kmax);
    % plot mask
    if (opt_plots), plot_2dgridded(flipud(go_k1),89.0,'',[str_dirout '/' str_nameout '.k1_out.USEREDITED'],['k1 out -- user edited version']); end
    % convert k-levels back to depth
    [go_topo] = fun_conv_k1(go_de,go_k1);
    %
    fprintf('       - User-editing complete.\n')
    %
end
%
if opt_maketopo
    % plot final k1
    if (opt_plots), plot_2dgridded(flipud(go_masknan.*go_k1),89.0,'',[str_dirout '/' str_nameout '.k1_out.FINAL'],['k1 out -- FINAL version']); end
    % plot final topo
    if (opt_plots), plot_2dgridded(flipud(go_masknan.*go_topo),99999.0,'',[[str_dirout '/' str_nameout] '.topo_out.FINAL'],['topo out -- FINAL version']); end
end
%
% *** CALCULATE RUNOFF & COMPLETE k1 FILE ******************************* %
%
n_step = n_step+1;
%
% NOTE: ordering is a little illogical becasue
%       make_grid_runoff_rnd requires the extended grid, while
%       make_grid_runoff_roof is easier done without ...
if opt_makeocean
    %
    disp(['>  ' num2str(n_step) '. CALCULATING RUN-OFF AND GENERATE .k1 FILE ...']);
    % (i) first, check for all ocean
    if (max(max(go_k1)) < 90), opt_makerunoff = false; end
    % (ii) create roof scheme (if selected)
    if (opt_makerunoff && par_runoffopt == 0),
        [go_k1] = make_grid_runoff_roof(go_mask,go_k1,str);
        loc_k1 = go_k1;
        loc_k1(find(loc_k1 < 91)) = 95;
        if (opt_plots), plot_2dgridded(flipud(loc_k1),95.0,'',[str_dirout '/' str_nameout '.k1_out.RUNOFF'],['k1 out -- RUNOFF']); end
    end
    % (iii) extend k1 grid
    % NOTE: mark first row: maxk+1, last as maxk+2
    %       (so, slightly different from standard/original GENIE format)
    goex_k1 = go_k1;
    goex_k1 = [goex_k1(1,:); goex_k1; goex_k1(end,:)];
    goex_k1(1,:) = kmax+1;
    goex_k1(end,:) = kmax+2;
    % add buffer columns for E-W wall
    goex_k1 = [goex_k1(:,end) goex_k1 goex_k1(:,1)];
    % (iv) create random runoff grid (if selected)
    if (opt_makerunoff && par_runoffopt == 1),
        [goex_k1] = make_grid_runoff_rnd(goex_k1,str,opt_debug);
        loc_k1 = goex_k1(2:end-1,2:end-1);
        loc_k1(find(loc_k1 < 91)) = 95;
        if (opt_plots), plot_2dgridded(flipud(loc_k1),95.0,'',[str_dirout '/' str_nameout '.k1_out.RUNOFF'],['k1 out -- RUNOFF']); end
    end
    % (v) save .k1 file
    fprint_2DM(goex_k1(:,:),[],[[str_dirout '/' str_nameout] '.k1'],'%3i','%3i',true,false);
    fprintf('       - .k1 file saved\n')
    %
end
%
% *** IDENTIFY ISLANDS ************************************************** %
%
n_step = n_step+1;
%
if opt_makeocean
    %
    disp(['>  ' num2str(n_step) '. IDENTIFY ISLANDS ...']);
    % initial islands count
    [go_islands,n_islands,i_islands] = find_grid_islands(go_mask);
    % plot islands
    if (opt_plots), plot_2dgridded(flipud(go_islands),999,'',[[str_dirout '/' str_nameout] '.islnd_out.INIT'],['island out -- INITIAL']); end
    %
end
%
% *** UPDATE ISLANDS AND ISLAND PATHS *********************************** %
%
n_step = n_step+1;
%
if opt_makeocean
    %
    disp(['>  ' num2str(n_step) '. UPDATING ISLANDS & PATHS ...']);
    % NOTE: generate all possible paths initially (and filter later)
    % (1) generate generic borders around all (initial) islands
    [go_borders] = find_grid_borders(go_mask);
    if (opt_plots), plot_2dgridded(flipud(go_borders),99999.0,'',[[str_dirout '/' str_nameout] '.brds_out.INIT'],['borders out -- INITIAL']); end
    % (2) update islands count
    %     identify islands that are insufficiently seperated (and combined)
    %     identify polar islands
    %     re-number all
    [go_islands,n_islands,i_islands,i_poles] = find_grid_islands_update(go_islands,n_islands,i_islands,go_borders,opt_makepoleswide);
    if (opt_plots), plot_2dgridded(flipud(go_islands),999,'',[[str_dirout '/' str_nameout] '.islnd_out.FINAL'],['islands out -- FINAL']); end
    % (3) update borders
    %     number borders as per bordering islands
    [go_borders] = find_grid_borders_update(go_borders,go_islands,go_mask,n_islands);
    if (opt_plots), plot_2dgridded(flipud(go_borders),999,'',[[str_dirout '/' str_nameout] '.brds_out.FILTERED'],['borders out -- FILTERED']); end
    % (4) border check
    [opt_user] = find_grid_borders_check(go_borders,opt_user);
    % (5) user editing of borders
    if opt_user
        % user-editing! what can go wrong?
        [go_borders] = fun_grid_edit_borders(go_borders,go_mask);
        % plot mask
        if (opt_plots), plot_2dgridded(flipud(go_borders),999,'',[str_dirout '/' str_nameout '.brds_out.USEREDITED'],['borders out -- user edited version']); end
    end
    % plot final borders
    if (opt_plots), plot_2dgridded(flipud(go_borders),999,'',[[str_dirout '/' str_nameout] '.brds_out.FINAL'],['borders out -- FINAL']); end
    %
end
%
% *** GENERATE ISLAND PATHS ********************************************* %
%
n_step = n_step+1;
%
if opt_makeocean
    %
    disp(['>  ' num2str(n_step) '. GENERATING .paths FILE ...']);
    % create paths
    [n_paths,v_paths,n_islands,go_paths] = find_grid_paths(go_borders,n_islands,i_poles);
    % plot paths data
    if (opt_plots), plot_2dgridded(flipud(go_paths),999,'',[[str_dirout '/' str_nameout] '.paths_out.FINAL'],['Paths file -- FINAL']); end
    % save .paths file
    fprint_paths(n_paths,v_paths,[[str_dirout '/' str_nameout] '.paths']);
    fprintf('       - .paths file saved\n')
    %
end
%
% *** GENERATE PSI ISLANDS ********************************************** %
%
n_step = n_step+1;
%
if opt_makeocean
    %
    disp(['>  ' num2str(n_step) '. GENERATING .psiles FILE ...']);
    % generate PSI islands data
    [go_psiles,n_islands_recnt] = make_grid_psiles(go_islands,i_poles);
    % plot PSI islands data
    if (opt_plots), plot_2dgridded(flipud(go_psiles),999,'',[[str_dirout '/' str_nameout] '.psiles_out.FINAL'],['PSI islands file -- FINAL']); end
    % save .psiles file
    fprint_2DM(go_psiles(:,:),[],[[str_dirout '/' str_nameout] '.psiles'],'%3i','%3i',true,false);
    fprintf('       - .psiles file saved\n')
    % carry out check on # islands
    if (n_islands ~= n_islands_recnt),
        disp(['       ! Something odd about the islands count or configuration ...']);
    end
    %
end
%
% *** GENERATE SEDIMENT GRID ******************************************** %
%
n_step = n_step+1;
%
if opt_makeseds
    %
    disp(['>  ' num2str(n_step) '. GENERATING SEDIMENT TOPO ...']);
    %
    % check sed topo re-gridding options
    switch par_sedsopt
        case 1
            % option %1 -- re-grid sediment topography
            switch str(1).gcm
                case {'hadcm3','hadcm3l','foam','cesm','rockee'}
                    % if 'high res' sed grid is requested => assume twice ocean resolution
                    % + generate new vectors of grid properties
                    if opt_highresseds
                        [gos_lonm,gos_lone,gos_latm,gos_late,gos_dm,gos_de] = make_genie_grid(2*imax,2*jmax,kmax,par_max_D,par_lon_off,opt_equalarea);
                    else
                        gos_lone = go_lone;
                        gos_late = go_late;
                    end
                    % re-grid
                    [gos_topo,gos_ftopo] = make_regrid_2d(gi_lonce,gi_latce,gi_topo',gos_lone,gos_late,opt_debug);
                    gos_topo  = gos_topo';
                    gos_ftopo = gos_ftopo';
                    disp(['       - Re-gridded sediment topo from GCM bathymetry.']);
                case {'k1','grid'}
                    % convert k1 to depth
                    [gos_topo] = fun_conv_k1(go_de,go_k1);
                    disp(['       - Converted k1 file data (nothing to re-grid).']);
                otherwise
                    %
                    gos_topo = par_max_D*go_mask;
                    disp(['       - Set uniform ocean depth (nothing to re-grid).']);
            end
        case 2
            % option %2 -- create random sediment bathymetry from mask
            [gos_topo] = make_grid_topo_sed_rnd(go_mask,opt_highresseds,go_de(kmax-(par_min_Dk-1)),par_max_D);
            disp(['       - Created randomized sediment topography (nothing to re-grid).']);
        otherwise
            % option %0 -- convert k1 to depth
            [gos_topo] = fun_conv_k1(go_de,go_k1);
            disp(['       - Converted k1 file data (nothing to re-grid).']);
    end
    %
    % set mask
    if opt_highresseds
        gos_mask = gos_topo;
        gos_mask(find(gos_mask > 0.0)) = 1.0;
    else
        gos_mask = go_mask;
    end
    % filter topo
    gos_topo(isnan(gos_topo))         = 0.0;
    gos_topo(find(gos_topo < -9.9E9)) = 0.0;
    % invert to depth (rather than height)
    gos_topo = -gos_topo;
    % apply mask
    gos_topo = gos_mask.*gos_topo;
    % plot sediment topo
    if (opt_plots), plot_2dgridded(flipud(gos_topo),9999,'',[[str_dirout '/' str_nameout] '.sedtopo_out.FINAL'],['Sediment topo -- FINAL']); end
    % save sediment topo
    fprint_2DM(gos_topo(:,:),[],[[str_dirout '/' str_nameout] '.depth.dat'],'%10.2f','%10.2f',true,false);
    fprintf('       - .depth.dat saved\n')
    % save other sediment files
    gos_sedc = 0.0*gos_mask;
    fprint_2DM(gos_sedc(:,:),gos_mask(:,:),[[str_dirout '/' str_nameout] '.sedcoremask.dat'],'%5.1f','%5i',true,false);
    fprintf('       - template file .sedcoremask.dat saved\n')
    gos_reef = 0.0*gos_mask;
    fprint_2DM(gos_reef(:,:),gos_mask(:,:),[[str_dirout '/' str_nameout] '.reefmask.dat'],'%5.1f','%5i',true,false);
    fprintf('       - template file .reefmask.dat saved\n')
end
%
% *** SWITCH GRIDS ****************************************************** %
%
n_step = n_step+1;
%
disp(['>  ' num2str(n_step) '. SWITCH GRIDS ...']);
% NOTE: only with HadCM3 do we need to switch from ocean to atm grid
switch str(1).gcm
    case {'hadcm3'}
        % re-read axes
        [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_hadcm3x(str);
        disp(['       - Axis info re-read.']);
        % re-read (atmopshere) mask
        [gi_mask] = fun_read_amask_hadcm3x(str);
        disp(['       - Mask info re-read.']);
    otherwise
        % DO NOTHING
        disp(['         (Nothing to re-read.)']);
end
%
% *** RE-GRID WIND SPEED/STRESS DATA ************************************ %
%
% n_step = n_step+1;
%
if (opt_makezonalwind)
    %
    disp(['>  ' num2str(n_step) '. CREATING ZONAL-ONLY WIND PRODUCTS ...']);
    % create GENIE grid wind products
    [wstr,wspd,g_wspd] = make_grid_winds_zonal(go_latm,go_late,go_mask,[str_dirout '/' str_nameout],par_tauopt);
    disp(['       - Generated zonal wind products.']);
elseif (opt_makewind)
    %
    disp(['>  ' num2str(n_step) '. CREATING WIND PRODUCTS ...']);
    % create GENIE grid wind products
    switch str(1).gcm
        case {'hadcm3','hadcm3l','foam','cesm','rockee'}
            % re-grid winds from GCM
            % NOTE: the sets of grids and their edges required differ
            %       between hadcm3/hadcm3l and foam
            if (strcmp(str(1).gcm,'hadcm3') || strcmp(str(1).gcm,'hadcm3l'))
                make_grid_winds_hadcm3x(gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str,opt_plots);
            elseif (strcmp(str(1).gcm,'foam'))
                make_grid_winds_foam(gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonam,gi_lonae,gi_latam,gi_latae,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str,opt_plots);
            elseif (strcmp(str(1).gcm,'cesm'))
                make_grid_winds_cesm(gi_lonce,gi_latce,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str,opt_plots);
            elseif (strcmp(str(1).gcm,'rockee'))
                make_grid_winds_rockee(gi_lonam,gi_lonae,gi_latam,gi_latae,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str,opt_plots);
            end
            disp(['       - Re-grided GCM wind products.']);
        otherwise
            make_grid_winds_zonal(go_latm,go_late,go_mask,[str_dirout '/' str_nameout],par_tauopt);
            disp(['       - Generated zonal wind products.']);
    end
end
%
% *** LOAD ALBEDO DATA ************************************************** %
%
n_step = n_step+1;
%
if opt_makealbedo
    %
    disp(['>  ' num2str(n_step) '. LOADING ALBEDO DATA ...']);
    %
    switch str(1).gcm
        case {'hadcm3','hadcm3l','foam','cesm','rockee'}
            % read albedo
            if (strcmp(str(1).gcm,'hadcm3') || strcmp(str(1).gcm,'hadcm3l'))
                [gi_albd] = fun_read_albd_hadcm3x(str);
            elseif (strcmp(str(1).gcm,'foam'))
                [gi_albd] = fun_read_albd_foam(str);
            elseif (strcmp(str(1).gcm,'cesm'))
                [gi_albd] = fun_read_albd_cesm(str);
            elseif (strcmp(str(1).gcm,'rockee'))
                [gi_albd] = fun_read_albd_rockee(gi_loncm,gi_latcm,str);
            end
            disp(['       - Read GCM albedo data.']);
            % plot input albedo
            if (opt_plots), plot_2dgridded(flipud(gi_albd),100.0,'',[[str_dirout '/' str_nameout] '.albd_in'],['albedo in']); end
        otherwise
            disp(['         (Nothing to load.)']);
    end
    %
end
%
% *** RE-GRID & PROCESS ALBEDO ****************************************** %
%
n_step = n_step+1;
%
if opt_makealbedo
    %
    disp(['>  ' num2str(n_step) '. CREATING ALBEDO DATA ...']);
    %
    switch str(1).gcm
        case {'hadcm3','hadcm3l','foam','cesm','rockee'}
            % re-grid
            [go_albd,go_falbd] = make_regrid_2d(gi_lonae,gi_latae,gi_albd',go_lone,go_late,opt_debug);
            go_albd  = go_albd';
            go_falbd = go_falbd';
            disp(['       - Re-gridded GCM albedo data.']);
            % plot output albedo
            if (opt_plots), plot_2dgridded(flipud(go_albd),100.0,'',[[str_dirout '/' str_nameout] '.albd_out'],['albedo out']); end
            % save 2D file
            fprint_2DM(go_albd(:,:),[],[[str_dirout '/' str_nameout] '.2Dalbd.dat'],'%8.4f','%8.4f',true,false);
            fprintf('       - 2D albedo file saved\n')
            % create zonal mean
            vo_albd = mean(go_albd');
            disp(['       - Generated zonal mean albedo profile.']);
        otherwise
            % NOTE: if age == 0, then default (modern) GENIE,
            %       otherwise for generic ice-free world
            vo_albd = make_grid_albd(go_latm,par_age);
            disp(['       - Created generic zonal mean albedo profile.'])
            % plot output albedo profile
            figure;
            scatter(go_latm,vo_albd);
            axis([-90 90 0.0 1.0]);
            xlabel('Latitude');
            ylabel('Albedo');
            title('Zonally averaged albedo profile');
            print('-dpsc2', [[str_dirout '/' str_nameout] '.zonalalbd.' str_date '.ps']);
            % reorientate albedo vector for saving
            vo_albd = fliplr(vo_albd);
    end
    % save albedo vector
    % NOTE: when the file is read in by GENIE, it counrs down in j value:
    %       such that the N pole is the first element in file;
    %       fprint_1Dn saves the 2st row at the top, hence the vector
    %       must be orientated as in a map orientation (N at top)
    fprint_1Dn(vo_albd(:),[[str_dirout '/' str_nameout] '.albd.dat'],'%8.4f','%8.4f',true,false);
    fprintf('       - .albd.dat file saved\n')
    %
end
%
% *** GENERATE CONFIG FILE PARAMETER LINES ****************************** %
%
n_step = n_step+1;
%
disp(['>  ' num2str(n_step) '. GENERATING CONFIG FILE PARAMETER LINES ...']);
%
fid = fopen([str_dirout '/' 'config_' str_date '.txt'],'w');
% START
fprintf(fid,'%s\n','##################################################################################');
fprintf(fid,'%s\n',['### cGENIE .config file parameter lines generated by muffingen v',num2str(par_muffingen_ver),' on: ',str_date,' ###']);
% data path
fprintf(fid,'%s\n','# INPUT FILE PATH');
fprintf(fid,'%s\n',['ea_1=''../../cgenie.muffin/genie-paleo/',par_wor_name,'''']);
fprintf(fid,'%s\n',['go_1=''../../cgenie.muffin/genie-paleo/',par_wor_name,'''']);
fprintf(fid,'%s\n',['gs_1=''../../cgenie.muffin/genie-paleo/',par_wor_name,'''']);
% Grid resolution
fprintf(fid,'%s\n','# Grid resolution');
fprintf(fid,'%s\n',['GENIENXOPTS=''$(DEFINE)GENIENX=',num2str(imax),'''']);
fprintf(fid,'%s\n',['GENIENYOPTS=''$(DEFINE)GENIENY=',num2str(jmax),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLONSOPTS=''$(DEFINE)GOLDSTEINNLONS=',num2str(imax),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLATSOPTS=''$(DEFINE)GOLDSTEINNLATS=',num2str(jmax),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLEVSOPTS=''$(DEFINE)GOLDSTEINNLEVS=',num2str(kmax),'''']);
% Topography
fprintf(fid,'%s\n','# Topography');
fprintf(fid,'%s\n',['ma_fname_topo=''',par_wor_name,'''']);
% Assumed longitudinal offset of the grid
fprintf(fid,'%s\n','# Assumed longitudinal offset of the grid');
fprintf(fid,'%s\n',['gm_par_grid_lon_offset=',num2str(par_lon_off)]);
% Equal area?
if ~opt_equalarea,
    fprintf(fid,'%s\n',['ea_grid=1']);
    fprintf(fid,'%s\n',['go_grid=1']);
    fprintf(fid,'%s\n',['gs_grid=1']);
end
% Ocean depth scalar (dsc)
fprintf(fid,'%s\n','# Ocean depth scalar (m) [internally, parameter: dsc]');
fprintf(fid,'%s\n',['go_par_dsc=',num2str(par_max_D)]);
% Boundary conditions: EMBM
fprintf(fid,'%s\n','# Boundary conditions: EMBM');
fprintf(fid,'%s\n',['ea_topo=''',par_wor_name,'''']);
fprintf(fid,'%s\n',['ea_taux_u=''',par_wor_name,'.taux_u.dat''']);
fprintf(fid,'%s\n',['ea_tauy_u=''',par_wor_name,'.tauy_u.dat''']);
fprintf(fid,'%s\n',['ea_taux_v=''',par_wor_name,'.taux_v.dat''']);
fprintf(fid,'%s\n',['ea_tauy_v=''',par_wor_name,'.tauy_v.dat''']);
fprintf(fid,'%s\n',['ea_adv_u=''',par_wor_name,'.wvelx.dat''']);
fprintf(fid,'%s\n',['ea_adv_v=''',par_wor_name,'.wvely.dat''']);
% Boundary conditions: GOLDSTEIN
fprintf(fid,'%s\n','# Boundary conditions: GOLDSTEIN');
fprintf(fid,'%s\n',['go_topo=''',par_wor_name,'''']);
% Boundary conditions: GOLDSTEIN sea-ice
fprintf(fid,'%s\n','# Boundary conditions: GOLDSTEIN sea-ice');
fprintf(fid,'%s\n',['gs_topo=''',par_wor_name,'''']);
% Boundary conditions: ALBEDO!
if opt_makealbedo
    fprintf(fid,'%s\n','# Boundary conditions: ALBEDO!');
    fprintf(fid,'%s\n',['ea_par_albedo1d_name=''',par_wor_name,'.albd.dat''']);
end
% Boundary conditions: BIOGEM
if (opt_makezonalwind)
    fprintf(fid,'%s\n',['bg_ctrl_force_windspeed=.false']);
    fprintf(fid,'%s\n','# gas transfer coeff');
    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.722)]);
elseif (opt_makewind)
    % windspeed
    % NOTE: bg_ctrl_force_windspeed is .true. by default
    % NOTE: par_wspeed_avstr is the averaging product code
    switch str(1).gcm
        case {'hadcm3','hadcm3l','foam','cesm'}
            fprintf(fid,'%s\n','# Boundary conditions: BIOGEM');
            fprintf(fid,'%s\n',['bg_par_pindir_name=''../../cgenie.muffin/genie-paleo/',par_wor_name,'/''']);
            fprintf(fid,'%s\n',['bg_par_windspeed_file=''',par_wor_name,'.windspeed_' str(1).wspd '.dat''']);
        otherwise
            fprintf(fid,'%s\n',['bg_ctrl_force_windspeed=.false']);
    end
    % air-sea gas exchange
    % NOTE: re-scale to give a modern global mean air-sea coefficient of
    %       ~0.058 mol m-2 yr-1 uatm-1
    %       (default is bg_par_gastransfer_a=0.310)
    fprintf(fid,'%s\n','# BIOGEM MISC');
    switch str(1).gcm
        case {'hadcm3','hadcm3l'}
            fprintf(fid,'%s\n','# gas transfer coeff');
            fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.904)]);
        case {'foam'}
            fprintf(fid,'%s\n','# gas transfer coeff');
            fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(1.044)]);
        case {'cesm'}
            fprintf(fid,'%s\n','# gas transfer coeff');
            switch str(1).wspd
                case {'uvaa'}
                    % NOTE: @ 0.310 --> 0.027482 mol m-2 yr-1 uatm-1
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.6542)]);
                case {'uvma'}
                    % NOTE: @ 0.310 --> 0.031064 mol m-2 yr-1 uatm-1
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.5788)]);
                case {'wsma'}
                    % NOTE: @ 0.310 --> 0.051051 mol m-2 yr-1 uatm-1
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.3522)]);
            end
        otherwise
            % zonal field
            % NOTE: for now: don't distinguish between different zonal
            %       wind stress assumptions (and associated scaling)
            %       => take a modern-like tau profile as corresponding
            %          to ~0.058 mol m-2 yr-1 uatm-1 (and 0.310)
            %       @ 0.310, drakeworld gives 0.024903 mol m-2 yr-1 uatm-1
            %       (approx the mean of waterworld and ridgeworld values)
            %       => a = 0.722
            switch par_tauopt
                case {1}
                    % (low) modern NH / paleo Eocene (both hemispheres)
                    % (0.0201 mol m-2 yr-1 uatm-1 @ a=0.310)
                    fprintf(fid,'%s\n','# gas transfer coeff');
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.722)]);
                case {2}
                    % (high) water world
                    % (0.0297 mol m-2 yr-1 uatm-1 @ a=0.310)
                    fprintf(fid,'%s\n','# gas transfer coeff');
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.722)]);
                otherwise
                    % intermediate case
                    % (0.024903 mol m-2 yr-1 uatm-1 @ a=0.310)
                    % NOTE: this is the 'automatically' determined case
                    fprintf(fid,'%s\n','# gas transfer coeff');
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.722)]);
            end
    end
end
% SEDGEM/ROKGEM
if opt_makeseds
    % Grid resolution of solid Earth components
    fprintf(fid,'%s\n','# Grid resolution of solid Earth components');
    fprintf(fid,'%s\n',['SEDGEMNLONSOPTS=''$(DEFINE)SEDGEMNLONS=',num2str(imax),'''']);
    fprintf(fid,'%s\n',['SEDGEMNLATSOPTS=''$(DEFINE)SEDGEMNLATS=',num2str(jmax),'''']);
    fprintf(fid,'%s\n',['ROKGEMNLONSOPTS=''$(DEFINE)ROKGEMNLONS=',num2str(imax),'''']);
    fprintf(fid,'%s\n',['ROKGEMNLATSOPTS=''$(DEFINE)ROKGEMNLATS=',num2str(jmax),'''']);
    % Topography for solid Earth components
    fprintf(fid,'%s\n','# Topography for solid Earth components');
    fprintf(fid,'%s\n',['sg_par_pindir_name=''../../cgenie.muffin/genie-paleo/',par_wor_name,'/''']);
    fprintf(fid,'%s\n',['sg_par_sed_topo_D_name=''',par_wor_name,'.depth.dat''']);
    fprintf(fid,'%s\n',['sg_par_sed_reef_mask_name=''',par_wor_name,'.reefmask.dat''']);
    fprintf(fid,'%s\n',['sg_par_sedcore_save_mask_name=''',par_wor_name,'.sedcoremask.dat''']);
    fprintf(fid,'%s\n',['rg_par_pindir_name=''../../cgenie.muffin/genie-paleo/',par_wor_name,'/''']);
    fprintf(fid,'%s\n',['rg_topo=''',par_wor_name,'.k1''']);
end
% GEOLOGIC AGE DEPENDENT PARAMETERS
fprintf(fid,'%s\n','# GEOLOGIC AGE DEPENDENT PARAMETERS');
if (par_age == 0.0),
    fprintf(fid,'%s\n','# Solar constant (W m-2) ... don''t forget to adjust it if not modern!!');
    fprintf(fid,'%s\n',['###ma_genie_solar_constant=','1368.0']);
    fprintf(fid,'%s\n','# ... also, salinity should be set 1 PSU lower if it an ice-free World');
    fprintf(fid,'%s\n',['###go_saln0=33.9']);
    fprintf(fid,'%s\n','# Orbital parameters (modern, defaults)');
    fprintf(fid,'%s\n',['###ea_par_orbit_osce=','0.0167',' # eccentricity']);
    fprintf(fid,'%s\n',['###ea_par_orbit_oscsob=','0.397789',' # sine of obliquity']);
    fprintf(fid,'%s\n',['###ea_par_orbit_oscgam=','102.92',' # longitude of perihelion']);
else
    loc_per = 100.0*(1-1/(1+(2/5)*(1-(4.570E+03-par_age)/4.570E+03)));
    loc_S0  = 1.368E+03*(100-loc_per)/100;
    fprintf(fid,'%s\n',['# Solar constant (W m-2) -- reduced (from 1368.0) by ', num2str(loc_per), '% corresponding to age ', num2str(par_age), ' Ma']);
    fprintf(fid,'%s\n',['ma_genie_solar_constant=',num2str(loc_S0)]);
    fprintf(fid,'%s\n','# Ocean salinity -- assuming an ice-free World (1 PSU lower than modern)');
    fprintf(fid,'%s\n',['go_saln0=33.9']);
    fprintf(fid,'%s\n','# Orbital parameters -- modern set => adjust as necessary');
    fprintf(fid,'%s\n',['ea_par_orbit_osce=','0.0167',' # eccentricity']);
    fprintf(fid,'%s\n',['ea_par_orbit_oscsob=','0.397789',' # sine of obliquity']);
    fprintf(fid,'%s\n',['ea_par_orbit_oscgam=','102.92',' # longitude of perihelion']);
end
if (par_age > 0.0) && (par_age <= 100.0),
    loc_Ca = 1E-3*(1.028E-02*1000 - 0.1966*(-par_age) - 0.001116*(-par_age)^2 - 0.000003374*(-par_age)^3 - 0.000000006584*(-par_age)^4);
    loc_Mg = 1E-3*(5.282E-02*1000 + 0.915*(-par_age) + 0.01308*(-par_age)^2 + 0.00008419*(-par_age)^3 + 0.000000201*(-par_age)^4);
    fprintf(fid,'%s\n','# Ocean Mg/Ca concentrations');
    fprintf(fid,'%s\n',['bg_ocn_init_35=',num2str(loc_Ca)]);
    fprintf(fid,'%s\n',['bg_ocn_init_50=',num2str(loc_Mg)]);
else
    fprintf(fid,'%s\n','# Ocean Mg/Ca concentrations (modern defaults, mol kg-1)');
    fprintf(fid,'%s\n',['bg_ocn_init_35=','1.028E-02']);
    fprintf(fid,'%s\n',['bg_ocn_init_50=','5.282E-02']);
end
% END
fprintf(fid,'%s\n','##################################################################################');
fclose(fid);
fprintf('       - .config file saved\n')
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
% final messages
disp([' ']);
disp(['------------------------------------------------------------']);
disp(['   Congratulations! SOMETHING was created ... ']);
disp(['   ... hope it was what you wished for! :o)']);
disp(['------------------------------------------------------------']);
disp([' ']);
% end logging
diary off;
% clean up
close all;
%
disp(['<<< END']);
disp([' ']);
%
% *********************************************************************** %
