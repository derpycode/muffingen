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
%
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
par_muffingen_ver = 0.61;
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
% *** check / filter options ******************************************** %
%
% zonal wind-stress generaton parameter
if ~exist('par_tauopt','var'),
    par_tauopt = 0;
end
% age parameter
if ~exist('par_age','var'),
    par_age = 0.0;
    par_age_emty = true;
else
    par_age_emty = false;
end
% process GCM name string
if strcmp(par_gcm,'hadcm3l'), par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HadCM3'),  par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HadCM3L'), par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HADCM3'),  par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'HADCM3L'), par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'um'),      par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'UM'),      par_gcm = 'hadcm3'; end
if strcmp(par_gcm,'FOAM'),    par_gcm = 'foam';   end
if strcmp(par_gcm,'K1'),      par_gcm = 'k1';     end
if strcmp(par_gcm,'.k1'),     par_gcm = 'k1';     end
if strcmp(par_gcm,'.K1'),     par_gcm = 'k1';     end
if strcmp(par_gcm,'MASK'),    par_gcm = 'mask';   end
if strcmp(par_gcm,'dat'),     par_gcm = 'mask';   end
if strcmp(par_gcm,'.dat'),    par_gcm = 'mask';   end
if strcmp(par_gcm,''),        par_gcm = 'blank';  end
if strcmp(par_gcm,'BLANK'),   par_gcm = 'blank';  end
if strcmp(par_gcm,'none'),    par_gcm = 'blank';  end
if strcmp(par_gcm,'NONE'),    par_gcm = 'blank';  end
% adjust options accroding to input (GCM) type
switch par_gcm
    case {'hadcm3','foam'}
    case {'k1','mask'}
    otherwise
        opt_makeall=false;
        opt_user=true;
end
% deal with meta selections
if opt_makeall,
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
if opt_makeseds,
    opt_maketopo=true;
end
% rename variables
% NOTE: for some earlier code consistency ... now redundant ...
imax = par_max_i;
jmax = par_max_j;
kmax = par_max_k;
str_nameout = par_wor_name;
% initialize optional netCDF filenames
% -> set default variable names
switch par_gcm
    case ('hadcm3')
        if ~exist('par_nc_topo_name','var'), par_nc_topo_name = [par_expid '.qrparm.omask']; end
        if ~exist('par_nc_axes_name','var'), par_nc_axes_name = [par_expid 'a.pdclann']; end
        if ~exist('par_nc_atmos_name','var'), par_nc_atmos_name = [par_expid '_sed']; end
        if ~exist('par_nc_ocean_name','var'), par_nc_ocean_name = ''; end
        if ~exist('par_nc_coupl_name','var'), par_nc_coupl_name = [par_expid 'a.pdclann']; end
    case ('foam')
        if ~exist('par_nc_topo_name','var'), par_nc_topo_name = 'topo'; end
        if ~exist('par_nc_axes_name','var'), par_nc_axes_name = par_nc_topo_name; end
        if ~exist('par_nc_atmos_name','var'), par_nc_atmos_name = 'atmos'; end
        if ~exist('par_nc_ocean_name','var'), par_nc_ocean_name = ''; end
        if ~exist('par_nc_coupl_name','var'), par_nc_coupl_name = ''; end
    otherwise
        par_nc_topo_name  = '';
        par_nc_axes_name  = '';
        par_nc_atmos_name = '';
        par_nc_ocean_name = '';
        par_nc_coupl_name = '';
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
str = setfield(str, {1}, 'nc', par_nc_axes_name);
str = setfield(str, {2}, 'nc', par_nc_topo_name);
str = setfield(str, {3}, 'nc', par_nc_atmos_name);
str = setfield(str, {4}, 'nc', par_nc_ocean_name);
str = setfield(str, {5}, 'nc', par_nc_coupl_name);
%
% *** initialize reporting ********************************************** %
%
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
% (0) (initialize)
% (1) CONFIRM OPTIONS
% (2) CREATING GENIE GRID
% (3) READ AXES INFORMATION
% (4) LOAD TOPO & MASK DATA
% (5) RE-GRID MASK
% (6) FILTER MASK
% (7) ADJUST MASK -- USER!
% (8) RE-GRID TOPO
% (9) RE-GRID VERTICALLY
% (10) ADJUST TOPO -- AUTOMATIC BATHYMETRY FILTERING
% (11) ADJUST TOPO -- USER!
% (12) CALCULATE RUNOFF & COMPLETE k1 FILE
% (13) IDENTIFY ISLANDS
% (14) UPDATE ISLANDS AND ISLAND PATHS
% (15) GENERATE ISLAND PATHS
% (16) GENERATE PSI ISLANDS
% (17) RE-GRID WIND SPEED/STRESS DATA
% (18) LOAD ALBEDO DATA
% (19) RE-GRID & PROCESS ALBEDO
% (20) GENERATE SEDIMENT GRID
% (21) GENERATE CONFIG FILE PARAMETER LINES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *** (1) CONFIRM OPTIONS *********************************************** %
%
disp(['>   1. CHECKING PRIMARY OPTIONS ...']);
% check world name
if (length(par_wor_name) ~= 8),
    disp(['       * ERROR: World name (par_wor_name) must be 8 characters long.']);
    disp(['--------------------------------------------------------']);
    disp([' ']);
    diary off;
    return;
end
% check GCM options
switch str(1).gcm
    case ('hadcm3')
        disp(['       * GCM == ' str(1).gcm ' (OK)']);
    case ('foam')
        disp(['       * GCM == ' str(1).gcm ' (OK)']);
    case {'k1','mask'}
        disp(['       * GENIE grid will be loaded directly from k1 or mask text file: ' str(1).exp]);
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
% *** (2) SET UP OUTPUT GRID ******************************************** %
%
disp(['>   2. CREATING GENIE GRID ...']);
% create GENIE grid
[go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea);
disp(['       - GENIE grid generated.']);
%
% *** (3) LOAD GRID (AXES) DATA ***************************************** %
%
disp(['>   3. READING AXES INFORMATION ...']);
%
switch par_gcm
    case {'hadcm3','foam'}
        % read axes
        if strcmp(par_gcm,'hadcm3')
            [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_hadcm3(str);
        else
            [gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_lonam,gi_lonae,gi_latam,gi_latae] = fun_read_axes_foam(str);
        end
        disp(['       - Axis info read.']);
    otherwise
        % DO NOTHING
        disp(['         (Nothing to load.)']);
end
%
% *** (4) LOAD TOPO & MASK DATA ***************************************** %
%
% NOTE: mask is defined with '1' for ocean
%
disp(['>   4. READING MASK & TOPO GRIDS ...']);
%
switch par_gcm
    case {'hadcm3','foam'}
        % read topo
        if strcmp(par_gcm,'hadcm3')
            [gi_topo,gi_mask] = fun_read_topo_hadcm3(str);
        else
            [gi_topo,gi_mask] = fun_read_topo_foam(str);
        end
        disp(['       - Mask & topo info read.']);
        % plot input mask & topo
        plot_2dgridded(flipud(gi_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_in'],['mask in']);
        plot_2dgridded(flipud(gi_topo),6000.0,'',[[str_dirout '/' str_nameout] '.topo_in'],['topo in']);
    case {'k1','mask'}
        % load topo directly
        [go_k1,go_mask,imax,jmax] = fun_read_k1(str);
        disp(['       - k1 read.']);
        % re-create GENIE grid with derived (/updated?) grid dimensions
        % NOTE: the value of kmax is taken from the config file
        %      (while imax and jmax are deduced from the file)
        [go_lonm,go_lone,go_latm,go_late,go_dm,go_de] = make_genie_grid(imax,jmax,kmax,par_max_D,par_lon_off,opt_equalarea);
        disp(['       - GENIE grid re-generated.']);
    otherwise
        % DO NOTHING
        disp(['         (Nothing to load.)']);
end
%
% *** (5) RE-GRID MASK ************************************************** %
%
disp(['>   5. RE-GRIDING MASK ...']);
%
switch par_gcm
    case {'hadcm3','foam'}
        % initial re-gridding of mask
        % NOTE: need to transpose around [gi_mask] to have correct input format
        %       to make_regrid_2d
        %       similarly, output needs to be transposed back again
        % NOTE: pass edges of c-grid
        [go_mask,go_fmask] = make_regrid_2d(gi_lonce,gi_latce,gi_mask',go_lone,go_late,false);
        disp(['       - Mask re-gridded.']);
        go_mask = go_mask';
        go_fmask = go_fmask';
        % create mask (<>= par_A_frac_threshold fractional area thresold)
        go_mask(find(go_mask>=par_A_frac_threshold)) = 1.0;
        go_mask(find(go_mask<par_A_frac_threshold))  = 0.0;
        % calculate respective fractional areas
        [si_farea,si_farearef] = fun_grid_calc_ftotarea(gi_mask,gi_lonce,gi_latce);
        [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
        disp(['       * Original land area fraction    = ', num2str(1.0-si_farea)]);
        disp(['       * Re-gridded land area fraction  = ', num2str(1.0-so_farea)]);
    case {'k1','mask'}
        disp(['         (Nothing to do ... k1/mask file already loaded.)']);
    otherwise
        go_mask = zeros(jmax,imax);
        go_mask = go_mask + 1;
        disp(['       - Blank mask created (nothing to re-grid).']);
end
% plot & save initial mask re-grid
plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.RAW'],['mask out -- RAW re-gridded']);
%
% *** (6) FILTER MASK *************************************************** %
%
% filter mask if requested
% NOTE: when loading in a default 'k1' file, best to skip this step
% set VERSION 0 (raw)
grid_ver = 0;
str_ver = num2str(grid_ver);
%
if opt_makemask && (opt_filtermask || (par_min_oceann > 0)),
    %
    disp(['>   6. FILTERING MASK ...']);
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
        plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]);
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
        plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]);
        %
    end
    %
    if (par_min_oceann > 0),
        %
        % SMALL WATER BODY FILTERING MASK FILTERING
        %
        [go_oceans,n_oceans,i_oceans] = find_grid_oceans(go_mask);
        % plot oceans!
        plot_2dgridded(flipud(go_oceans),999,'',[[str_dirout '/' str_nameout] '.ocean_out.INIT'],['oceans out -- INITIAL']);
        % increment VERSION
        grid_ver = grid_ver + 1;
        str_ver = num2str(grid_ver);
        % clean up small water bodies
        [go_mask,go_oceans,n_oceans] = find_grid_oceans_update(go_mask,go_oceans,n_oceans,par_min_oceann);
        fprintf('       - Small water bodies cleaned up.\n')
        % plot mask
        plot_2dgridded(flipud(go_mask),2.0,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]);
        %
    end
    %
    % calculate new fractional area
    [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
    disp(['       * Revised land area fraction = ', num2str(1.0-so_farea)]);
    %
end
%
% *** (7) ADJUST MASK -- USER! ****************************************** %
%
if opt_user
    %
    disp(['>   7. USER EDITING OF MASK ...']);
    %
    [go_mask]  = fun_grid_edit_mask(go_mask);
    % increment VERSION
    grid_ver = grid_ver + 1;
    str_ver = num2str(grid_ver);
    % plot mask
    plot_2dgridded(flipud(go_mask),2,'',[[str_dirout '/' str_nameout] '.mask_out.v' str_ver],['mask out -- version ' str_ver]);
    % calculate new fractional area
    [so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
    disp(['       * Revised land area fraction = ', num2str(1.0-so_farea)]);
    %
    fprintf('       - User-editing complete.\n')
    %
end
%
% create GENIE NaN mask
go_masknan = go_mask;
go_masknan(find(go_masknan == 0)) = NaN;
%
% plot final mask
plot_2dgridded(flipud(go_mask),99999.0,'',[[str_dirout '/' str_nameout] '.mask_out.FINAL'],['mask out -- FINAL version']);
%
% calculate new fractional area
[so_farea,so_farearef] = fun_grid_calc_ftotarea(go_mask,go_lone,go_late);
disp(['       * Final land area fraction   = ', num2str(1.0-so_farea)]);
%
% *** (8) RE-GRID TOPO ************************************************** %
%
if opt_maketopo
    %
    disp(['>   8. RE-GRIDING TOPOGRAPHY ...']);
    %
    switch par_gcm
        case {'hadcm3','foam'}
            % initial re-gridding of topo
            % NOTE: need to transpose around [gi_topo] to have correct input format
            %       to make_regrid_2d
            %       similarly, output needs to be transposed back again
            % NOTE: pass edges of c-grid
            [go_topo,go_ftopo] = make_regrid_2d(gi_lonce,gi_latce,gi_topo',go_lone,go_late,false);
            disp(['       - Topography re-gridded.']);
            go_topo = go_topo';
            go_ftopo = go_ftopo';
        case {'k1'}
            % convert k1 to depth
            [go_topo] = fun_conv_k1(go_de,go_k1);
            disp(['         (Nothing to re-grid -- convert k1 file data.)']);
        otherwise
            go_topo = par_max_D*go_mask;
            disp(['         (Nothing to re-grid -- set uniform ocean depth.)']);
    end
    % plot & save initial topo re-grid
    if ~strcmp(par_gcm,'k1'),
        plot_2dgridded(flipud(go_topo),99999.0,'',[[str_dirout '/' str_nameout] '.topo_out.RAW'],['topo out -- RAW']);
    end
    %
end
%
% *** (9) RE-GRID VERTICALLY ******************************************** %
%
if opt_maketopo,
    %
    disp(['>   9. RE-GRIDING OCEAN BATHYMETRY ...']);
    %
    switch par_gcm
        case {'k1'}
            disp(['         (Nothing to re-grid as k1 file already loaded.)']);
        otherwise
            % convert depth into k levels (and create k1 grid)
            [go_k1] = find_grid_k(par_min_Dk,go_dm,go_de,go_mask,go_topo);
            fprintf('       - Bathymetry re-gridding complete.\n')
    end
    % filter min k value
    go_k1(find(go_k1 < par_min_k)) = par_min_k;
    % plot initial k1 re-grid
    plot_2dgridded(flipud(go_k1),89.0,'',[[str_dirout '/' str_nameout] '.k1_out.RAW'],['k1 out -- RAW re-gridded']);
    %
end
%
% *** (10) ADJUST TOPO -- AUTOMATIC BATHYMETRY FILTERING **************** %
%
% carry out basic automatic topo filtering
if opt_maketopo && opt_filtertopo,
    %
    disp(['>  10. FILTERING BATHYMETRY ...']);
    %
    [go_k1] = fun_grid_topo_filter(go_k1);
    fprintf('       - Topography filtered.\n')
    % plot adjusted k1 re-grid
    plot_2dgridded(flipud(go_k1),89.0,'',[[str_dirout '/' str_nameout] '.k1_out.FILTERED'],['k1 out -- auto filtered']);
end
%
% *** (11) ADJUST TOPO -- USER! ***************************************** %
%
if opt_maketopo && opt_user
    %
    disp(['>  11. USER EDITING OF TOPO ...']);
    % user-editing! what can go wrong?
    [go_k1] = fun_grid_edit_k1(go_k1,kmax);
    % plot mask
    plot_2dgridded(flipud(go_k1),89.0,'',[str_dirout '/' str_nameout '.k1_out.USEREDITED'],['k1 out -- user edited version']);
    % convert k-levels back to depth
    [go_topo] = fun_conv_k1(go_de,go_k1);
    %
    fprintf('       - User-editing complete.\n')
    %
end
%
if opt_maketopo,
    % plot final k1
    plot_2dgridded(flipud(go_k1),89.0,'',[str_dirout '/' str_nameout '.k1_out.FINAL'],['k1 out -- FINAL version']);
    % plot final topo
    plot_2dgridded(flipud(go_masknan.*go_topo),99999.0,'',[[str_dirout '/' str_nameout] '.topo_out.FINAL'],['topo out -- FINAL version']);
end
%
% *** (12) CALCULATE RUNOFF & COMPLETE k1 FILE ************************** %
%
% NOTE: ordering is a little illogical becasue
%       make_grid_runoff_rnd requires the extended grid, while
%       make_grid_runoff_roof is easier done without ...
if opt_makeocean
    %
    disp(['>  12. CALCULATING RUN-OFF AND GENERATE .k1 FILE ...']);
    % (i) first, check for all ocean
    if (max(max(go_k1)) < 90), opt_makerunoff = false; end
    % (ii) create roof scheme (if selected)
    if (opt_makerunoff && par_runoffopt == 0),
        [go_k1] = make_grid_runoff_roof(go_mask,go_k1,str);
        loc_k1 = go_k1;
        loc_k1(find(loc_k1 < 91)) = 95;
        plot_2dgridded(flipud(loc_k1),95.0,'',[str_dirout '/' str_nameout '.k1_out.RUNOFF'],['k1 out -- RUNOFF']);
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
        plot_2dgridded(flipud(loc_k1),95.0,'',[str_dirout '/' str_nameout '.k1_out.RUNOFF'],['k1 out -- RUNOFF']);
    end
    % (v) save .k1 file
    fprint_2DM(goex_k1(:,:),[],[[str_dirout '/' str_nameout] '.k1'],'%3i','%3i',true,false);
    fprintf('       - .k1 file saved\n')
    %
end
%
% *** (13) IDENTIFY ISLANDS ********************************************* %
%
if opt_makeocean
    %
    disp(['>  13. IDENTIFY ISLANDS ...']);
    % initial islands count
    [go_islands,n_islands,i_islands] = find_grid_islands(go_mask);
    % plot islands
    plot_2dgridded(flipud(go_islands),999,'',[[str_dirout '/' str_nameout] '.islnd_out.INIT'],['island out -- INITIAL']);
    %
end
%
% *** (14) UPDATE ISLANDS AND ISLAND PATHS ****************************** %
%
if opt_makeocean
    %
    disp(['>  14. UPDATING ISLANDS & PATHS ...']);
    % NOTE: generate all possible paths initially (and filter later)
    % (1) generate generic borders around all (initial) islands
    [go_borders] = find_grid_borders(go_mask);
    plot_2dgridded(flipud(go_borders),99999.0,'',[[str_dirout '/' str_nameout] '.brds_out.INIT'],['borders out -- INITIAL']);
    % (2) update islands count
    %     identify islands that are insufficiently seperated (and combined)
    %     identify polar islands
    %     re-number all
    [go_islands,n_islands,i_islands,i_poles] = find_grid_islands_update(go_islands,n_islands,i_islands,go_borders,opt_makepoleswide);
    plot_2dgridded(flipud(go_islands),999,'',[[str_dirout '/' str_nameout] '.islnd_out.FINAL'],['islands out -- FINAL']);
    % (3) update borders
    %     number borders as per bordering islands
    [go_borders] = find_grid_borders_update(go_borders,go_islands,go_mask,n_islands);
    plot_2dgridded(flipud(go_borders),999,'',[[str_dirout '/' str_nameout] '.brds_out.FILTERED'],['borders out -- FILTERED']);
    % (4) border check
    [opt_user] = find_grid_borders_check(go_borders,opt_user);
    % (5) user editing of borders
    if opt_user
        % user-editing! what can go wrong?
        [go_borders] = fun_grid_edit_borders(go_borders,go_mask);
        % plot mask
        plot_2dgridded(flipud(go_borders),999,'',[str_dirout '/' str_nameout '.brds_out.USEREDITED'],['borders out -- user edited version']);
    end
    % plot final borders
    plot_2dgridded(flipud(go_borders),999,'',[[str_dirout '/' str_nameout] '.brds_out.FINAL'],['borders out -- FINAL']);
    %
end
%
% *** (15) GENERATE ISLAND PATHS **************************************** %
%
if opt_makeocean
    %
    disp(['>  15. GENERATING .paths FILE ...']);
    % create paths
    [n_paths,v_paths,n_islands,go_paths] = find_grid_paths(go_borders,n_islands,i_poles);
    % plot paths data
    plot_2dgridded(flipud(go_paths),999,'',[[str_dirout '/' str_nameout] '.paths_out.FINAL'],['Paths file -- FINAL']);
    % save .paths file
    fprint_paths(n_paths,v_paths,[[str_dirout '/' str_nameout] '.paths']);
    fprintf('       - .paths file saved\n')
    %
end
%
% *** (16) GENERATE PSI ISLANDS ***************************************** %
%
if opt_makeocean
    %
    disp(['>  16. GENERATING .psiles FILE ...']);
    % generate PSI islands data
    [go_psiles,n_islands_recnt] = make_grid_psiles(go_islands,i_poles);
    % plot PSI islands data
    plot_2dgridded(flipud(go_psiles),999,'',[[str_dirout '/' str_nameout] '.psiles_out.FINAL'],['PSI islands file -- FINAL']);
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
% *** (17) RE-GRID WIND SPEED/STRESS DATA ******************************* %
%
if opt_makewind
    %
    disp(['>  17. CREATING WIND PRODUCTS ...']);
    % create GENIE grid wind products
    switch par_gcm
        case {'hadcm3','foam'}
            % re-grid winds from GCM
            % NOTE: the sets of grids and their edges required differ
            %       between hadcm3 and foam
            if strcmp(par_gcm,'hadcm3')
                [wstr,wspd,g_wspd] = make_grid_winds_hadcm3(gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonpm,gi_lonpe,gi_latpm,gi_latpe,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str);
            else
                [wstr,wspd,g_wspd] = make_grid_winds_foam(gi_loncm,gi_lonce,gi_latcm,gi_latce,gi_lonam,gi_lonae,gi_latam,gi_latae,gi_mask,go_lonm,go_lone,go_latm,go_late,go_mask,str);
            end
            disp(['       - Re-grided GCM wind products.']);
        otherwise
            [wstr,wspd,g_wspd] = make_grid_winds_zonal(go_latm,go_late,go_mask,[str_dirout '/' str_nameout],par_tauopt);
            disp(['       - Generated zonal wind products.']);
    end
end
%
% *** (18) LOAD ALBEDO DATA ********************************************* %
%
if opt_makealbedo
    %
    disp(['>  18. LOADING ALBEDO DATA ...']);
    %
    switch par_gcm
        case {'hadcm3','foam'}
            % read albedo
            if strcmp(par_gcm,'hadcm3')
                [gi_albd] = fun_read_albd_hadcm3(str);
            else
                [gi_albd] = fun_read_albd_foam(str);
            end
            disp(['       - Read GCM albedo data.']);
            % plot input albedo
            plot_2dgridded(flipud(gi_albd),100.0,'',[[str_dirout '/' str_nameout] '.albd_in'],['albedo in']);
        otherwise
            disp(['         (Nothing to load.)']);
    end
    %
end
%
% *** (19) RE-GRID & PROCESS ALBEDO ************************************* %
%
if opt_makealbedo
    %
    disp(['>  19. CREATING ALBEDO DATA ...']);
    %
    switch par_gcm
        case {'hadcm3','foam'}
            % re-grid
            [go_albd,go_falbd] = make_regrid_2d(gi_lonae,gi_latae,gi_albd',go_lone,go_late,false);
            go_albd  = go_albd';
            go_falbd = go_falbd';
            disp(['       - Re-gridded GCM albedo data.']);
            % plot output albedo
            plot_2dgridded(flipud(go_albd),100.0,'',[[str_dirout '/' str_nameout] '.albd_out'],['albedo out']);
            % save 2D file
            fprint_2DM(go_albd(:,:),[],[[str_dirout '/' str_nameout] '.2Dalbd.dat'],'%8.4f','%8.4f',true,false);
            fprintf('       - 2D albedo file saved\n')
            % create zonal mean
            vo_albd = mean(go_albd');
            disp(['       - Generated zonal mean albedo profile.']);
        otherwise;
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
% *** (20) GENERATE SEDIMENT GRID *************************************** %
%
if opt_makeseds
    %
    disp(['>  20. GENERATING SEDIMENT TOPO ...']);
    %
    % check sed topo re-gridding options
    switch par_sedsopt
        case 1
            % option %1 -- re-grid sediment topography
            switch par_gcm
                case {'hadcm3','foam'}
                    % if 'high res' sed grid is requested => assume twice ocean resolution
                    % + generate new vectors of grid properties
                    if opt_highresseds,
                        [gos_lonm,gos_lone,gos_latm,gos_late,gos_dm,gos_de] = make_genie_grid(2*imax,2*jmax,kmax,par_max_D,par_lon_off,opt_equalarea);
                    else
                        gos_lone = go_lone;
                        gos_late = go_late;
                    end
                    % re-grid
                    [gos_topo,gos_ftopo] = make_regrid_2d(gi_lonce,gi_latce,gi_topo',gos_lone,gos_late,false);
                    gos_topo  = gos_topo';
                    gos_ftopo = gos_ftopo';
                    disp(['       - Re-gridded sediment topo from GCM bathymetry.']);
                case {'k1'}
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
    % plot sediment topo
    plot_2dgridded(flipud(gos_topo),9999,'',[[str_dirout '/' str_nameout] '.sedtopo_out.FINAL'],['Sediment topo -- FINAL']);
    % save sediment topo
    fprint_2DM(gos_topo(:,:),[],[[str_dirout '/' str_nameout] '.depth.dat'],'%8.2f','%8.2f',true,false);
    fprintf('       - .depth.dat saved\n')
    % save other sediment files
    gos_mask = gos_topo;
    gos_mask(find(gos_mask > 0)) = 1;
    gos_sedc = 0.0*gos_mask;
    fprint_2DM(gos_sedc(:,:),gos_mask(:,:),[[str_dirout '/' str_nameout] '.sedcoremask.dat'],'%4.1f','%4i',true,false);
    fprintf('       - template file .sedcoremask.dat saved\n')
    gos_reef = 0.0*gos_mask;
    fprint_2DM(gos_reef(:,:),gos_mask(:,:),[[str_dirout '/' str_nameout] '.reefmask.dat'],'%4.1f','%4i',true,false);
    fprintf('       - template file .reefmask.dat saved\n')
end
%
% *** (21) GENERATE CONFIG FILE PARAMETER LINES ************************* %
%
disp(['>  21. GENERATING CONFIG FILE PARAMETER LINES ...']);
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
fprintf(fid,'%s\n',['GENIENXOPTS=''$(DEFINE)GENIENX=',num2str(par_max_i),'''']);
fprintf(fid,'%s\n',['GENIENYOPTS=''$(DEFINE)GENIENY=',num2str(par_max_j),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLONSOPTS=''$(DEFINE)GOLDSTEINNLONS=',num2str(par_max_i),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLATSOPTS=''$(DEFINE)GOLDSTEINNLATS=',num2str(par_max_j),'''']);
fprintf(fid,'%s\n',['GOLDSTEINNLEVSOPTS=''$(DEFINE)GOLDSTEINNLEVS=',num2str(par_max_k),'''']);
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
% Boundary conditions: EMBM
fprintf(fid,'%s\n','# Boundary conditions: EMBM');
fprintf(fid,'%s\n',['ea_topo=''',par_wor_name,'''']);
fprintf(fid,'%s\n',['ea_taux_u=''',par_wor_name,'_taux_u.dat''']);
fprintf(fid,'%s\n',['ea_tauy_u=''',par_wor_name,'_tauy_u.dat''']);
fprintf(fid,'%s\n',['ea_taux_v=''',par_wor_name,'_taux_v.dat''']);
fprintf(fid,'%s\n',['ea_tauy_v=''',par_wor_name,'_tauy_v.dat''']);
fprintf(fid,'%s\n',['ea_adv_u=''',par_wor_name,'_wvelx.dat''']);
fprintf(fid,'%s\n',['ea_adv_v=''',par_wor_name,'_wvely.dat''']);
% Boundary conditions: GOLDSTEIN
fprintf(fid,'%s\n','# Boundary conditions: GOLDSTEIN');
fprintf(fid,'%s\n',['go_topo=''',par_wor_name,'''']);
% Boundary conditions: GOLDSTEIN sea-ice
fprintf(fid,'%s\n','# Boundary conditions: GOLDSTEIN sea-ice');
fprintf(fid,'%s\n',['gs_topo=''',par_wor_name,'''']);
% Boundary conditions: ALBEDO!
if opt_makealbedo,
    fprintf(fid,'%s\n','# Boundary conditions: ALBEDO!');
    fprintf(fid,'%s\n',['ea_par_albedo1d_name=''',par_wor_name,'.albd.dat''']);
end
% Boundary conditions: BIOGEM
if opt_makewind
    % windspeed
    % NOTE: bg_ctrl_force_windspeed is .true. by default
    switch par_gcm
        case {'hadcm3','foam'}
            fprintf(fid,'%s\n','# Boundary conditions: BIOGEM');
            fprintf(fid,'%s\n',['bg_par_pindir_name=''../../cgenie.muffin/genie-paleo/',par_wor_name,'/''']);
            fprintf(fid,'%s\n',['bg_par_windspeed_file=''',par_wor_name,'_windspeed.dat''']);
        otherwise
            fprintf(fid,'%s\n',['bg_ctrl_force_windspeed=.false']);
    end
    % air-sea gas exchange
    % NOTE: re-scale to give a modern global mean air-sea coefficient of
    %       ~0.058 mol m-2 yr-1 uatm-1
    %       (default is bg_par_gastransfer_a=0.310)
    fprintf(fid,'%s\n','# BIOGEM MISC');
    switch par_gcm
        case {'hadcm3'}
            fprintf(fid,'%s\n','# gas transfer coeff');
            fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.904)]);
        case {'foam'}
            fprintf(fid,'%s\n','# gas transfer coeff');
            fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(1.044)]);
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
                    fprintf(fid,'%s\n','# gas transfer coeff');
                    fprintf(fid,'%s\n',['bg_par_gastransfer_a=',num2str(0.722)]);
            end
    end
end
% SEDGEM/ROKGEM
if opt_makeseds
    % Grid resolution of solid Earth components
    fprintf(fid,'%s\n','# Grid resolution of solid Earth components');
    fprintf(fid,'%s\n',['SEDGEMNLONSOPTS=''$(DEFINE)SEDGEMNLONS=',num2str(par_max_i),'''']);
    fprintf(fid,'%s\n',['SEDGEMNLATSOPTS=''$(DEFINE)SEDGEMNLATS=',num2str(par_max_j),'''']);
    fprintf(fid,'%s\n',['ROKGEMNLONSOPTS=''$(DEFINE)ROKGEMNLONS=',num2str(par_max_i),'''']);
    fprintf(fid,'%s\n',['ROKGEMNLATSOPTS=''$(DEFINE)ROKGEMNLATS=',num2str(par_max_j),'''']);
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
