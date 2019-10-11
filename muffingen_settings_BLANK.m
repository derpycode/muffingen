% muffingen_settings
%
%   ***********************************************************************
%   *** PARAMETER SETTINGS FOR muffingen CONFIG GENERATOR *****************
%   ***********************************************************************
%
%   Edit this file directly for additional user settings
%
%   ***********************************************************************

% *********************************************************************** %
% *** USER SETTINGS ***************************************************** %
% *********************************************************************** %
%
% PARAM NAME & DEFAULT VALUE   % [FORMAT] BRIEF DESCRIPTION
%
% *** CONFIG NAME AND MAIN INPUT SETTINGS ******************************* %
%
par_wor_name='worblank';       % ['worblank'] 8-char (output) config name
par_gcm='';                    % [''] input format
par_expid='';                  % [''] input experiment/data name
%
% *** FILE PATHS ******************************************************** %
%
par_pathin='INPUT.EXAMPLES';   % ['EXAMPLES.INPUT'] path to input dir
par_pathout='OUTPUT.EXAMPLES'; % ['EXAMPLES.OUTPUT'] path to output dir
opt_outputdir=false;           % [false/true] ask for output directory?
%
% *** GCM netCDF FILENAMES ********************************************** %
%
par_nc_topo_name  = '';        % [''] optional .nc name
par_nc_mask_name  = '';        % [''] optional .nc name
par_nc_axes_name  = '';        % [''] optional .nc name
par_nc_atmos_name = '';        % [''] optional .nc name
par_nc_ocean_name = '';        % [''] optional .nc name
par_nc_coupl_name = '';        % [''] optional .nc name
%
% *** GRID RESOLUTION *************************************************** %
%
par_max_i=36;                  % [36] # grid cells in longitude dir (i)
par_max_j=36;                  % [36] # grid cells in latitude  dir (j)
par_max_k=16;                  % [16] # depth levels in ocean
opt_equalarea=true;            % [false/true] equal area grid?
%
% *** REGRIDDING SETTINGS *********************************************** %
%
par_max_D=5000.0;              % [5000.0] max grid depth (m)
par_sur_D=0.0;                 % [0.0] reference surface depth
par_min_Dk=2;                  % [2] minimum ocean depth (as # levels)
par_min_k=3;                   % [1] maximum ocean depth (k value)
par_lon_off=-180.0;            % [-180.0] longitude offset of grid start
par_A_frac_threshold=0.50;     % [0.5] fractional area threshold for 'land'
par_mask_mask_name = '';       % [''] mask to change land/ocean features
par_sedsopt=0;                 % [0/1/2] sediment re-gridding option
opt_highresseds=false;         % [false/true] create 2x res sediment grid
%
% *** BOUNDARY CONDITION SETTINGS *************************************** %
%
par_runoffopt=0;               % [0/1] run-off generation option
par_tauopt=1;                  % [0/1/2] zonal windstress generation option
par_age=0.0;                   % [0.0] optional age of paleo configuration
%
% *** OPTIONS -- MAIN *************************************************** %
%
opt_makeall=false;             % [false/true] apply all common options?
opt_user=false;                % [false/true] enable user input to grid
opt_plots=true;                % plot all input and output?
%
% *** OPTIONS -- DATA GENERATION **************************************** %
%
opt_makemask=true;             % [false/true] re-grid mask?
opt_maketopo=true;             % [false/true] re-grid bathymetry?
opt_makeocean=true;            % [false/true] create ocean files?
opt_makerunoff=true;           % [false/true] create runoff pattern?
opt_makewind=true;             % [false/true] re-grid wind products?
opt_makezonalwind=false;       % [false/true] force zonal wind generation
opt_makealbedo=true;           % [false/true] make albedo file
opt_makeseds=false;            % [false/true] make sediment files
%
% *** OPTIONS -- GRID FILTERING ***************************************** %
%
opt_filtermask=false;          % [false/true] filter land-sea mask?
opt_filtertopo=false;          % [false/true] filter topography?
opt_makepoleswide=false;       % [false/true] force wide polar island zone
par_min_oceann=20;             % [20] minimum allowed lake size (# cells)
%
% *** ENVIRONMENT/OTHER SETTINGS **************************************** %
%
par_dpath_source='source';     % ['source'] relative path to muffingen code
opt_debug=false;               % [false/true] debug output?
%
% *********************************************************************** %
