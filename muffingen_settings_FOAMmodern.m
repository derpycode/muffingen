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
% *** EXPERIMENT INPUT AND OUTPUT *************************************** %
%
par_wor_name='FM000a';         % ['wworld'] 6-char (output) config name
par_gcm='foam';             % ['hadcm3l'] input format/GCM name
par_expid='ctrlAC_1370W';             % ['xbomg'] experiment name
par_pathin='../DATA.FOAM';   % ['EXAMPLES.INPUT'] path to input dir
par_pathout='../OUTPUT'; % ['EXAMPLES.OUTPUT'] path to output dir
opt_outputdir=false;           % [false/true] ask for output directory?
%par_nc_topo_name='';  % [''] -- optional/specific .nc file name
%par_nc_atmos_name=''; % [''] -- optional/specific .nc file name
%par_nc_coupl_name=''; % [''] -- optional/specific .nc file name
%par_nc_ocean_name=''; % [''] -- optional/specific .nc file name
%
% *** GRID RESOLUTION & REGRIDDING CONTROLS ***************************** %
%
par_max_i=36;                  % [36] # grid cells in longitude dir (i)
par_max_j=36;                  % [36] # grid cells in latitude  dir (j)
par_max_k=16;                  % [36] # depth leves in ocean
par_max_D=5000.0;              % [5000.0] max grid depth (m)
par_lon_off=-180.0;            % [-180.0] longitude offset of grid start
par_min_Dk=2;                  % [2] minimum # depth levels of ocean cells
par_A_frac_threshold=0.5;     % [0.5] fractional area threshold for 'land'
opt_equalarea=true;            % [false/true] equal area grid?
opt_highresseds=false;         % [false/true] create 2x res sediment grid
par_sedsopt=0;                 % [0/1/2] sediment re-gridding option
%
% *** OPTIONS -- MAIN *************************************************** %
%
opt_makeall=false;              % [false/true] apply all common options?
opt_user=true;                 % [false/true] enable user input to grid
%
% *** OPTIONS -- OTHER ************************************************** %
%
opt_makemask=true;             % [false/true] re-grid mask?
opt_maketopo=true;             % [false/true] re-grid bathymetry?
opt_makeocean=true;            % [false/true] creat ocean files?
opt_makerunoff=true;           % [false/true] create runoff pattern?
opt_makewind=true;             % [false/true] re-grid wind products?
opt_makealbedo=true;           % [false/true] make albedo file
opt_makeseds=false;             % [false/true] make sediment files
%
% *** GRID FILTERING **************************************************** %
%
opt_filtermask=true;           % [false/true] filter land-sea mask?
opt_filtertopo=true;           % [false/true] filter topography?
opt_makepoleswide=true;        % [false/true] force wide polar island zone
par_min_oceann=20;             % [20] minimum allowed lake size (# cells)
%
% *** ENVIRONMENT SETTINGS ********************************************** %
%
par_age=0.0;                   % [0.0] optional age of paleo configuration
par_dpath_source='source';     % ['source'] relative path to muffingen code
opt_debug=false;               % [false/true] debug output?
%
% *********************************************************************** %
