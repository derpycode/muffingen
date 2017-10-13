function [vo_albd]  = make_grid_albd(go_latm,par_age)
%
%%

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% determine output vector size
[jmax] = length(go_latm);
% create otput vector
vo_albd = 0.0*go_latm;
%
loc_latm = (pi/180.0)*go_latm;
%
% *********************************************************************** %

% *********************************************************************** %
% *** CREATE GENERIC ALBEDO ********************************************* %
% *********************************************************************** %
%
% test for modern and set albedo profile parameters
if (par_age == 0.0),
    % assume modern
    % NOTE: this is the GENIE modern default
    climalbedo_offset = 0.20;
    climalbedo_amp    = 0.36;
    climalbedo_mod    = 0.0;
else
    % assume ice-free
    % NOTE: this is the Ridgwell and Schmidt [2010] early Eocene profile
    climalbedo_offset = 0.20;
    climalbedo_amp    = 0.26;
    climalbedo_mod    = 0.25;
end
% calculate the albedo profile itself
for j = 1:jmax,
    vo_albd(j) = climalbedo_offset + climalbedo_amp*0.5*(1.0 - cos(2.0*loc_latm(j)) + climalbedo_mod*cos(6.0*loc_latm(j)));
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
