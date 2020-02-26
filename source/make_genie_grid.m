function [axis_imid,axis_iedge,axis_jmid,axis_jedge,axis_kmid,axis_kedge,par_max_D] = make_genie_grid(DUM_NI,DUM_NJ,DUM_NK,DUM_MAXD,DUM_OFF,DUM_EA,DUM_DK);
% MAKE_GENIE_GRID
%
%   ***********************************************************************
%   *** Generate GENIE grid ***********************************************
%   ***********************************************************************
%
% NOTE: GOLDSTEIn variable values calculated according to gseto.f
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   16/09/09: started ...
%             mostly copied from make_regrid_WOA2GENIE.m
%   17/01/21: renamed function, renamed parameters
%             tarted up
%   17/02/23: adjsted formatting of messages
%   17/03/07: added option for non-equal area
%             adjsted formatting of messages ...
%   20/02/25: improved the deeper ocean generation capability
%             -> now additional deeper layers can be added to
%                a known (e.g. 16-level) initial ocean layer structure
%
%   ***********************************************************************
%%

% *********************************************************************** %
% *** GENERATE GENIE GRID *********************************************** %
% *********************************************************************** %

% START
%%%disp(['       START [make_genie_grid] >>>'])
close all;

% INITIALIZE
% dummy variables
n_i = DUM_NI;
n_j = DUM_NJ;
n_k = DUM_NK;
par_max_D            = DUM_MAXD;
par_grid_i_offset_in = DUM_OFF;
par_grid_equalarea   = DUM_EA;
dk                   = DUM_DK; % additional k levels to be added to base
% local constants
par_ez0 = 0.1;
% adjust n_k
%
% *** create GENIE grid ************************************************* %
%
% NOTE: par_grid_i_offset_in is to enable the GENIE grid to be matched
%       to the input grid (i.e. align the Prime Meridian)
% NOTE: at this point, no account is taken of whether the final GENIE
%       grid should start at e.g. -260E (set by par_grid_i_offset_out)
%
% lon (west boundary)
for i=1:n_i
    axis_iedge(i) = (i-1)*(360.0/n_i) + par_grid_i_offset_in;
    axis_di(i)    = (360.0/n_i);
end
axis_iedge(n_i+1)   = (n_i)*(360.0/n_i) + par_grid_i_offset_in;
axis_imid           = axis_iedge(1:n_i) + 0.5*axis_di;
axis_ibnds(1:n_i,1) = axis_iedge(1:n_i);
axis_ibnds(1:n_i,2) = axis_iedge(2:n_i+1);
axis_ibnds          = axis_ibnds';
% lat (south boundary)
if par_grid_equalarea
    for j=1:n_j
        axis_jedge(j) = (180.0/pi)*asin((2*(j-1)/n_j) - 1.0);
        axis_jmid(j)  = (180.0/pi)*asin(((1 + 2*(j-1))/n_j) - 1.0);
    end
    axis_jedge(n_j+1)   = (180.0/pi)*asin((2*n_j/n_j) - 1.0);
    axis_jbnds(1:n_j,1) = axis_jedge(1:n_j);
    axis_jbnds(1:n_j,2) = axis_jedge(2:n_j+1);
    axis_jbnds          = axis_jbnds';
else
    for j=1:n_j
        axis_jedge(j) = 90.0*(2*(j-1) - n_j)/n_j;
        axis_jmid(j)  = 90.0*(1 + 2*(j-1) - n_j)/n_j;
    end
    axis_jedge(n_j+1)   = 90.0;
    axis_jbnds(1:n_j,1) = axis_jedge(1:n_j);
    axis_jbnds(1:n_j,2) = axis_jedge(2:n_j+1);
    axis_jbnds          = axis_jbnds';    
end
% depth (bottom boundary)
% NOTE: set scale depth using adjusted value of n_k
z1 = par_ez0*((1.0 + 1/par_ez0)^(1.0/(n_k-dk)) - 1.0);
tv4 = par_ez0*((z1/par_ez0+1)^0.5-1);
tv2 = 0;
tv1 = 0;
% initialize arrays
zro(n_k) = -tv4;
zw(n_k+1) = tv2;
% calculate (scaled) layer depths
for k=1:1:n_k
    if par_ez0 > 0
        tv3 = par_ez0*((z1/par_ez0+1)^k-1);
        dz(n_k-k+1) = tv3 - tv2;
        tv2 = tv3;
        tv5 = par_ez0*((z1/par_ez0+1)^(k+0.5)-1);
        if k < n_k
            dza(n_k-k) = tv5 - tv4;
        end
        tv4 = tv5;
        tv1 = tv1 + dz(n_k-k+1);
    else
        dz(k) = 1d0/n_k;
        dza(k) = 1d0/n_k;
    end
end
for k=n_k:-1:1
    if k > 1
        zro(k-1) = zro(k) - dza(k-1);
    end
    zw(k) = zw(k+1) - dz(k);
end
% set depth grid bounds
% NOTE: k counts from TOP to BOTTOM; 
%       bnd #1 is top
axis_kmid(1:n_k)    = -par_max_D*zro(:);
axis_kbnds(1:n_k,1) = -par_max_D*zw(2:n_k+1);
axis_kbnds(1:n_k,2) = -par_max_D*zw(1:n_k);
axis_kedge(1:n_k+1) = -par_max_D*zw(1:n_k+1);
axis_kth(1:n_k)     = -par_max_D*dz(1:n_k);
axis_kbnds = axis_kbnds';
% update maximum ocean depth
par_max_D = axis_kbnds(2,1);
%
% END
%%%disp(['       <<< END [make_genie_grid]'])
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
