function [grid_topos] = make_grid_topo_sed_rnd(grid_mask,opt_highresseds,bath_min,bath_max)
% make_grid_topo_sed_rnd
%
%   ***********************************************************************
%   *** Create randomized SEDGEM topo *************************************
%   ***********************************************************************
%
%   make_grid_topo_sed_rnd(grid_mask,opt_highresseds)
%   BLAH
%
%   ***********************************************************************
%   *** HISTORY ***********************************************************
%   ***********************************************************************
%
%   14/11/22: CREATED [derived from Goldschmidt 2013 era process_hyps.m]
%   17/03/13: re-write for muffingen
%
%   ***********************************************************************

% *********************************************************************** %
% *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
% *********************************************************************** %
%
% determine output grid size (remember: [rows columns])
[jmax, imax] = size(grid_mask);
if opt_highresseds,
    imaxs = 2*imax;
    jmaxs = 2*jmax;
else
    imaxs = imax;
    jmaxs = jmax;
end
% create output array
grid_topos = zeros(jmaxs,imaxs);
%
% *********************************************************************** %

% *********************************************************************** %
% *** LOAD HYPSOGRAPHIC DATA ******************************************** %
% *********************************************************************** %
%
% load data (processed with 10 m intervals bins)
hyps=load('elev0010_area_int_norm.dat','-ascii');
% invert depth
hyps(:,1)=-hyps(:,1);
% filter out too shallow/deep data
hyps(find(hyps(:,1)<bath_min),:)=[];
hyps(find(hyps(:,1)>bath_max),:)=[];
%
% *** process data ****************************************************** %
%
% calculate total fractional area remaining (after filtering)
loc_total_frac_seafloor = hyps(end,2)-hyps(1,2);
% re-normalize fractional area
hyps(:,2) = (hyps(:,2)-hyps(1,2))/loc_total_frac_seafloor;
%
% *********************************************************************** %

% *********************************************************************** %
% *** CREATE TOPO ******************************************************* %
% *********************************************************************** %
%
for i = 1:imaxs,
    for j = 1:jmaxs,
        if grid_mask(i,j),
            if opt_highresseds
                for ii = 0:1
                    for jj = 0:1
                        roll = rand(1);
                        grid_topos(i+ii,j+jj) = double(hyps(find(abs(roll - hyps(:,2)) == min(abs(roll - hyps(:,2)))),1));
                    end
                end
            else
                roll = rand(1);
                grid_topos(i,j) = double(hyps(find(abs(roll - hyps(:,2)) == min(abs(roll - hyps(:,2)))),1));
            end
        end
    end
end%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
