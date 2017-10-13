function [farea,farearef] = fun_grid_calc_ftotarea(mask,lone,late)
%
%%

% *********************************************************************** %
% *** CALCULATE TOTAL AREA ********************************************** %
% *********************************************************************** %
%
% determine grid size
lonmax = length(lone)-1;
latmax = length(late)-1;
%
loc_farea    = 0.0;
loc_farearef = 0.0;
%
% NOTE: lon then lat in loop order
for lon = 1:lonmax
% area fraction component in longitude
    loc_flon = (lone(lon+1) - lone(lon))/360.0;
    for lat = 1:latmax
% area fraction component in latitude
        loc_flat = (sin(pi*late(lat+1)/180.0) - sin(pi*late(lat)/180.0))/2.0;
        % sum masked and total area fractions
        loc_farea    = loc_farea + mask(lat,lon)*loc_flon*loc_flat;
        loc_farearef = loc_farearef + loc_flon*loc_flat;
    end
end
% return function value
farea    = loc_farea;
farearef = loc_farearef;
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
