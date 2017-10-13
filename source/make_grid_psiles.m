function [grid_psiles,n_islands] = make_grid_psiles(grid_islands,i_poles)
%
%%

% *********************************************************************** %
% *** CREATE PSILES GRID ************************************************ %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_islands);
% make local copies
gi   = grid_islands;
gpsi = gi;
% create PSI islands grid
% NOTE: format has the effect of:
%       (1) shaddow N-wards (creating a new jmax+1 row)
%       (2) shaddow W-wards (wrap-around ... do no create a new column)
%       the reason for this is that the PSI grid is at the SE corner of the
%       c-grid cell, in the u and -v direction (E and S)
%       hence a given land point invalidaes the PSI point to the N, and W
%       the additional row of cells at the N pole, completes the grid
%       (such that each c-point is bounded by PSI points
% (1) deal with u direction
for j = 1:jmax
    % first deal with E-W wall
    % (marking cells on the far E boundary if the corresponding cell is an
    % island on the W boundary)
    i = 1;
    if (gi(j,i) > 0) || (gi(j,i) == -3),
        gpsi(j,imax) = gi(j,i);
    end
    % now for the rest of the grid
    for i = 2:imax
        if (gi(j,i) > 0) || (gi(j,i) == -3),
            gpsi(j,i-1) = gi(j,i);
        end
    end
end
% (2) deal with -v direction
for j = 2:jmax
    for i = 1:imax
        if (gi(j,i) > 0) || (gi(j,i) == -3),
            gpsi(j-1,i) = gi(j,i);
        end
    end
end
% (3) deal with u AND -v direction (W-N diagianol cells)
for j = 2:jmax
    % first deal with E-W wall
    % (marking cells on the far E boundary if the corresponding cell is an
    % island on the W boundary)
    i = 1;
    if (gi(j,i) > 0) || (gi(j,i) == -3),
        gpsi(j-1,imax) = gi(j,i);
    end
    % now for the rest of the grid
    for i = 2:imax
        if (gi(j,i) > 0) || (gi(j,i) == -3),
            gpsi(j-1,i-1) = gi(j,i);
        end
    end
end
% (4) extend j grid (at top)
gpsi = [gpsi(1,:); gpsi];
% (5) add a N pole row
if ~isempty(find(i_poles == -1)),
    loc_isnl = max(max(gpsi)) + 1;
else
    loc_isnl = max(gpsi(1,:));
end
gpsi(1,:) = loc_isnl;
% (6) fix S pole row
if ~isempty(find(i_poles == -2)),
    loc_isnl = max(max(gpsi)) + 1;
else
    loc_isnl = max(gpsi(end,:));
end
gpsi(end,:) = loc_isnl;
% (7) N-S connected island ...
if ~isempty(find(i_poles == -3)),
    loc_isnl = max(max(gpsi)) + 1;
    gpsi(find(gpsi == -3)) = loc_isnl;
    gpsi(1,:) = loc_isnl;
    gpsi(end,:) = loc_isnl;
end
% (8) return psiles grid
grid_psiles = gpsi;
n_islands   = loc_isnl;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
