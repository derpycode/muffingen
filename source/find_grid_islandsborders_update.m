function [grid_islands,grid_borders,n_islands,i_islands]  = find_grid_islandsborders_update(grid_islands,grid_borders,n_islands,i_islands,grid_mask)
%
%%

% *********************************************************************** %
% *** UPDATE THE ISLANDS COUNT ****************************************** %
% *********************************************************************** %
%
% initialize
% determine mask size (remember: [rows columns] == [j i])
[jmax, imax] = size(grid_islands);
% create search array
vdsrch = [1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1; 0 1]; 
%
% *********************************************************************** %

% *********************************************************************** %
% *** UPDATE THE ISLANDS COUNT ****************************************** %
% *********************************************************************** %
%
% copy & expand grids
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0.0;
gb_ex(end,:) = 0.0;
gi_ex = grid_islands;
gi_ex = [gi_ex(:,end) gi_ex gi_ex(:,1)];
gi_ex = [gi_ex(1,:); gi_ex; gi_ex(end,:)];
gi_ex(1,:)   = 0.0;
gi_ex(end,:) = 0.0;
% ADJACENT ISLANDS AWAY FROM POLES
isearch  = true;
ireplace = false;
while isearch
    % search across grid
    for i = 2:imax+1
        for j = 2:jmax+1
            % set results vector empty
            isrch = [];
            % test for cell being a border
            if gb_ex(j,i)
                % search surrounding cells
                for s = 1:8
                    loc_j = j + vdsrch(s,1);
                    loc_i = i + vdsrch(s,2);
                    if gi_ex(loc_j,loc_i) > 0,
                        isrch = [isrch; gi_ex(loc_j,loc_i)];
                    end
                end
            end
            % test for > 1 different island adjacent to border cell
            if (length(isrch) > 0) && (min(isrch) ~= max(isrch))
                ireplace = true;
            end
            if ireplace, break; end
        end
        if ireplace, break; end
    end
    % test for:
    % island re-numbering required
    % or no further changes after grid search complete
    if ireplace,
        % update island numbering in  original and expanded islands grids
        gi_ex(find(gi_ex == max(isrch))) = min(isrch);
        grid_islands(find(grid_islands == max(isrch))) = min(isrch);
        % update island numbers list
        i_islands(find(i_islands == max(isrch))) = [];
        % update islands count
        n_islands = n_islands - 1;
        % reset islands replacement flag
        ireplace = false;
    else
        isearch = false;
    end
end
% INDENTIFY ISLANDS CONNECTED BY POLES
% NOTE: at this point, polar-only islands are not identified
% search N polar row
for i = 1:imax
    % test for cell being a N island cell
    if (grid_islands(1,i) > 0)
        % update island numbers list
        % (do this before replacing the critical number in question  ...)
        i_islands(find(i_islands == grid_islands(1,i))) = [];
        % update islands count
        n_islands = n_islands - 1;
        % replace island number
        grid_islands(find(grid_islands == grid_islands(1,i))) = -1;
    end
end
% search S polar row
for i = 1:imax
    % test for cell being a N island cell
    if (grid_islands(jmax,i) > 0)
        % update island numbers list
        % (do this before replacing the critical number in question  ...)
        i_islands(find(i_islands == grid_islands(jmax,i))) = [];
        % update islands count
        n_islands = n_islands - 1;
        % replace island number
        grid_islands(find(grid_islands == grid_islands(jmax,i))) = -2;
    end
end
% renumber islands -- non-polar
if (n_islands > 0)
    for islnd = 1:n_islands,
        grid_islands(find(grid_islands == i_islands(islnd))) = islnd;
        i_islands(islnd) = islnd;
        disp(['   * updated find: ' num2str(islnd) ' islands']);
    end
end
% renumber islands -- polar
if ~isempty(find(grid_islands == -1))
    n_islands = n_islands + 1;
    grid_islands(find(grid_islands == -1)) = n_islands;
end
if ~isempty(find(grid_islands == -2))
    n_islands = n_islands + 1;
    grid_islands(find(grid_islands == -2)) = n_islands;
end
%
% *********************************************************************** %


% *********************************************************************** %
% *** FURTHER REFINE ISLANDS AND REFINE ISLANDS BORDERS ***************** %
% *********************************************************************** %
%
% copy & expand grids
% NOTE: mask is defined with '1' for ocean
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0;
gb_ex(end,:) = 0;
gi_ex = grid_islands;
gi_ex = [gi_ex(:,end) gi_ex gi_ex(:,1)];
gi_ex = [gi_ex(1,:); gi_ex; gi_ex(end,:)];
gi_ex(1,:)   = 0;
gi_ex(end,:) = 0;
gm_ex = grid_mask;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 0;
gm_ex(end,:) = 0;
% find any non zero island number cell bordering N pole
% and set extended array accordingly
if (max(grid_islands(end,:)) > 0)
    gi_ex(1,:) = max(grid_islands(end,:));
else
    gi_ex(end,:) = n_islands + 1;
    disp(' *** island/border situation excepton *** ');
    return;
end
% find any non zero island number cell bordering S pole
% and set extended array accordingly
if (max(grid_islands(end,:)) > 0)
    gi_ex(end,:) = max(grid_islands(end,:));
else
    gi_ex(end,:) = n_islands + 1;
    disp(' *** island/border situation excepton *** ');
    return;
end
% search across grid
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being a border
        if gb_ex(j,i)
            % set results index empty
            isrch = 0;
            % search surrounding cells
            for s = 1:8
                loc_j = j + vdsrch(s,1);
                loc_i = i + vdsrch(s,2);
                % test for adjacent island cell
                if (gi_ex(loc_j,loc_i) ~= 0),
                    isrch = gi_ex(loc_j,loc_i);
                end
            end
            % set border number
            % NOTE: remember to convert back array indices
            grid_borders(j-1,i-1) = isrch;
        end
    end
end
% update extended borders array
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0.0;
gb_ex(end,:) = 0.0;
% search across grid -- remove buried border cells
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being a border
        if gb_ex(j,i) ~= 0
            % set island index
            iislnd = gb_ex(j,i);
            isrch  = 0;
            % search surrounding cells
            for s = 1:8
                loc_j = j + vdsrch(s,1);
                loc_i = i + vdsrch(s,2);
                % test adjacent cells
                % (if ocean and not the same island border)
                if (gm_ex(loc_j,loc_i) && (gb_ex(loc_j,loc_i) ~= iislnd)),
                    isrch = 1;
                end
            end
            % reset border if no adjoining ocean grid point
            % of different border or no border exists
            % NOTE: remember to convert back array indices
            if ~isrch,
                grid_borders(j-1,i-1) = 0;
            end
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
