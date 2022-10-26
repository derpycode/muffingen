function [grid_k1]  = make_grid_runoff_roof(grid_mask,grid_k1,str)
%
%%

% *********************************************************************** %
% *** GENERATE ROOF RUNOFF PATTERN ************************************** %
% *********************************************************************** %
%
% KEY:
%      91 == E == [0 1]
%      92 == S == [1 0]
%      93 == W == [0 -1]
%      94 == N == [-1 0]
% NOTE: the grid is flipped up-down re. row indexing ...
%
% *** INITIALIZE ******************************************************** %
%
% determine grid size (remember: [rows columns])
[jmax imax] = size(grid_mask);
% filter input grid
% NOTE: input grid may have:
%       ocean > 0, land 0
%       ocean > 0, land >90
% NOTE: retain existing runoff information
gr = grid_k1;
gr(find(gr > 94)) = 90; % gr(find(gr > 90)) = 90;
gr(find(gr == 0)) = 90;
gr(find(gr < 90)) = 0;
% create search mask
gm = grid_mask;
% update mask to remove existing run-off cells from search
gm(find(gr > 90)) = 1;
%
% *** ITERATIVELY POPULATE GRID ***************************************** %
%
search = true;
while search
    % search across grid
    % NOTE: '1' in the grid mask == ocean or already connected or ocean
    for i = 1:imax
        for j = 1:jmax
            % test for unassigned (land) cell
            if ~gm(j,i)
                % search for cell adjacent to watershed: assign runoff dir
                [gr] = fun_grid_cell_runoff_search(j,i,gm,gr);
            end
        end
    end
    % update mask
    gm(find(gr > 90)) = 1;
    % check for grid being fully populated
    if min(min(gm)) == 1, search = false; end
    %
end
%
% *** MERGE ARRAYS ****************************************************** %
%
% create and return final k1 array
% NOTE: runoff array has zero for ocean and >90 for land
%       original grid_k1 adjusted to have zero for land, >0 for ocean
%       (so they'll happily sum/combine)
grid_k1(find(grid_k1 >= 90)) = 0;
grid_k1 = grid_k1 + gr;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
