function [grid_runoff] = fun_grid_cell_runoff_search(j,i,grid_mask,grid_runoff)
%
%%

% *********************************************************************** %
% *** GRID CELL SEARCH ************************************************** %
% *********************************************************************** %
%
% KEY:
%      91 == E == [0 1]
%      92 == S == [1 0]
%      93 == W == [0 -1]
%      94 == N == [-1 0]
% NOTE: the grid is flipped up-down re. row indexing ...
%
% determine grid size
[jmax, imax] = size(grid_mask);
% create search array
vdsrch = [0 1 91; 0 -1 93; -1 0 94; 1 0 92];
n_max = length(vdsrch);
% expand grids
% NOTE: set expended N boundary as void (dry cells)
%       (also S boundary as void)
%       (remember wet == 1, searched == 1)
gm_ex = grid_mask;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 0;
gm_ex(end,:) = 0; 
% set local extended grid (j,i) location indices
j = j+1;
i = i+1;
% create random order of search
% (to avoid there being a systematic basi in runoff direction)
% NOTE: only randomize E vs. W, and N vs. S (seperately)
vdsrchdir = [randperm(n_max/2), randperm(n_max/2)+n_max/2];
% search surrounding cells
for n = 1:n_max
    s = vdsrchdir(n);
    loc_j = j + vdsrch(s,1);
    loc_i = i + vdsrch(s,2);
    if gm_ex(loc_j,loc_i),
        % NOTE: remember to convert back array indices
        %       but no East or West boundary issue
        %       (grid_runoff is jmax x imax)
        grid_runoff(j-1,i-1) = vdsrch(s,3);
        % stop searching!
        break 
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
