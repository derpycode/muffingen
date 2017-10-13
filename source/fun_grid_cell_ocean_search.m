function [gsrch,vsrch]  = fun_grid_cell_ocean_search(gm,gsrch,vsrch)
%
%%

% *********************************************************************** %
% *** GRID CELL SEARCH ************************************************** %
% *********************************************************************** %
%
% determine grid size
[jmax, imax] = size(gm);
% create search array
vdsrch_nsew = [1 0; 0 1; -1 0; 0 -1]; 
% set local extended grid (j,i) location indices
j = vsrch(1,1)+1;
i = vsrch(1,2)+1;
% expand grids
% NOTE: set expended N and S boundaries as void (searched/wet cells)
%       (remember wet == 1, searched == 1)
gm_ex = gm;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 1.0;
gm_ex(end,:) = 1.0;
gsrch_ex = gsrch;
gsrch_ex = [gsrch_ex(:,end) gsrch_ex gsrch_ex(:,1)];
gsrch_ex = [gsrch_ex(1,:); gsrch_ex; gsrch_ex(end,:)];
gsrch_ex(1,:)   = 1.0;
gsrch_ex(end,:) = 1.0;
% search surrounding cells
% NOTE: only count N-S-E-W connections when determining ocean connectivity
for s = 1:length(vdsrch_nsew)
    loc_j = j + vdsrch_nsew(s,1);
    loc_i = i + vdsrch_nsew(s,2);
    if gm_ex(loc_j,loc_i) && ~gsrch_ex(loc_j,loc_i),
        % test for hitting East or West boundaries,
        % and add new entry to vector
        % NOTE: remember to convert back array indices
        % mark cell as searched (even if not yet assigned an island num)
        % so that it does not get added to vector again
        if (loc_i == 1)
            vsrch = [vsrch; loc_j-1 imax];
            gsrch(loc_j-1,imax) = 1;
        elseif (loc_i == imax+2)
            vsrch = [vsrch; loc_j-1 1];
            gsrch(loc_j-1,1) = 1;
        else
            vsrch = [vsrch; loc_j-1 loc_i-1];
            gsrch(loc_j-1,loc_i-1) = 1;
        end
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
