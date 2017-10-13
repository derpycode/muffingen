function [gbrdrs] = fun_grid_cell_borders_search(j,i,gm,gbrdrs)
%
%%

% *********************************************************************** %
% *** GRID CELL SEARCH ************************************************** %
% *********************************************************************** %
%
% determine grid size
[jmax, imax] = size(gm);
% create search array
vdsrch = [1 1; 1 0; 1 -1; 0 -1; -1 -1; -1 0; -1 1; 0 1]; 
% set local extended grid (j,i) location indices
j = j+1;
i = i+1;
% expand grids
% NOTE: set expended N boundary as void (dry cells)
%       (also S boundary as void)
%       (remember wet == 1, searched == 1)
gm_ex = gm;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 0;
gm_ex(end,:) = 0;
% search surrounding cells
for s = 1:length(vdsrch)
    loc_j = j + vdsrch(s,1);
    loc_i = i + vdsrch(s,2);
    if gm_ex(loc_j,loc_i),
        % test for hitting East or West boundaries
        % NOTE: remember to convert back array indices
        if (loc_i == 1)
            gbrdrs(loc_j-1,imax) = 1;
        elseif (loc_i == imax+2)
            gbrdrs(loc_j-1,1) = 1;
        else
            gbrdrs(loc_j-1,loc_i-1) = 1;
        end
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
