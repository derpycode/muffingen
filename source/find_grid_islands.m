function [grid_islands,n_islands,i_islands]  = find_grid_islands(grid_mask)
%
%%

% *********************************************************************** %
% *** COUNT THE ISLANDS! ************************************************ %
% *********************************************************************** %
%
% determine output grid size (remember: [rows columns])
[jmax, imax] = size(grid_mask);
% copy mask
gm = grid_mask;
% initialzie island array, searched array
gislnd = zeros(jmax,imax);
gsrch  = zeros(jmax,imax);
% initialize island count
n_islnd = 0;
% 
i_islands = [];
% search across grid
for i = 1:imax
    for j = 1:jmax
        % test for cell being land and not yet searched
        % NOTE: mask is defined with '1' for ocean
        if ~gm(j,i) && ~gsrch(j,i)
            % update island count
            n_islnd = n_islnd + 1;
            disp(['       * found ' num2str(n_islnd) ' land masses']);
            i_islands = [i_islands n_islnd];
            % seed the search vector with just the current point
            vsrch = [j i];
            % mark cell as searched
            gsrch(vsrch(1,1),vsrch(1,2)) = 1;
            % set initial search vector length
            n = 1;
            % now conduct recursive search
            while n > 0
                % carry out surrounding cell search and update
                % vector of locations to be further searched from
                [gsrch vsrch] = fun_grid_cell_island_search(gm,gsrch,vsrch);
                % assign island number
                gislnd(vsrch(1,1),vsrch(1,2)) = n_islnd;
                % remove used search vector entry
                vsrch(1,:) = [];
                % update search vector length
                n = length(vsrch);
            end
        end
    end
end
% return variables
grid_islands = gislnd;
n_islands = n_islnd;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
