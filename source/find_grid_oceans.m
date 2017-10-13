function [grid_oceans,n_oceans,i_oceans]  = find_grid_oceans(grid_mask)
%
%%

% *********************************************************************** %
% *** COUNT THE OCEANS! ************************************************* %
% *********************************************************************** %
%
% determine output grid size (remember: [rows columns])
[jmax, imax] = size(grid_mask);
% copy mask
gm = grid_mask;
% initialzie ocean array, searched array
gocean = zeros(jmax,imax);
gsrch  = zeros(jmax,imax);
% initialize island count
n_ocean = 0;
% 
i_oceans = [];
% search across grid
for i = 1:imax
    for j = 1:jmax
        % test for cell being ocean and not yet searched
        % NOTE: mask is defined with '1' for ocean
        if gm(j,i) && ~gsrch(j,i)
            % update ocean count
            n_ocean = n_ocean + 1;
            disp(['       * found ' num2str(n_ocean) ' water masses']);
            i_oceans = [i_oceans n_ocean];
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
                [gsrch vsrch] = fun_grid_cell_ocean_search(gm,gsrch,vsrch);
                % assign island number
                gocean(vsrch(1,1),vsrch(1,2)) = n_ocean;
                % remove used search vector entry
                vsrch(1,:) = [];
                % update search vector length
                n = length(vsrch);
            end
        end
    end
end
% return variables
grid_oceans = gocean;
n_oceans = n_ocean;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
