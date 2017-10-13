function [grid_mask,grid_oceans,n_oceans]  = find_grid_oceans_update(grid_mask,grid_oceans,n_oceans,par_min_oceann)
%
%%

% *********************************************************************** %
% *** UPDATE THE OCEANS COUNT ******************************************* %
% *********************************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax, imax] = size(grid_mask);
gm = grid_mask;
go = grid_oceans;
g = zeros(jmax,imax);
g = g + 1;
%
% *** REMOVE TOO-SMALL WATER BODIES ************************************* %
%
loc_n_oceans = n_oceans;
for ocean = 1:n_oceans,
    v_ocean = find(go == ocean);
    if length(v_ocean) <= par_min_oceann,
        g(v_ocean) = 0;
        loc_n_oceans = loc_n_oceans-1;
    end
end
%
% *** RENUMBER OCEANS *************************************************** %
%
% renumber oceans
%
% %%% DON'T BOTHER FOR NOW!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% *** UDATE AND RETURN VARIABLES **************************************** %
%
n_oceans = loc_n_oceans;
grid_mask = g.*gm;
grid_oceans = go;
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%
% *********************************************************************** %
