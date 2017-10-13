function [grid_borders]  = find_grid_borders(grid_mask)
%
%%

% *********************************************************************** %
% *** IDENTIFY BORDERS AROUND ISLANDS *********************************** %
% *********************************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax imax] = size(grid_mask);
% copy mask
gm = grid_mask;
% initialzie borders array
gbrdrs = zeros(jmax,imax);
% search across grid
for i = 1:imax
    for j = 1:jmax
        % test for cell being land
        if ~gm(j,i)
            % look for ocean cells and mark borders
            [gbrdrs] = fun_grid_cell_borders_search(j,i,gm,gbrdrs);
        end
    end
end
% add a border to any remaining S pole bordering ocean cells 
% (and not already adjacent to land)
for i = 1:imax
    if gm(jmax,i) && gm(jmax-1,i),
        gbrdrs(jmax,i) = 1;
    end
end
% add a border to any remaining N pole bordering ocean cells
% (and not already adjacent to land)
for i = 1:imax
    if gm(1,i) && gm(2,i),
        gbrdrs(1,i) = 1;
    end
end
% return array
grid_borders = gbrdrs;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
