function [grid_cells]  = find_grid_4cells(grid_mask)
%
%%

% *********************************************************************** %
% *** FIND 4 SURROUNDING LAND CELLS ************************************* %
% *********************************************************************** %
%
%
% determine mask size (remember: [rows columns])
[jmax imax] = size(grid_mask);
gc = zeros(jmax,imax);
% add boundaries to mask
gm = grid_mask;
gm = [gm(:,end) gm gm(:,1)];
gm = [gm(1,:); gm; gm(end,:)];
gm(1,:)   = 1.0;
gm(end,:) = 1.0;
% search across inner, *original* grid
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being ocean
        if gm(j,i)
            cnt = gm(j+1,i) + gm(j-1,i) + gm(j,i+1) + gm(j,i-1);
            % test for cell surrounded by land (in 4 edge directions)
            if (cnt == 0)
                gc(j-1,i-1) = 1.0;
            end
        end
    end
end
grid_cells = gc;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
