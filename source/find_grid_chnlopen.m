function [grid_cells]  = find_grid_chnlopen(grid_mask)
%
% OPEN:
%
% ---------------------
% ! X ! X ! X ! X ! X !
% ---------------------
% !   !   ! X !   !   !
% ---------------------
% !   !   !   !   !   !
% ---------------------
% ! X ! X ! X ! X ! X !
% ---------------------
%
% ---------------------
% ! X ! X ! X ! X ! X !
% ---------------------
% !   !   !   ! X !   !
% ---------------------
% !   !   !   !   !   !
% ---------------------
% ! X ! X ! X !   !   !
% ---------------------
%
%%

% *********************************************************************** %
% *** FIND NARROW CHANNEL *********************************************** %
% *********************************************************************** %
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
            % test for all combinations of:
            % a cell surrounded by land on opposite edges
            if gm(j,i+1)&&gm(j,i-1), gc(j-1,i-1) = 1.0; end
            if gm(j,i-1)&&gm(j,i+1), gc(j-1,i-1) = 1.0; end
            if gm(j+1,i)&&gm(j-1,i), gc(j-1,i-1) = 1.0; end
            if gm(j-1,i)&&gm(j+1,i), gc(j-1,i-1) = 1.0; end
        end
    end
end
grid_cells = gc;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
