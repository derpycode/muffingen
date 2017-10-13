function [grid_cells]  = find_grid_poleclsd(grid_mask)
%
% CHANGE:
%
% =====================
% !   !   ! X ! X !   ! j == 1
% ---------------------
% ! X ! X ! X ! X ! X ! j == 2
% ---------------------
% ! X ! X ! X ! X ! X ! j == 3
% ---------------------
% ! X ! X ! X ! X ! X ! ...
% ---------------------
%
% ---------------------
% ! X ! X ! X ! X ! X ! ...
% ---------------------
% ! X ! X ! X ! X ! X ! j == jmax-2
% ---------------------
% ! X ! X ! X ! X ! X ! j == jmax-1
% ---------------------
% !   !   ! X ! X !   ! j == jmax
% =====================
%
% TO:
%
% =====================
% ! X ! X ! X ! X ! X ! j == 1
% ---------------------
% ! X ! X ! X ! X ! X ! j == 2
% ---------------------
% ! X ! X ! X ! X ! X ! j == 3
% ---------------------
% ! X ! X ! X ! X ! X ! ...
% ---------------------
%
% ---------------------
% ! X ! X ! X ! X ! X ! ...
% ---------------------
% ! X ! X ! X ! X ! X ! j == jmax-2
% ---------------------
% ! X ! X ! X ! X ! X ! j == jmax-1
% ---------------------
% ! X ! X ! X ! X ! X ! j == jmax
% =====================
%
%%

% *********************************************************************** %
% *** FIND SINGLE POLE CONNECTED CELLS ********************************** %
% *********************************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax imax] = size(grid_mask);
gc = zeros(jmax,imax);
% add boundaries to mask
gm = grid_mask;
gm = [gm(:,end) gm gm(:,1)];
gm = [gm(1,:); gm; gm(end,:)];
gm(1,:)   = 0.0;
gm(end,:) = 0.0;
% search across borders of poles
for i = 2:imax+1
    j = 2;
    % test for cell being ocean
    if gm(j,i)
        % adjacent land to N and S (counting pole as land)
        if (~gm(j+1,i)&&~gm(j-1,i)), gc(j-1,i-1) = -1; end
    end
    j = jmax+1;
    % test for cell being ocean
    if gm(j,i)
        % adjacent land to N and S (counting pole as land)
        if (~gm(j+1,i)&&~gm(j-1,i)), gc(j-1,i-1) = -1; end
    end
end
grid_cells = gc;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
