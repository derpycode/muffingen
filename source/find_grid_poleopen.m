function [grid_cells]  = find_grid_poleopen(grid_mask)
%
% CLEAN:
%
% =====================
% !   !   ! X !   !   ! j == 1
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
% !   !   ! X !   !   ! j == jmax
% =====================
%
% DON'T CLEAN:
%
% =====================
% !   ! X ! X !   !   ! j == 1
% ---------------------
% ! X ! X ! X ! X ! X ! j == 2
% ---------------------
% ! X ! X ! X ! X ! X ! j == 3
% ---------------------
% ! X ! X ! X ! X ! X ! ...
% ---------------------
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
    % test for cell being land
    if ~gm(j,i)
        % adjacent land to N and S (counting pole as land)
        % AND ocean to W and E
        if (~(gm(j+1,i)||gm(j-1,i)))&&gm(j,i-1)&&gm(j,i+1), gc(j-1,i-1) = 1.0; end
    end
    j = jmax+1;
    % test for cell being land
    if ~gm(j,i)
        % adjacent land to N and S (counting pole as land)
        % AND ocean to W and E
        if (~(gm(j+1,i)||gm(j-1,i)))&&gm(j,i-1)&&gm(j,i+1), gc(j-1,i-1) = 1.0; end
    end
end
grid_cells = gc;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
