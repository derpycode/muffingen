function [grid_n]  = find_grid_ncells(grid,pole)
%
%%

% *********************************************************************** %
% *** FIND # OF SURROUNDING (N-S-E-W) CELLS ***************************** %
% *********************************************************************** %
%
%
% determine mask size (remember: [rows columns])
[jmax, imax] = size(grid);
gn = zeros(jmax,imax);
% add boundaries to mask
% NOTE: poles are counted as whatever is passed as dummy variable 'pole'
g = grid;
g = [g(:,end) g g(:,1)];
g = [g(1,:); g; g(end,:)];
g(1,:)   = pole;
g(end,:) = pole;
% search across inner, *original* grid
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being true
        if g(j,i)
            cnt = g(j+1,i) + g(j-1,i) + g(j,i+1) + g(j,i-1);
        else
            cnt = NaN;
        end
        % populate n (number of N-S-E-W connecting cells) in output array
        gn(j-1,i-1) = cnt;
    end
end
grid_n = gn;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
