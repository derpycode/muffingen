function [grid_mask]  = fun_grid_edit_mask(grid_mask)
%
% basic mask alteration code after:
% Andrew Yool (axy@soc.soton.ac.uk), October 2004.
%
%%

% *********************************************************************** %
% *** USER EDITING OF MASK ********************************************** %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_mask);
% local copy
% NOTE: extend to the N and E for successful plotting using pcolor ...
%       (duplicate E column and N row)
gm_ex = grid_mask;
gm_ex = [gm_ex gm_ex(:,end)];
gm_ex = [gm_ex(1,:); gm_ex];
% Do the topography alteration
fprintf('       * Mask alteration procedure :\n');
fprintf('         (1) left button to turn grid cell to ocean\n');
fprintf('         (2) right button to turn grid cell to land\n');
fprintf('         (3) middle button, or clickoutside the grid, to finish\n');
cmap = [0.5 1 0.5; 0.5 0.5 1]; 
%
flag = 0;
while flag == 0
    figure(1); clf
    colormap(cmap);
    pcolor (flipud(gm_ex)); axis image;
    caxis ([-0.5 1.5]); h = colorbar ('horiz');
    %set(h,'XTick',0:1,'XTickLabel',cbartxt);
    set(h,'XTick',0:1);
    title ('Ocean land/sea mask');
    [x,y,button] = ginput(1);
    ix = floor(x); iy = floor(y);
    if ix < 1 | ix > imax | iy < 1 | iy > jmax
        fprintf('       - Out of grid range\n');
        flag = 1;
        fprintf('       * Mask alteration complete\n');
    elseif button == 2
        flag = 1;
        fprintf('       * Mask alteration complete\n');
    else
        loc_j = jmax-iy+2;
        if button == 1
            fprintf('         -> Cell at (%d, %d) now ocean\n', ix, iy);
            gm_ex(loc_j, ix) = 1;
        elseif button == 3
            fprintf('         -> Cell at (%d, %d) now land\n', ix, iy);
            gm_ex(loc_j, ix) = 0;
        else
            fprintf('       - Cannot switch mask value of cell at (%d, %d)\n',ix,iy);
        end
    end
    clf
end
% return new mask
grid_mask = gm_ex(2:end,1:end-1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
