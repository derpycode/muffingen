function [grid_borders]  = fun_grid_edit_borders(grid_borders,grid_mask)
%
%
%%

% *********************************************************************** %
% *** USER EDITING OF BORDERS ******************************************* %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_borders);
% local copy
% NOTE: extend to the N and E for successful plotting using pcolor ...
%       (duplicate E column and N row)
gb_ex = grid_borders;
gb_ex = [gb_ex gb_ex(:,end)];
gb_ex = [gb_ex(1,:); gb_ex];
gm_ex = grid_mask;
gm_ex = [gm_ex gm_ex(:,end)];
gm_ex = [gm_ex(1,:); gm_ex];
% Do the topography alteration
fprintf('       * Borders alteration procedure :\n');
fprintf('         (1) left button to increase border #\n');
fprintf('         (2) right button to decrease border # (0 == no border)\n');
fprintf('         (3) middle button, or clickoutside the grid, to finish\n');
fprintf('         [exit (middle button) if unsure WTF]\n');
%
flag = 0;
while flag == 0
    figure(1); clf
    pcolor (flipud(gb_ex)); % axis image;
    caxis ([0.5 (max(max(grid_borders))+1.5)]); h = colorbar ('horiz');
    set(h,'XTick',1:(max(max(grid_borders))+1));
    title ('Borders');
    [x,y,button] = ginput(1);
    ix = floor(x); iy = floor(y);
    if ix < 1 | ix > imax | iy < 1 | iy > jmax
        fprintf('       - Out of grid range\n');
        flag = 1;
        fprintf('       * Borders alteration complete\n');
    elseif button == 2
        flag = 1;
        fprintf('       * Borders alteration complete\n');
    else
        loc_j = jmax-iy+2;
        if (gm_ex(loc_j,ix) && (button == 1)),
            fprintf('         -> Increasing border # at (%d, %d) to %d\n',ix,iy,gb_ex(loc_j,ix)+1);
            gb_ex(loc_j,ix) = gb_ex(loc_j,ix)+1; 
        elseif ((gb_ex(loc_j,ix) >= 1) && (button == 3)),
            fprintf('         -> Decreasing border # at (%d, %d) to %d\n',ix,iy,gb_ex(loc_j,ix)-1);
            gb_ex(loc_j,ix) = gb_ex(loc_j,ix)-1; 
        else
            fprintf('       - Cannot change cell at (%d, %d)\n',ix,iy);
        end
    end
    clf 
end
% return new k1
grid_borders = gb_ex(2:end,1:end-1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
