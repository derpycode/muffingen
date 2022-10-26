function [grid_k1]  = fun_grid_edit_k1(grid_k1,k_max)
%
% basic mask alteration code after:
% Andrew Yool (axy@soc.soton.ac.uk), October 2004.
%
%%

% *********************************************************************** %
% *** USER EDITING OF k1 (TOPOGRAPY) ************************************ %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_k1);
kmax = k_max;
% local copy
% NOTE: extend to the N and E for successful plotting using pcolor ...
%       (duplicate E column and N row)
gk_ex = grid_k1;
gk_ex = [gk_ex gk_ex(:,end)];
gk_ex = [gk_ex(1,:); gk_ex];
% Do the topography alteration
fprintf('       * Topography alteration procedure :\n');
fprintf('         (1) left button to deepen grid cell\n');
fprintf('         (2) right button to shallow grid cell\n');
fprintf('         (3) middle button, or clickoutside the grid, to finish\n');
%cmap = [0.5 1 0.5; 0.5 0.5 1]; colormap (cmap);
%
colormap default; 
%
flag = 0;
while flag == 0
%     figure(1); clf
%     pcolor(flipud(gk_ex)); % axis image;
%     caxis([0.5 (kmax+1.5)]); h = colorbar ('horiz');
%     set(h,'XTick',1:kmax+1);
%     title ('Ocean land/sea mask and bathymetry');
    figure(1); clf
    loc_gk_ex = flipud(gk_ex); loc_gk_ex(find(loc_gk_ex>kmax)) = NaN; % create array for plotting only
    pcolor(loc_gk_ex); % axis image;
    h = colorbar ('horiz');
    set(h,'XTick',1:kmax);
    title ('Ocean bathymetry');
    [x,y,button] = ginput(1);
    ix = floor(x); iy = floor(y);
    if ix < 1 | ix > imax | iy < 1 | iy > jmax
        fprintf('       - Out of grid range\n');
        flag = 1;
        fprintf('       * Mask alteration complete\n');
    elseif button == 2
        flag = 1;
        fprintf('       * Topography alteration complete\n');
    else
        loc_j = jmax-iy+2;
        if ((gk_ex(loc_j,ix) > 1) && (gk_ex(loc_j,ix) <= kmax) && (button == 1))
            fprintf('         -> Deepening cell at (%d, %d) to k=%d\n',ix,iy,gk_ex(loc_j,ix)-1);
            gk_ex(loc_j,ix) = gk_ex(loc_j,ix)-1; 
        elseif ((gk_ex(loc_j,ix) < kmax) && (button == 3))
            fprintf('         -> Shallowing cell at (%d, %d) to k=%d\n',ix,iy,gk_ex(loc_j,ix)+1);
            gk_ex(loc_j,ix) = gk_ex(loc_j,ix)+1; 
        else
            fprintf('       - Cannot deepen/shallow cell at (%d, %d)\n',ix,iy);
        end
    end
    clf 
end
% return new k1
grid_k1 = gk_ex(2:end,1:end-1);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
