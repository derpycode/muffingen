function [grid_topo]  = fun_conv_k1(v_de,grid_k1)
%
%%

% *********************************************************************** %
% *** CONVERT DEPTH LEVELS INTO TOPOGRAPHY ****************************** %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_k1);
kmax = length(v_de)-1;
% local copy
gk = grid_k1;
% initialzie output grid
grid_topo = zeros(jmax,imax);
% de-grid k1
for j = 1:jmax
    for i = 1:imax
        if (gk(j,i) < 90),
            grid_topo(j,i) = v_de(gk(j,i));
        end
    end
end
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
