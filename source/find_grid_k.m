function [grid_k1]  = find_grid_k(par_min_Dk,v_dm,v_de,grid_mask,grid_topo)
%
%%

% *********************************************************************** %
% *** CONVERT TOPOGRAPHY INTO DEPTH LEVELS ****************************** %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_topo);
kmax = length(v_dm);
% create top and bottom boundary vectors
v_dbot = v_de(1:kmax);
v_dtop = v_de(2:kmax+1);
% first process topo:
% (i) apply mask
% (ii) turn negative heights into positive depths
% (iii) filter any over-deep cells
gd = grid_mask.*grid_topo;
gd = abs(gd);
gd(find(gd>v_dbot(1))) = v_dbot(1);
% initialzie output grid
grid_k1 = zeros(jmax,imax);
% re-grid depth
for j = 1:jmax
    for i = 1:imax
        if (gd(j,i) > 0.0),
            loc_k = intersect( find(gd(j,i)<=v_dbot), find(gd(j,i)>v_dtop) );
            if loc_k > (kmax - (par_min_Dk - 1))
                grid_k1(j,i) = kmax - (par_min_Dk - 1);
            else
                grid_k1(j,i) = loc_k;
            end
        end
    end
end
grid_k1(find(grid_k1==0))=90;
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
