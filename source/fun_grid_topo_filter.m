function [grid_k1] = fun_grid_topo_filter(grid_k1)
% fun_grid_topo_filter
%   fun_grid_topo_filter(K1) filter the ocean floor topography
%   to remove single cell 'holes'.
%%

% *********************************************************************** %
% *** FILTER k1 DATA **************************************************** %
% *********************************************************************** %
%
% determine grid size
[jmax imax] = size(grid_k1);
% adjust for extended grid
jmax = jmax + 2;
imax = imax + 2;
% determine kmax
kmax = max(max(grid_k1(find(grid_k1 < 90))));
% add boundaries to k1 and set as 'old' array
% NOTE: set N and S pole boundaries to a value (topo index) of 90
%       to ensure that pole-bordering single cells are caught
k1_old = grid_k1;
k1_old = [k1_old(:,end) k1_old k1_old(:,1)];
k1_old = [k1_old(1,:); k1_old; k1_old(end,:)];
k1_old(1,:)   = 90;
k1_old(end,:) = 90;
% initialize parameters
par_opt1 = true;
par_opt2 = true;
loop_count_max = 1000;
%  filter topo
% NOTE: as compared to original, now no need to check for shallow cells
%       (vertical re-gridding has already removed them)
C = true;
loop_count = 0;
k1_new = k1_old;
while C,
    C = false;
    loop_count = loop_count + 1;
    for i = 2:imax-1,
        for j = 2:jmax-1,
            % check for ocean cell
            if k1_new(j,i) < 90
                % first look around - check the cell values in all 4 directions
                k1_neigh(1) = k1_new(j-1,i);
                k1_neigh(2) = k1_new(j,i+1);
                k1_neigh(3) = k1_new(j+1,i);
                k1_neigh(4) = k1_new(j,i-1);
                % find land neighbours
                n_l = length(k1_neigh(find(k1_neigh >= 90)));
                % find shallower and deeper neighbours
                k1_neigh_o = k1_neigh(find(k1_neigh < 90));
                n_s = length(k1_neigh_o(find(k1_neigh_o > k1_new(j,i))));
                n_d = length(k1_neigh_o(find(k1_neigh_o < k1_new(j,i))));
                % STEP #1 -- fill in holes
                % find shallower neighbours + land
                if ( (n_s >= 3) || ((n_l == 2) && (n_s == 2) && par_opt1) || ((n_l == 1) && (n_s >= 2) && par_opt2) ),
                    disp(['         -> MAKE SHALLOWER :: ' 'Loop count = ' num2str(loop_count) ' @ (' num2str(i-1) ',' num2str(j-1) ',' num2str(k1_new(j,i)) ')' '; s = ' num2str(n_s) '; l = ' num2str(n_l)])
                    if (k1_new(j,i) < kmax),
                        k1_new(j,i) = k1_new(j,i) + 1;
                        C = true;
                    end
                end
                % STEP #2 -- smooth down peaks
                % find deeper neighbours + land
                if ( (n_d > 3) || ((n_l >= 2) && (n_d > 1)  && par_opt1) || ((n_l == 1) && (n_d == 3) && par_opt2) ),
                    disp(['         -> MAKE DEEPER :: ' 'Loop count = ' num2str(loop_count) ' @ (' num2str(i-1) ',' num2str(j-1) ',' num2str(k1_new(j,i)) ')' '; d = ' num2str(n_d) '; l = ' num2str(n_l)])
                    if (k1_new(j,i) > 1), 
                        k1_new(j,i) = k1_new(j,i) - 1; 
                        C = true;
                    end
                end
            end
        end
    end
    % update East and West grid pseudo-data
    k1_new(:,1) = k1_new(:,imax-1);
    k1_new(:,imax) = k1_new(:,2);
    %
    if (loop_count == loop_count_max), C = false; end
end
% return revised core k1 grid
grid_k1 = k1_new(2:end-1,2:end-1);
%
% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
