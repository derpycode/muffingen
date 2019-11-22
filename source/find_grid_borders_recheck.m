function [opt_user]  = find_grid_borders_recheck(grid_borders,opt_user)
%
%%

% *********************************************************************** %
% *** CHECK ISLANDS BORDERS ********************************************* %
% *********************************************************************** %
%
% determine mask size (remember: [rows columns])
[jmax, imax] = size(grid_borders);
% create search array
vdsrch_nsew = [1 0; 0 1; -1 0; 0 -1];
% copy & expand grid
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0;
gb_ex(end,:) = 0;
%
% *** SEARCH FOR TRIPPLE POINTS! **************************************** %
%
% search across grid -- find any tripple points
for i = 2:imax+1
    for j = 2:jmax+1
        % test for cell being a border
        if gb_ex(j,i) > 0
            isrch  = 0;
            % search surrounding cells
            for s = 1:length(vdsrch_nsew)
                loc_j = j + vdsrch_nsew(s,1);
                loc_i = i + vdsrch_nsew(s,2);
                % test adjacent cells
                % (if the same border #)
                if (gb_ex(j,i) == gb_ex(loc_j,loc_i))
                    isrch  = isrch + 1;
                end
            end
            % report issue
            % NOTE: IRNORE the first border as this will not become a path
            if ((isrch >= 3) && (gb_ex(j,i) > 1))
                opt_user = true;
            end
        end
    end
end
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% *********************************************************************** %
%
%%%
%
% *********************************************************************** %
