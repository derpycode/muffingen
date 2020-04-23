function [grid_borders,opt_user]  = find_grid_borders_check(grid_borders,grid_mask,opt_user)
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
vdsrch_diag = [1 1; 1 -1; -1 -1; -1 1];
% copy & expand grid
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0;
gb_ex(end,:) = 0;
% copy & expand grid
gm_ex = grid_mask;
gm_ex = [gm_ex(:,end) gm_ex gm_ex(:,1)];
gm_ex = [gm_ex(1,:); gm_ex; gm_ex(end,:)];
gm_ex(1,:)   = 0;
gm_ex(end,:) = 0;
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
            % NOTE: IGNORE the first border as this will not become a path
            if ((isrch >= 3) && (gb_ex(j,i) > 1))
                disp([' *** Problem with island path @ (',num2str(loc_i),',',num2str(jmax-loc_j+2),') ...']);
                opt_user = true;
                if (isrch == 3)
                    isrch  = 0;
                    % search surrounding cells
                    for s = 1:length(vdsrch_diag)
                        loc_j = j + vdsrch_diag(s,1);
                        loc_i = i + vdsrch_diag(s,2);
                        % test adjacent cells for unassigned border and ocean
                        if (gb_ex(loc_j,loc_i)==0 && gm_ex(loc_j,loc_i)==1 )
                            loc_gb_ex = gb_ex;
                            loc_gb_ex(loc_j,loc_i) = loc_gb_ex(j,i);
                            loc_gb_ex(j,i) = 0;
                            [loc_user] = find_grid_borders_recheck(loc_gb_ex,false);
                            if ~loc_user
                                gb_ex(loc_j,loc_i) = gb_ex(j,i);
                                gb_ex(j,i) = 0;
                                grid_borders(loc_j-1,loc_i-1) = grid_borders(j-1,i-1);
                                grid_borders(j-1,i-1) = 0;
                                disp(['     ... junction was only a tripple junction and has now been fixed ... hopefully!']);
                                opt_user = false;
                                break
                            else
                                isrch = isrch + 1;
                            end
                        end
                    end
                else
                    disp([' *** User intevention required @ (',num2str(loc_i),',',num2str(jmax-loc_j+2),') to resolve problem with island path.']);
                    opt_user = true;
                end
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
