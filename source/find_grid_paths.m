function [n_paths,v_paths,n_islands,grid_paths] = find_grid_paths(grid_borders,n_islands,i_poles);
%
%%

% *********************************************************************** %
% *** ELUCIDATE PATHS *************************************************** %
% *********************************************************************** %
%
% NOTE: polar borders can be incomplete, hence the pole must be followed
%       in the same was as an island is
%       (special borders are created to enable this to occur automatically)
% NOTE: don't worry initially about the rows being counted from
%       top-top-bottom (opposite of GENIE grid)
%       => this is fixed at the very end
% NOTE: format is of first parameter being the direction to the *next*
%       cell
% NOTE: GOLDSTEIN directions:
%        2 == North (+ve v)
%       -2 == South (-ve v)
%        1 == East  (+ve u)
%       -1 == West  (-ve u)
% determine mask size (remember: [rows columns])
[jmax imax] = size(grid_borders);
% create search array
% directions: 1==S, 2==E, 3==N, 4==W
vdsrch = [1 0; 0 1; -1 0; 0 -1];
% copy & expand border grid
gb_ex = grid_borders;
gb_ex = [gb_ex(:,end) gb_ex gb_ex(:,1)];
gb_ex = [gb_ex(1,:); gb_ex; gb_ex(end,:)];
gb_ex(1,:)   = 0.0;
gb_ex(end,:) = 0.0;
% created 'searched' grid
gs_ex = zeros(jmax+2,imax+2);
% initialize paths arrays
v_paths = [];
n_paths = [];
% create naked polar island borders
if (~isempty(i_poles)),
    % N pole
    if ~isempty(find(i_poles == -1)),
        n_islands = n_islands + 1;
        gb_ex(2,:) = n_islands;
    end
    % S pole
    if ~isempty(find(i_poles == -2)),
        n_islands = n_islands + 1;
        gb_ex(end-1,:) = n_islands;
    end
end
%
if (n_islands >= 2)
    %
    disp(['       * Ignoring border #1']);
    % set dummy first path data
    v_paths = [v_paths; 0 0 0];
    n_paths = [n_paths 1];
    % LOOP >>>
    for islnd = 2:n_islands
        %
        disp(['       * Creating raw path #' num2str(islnd) ' ...']);
        % find top LH border cell of the current island border
        % raster left-to-right, top-to-bottom
        corner = false;
        for j = 2:jmax+1
            for i = 2:imax+1
                if (gb_ex(j,i) == islnd), corner = true; end
                if corner, break; end
            end
            if corner, break; end
        end
        %
        loc_j = j;
        loc_i = i;
        % mark initial cell as searched
        gs_ex(loc_j,loc_i) = 1;
        % test for E-W wall and also mark wrap-around cell as searched
        if (loc_i == 2), gs_ex(loc_j,imax+2) = 1; end
        if (loc_i == imax+1), gs_ex(loc_j,1) = 1; end
        % check assumption of East being A-OK ...
        if (gb_ex(loc_j,loc_i+1) == islnd),
            % initial direction to the next cell *must* be East [1] [assumption!]
            % => take first step in that direction
            % record direction and current location (in core grid indices)
            v_paths = [v_paths; 1 loc_i-1 loc_j-1];
            % now move 1 East
            loc_i = loc_i + 1;
            % test for E-W wall ...
            % mark i==1 cell as implicitly, already searched
            if (loc_i == imax+2),
                loc_i = 2;
                gs_ex(loc_j,loc_i)  = 1;
                gs_ex(loc_j,imax+2) = 1;
            end
        elseif (gb_ex(loc_j+1,loc_i) == islnd),
            disp('       ! Initial E direction path follow step failed ...');
            disp('         ... trying the S direction ...');
            % try South ...
            % record direction and current location (in core grid indices)
            v_paths = [v_paths; -2 loc_i-1 loc_j-1];
            % now move 1 South
            loc_j = loc_j + 1;
        else
            disp(' *** Path follow failed (both E and S directions) :(');
            diary off;
            return;
        end
        % mark as searched
        gs_ex(loc_j,loc_i) = 1;
        % initialize vector length at 1
        % (as the vector has already been populaed with its first line)
        n_path = 1;
        % now follow path around island -- clockwise
        follow = true;
        while follow
            % find adajacent, unmarked border cell
            % NOTE: if none exist, finish!
            % search surrounding cells
            follow = false;
            for s = 1:4
                loc_jj = loc_j + vdsrch(s,1);
                loc_ii = loc_i + vdsrch(s,2);
                % test for adjacent border cell
                if (gb_ex(loc_jj,loc_ii) == islnd) && ~gs_ex(loc_jj,loc_ii),
                    % record current location and direction to next cell
                    % directions: 1==S, 2==E, 3==N, 4==W
                    % reminder:
                    %        2 == North
                    %       -2 == South
                    %        1 == East
                    %       -1 == West
                    % NOTE: take into account expanded i,j indices and
                    %       record location in core grid coordinates
                    switch s
                        case (1)
                            v_paths = [v_paths; -2 loc_i-1 loc_j-1];
                        case (2)
                            v_paths = [v_paths;  1 loc_i-1 loc_j-1];
                        case (3)
                            v_paths = [v_paths;  2 loc_i-1 loc_j-1];
                        case (4)
                            v_paths = [v_paths; -1 loc_i-1 loc_j-1];
                    end
                    % update path length count
                    n_path = n_path + 1;
                    % mark tested cell as searched
                    loc_j = loc_jj;
                    loc_i = loc_ii;
                    gs_ex(loc_j,loc_i) = 1;
                    % test for E-W wall:
                    % adjust (j,i) location if necessary and also
                    % mark as searched
                    if (loc_ii == 1)
                        gs_ex(loc_j,imax+1) = 1; % wrap-around location
                        loc_i = imax+1;
                        gs_ex(loc_j,imax+2) = 1; % implicitly, already searched
                    elseif (loc_ii == 2)
                        gs_ex(loc_j,imax+2) = 1; % wrap-around location
                    elseif (loc_ii == imax+2)
                        gs_ex(loc_j,2) = 1;      % wrap-around location
                        loc_i = 2;
                        gs_ex(loc_j,1) = 1;      % implicitly, already searched
                    elseif (loc_ii == imax+1)
                        gs_ex(loc_jj,1) = 1;     % wrap-around location
                    end
                    % continue ...
                    follow = true;
                    % exit (s) loop
                    break;
                end
            end
        end
        % add final location, calculate direction to start, update count
        % NOTE: take into account upside-down GENIE array indexing in MATLAB
        % reminder:
        %        2 == North
        %       -2 == South
        %        1 == East
        %       -1 == West
        if ((loc_i-i+1) == imax) || ((loc_i-i) == -1), %East
            v_paths = [v_paths;  1 loc_i-1 loc_j-1];
        elseif ((loc_i-i-1) == -imax) || ((loc_i-i) == 1), %West
            v_paths = [v_paths; -1 loc_i-1 loc_j-1];
        elseif (loc_j-j) == 1, %North
            v_paths = [v_paths;  2 loc_i-1 loc_j-1];
        elseif (loc_j-j) == -1, %South
            v_paths = [v_paths; -2 loc_i-1 loc_j-1];
        else
            % NOTE: The requirement for 2-cell seperation could be fixed
            %       (removed) in a future release.
            disp([' *** You may have insufficient seperation between the South Pole and a land mass:']);
            disp(['     => create a 2 cell seperation of ocean, or join landmass to S. Pole.']);
            disp(['       (but there could be other issues ...)']);
            disp([' ']);
            diary off;
            error(['Error. \nFailed to complete path loop @ (',num2str(loc_i),',',num2str(loc_j),'): %s'],'Exiting ...');
            return;
        end
        n_path = n_path + 1;
        % write out path length
        n_paths = [n_paths n_path];
        %
    end
    % <<< LOOP
end
% now fix the up-side-down World
for p = 1:sum(n_paths),
    v_paths(p,3) = jmax-v_paths(p,3)+1;
end
%
% *** PATH POST-PROCESSING ********************************************** %
%
% Included in the automatically generated paths, are 2 illegal moves that
% result in the c-grid being transversed twice
% (and hence the path direction ambigeous).
% There are:
% The top right corner cell of a (bend in a) path,
% whether an outer or innter bend, and characterized by:
% path cells to the E and the S
% (and ocean and/or land to the W and the N)
% Both other corner configurations are valid.
%
% The search is basically for:
% (a) a '+1' followed by a '-2'
% (b) a '+2' followed by a '-1'
% (and remembering to deal with the wrap-around E-W grid)
%
% NOTE: don't bother replacing the 1st path (not written out!)
%
% REMEMBER:
%
%   ----v----
%   |       |
%   |   c   u
%   |       |
%   ---------
%
%        2 == North (+ve v)
%       -2 == South (-ve v)
%        1 == East  (+ve u)
%       -1 == West  (-ve u)
%
if n_islands >= 2
    %
    for islnd = 2:n_islands
        %
        disp(['       * Building path #' num2str(islnd) ' ...']);
        %
        pmin = sum(n_paths(1:islnd-1))+1;
        pmax = sum(n_paths(1:islnd));
        dp = pmax-pmin+1;
        % create tempoary wrap-around array
        % (in terms of following the path loop)
        vp_ex = v_paths(pmin:pmax,:);
        vp_ex = [vp_ex(end,:); vp_ex; vp_ex(1,:)];
        %
        vpp = [];
        %
        p = 2;
        pp = 1;
        while p <= dp+1,
            if (vp_ex(p-1,1)==1 && vp_ex(p,1)==-2) || (vp_ex(p-1,1)==2 && vp_ex(p,1)==-1),
                % top RH corner, clockwise
                % [p-1] == 1; [p] == -2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ----------v------
                % !-2 >-1 > p ! o !   !-2 >-1 >   ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 1 ! o !   ! X ! X ! 1 ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 2 ! o !   ! X ! X ! 2 ! o !
                % ----------v------   -----------------
                % top RH corner, anticlockwise
                % [p-1] == 2; [p] == -1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % < 2 < 1 < p ! o !   ! 2 < 1 <   < o !
                % ----------^------   ----------^------
                % ! X ! X !-1 ! o !   ! X ! X !-1 ! o !
                % ----------^------   ----------^------
                % ! X ! X !-2 ! o !   ! X ! X !-2 ! o !
                % -----------------   -----------------
                disp(['         -> NE corner :: ' 'Skip path entry @ (' num2str(vp_ex(p,2)) ',' num2str(vp_ex(p,3)) ')']);
                % NOTE: add no path component (or pp count update)
                % update p count
                p = p+1;
            elseif (vp_ex(p-1,1)==-2 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==-1 && vp_ex(p,1)==2),
                % bottom LH corner, anticlockwise
                % [p-1] == -2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ------v----------
                % ! o !-2 ! X ! X !   ! o !-2 ! X ! X !
                % ------v----------   ------v----------
                % ! o !-2 ! X ! X !   ! o !-1 ! X ! X !
                % ------v----------   ------v----------
                % ! o ! p > 1 > 2 >   ! o !pp > 1 > 2 >
                % -----------------   -----------------
                % bottom LH corner, clockwise
                % [p-1] == -1; [p] == 2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ------^----------   ------^----------
                % ! o ! 2 ! X ! X !   ! o ! 2 ! X ! X !
                % ------^----------   ------^----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! X ! X !
                % ------^----------   ------^----------
                % ! o ! p <-1 <-2 !   ! o !pp <-1 <-2 <
                % -----------------   -----------------
                disp(['         -> SE corner :: ' 'Add additional path entry @ (' num2str(vp_ex(p,2)) ',' num2str(vp_ex(p,3)) ')']);
                % duplicate path component
                % NOTE: adjust vector in first entry to complete path turn
                vpp(pp,:) = vp_ex(p,:);
                vpp(pp,1) = vp_ex(p-1,1);
                pp = pp+1;
                vpp(pp,:) = vp_ex(p,:);
                pp = pp+1;
                % update p count
                p = p+1;
            elseif (vp_ex(p-1,1)==2 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==1 && vp_ex(p,1)==2) || (vp_ex(p-1,1)==1 && vp_ex(p,1)==1) || (vp_ex(p-1,1)==2 && vp_ex(p,1)==2),
                % top LH corner, clockwise
                % [p-1] == 2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % ! o ! p > 1 > 2 >   ! o ! p > 1 > 2 >
                % ------^----------   ------^----------
                % ! o !-1 ! X ! X !   ! o !-1 ! X ! X !
                % ------^----------   ------^----------
                % ! o !-2 ! X ! X !   ! o !-2 ! X ! X !
                % -----------------   -----------------
                % bottom RH corner, anticlockwise
                % [p-1] == 2; [p] == 1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ----------^------   ----------^------
                % ! X ! X !-2 ! o !   ! X ! X !-2 ! o !
                % ----------^------   ----------^------
                % ! X ! X !-1 ! o !   ! X ! X !-1 ! o !
                % ----------^------   ----------^------
                % ! 2 > 1 > p ! o !   ! 2 > 1 > p ! o !
                % -----------------   -----------------
                % East (right), North (up)
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % ------^----------   ------^----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! o ! o !
                % ------^----------   ------^----------
                % !-1 > p > 1 > 1 >   !-1 > p > 1 > 2 >
                % ------^----------   ------^----------
                % ! o !-1 ! o ! o !   ! o !-1 ! o ! o !
                % -----------------   -----------------
                vpp(pp,:) = vp_ex(p,:);
                pp = pp+1;
                p = p+1;
            elseif (vp_ex(p-1,1)==-1 && vp_ex(p,1)==-2) || (vp_ex(p-1,1)==-2 && vp_ex(p,1)==-1) || (vp_ex(p-1,1)==-1 && vp_ex(p,1)==-1) || (vp_ex(p-1,1)==-2 && vp_ex(p,1)==-2),
                % top LH corner, anticlockwise
                % [p-1] == -1; [p] == -2
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   -----------------
                % ! o ! p <-1 <-2 !   ! o ! p <-1 <-2 <
                % ------v----------   ------v----------
                % ! o ! 1 ! X ! X !   ! o ! 1 ! X ! X !
                % ------v----------   ------v----------
                % ! o ! 2 ! X ! X !   ! o ! 2 ! X ! X !
                % ------v----------   -----------------
                % bottom RH corner, clockwise
                % [p-1] == -2; [p] == -1
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ----------v------
                % ! X ! X ! 2 ! o !   ! X ! X ! 2 ! o !
                % ----------v------   ----------v------
                % ! X ! X ! 1 ! o !   ! X ! X ! 1 ! o !
                % ----------v------   ----------v------
                % <-2 <-1 < p ! o !   !-2 <-1 < p ! o !
                % -----------------   -----------------
                % West (left), South (down)
                % BORDER DIRECTIONS   GOLDSTEIN VELOCITY
                % -----------------   ------v----------
                % ! o !-1 ! X ! X !   ! o !-1 ! o ! o !
                % ------V----------   ------v----------
                % < 1 < p <-1 <-1 !   ! 1 < p <-1 <-2 <
                % ------V----------   ------v----------
                % ! o ! 1 ! o ! o !   ! o ! 1 ! o ! o !
                % ------V----------   -----------------
                vpp(pp,:) = vp_ex(p,:);
                vpp(pp,1) = vp_ex(p-1,1);
                pp = pp+1;
                p = p+1;
            else
                disp(' *** IMPOSSIBLE!');
                return;
            end
            %
        end
        % add new paths back in array
        v_paths(pmin:pmax,:) = vpp(:,:);
    end
    %
end
%
%
% *** CREATE 2D PATHS MAP *********************************************** %
%
% create empty array
grid_paths = zeros(jmax,imax);
% populate array
% NOTE: don't bother outputting 1st path
if n_islands >= 2,
    for p = (n_paths(1)+1):sum(n_paths),
        grid_paths(v_paths(p,3),v_paths(p,2)) = v_paths(p,1);
    end
end
% re-orientate
grid_paths = flipud(grid_paths);
%
% *********************************************************************** %

% *********************************************************************** %
% *** END *************************************************************** %
% ***********************************************************************
%
%%%
%
% *********************************************************************** %
