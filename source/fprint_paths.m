function [] = fprint_paths(PNPATHS,PVPATHS,PNAME)
% fprint_paths()
%%

% initialize local variables
n_paths  = PNPATHS;
v_paths  = PVPATHS;
filename = PNAME;
% determine total number of paths
n_max = length(n_paths);
% open file for output
fid = fopen(filename,'w');
% write data
% NOTE: only write output if there is more than 1 path
if n_max >= 2,
    % write paths length data
    fprintf(fid,'%4i',n_paths(2:end));
    fprintf(fid,'\n');
    % initialize paths start & end
    m_start = 1;
    m_end   = n_paths(1);
    for n = 2:n_max
        % update paths start & end
        m_start = m_start + n_paths(n-1);
        m_end   = m_end + n_paths(n);
        % add blank line
        fprintf(fid,'\n');
        % write paths data
        for m = m_start:m_end
            fprintf(fid,'%4i',v_paths(m,:));
            fprintf(fid,'\n');
        end
    end
else
    % PRINT NOTHING!
    fprintf(fid,'\n');
end
% close file
fclose(fid);
