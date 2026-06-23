function [idvx, idx_set, ide_const] = read_landmarks(filename, Src)
%READINTEGERFILE Reads a custom-formatted file containing integers.
%   [idx_line, idx_cycle, idvx] = readIntegerFile(filename)
%
%   - Lines starting with '#' are comments → ignored.
%   - Lines starting with 'l' contain integers → stored in cell array idx_line.
%   - Lines starting with 'c' contain integers → stored in cell array idx_cycle.
%   - Lines starting directly with an integer → stored in numeric array idvx.

    % Initialize outputs
    idx_line = {};
    idx_cycle = {};
    idx_set = {};
    idvx = [];
    ide_const = [];

    % Open file
    fid = fopen(filename, 'r');
    if fid == -1
        disp(['Could not open file: %s', filename]);
        disp('Assume that there is no landmarks.');
        return;
    end

    % Read file line-by-line
    line = fgetl(fid);
    while ischar(line)
        line = strtrim(line);

        % Skip empty lines
        if isempty(line)
            line = fgetl(fid);
            continue;
        end

        % Case 1: comment
        if startsWith(line, '#')
            line = fgetl(fid);
            continue;
        end

        % Case 2: starts with 'l'
        if startsWith(line, 'l')
            nums = sscanf(line(2:end), '%d')'; % read integers after 'l'
            idx_line{end+1} = nums' + 1;
        
        % Case 3: starts with 'c'
        elseif startsWith(line, 'c')
            nums = sscanf(line(2:end), '%d')'; % read integers after 'c'
            idx_cycle{end+1} = nums' + 1;

        % Case 4: starts with an integer → idvx
        else
            % verify line starts with a digit or minus sign
            if ~isempty(regexp(line(1), '[0-9-]', 'once'))
                nums = sscanf(line, '%d')'; 
                idvx = [idvx; nums' + 1];  % append to numeric array
            else
                warning(['Unrecognized line format: ', line]);
            end
        end

        line = fgetl(fid);
    end

    fclose(fid);

    idx_set = [idx_line; idx_cycle];

    if ~isempty(idx_set) && nargout > 2
        % Compute edge IDs
        n = length(idx_line) + length(idx_cycle);
        ide_const = cell(n,1);
        for i = 1:length(idx_line)
            P = sort([idx_line{i}(1:end-1), idx_line{i}(2:end)], 2);
            [~,~,ide] = intersect(P, Src.E2V, 'rows', 'stable');
            ide_const{i} = ide;
        end
        for i = 1:length(idx_cycle)
            P = sort([idx_line{i}, circshift(idx_line{i}, [1,0])], 2);
            [~,~,ide] = intersect(P, Src.E2V, 'rows', 'stable');
            ide_const{length(idx_cycle)+i} = ide;
        end
        ide_const = unique(cell2mat(ide_const));
    else
        ide_const = [];
    end
end
