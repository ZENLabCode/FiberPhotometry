% findFiles
%-------------------------------------------------------------------------
% Checks if file(s) exist or automatic searches for files.
%
% SYNTAX
%    [list,removed] = findFiles(rFiles)
%       If rFiles is of type cell,
%           output list are the existing files from rFiles
%           output removed are the in-existing files from rFiles
%       If rFiles is of type char (output removed will be empyt)
%           output list = rFiles if rFiles is an existing file
%           output list will be a file list of automatically searched
%             files (using with wildcard, e.g. 'C:\data\*.mat)
%       PS: all outputs are of type cell.
%
%
% Thomas Rusterholz, 05 Nov 2020
%-------------------------------------------------------------------------
function [list,removed] = findFiles(rFiles)

removed = {}; %init
if iscell(rFiles)
    ind = cellfun(@(x)exist(x,'file')==2,rFiles);
    list    = rFiles(ind);
    removed = rFiles(~ind);
elseif ischar(rFiles)
    if exist(rFiles,'file')==2
        list = {rFiles};
    else
        [stat,list] = dos(sprintf('dir "%s" /S/B/A:-H-D',rFiles));
        if stat==1
            list = {};
        else
            %split list
            list = regexprep(list,sprintf('\r'),'');
            list = regexp(list,newline,'split');
            list(cellfun(@isempty,list)) = [];
            %add path if missing (just in case is missing)
            if isempty(fileparts(list{1}))
                list = fullfile(fileparts(rFiles),list);
            end
            %exclude dot-files (hidden mac-files)
            ind = cellfun(@(x)x(find(x==filesep,1,'last')+1)=='.',list);
            list(ind) = [];
            list = list(:);
        end
    end
else
    error('Class ''%s'' not supported for variable rFiles',class(rFiles))
end
end

