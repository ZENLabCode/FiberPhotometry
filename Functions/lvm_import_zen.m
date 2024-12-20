% lvm_import_zen
%------------------------------------------------------------------------
% Imports data from a LabView lvm-file (zen-lab version).
%
% SYNTAX
%   [data,info,bad] = lvm_import_zen(file);
%     Input:
%       - file : lvm-filename
%     Output
%       - data : read data as a structure
%         info : reading info (if some special cases happened)
%         bad  : cell with sample correction infos
%                 - 1st column, sample index
%                 - 2nd column, text from where the script tried to
%                   evaluate missing data samples.
%
% NOTES:
%  - This script is an adaptation from lvm_import.m
%    Some lvm-files failed to read with lvm_import, because files were
%    corrupted (basically LabView script needs to be corrected).
%        The reason was, that some numeric values had some non-numeric
%    charachters (e.g: comma, double point, semicolon, vertical tab, ...).
%    This made MATLAB's textscan stop reading and script continued with a
%    2nd segment but missing header infos (files had only one segment).
%        This script will care for that (hopefully correct). It tries to
%    evaluate the correct values and continues reading the remaining data.
%    Fortunately, this happend very seldom and only affected single
%    samples (I once found two samples in a 512 Hz, 10-h recording). So if
%    correction is wrong, I can live with that.
%    NOTE: correction will fail in some 
%  - Additionally, this script will decline comma separated digit values.
%    Maybe I will implement it later. However, it's anyway bad style using 
%    comma separated digits working with data on computers.
%  - Additional infos for lvm-files, see:
%    https://www.ni.com/de-ch/support/documentation/supplemental/06/specification-for-the-labview-measurement-file---lvm-.html
%  - Output data.clock is a datevec (seconds are rounded down)
%
%
% Thomas Rusterholz, 7 Oct 2020
%-------------------------------------------------------------------------

function [data,info,bad] = lvm_import_zen(file)

%INIT
%defaults
dlm  = '\t'; %delimiter (hopefully)
info = '';   %info string, for failures and other informations
bad = cell(0,2); %{sample, bad strings for numeric values}
formDate = 'yyyy/mm/dd'; %date format (for clock)
formTime = 'HH:MM:SS';   %time format (for clock), no [ms]
%open file
fid = fopen(file,'r');
if fid==-1
    info = 'File not found';
    data = []; fclose(fid); return
end
%find first line, ignoring empty lines
str = '';
while isempty(str) && ~feof(fid)
    str = strtrim(fgetl(fid));
end
if isempty(str) && feof(fid)
    info = 'File has no data.';
    data = []; fclose(fid); return
end

%FILE WITHOUT HEADER
if ismember(str(1),'-+.0123456789')
    info = 'File has no header information, better check output!';
    fclose(fid);
    data.Segment1.data = dlmread(file,dlm); return
elseif ~strcmpi(str,'LabVIEW Measurement')
    info = 'This is not a correct lvm-file.';
    data = []; fclose(fid); return
end

%MAIN HEADER
hdr = {}; %header lines
while ~strcmpi(str,'***End_of_Header***') && ~feof(fid)
    hdr{end+1,1} = str;
    str = strtrim(fgetl(fid));
end
str = strtrim(fgetl(fid)); %next line (after '***End_of_Header***')
hdr(cellfun(@isempty,hdr)) = []; %remove empty lines
%delimiter
tmp = hdr{strncmpi(hdr,'Separator',numel('Separator'))};
if contains(lower(tmp),'tab')
    dlm = '\t';
elseif contains(lower(tmp),'comma')
    dlm = ',';
end
%header fields
for k = 1:numel(hdr)
    tmp = strtrim(regexp(hdr{k},dlm,'split'));
    data.(strrep(tmp{1},' ','_')) = strjoin(tmp(2:end),dlm);
end
%decimal separator
if isfield(data,'Decimal_Separator') && data.Decimal_Separator~='.'
    info = sprintf(['%s.m only works with dot as decimal separator, ',...
        'not with ''%s''.'],mfilename,data.Decimal_Separator);
    data = []; fclose(fid); return
end
%clock field (good to have), rounded down
if all(isfield(data,{'Time','Date'}))
    tmp = sprintf('%s %s',data.Date,data.Time);
    data.clock = datevec(tmp,[formDate,' ',formTime]);
end

%SEGMENT LOOP
seg = 0; %count segments
while ~feof(fid)
    clear dat; %segment data
    seg = seg+1;
    
    %SEGMENT HEADER
    hdr = {}; %header lines
    while ~strcmpi(str,'***End_of_Header***') && ~feof(fid)
        hdr{end+1,1} = str;
        str = strtrim(fgetl(fid));
    end
    hdr(cellfun(@isempty,hdr)) = []; %remove empty lines
    %ignore Specials
    ind1 = find(strcmpi(str,'***Start_Special***'));
    ind2 = find(strcmpi(str,'***End_Special***'));
    for k = numel(ind1):-1:1
        hdr(ind1(k):ind2(k)) = [];
    end
    %append
    clear dat
    for k = 1:numel(hdr)
        vals = strtrim(regexp(hdr{k},dlm,'split')); %label and values
        label = vals{1}; vals(1) = [];
        vals(cellfun(@isempty,vals)) = []; %double tabs
        nums = cellfun(@str2double,vals);
        if all(~isnan(nums))
            dat.(label) = nums;
        else
            dat.(label) = vals;
        end
    end
    if isfield(dat,'Delta_X') && ~isfield(dat,'fs')
        dat.fs = round(1./dat.Delta_X); %add sampling rates
    end
    
    %DATA
    %column labels
    while (strcmpi(str,'***End_of_Header***')||isempty(str)) && ~feof(fid)
        str = strtrim(fgetl(fid));
    end
    labels = strtrim(regexp(str,dlm,'split'));
    indC = strcmpi(labels,'comment');
    indX = strcmpi(labels,'X_Value');
    indY = ~indC & ~indX;
    %values
    noCHA = numel(labels);
    form = cell(noCHA,1); %init, read format
    form(:) = {'%f'}; form(indC) = {'%s'};
    form = strjoin(form);
    values = textscan(fid,form,'delimiter',dlm);
    
    %CORRUPT FILE CORRECTION (the big magic)
    noSAMS = cellfun(@numel,values); %count samples per channel
    minSAM = min(noSAMS);
    maxSAM = max(noSAMS);
    while maxSAM~=minSAM
        if feof(fid)
            %if diff==1, probably stopped  before saving was finishid
            info = strtrim(sprintf('%s %i samples removed',info,...
                maxSAM-minSAM));
            values = cellfun(@(x)x(1:minSAM),values,'uniformoutput',false);
            if maxSAM-minSAM>1 %should never happen
                warning(['Check data file, channels had different ',...
                    'amount of samples (cut to min). Max difference ',...
                    'is %i'],...
                    maxSAM-minSAM)
            end
        else %works in the files tested, but no guarantee
            if maxSAM-minSAM>1
                error('grrrhhhh')
            end
            info = strtrim(sprintf('%s Corrected characters.',info));            
            %read missing values
            str = fgetl(fid);
            bad(end+1,:) = {maxSAM,str};
            str = strtrim(regexp(str,dlm,'split'));
            %correct str size (hopefully)
            ind0 = find(noSAMS~=maxSAM,1,'first');
            d = numel(str)+ind0-1-noCHA;
            if d<0
                str = [str,repmat({''},1,abs(d))];
            elseif d>0
                str(1:d) = [];
            end
            %class
            indNum = find(indX|indY); %channels with numeric values
            for k = 1:numel(str)
                if ismember(k+ind0-1,indNum)
                    a = str{k};
                    b = str2double(a(ismember(a,'-+.0123456789')));
                    if isnan(b) %2nd try
                        b = str2double(a(ismember(a,'-+.0123456789eE')));
                    end
                    str{k} = b;
                end
            end
            %correct value size
            for cha = ind0:noCHA
                values{cha} = [values{cha};str{cha-ind0+1}];
            end
            %append next data
            tmp = textscan(fid,form,'delimiter',dlm);
            for cha = 1:noCHA
                values{cha} = [values{cha};tmp{cha}];
            end
        end
        %re-check
        noSAMS = cellfun(@numel,values);
        minSAM = min(noSAMS);
        maxSAM = max(noSAMS);
    end
    
    %APPEND
    dat.column_labels = labels(indY);
    tmp = cell2mat(values(indX));
    if all(isnan(tmp(:)))
        dat.dataX = [];
    else
        dat.dataX = tmp;
    end
    dat.dataY = cell2mat(values(indY));
    tmp = [values{indC}];
    if all(cellfun(@isempty,tmp(:)))
        dat.comments = {};
    else
        dat.comments = tmp;
    end
    data.(sprintf('Segment%i',seg)) = dat;
    
    %CHECK (number of columns, works so far)
    f = fieldnames(dat); %all fields
    f(ismember(lower(f),... remove fields that are/might be different
        {'channels','comments','datax'})) = [];
    for k = 1:numel(f)
        if size(dat.(f{k}),2)~=dat.Channels
            tmp = {... {label, label format, data, data format}
                'File'    ,' %-s : ' ,file ,'%s';...
                'Segment' ,' %-s : ' ,seg  ,'%i';...
                'Field'   ,' %-s : ' ,f{k} ,'%s';...
                };
            n = num2str(-max(cellfun(@numel,tmp(:,1))));
            tmp(:,2) = cellfun(@(x)strrep(x,'-',n),tmp(:,2),...
                'uniformoutput',false);
            tmp = tmp';
            tmp = sprintf(sprintf('%s%s\n',tmp{2:2:end,:}),tmp{1:2:end,:});
            warning(['Check, number of columns should be %i, ',...
                'but is %i\n%s'],...
                dat.Channels,size(dat.(f{k}),2),tmp)
        end
    end
end %segement loop
fclose(fid);
end