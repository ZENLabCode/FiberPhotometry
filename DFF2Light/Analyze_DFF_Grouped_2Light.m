% Analyze_DFF_Grouped_2Light
%-------------------------------------------------------------------------
% Analyze mean DFF data of each animals and condition (control + exp)
% MB - Last update: 03/03/2021 - Matlab R2019b

clc; clear; close all;
addpath('C:\FiberPhotometry\functions')


%% PARAMETERS
%---------------
%READ FILE LIST
[~,~,DATA] = xlsread('C:\FiberPhotometry\filelist.xlsx');
%corresponding column labels (first row of DATA)
Columns = { ... {variable name (do not change), label in xls-file}
    ... PS: - all labels must be unique (case in-sensitive)
    ...     - only the ones needed (others will be ignored)
    ...     - order is not important
    'path'      ,'Path';...
    'group'     ,'Group';...
    'trial'     ,'Trial';...
    'condition' ,'Condition';...
    'channel'   , 'Channel';...
    };

%Files (data file to read, same in all data paths)
rFile = 'Ch?.mat'; %with extention!

%ANALYSIS LABELS/TITLES
Labels = {... {fieldname (do not change), title string} !!
    % choose only one "before" at the time !!
%     'dff_before15' ,'15min before light pulse';...
    'dff_before5' ,'5min before light pulse';...
    'dff_stim'   ,'ZT14 - 15 min light pulse';...
    'dff_after'  ,'15min after light pulse';...
    };
noLAB = size(Labels,1);

Profiles = {... {filedname mean, duration [min], fieldname data }, this order
    % choose every profile for each selected Labels
%     'profile1',  15, 'dff_before15';...
    'profile1',  5, 'dff_before5';...
    'profile2', 14, 'dff_stim';...
    'profile3', 15, 'dff_after';...
    };

%ANALYSIS OPTIONS
boutDur = 60; %[s] for bouts of 1 min = 60 // 5 min = 60*5
delta_cut = 5; %[min]

%STATISTIC
opt.stat = 'SEM'; %'STD' or 'SEM'

%SAVING OPTIONS
saveData = false; %false for testing
sPath = 'C:\FiberPhotometry\Analysis';
sFile = 'DFF_Grouped'; %save filename



%% MAIN SCRIPT
%---------------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)));

%PREPARE DATA
%string trim (just in case)
ind = cellfun(@ischar,DATA);
DATA(ind) = cellfun(@strtrim,DATA(ind),'UniformOutput',false);
Columns   = cellfun(@strtrim,Columns,'UniformOutput',false);
labels = DATA(1,:);
DATA(1,:) = [];
%get column number
for k = 1:size(Columns)
    [timing,label] = Columns{k,:};
    col = find(strcmpi(labels,label));
    n   = numel(col);
    switch n
        case 0
            error('Label ''%s'' not found in DATA',label)
        case 1
            columns.(timing) = col;
        otherwise
            error('Label ''%s'' must be unique (found N = %i)',label,n)
    end
end

%SELECT DATA
tmp = selectionList(DATA(:,columns.path));
ind = ismember(DATA(:,columns.path),tmp);
DATA(~ind,:) = [];

%conditions,groups,trials
conditions = natsort(unique(DATA(:,columns.condition)));
groups = natsort(unique(DATA(:,columns.group)));
trials = natsort(unique(DATA(:,columns.trial)));

%number of ...
noDATA = numel(DATA(:,1));
noTRI = numel(trials);
noCOND = numel(conditions);
noGRO = numel(groups);

%PRINT OUT
n = numel(num2str(max([noDATA,noCOND])));

noTC = NaN(noGRO,noCOND);   %number of trials per condition/group
noMG = NaN(noGRO,noCOND);   %number of mouse per condition/group

for i = 1:noGRO
    fprintf('------------- %s -------------\n',strjoin(groups(i),' '))
    for q = 1:noCOND
        ind1 = strcmpi(DATA(:,columns.group),groups{i});
        ind2 = strcmpi(DATA(:,columns.condition),conditions{q});
        lastInd = ind1.*ind2;
        noTC(i,q) = sum(lastInd);
        
        lastInd = unique(lastInd);
        noMG(i,q) = sum(lastInd);
        fprintf('%s: n = %i with N = %i trials \n',strjoin(conditions(q),' '),noMG(i,q),noTC(i,q))
    end
end
fprintf('\n------------------------------------------\n')
fprintf('Data       (N = %*i)\n',n,noDATA)
fprintf('Animals    (N = %*i)\n',n,sum(sum(noMG)))
fprintf('Conditions (N = %*i): %s\n',n,noCOND,strjoin(conditions,' '))
fprintf('Groups     (N = %*i): %s\n',n,noGRO,strjoin(groups,' '))
fprintf('Trials     (N = %*i): %s\n',n,noTRI,strjoin(trials,' '))


legSTRg = cell(noCOND,1);   %for legend
for i = 1:noCOND
    tmp1 = sum(noMG(:,i));
    tmp2 = sum(noTC(:,i));
    legSTRg{i} = sprintf('%s: N = %i trials (n = %i mice)',conditions{i},tmp2,tmp1);
end

%Condition colors
colsCOND = [0 0 0;1 0 0];

%% READ DATA
fprintf('\nREAD DATA\n')

%CONDITIONS LOOP (control/experimental)
clear res
NoDat = NaN(noCOND,1);
for cond = 1:noCOND
    condition = conditions{cond};
    ind = strcmpi(DATA(:,columns.condition),condition);
    Data = DATA(ind,:);
    noDATA = size(Data,1);
    NoDat(cond) = noDATA;
    %init
    res(cond).mean = NaN(noDATA,noLAB-1);
    for d = 1:noDATA
        [rPath,~,~,~,channel] = Data{d,:};
        file = fullfile(rPath,[channel,'.mat']);
        data = load(file);
        
        % FILL
        for k = 1:noLAB
            label = Labels{k,1};
            dat = data.(label);
            %total mean
            res(cond).mean(d,k) = mean(dat);
            
            %mean per 1 min bins
            ind = find(strcmpi(Profiles(:,3),label));
            if ~isempty(ind)
                [lab2,cols] = Profiles{ind,1:2};
                rows = data.SampRate*boutDur;
                tmp  = reshape(dat(1:rows*cols),[rows,cols]);
                %init
                if ~isfield(res(cond),lab2) || isempty(res(cond).(lab2))
                    res(cond).(lab2) = NaN(noDATA,cols);
                end
                %append
                res(cond).(lab2)(d,:) = nanmean(tmp,1);
                res(cond).raw.(lab2)(d,:) = dat(1:rows*cols)';
            end
        end
    end
end %condition loop

%% MEAN DATA
clear mRes sRes;
for cond = 1:noCOND
    tmp = res(cond).mean;
    dRes(cond).mouse = tmp;
    mRes(cond).mean = mean(tmp,1);
    sRes(cond).mean = std(tmp,[],1);
    
    for k = 1:size(Profiles,1)
        label = Profiles{k,1};
        tmp = res(cond).(label);
        mRes(cond).(label) = mean(tmp,1);
        sRes(cond).(label) = std(tmp,[],1);
        dRes(cond).(label) = tmp;
        dRes(cond).raw.(label) = res(cond).raw.(label);
        
        meaCut_first(cond).(label) = mean(mean(tmp(:,1:delta_cut),2));
        stdCut_first(cond).(label) = std(mean(tmp(:,1:delta_cut),2));
        meaCut_last(cond).(label) = mean(mean(tmp(:,end-(delta_cut-1):end)));
        stdCut_last(cond).(label) = std(mean(tmp(:,end-(delta_cut-1):end),2));
        
        
        deltRes(cond).(label) = [meaCut_first(cond).(label),...
            meaCut_last(cond).(label)];
        stRes(cond).(label) = [stdCut_first(cond).(label),...
            stdCut_last(cond).(label)];
        
    end
end

%% FIGURE I
% Plot histogram DFF control vs experimental
fprintf('PLOT FIGURE DFF GROUP MEAN\n')
hf(1) = figure('WindowState','maximize');
hold on

x0 = 1:noLAB;
dx = 1/(noCOND+1);
hp = NaN(noCOND,1);
legSTR = cell(noCOND,1);

%CONDITION LOOP
for cond = 1:noCOND
    condition = conditions{cond};
    x = x0+(cond-(noCOND+1)/2)*dx;
    
    %DATA
    dat = dRes(cond).mouse;
    mea = mRes(cond).mean;
    n   = NoDat(cond);
    switch lower(opt.stat)
        case 'std'
            err = sRes(cond).mean;
        case 'sem'
            err = (sRes(cond).mean)./sqrt(n);
        otherwise
            error('opt.stat = ''%s'' is not implemented',opt.stat)
    end
    
    %Plot
    hp(cond) = bar(x,mea,'barwidth',0.8*dx,'facecolor',colsCOND(cond,:));
    hold on;
    errorbar(x,mea,NaN(size(err)),err ,'linestyle','none','color','k');
    
    for k = 1:numel(x)
        tmp = dat(:,k);
        h = plot(repmat(x(k),size(tmp)),tmp,'o',...
            'markersize',5,'linewidth',1,'color','c');
    end
    
    legSTR{cond} = sprintf('%s (n = %i)',condition,n);
    
end

maxY = max(max(dRes(cond).mouse));

%text
sgtitle((sprintf('%s vs %s',conditions{1},conditions{2})),'FontWeight','bold','FontSize',20)
legend(hp,legSTRg);
ylabel(sprintf('\\DeltaF/F Mean + %s',opt.stat))

%settings
set(gca,'xtick',x0,'xticklabel',Labels(:,2),'tickdir','both',...
    'ticklength',[0.005, 0.1],'ylim',[0,1.1*maxY],'FontSize',15)
vline(x0(1:2)+.5,':k')


%% FIGURE II
% Plot DFF profile
fprintf('PLOT FIGURE DFF PROFILE\n')
hf(2) = figure('WindowState','maximize');
hold on
noPRO = size(Profiles,1);
%init
x = 0;
mima = [inf,-inf];
xticklabel = [];
N = NaN(1,noPRO);
for pro = 1:noPRO
    [field,~,label] = Profiles{pro,:};
    n = numel(mRes.(field));
    xticklabel = [xticklabel,1:n];
    x = (1:n)+x(end);
    N(pro) = x(end);
    hp = NaN(n,1);
    legSTR = cell(noCOND,1);
    
    %CONDITION LOOP
    for cond = 1:noCOND
        condition = conditions{cond};
        
        %DATA
        mea = mRes(cond).(field);
        n   = NoDat(cond);
        switch lower(opt.stat)
            case 'std'
                err = sRes(cond).(field);
            case 'sem'
                err = (sRes(cond).(field))./sqrt(n);
            otherwise
                error('opt.stat = ''%s'' is not implemented',opt.stat)
        end
        
        %Plot
        hp(cond) = plot(x,mea,'-o','color',colsCOND(cond,:));
        curve1 = mea + err;
        curve2 = mea - err;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween,colsCOND(cond,:),'facealpha',.2,'EdgeColor','none');
        mima(1) = min([mima(1),min(curve2)]);
        mima(2) = max([mima(2),max(curve1)]);
        legSTR{cond} = sprintf('%s (n = %i)',condition,n);
    end
end
%text
sgtitle((sprintf('%s vs %s',conditions{1},conditions{2})),'FontWeight','bold','FontSize',20)
xlabel('Time [min]');
ylabel(sprintf('\\DeltaF/F Mean + %s',opt.stat))
legend(hp(1:noCOND),legSTRg)
%setting
yl = mima+[-1,1]*diff(mima)*0.05;
set(gca,'xtick',1:1:x(end),'xticklabel',xticklabel,'tickdir','both',...
    'ticklength',[0.005, 0.1],'xlim',[.5,x(end)+.5],...
    'ylim',yl,'FontSize',15)
N2 = [0,N];
for k = 1:noPRO
    ind = strcmpi(Labels(:,1),Profiles{k,3});
    str = Labels{ind,2};
    text((sum(N2(k:k+1))+1)/2,mima(1),str,'fontsize',12,...
        ...text((sum(N2(k:k+1))+1)/2,yl(1),{str,''},...
        'horizontalalignment','center',...
        'verticalalignment','middle')
    if k<noPRO
        vline(N(k)+0.5,':k')
    end
end


%% FIGURE III
% Plot DFF profile
fprintf('PLOT FIGURE DELTA DFF PROFILE\n')
hf(3) = figure('WindowState','maximize');
hold on

x0 = 1:noLAB;
dx = 1/(noCOND+1);
hp = NaN(noCOND,1);
legSTR = cell(noCOND,1);


for pro = 1:noPRO
    [field,~,label] = Profiles{pro,:};
    n = numel(deltRes.(field));
    x = [x0(pro)-dx/2 x0(pro)+dx/2];
    
    hp = NaN(n,1);
    legSTR = cell(noCOND,1);
    
    %CONDITION LOOP
    for cond = 1:noCOND
        condition = conditions{cond};
        
        %DATA
        mea = deltRes(cond).(field);
        n   = NoDat(cond);
        switch lower(opt.stat)
            case 'std'
                err = stRes(cond).(field);
            case 'sem'
                err = (stRes(cond).(field))./sqrt(n);
            otherwise
                error('opt.stat = ''%s'' is not implemented',opt.stat)
        end
        
        %Plot
        hp(cond) = plot(x,mea,'-o','color',colsCOND(cond,:));
        hold on;
        errorbar(x,mea,err,err ,'linestyle','none','color','k');
        
        mima(1) = min([mima(1),min(curve2)]);
        mima(2) = max([mima(2),max(curve1)]);
        legSTR{cond} = sprintf('%s (n = %i)',condition,n);
        
    end
end
%text
sgtitle((sprintf('%s vs %s',conditions{1},conditions{2})),'FontWeight','bold','FontSize',20)
ylabel(sprintf('%imin Delta \\DeltaF/F Mean + %s',delta_cut,opt.stat))
legend(hp(1:noCOND),legSTRg)

%settings
yl = mima+[-1,1]*diff(mima)*0.05;
set(gca,'xtick',x0,'xticklabel',Labels(:,2),'tickdir','both',...
    'ticklength',[0.005, 0.1],'xlim',[.5,x(end)+.5],...
    'ylim',yl,'FontSize',15)

vline(x0(1:2)+.5,':k')

%% FIGURE IV
% Plot DFF profile all
fprintf('PLOT FIGURE DFF PROFILE RAW ALL TRACES\n')
hf(4) = figure('WindowState','maximize');
hold on
noPRO = size(Profiles,1);
%init
x = 0;
mima = [inf,-inf];
xticklabel = [];
N = NaN(1,noPRO);
for pro = 1:noPRO
    [field,~,label] = Profiles{pro,:};
    n = size(dRes(cond).raw.(field),2);
    x = (1:n)+x(end);
    N(pro) = x(end);
    xticklabel = [xticklabel,1:n/(data.SampRate*boutDur)];
    hp = NaN(n,1);
    legSTR = cell(noCOND,1);
    
    %CONDITION LOOP
    for cond = 1:noCOND
        condition = conditions{cond};
        
        %DATA
        mea = nanmean(dRes(cond).raw.(field),1);
        n   = NoDat(cond);
        err = nanstd(dRes(cond).raw.(field),[],1);
        switch lower(opt.stat)
            case 'std'
                err = err;
            case 'sem'
                err = (err)./sqrt(n);
            otherwise
                error('opt.stat = ''%s'' is not implemented',opt.stat)
        end
        
        %Plot
        hp(cond) = plot(x,mea,'-','color',colsCOND(cond,:));
        curve1 = mea + err;
        curve2 = mea - err;
        x2 = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(x2, inBetween,colsCOND(cond,:),'facealpha',.2,'EdgeColor','none');
        
        
        %         for i = 8
        %             hp(cond) = plot(x,mea(i,:),'-','color',colsCOND(cond,:));
        %         end
        mima(1) = min([mima(1),min(min(mea-err))]);
        mima(2) = max([mima(2),max(max(mea+err))]);
        legSTR{cond} = sprintf('%s (n = %i)',condition,n);
    end
end
%text
sgtitle((sprintf('%s vs %s - Mice #%i',conditions{1},conditions{2},i)),'FontWeight','bold','FontSize',20)
xlabel('Time [min]');
ylabel(sprintf('\\DeltaF/F Mean + %s',opt.stat))
legend(hp(1:noCOND),legSTRg)
%setting
yl = mima+[-1,1]*diff(mima)*0.05;
set(gca,'xtick',data.SampRate*boutDur:data.SampRate*boutDur:x(end),'xticklabel',xticklabel,'tickdir','both',...
    'ticklength',[0.005, 0.1],'xlim',[.5,x(end)+.5],...
    'ylim',yl,'FontSize',15)
N2 = [0,N];
for k = 1:noPRO
    ind = strcmpi(Labels(:,1),Profiles{k,3});
    str = Labels{ind,2};
    text((sum(N2(k:k+1))+1)/2,mima(1),str,'fontsize',12,...
        ...text((sum(N2(k:k+1))+1)/2,yl(1),{str,''},...
        'horizontalalignment','center',...
        'verticalalignment','middle')
    if k<noPRO
        vline(N(k)+0.5,':k')
    end
end



%% FIGURE V
% Plot DFF profile all - single
fprintf('PLOT FIGURE DFF PROFILE RAW ALL TRACES\n')
hf(5) = figure('WindowState','maximize');
hold on
noPRO = size(Profiles,1);
%init
x = 0;
mima = [inf,-inf];
xticklabel = [];
N = NaN(1,noPRO);

for pro = 1:noPRO
    [field,~,label] = Profiles{pro,:};
    n = size(dRes(cond).raw.(field),2);
    x = (1:n)+x(end);
    N(pro) = x(end);
    xticklabel = [xticklabel,1:n/(data.SampRate*boutDur)];
    hp = NaN(n,1);
    legSTR = cell(noCOND,1);
    
    %CONDITION LOOP
    cond = 1;
    condition = conditions{cond};
    mea = dRes(cond).raw.(field);
    ha(cond) = subplot(noCOND,1,cond);
    
    %DATA
    n   = NoDat(cond);
    for i = 1:size(mea,1)
        %Plot
        hp = plot(x,mea(i,:),'-','color',colsCOND(cond,:));
        hold on
        mima(1) = min([mima(1),min(min(mea(i,:)))]);
        mima(2) = max([mima(2),max(max(mea(i,:)))]);
    end
end

%setting
legend(sprintf('%s (n = %i)',condition,n));
xlabel('Time [min]');
ylabel(sprintf('\\DeltaF/F Mean + %s',opt.stat));
yl = mima+[-1,1]*diff(mima)*0.05;
set(gca,'xtick',data.SampRate*boutDur:data.SampRate*boutDur:x(end),'xticklabel',xticklabel,'tickdir','both',...
    'ticklength',[0.005, 0.1],'xlim',[.5,x(end)+.5],...
    'ylim',[-0.2,1.2],'FontSize',15)
box('off')
N2 = [0,N];
for k = 1:noPRO
    ind = strcmpi(Labels(:,1),Profiles{k,3});
    str = Labels{ind,2};
    text((sum(N2(k:k+1))+1)/2,mima(1),str,'fontsize',12,...
        ...text((sum(N2(k:k+1))+1)/2,yl(1),{str,''},...
        'horizontalalignment','center',...
        'verticalalignment','middle')
    if k<noPRO
        vline(N(k)+0.5,':k')
    end
end

%init
x = 0;
mima = [inf,-inf];
xticklabel = [];
N = NaN(1,noPRO);
for pro = 1:noPRO
    [field,~,label] = Profiles{pro,:};
    n = size(dRes(cond).raw.(field),2);
    x = (1:n)+x(end);
    N(pro) = x(end);
    xticklabel = [xticklabel,1:n/(data.SampRate*boutDur)];
    hp = NaN(n,1);
    legSTR = cell(noCOND,1);
    
    %CONDITION LOOP
    cond = 2;
    condition = conditions{cond};
    mea = dRes(cond).raw.(field);
    ha(cond) = subplot(noCOND,1,cond);
    %DATA
    n   = NoDat(cond);
    %Plot
    for i = 1:size(mea,1)
        
        hp = plot(x,mea(i,:),'-','color',colsCOND(cond,:));
        hold on
        mima(1) = min([mima(1),min(min(mea(i,:)))]);
        mima(2) = max([mima(2),max(max(mea(i,:)))]);
    end
    
    
end

%setting
legend(sprintf('%s (n = %i)',condition,n));
xlabel('Time [min]');
ylabel(sprintf('\\DeltaF/F Mean + %s',opt.stat));
yl = mima+[-1,1]*diff(mima)*0.05;
set(gca,'xtick',data.SampRate*boutDur:data.SampRate*boutDur:x(end),'xticklabel',xticklabel,'tickdir','both',...
    'ticklength',[0.005, 0.1],'xlim',[.5,x(end)+.5],...
    'ylim',[-0.2,1.2],'FontSize',15)
box('off')
N2 = [0,N];
for k = 1:noPRO
    ind = strcmpi(Labels(:,1),Profiles{k,3});
    str = Labels{ind,2};
    text((sum(N2(k:k+1))+1)/2,mima(1),str,'fontsize',12,...
        ...text((sum(N2(k:k+1))+1)/2,yl(1),{str,''},...
        'horizontalalignment','center',...
        'verticalalignment','middle')
    if k<noPRO
        vline(N(k)+0.5,':k')
    end
end

%text
sgtitle((sprintf('%s vs %s',conditions{1},conditions{2})),'FontWeight','bold','FontSize',20)



%% SAVING
%---------------
%Figure
fprintf('SAVING\n')

if saveData
    fprintf(' Path : %s\n',sPath)
    %Figures
    for fig = 1:numel(hf)
        tmp = sprintf('%s_%i',sFile,fig);
        print(hf(fig),fullfile(sPath,tmp),'-dpng','-r300')
        print(hf(fig),fullfile(sPath,tmp),'-depsc','-r300')
        fprintf(' Fig%i : %s .png and .eps\n',fig,tmp)
    end
    
    %Data
    clear info; %just in case
    info.dimensions        = 'control vs experimental,before15-5-1min/during stim/after';
    info.control.n            = NoDat(1);
    info.experimental.n       = NoDat(2);
    
    info.control.mean         = mRes(1).mean;
    info.control.mstd         = sRes(1).mean;
    info.control.deltaCut     = sprintf('%i min cut',delta_cut);
    
    for pro = 1:noPRO
        field = Profiles{pro,:};
        info.control.(field).mean  = mRes(1).(field);
        info.control.(field).std   = sRes(1).(field);
        info.control.(field).delta = deltRes(1).(field);
        info.control.(field).dstd  = stRes(1).(field);
        
        info.experimental.(field).mean   = mRes(2).(field);
        info.experimental.(field).std    = sRes(2).(field);
        info.experimental.(field).delta  = deltRes(2).(field);
        info.experimental.(field).dstd   = stRes(2).(field);
        
    end
    info.experimental.mean    = mRes(2).mean;
    info.experimental.std     = sRes(2).mean;
    
    info.trials               = trials(:);
    info.timing               = Labels(:,2);
    save(fullfile(sPath,[sFile,'.mat']),'info')
    fprintf(' Data : %s.mat \n\n',sFile)
else
    fprintf(2,' Nothing Saved!\n')
end
