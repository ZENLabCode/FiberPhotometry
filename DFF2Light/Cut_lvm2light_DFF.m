% Cut_lvm2light_DFF
%-------------------------------------------------------------------------
% Define a refererence channel (refCHA) and a the signal drop (RefStep) of the light ON
% Calculate the DFF of the others channels
% Crop the DFF according to before, during and after the light
% MB - Last update: 03/03/2021

clc; clear; close all;
addpath('C:\Photometry\MATLAB\LabView\functions')
addpath('C:\Photometry\MATLAB\OptoLab_v4.1\function\misc')

%% PARAMETERS
%---------------
%READ FILES
rFiles = 'C:\Photometry\Brenna_shCDK5\*.lvm';

%SAVING OPTIONS
saving = false;

%OPTION FOR TIME
units = '[min]';
q = 60; %60 for minutes, 3600 for hours

%OPTION FOR REF CHANNEL
RefStep = 0.02;  % set the drop in the ref channel
refCHA = 1; %change according to the reference channel (1-2-3) 

%% MAIN SCRIPT
%---------------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FILE-LIST
rFiles = selectionList(findFiles(rFiles));

%FILE LOOP
noFIL = numel(rFiles);
if noFIL==0
    fprintf(2,'No File to convert!\n')
    return
end
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    close all;
    rFile = rFiles{fil};
    fprintf('\n%*i/%i: %s\n',nnFIL,fil,noFIL,rFile)
    [~,label] = fileparts(rFile);
    
    %% READ DATA
    data = lvm_import_zen(rFile);
    
    %% CHANNELS
    
    %Reference channel
    %demodulate signal
    signalR = data.Segment1.dataY(:,refCHA); %cropped, in [V]
    fsR = unique(data.Segment1.fs);
    SampRate = 10; %demod signal
    n1 = fsR/SampRate;
    n2 = floor(numel(signalR)/n1); %o samples
    signal = reshape(signalR(1:n1*n2),[n1,n2]);
    signal = bandpower(signal,fsR,480+[-10,10]);
    signal = signal(:);
    noSAM  = numel(signal);
    
    %Find stim start-end to use for cutting
    x = [...
        max([min(find(diff(signal)<-RefStep))+1,1]),...
        min([max(find(diff(signal)>RefStep))+1,length(signal)])];
    StimSTA = min(x);
    StimEND = max(x);
    
    %Check stim points
    plot(signal,'marker','.')
    set(gcf,'WindowState','maximize');
    hold on
    hp1 = plot(StimSTA,signal(StimSTA),'ro');
    hp2 = plot(StimEND,signal(StimEND),'ro');
    title('Checking Stimulation START and END (in red)');
    %Recheck/reset data points
    p = get(gca,'position');
    dx = 1-sum(p([1,3]));
    uicontrol(gcf,'style','pushbutton',...
        'unit','normalized',...
        'position',[sum(p([1,3]))+dx/10,sum(p([2,4]))-0.1,dx-2*dx/10,0.1],...
        'string','OK','FontWeight','bold','FontSize',12,...
        'BackgroundColor','g',...
        'callback','close(gcf);');
    uicontrol(gcf,'style','pushbutton',...
        'unit','normalized',...
        'position',[sum(p([1,3]))+dx/10,sum(p([2,4]))-0.21,dx-2*dx/10,0.1],...
        'string','Change Points','FontWeight','bold','FontSize',12,...
        'BackgroundColor','r',...
        'callback',['x=inputdlg({''Stim Start'',''Stim End''},''Change'',',...
        '[1 50; 1 50],{num2str(StimSTA),num2str(StimEND)});',...
        'x = sort(cellfun(@str2double,x));',...
        'StimSTA = min([ max( [1,x(1)] ),numel(signal)-1 ]);',...
        'StimEND = min([ max( [StimSTA+1,x(2)] ),numel(signal) ]);',...
        'set(hp1,''xdata'',StimSTA,''ydata'',signal(StimSTA));',...
        'set(hp2,''xdata'',StimEND,''ydata'',signal(StimEND));',...
        ]);
   
    zoom xon
    waitfor(gcf)
    StimDur = (StimEND - StimSTA)/q/SampRate;
    fprintf('\n%s Stim duration = %f %s\n',indent,StimDur,units)
    
    %Channels calcium
    channels = 1:numel(data.Segment1.column_labels);
    channels(channels==refCHA) = [];
    for channel = channels
        
        %demodulate signal
        signalR = data.Segment1.dataY(:,channel); %cropped, in [V]
        signal = reshape(signalR(1:n1*n2),[n1,n2]);
        signal = bandpower(signal,fsR,480+[-10,10]);
        signal = signal(:);
        
        %HUG-LINE
        %filter signal
        %polynomal filter across rougly 100-s moving windows, order 2
        signalF = sgolayfilt(signal,3,10*SampRate+1); %Savitzky-Golay filter
        %convex hull index (and removing from 1st increasing signalF(indH))
        indH = convhull((1:noSAM)',signalF,'simplify',true); % hug points
        tmp = find([diff(signalF(indH));inf]>0,1,'first'); %removed +1 !!!
        indH(tmp+1:numel(indH)) = [];
        %create hug-line
        if numel(indH)>1
            sigF = signalF(indH);
            slope = (sigF(end)-sigF(end-1))/(indH(end)-indH(end-1));
            indH(end+1) = noSAM+1; %add for end point
            sigF(end+1) = sigF(end)+(noSAM-indH(end))*slope; %check again
            signalH = NaN(size(signal));
            ind = 0; %init
            for k = 1:numel(indH)-1
                n = diff(indH(k:k+1));
                lin = linspace(sigF(k),sigF(k+1),n+1);
                ind = (1:n)+ind(end);
                signalH(ind) = lin(1:end-1);
                if ~isequal(size(signalH),size(signal))
                    disp(k)
                    return
                end
            end
        else
            signalH = NaN(size(signalF));
        end
        
        %detrended signal
        signalD = signalF-signalH;
        
        %DFF
        mi = prctile(signalD,1);
        ma = prctile(signalD,99);
        dff = (signalD-mi)/(ma-mi);
        
        %Crop
        dff_before15 = dff(StimSTA-(15*q*SampRate):StimSTA-1);
        dff_before1 = dff(StimSTA-(1*q*SampRate):StimSTA-1);
        dff_before5 = dff(StimSTA-(5*q*SampRate):StimSTA-1);
        dff_stim   = dff(StimSTA:StimEND-1);
        dff_after  = dff(StimEND:StimEND+(15*q*SampRate));
        
        if channel == channels(1);
            %Print lengths
            fprintf('%s dff_before15: %.2f %s\n',indent,length(dff_before15)/q/SampRate,units)
            fprintf('%s dff_before5:  %.2f %s\n',indent,length(dff_before1)/q/SampRate,units)
            fprintf('%s dff_before1:  %.2f %s\n',indent,length(dff_before5)/q/SampRate,units)
            fprintf('%s dff_stim:     %.2f %s\n',indent,length(dff_stim)/q/SampRate,units)
            fprintf('%s dff_after:    %.2f %s\n',indent,length(dff_after)/q/SampRate,units)
        end
        
        %% SAVE cropped DFF
        sname = fullfile(fileparts(rFile),sprintf('Ch%i.mat',channel));
        if saving
            fprintf('%s SAVING DATA: %s\n',indent,sname)
            save(sname,'StimSTA','StimEND','dff_before15','dff_before1','dff_before5','dff_stim','dff_after','SampRate')
        else
            fprintf(2,'%s%s NOT saved!\n',indent,sname)
        end
    end
end
fprintf(2,'ANALYSIS OF %i FILES DONE !\n',noFIL)
