function [TM, MC, T] = check_srate(TM,MC,trialTypes,savepath)

% checks the sampling rate of the data files
% why? saving rates don't always match stated expected rates
% sampling rate of 1000 Hz may really be 999 Hz
%
% observation: looks like there are ~240 samples per second but the
% intervals between them are not as consistent as I thought they'd be.
% Implies that we should always interp your data to have the exact interval
% that matches the sampling rate?
%
% INPUT VARIABLES:
%   TM              Treadmill data structure from import_dflow
%   MC              Mocap, marker data structure from import_dflow
%   trialsTypes     String array of trialTypes in TM
%   savepath        Folder of where to save table and figures. if it
%                   does not exist or is empty, nothing is saved.
%
% OUTPUT VARIABLES:
%   TM              TM data structure with speed additions to checks field
%   MC              MC data structure with speed additions to checks field
%   T               Table of the srate mean+/-std for TM and MC
%
% written by Helen J. Huang (July 2023)
% updates tracked in git

hf=300;
figure(hf)
ntrials = length(trialTypes);
ct=1;
for t = trialTypes

    strparts = split(TM.(t).file,"/");
    filename=strparts(end);

    strparts = split(filename,"_");
    subj=strparts(1);

    TM.(t).checks.interval(1) = mean(diff(TM.(t).data.Time));
    TM.(t).checks.interval(2) = std(diff(TM.(t).data.Time));
    TM.(t).checks.srate = 1/TM.(t).checks.interval(1);

    MC.(t).checks.interval(1) = mean(diff(MC.(t).data.TimeStamp));
    MC.(t).checks.interval(2) = std(diff(MC.(t).data.TimeStamp));
    MC.(t).checks.srate = 1/MC.(t).checks.interval(1);

    % there are not the same... 
%     mean(1./diff(TM.(t).data(:,1)))
%     1/mean(diff(TM.(t).data(:,1)))

    % these times should be nearly the same.
    TM.(t).checks.trialtime = diff(TM.(t).data.Time([1 end]));
    MC.(t).checks.trialtime = diff(MC.(t).data.TimeStamp([1 end]));

    %% PLOTS
    %
    % first subplot shows the time vectors for the marker dflow file and
    % the treadmill dflow file. This will show that both files have the same
    % start and end times, but the treadmill file has more points.
    % This can be seen at 2 lines with different slopes.
    % Other 2 subplots will show the "real-time" sampling rate value.

    subplot(ntrials,3,(ct-1)*3+1);
    plot(TM.(t).data.Time,"m");
    hold on
    plot(MC.(t).data.TimeStamp,"b");
    legend("TM time", "MC time","location","southeast");
    title(subj+" "+t)
    ylabel("time")
    xlabel("frames")
    %%
    % A plot of the difference in the TimeStamp data in the mocap dflowfile

    subplot(ntrials,3,(ct-1)*3+2);
    plot(diff(TM.(t).data.Time),"m", "Marker",".","LineStyle","none");
    hold on
    plot(xlim(gca),TM.(t).checks.interval(1)*[1 1],'k','LineWidth',1);
    text(0.1*max(xlim(gca)),0.85*max(ylim(gca)),num2str(TM.(t).checks.interval(1))+ " s");
    title("diff in time of TM");
    ylabel("interval (s)")
    xlabel("frames")

    subplot(ntrials,3,(ct-1)*3+3);
    plot(diff(MC.(t).data.TimeStamp),"b","Marker",".","LineStyle","none");
    hold on
    plot(xlim(gca),MC.(t).checks.interval(1)*[1 1],'k','LineWidth',1);
    text(0.1*max(xlim(gca)),0.85*max(ylim(gca)),num2str(MC.(t).checks.interval(1))+ " s");
    title("diff in time of MC");
    ylabel("interval (s)")
    xlabel("frames")

    srateinfo(ct,1) = TM.(t).checks.interval(1);
    srateinfo(ct,2) = TM.(t).checks.interval(2);
    srateinfo(ct,3) = TM.(t).checks.srate;

    srateinfo(ct,4) = MC.(t).checks.interval(1);
    srateinfo(ct,5) = MC.(t).checks.interval(2);
    srateinfo(ct,6) = MC.(t).checks.srate;

    srateinfo(ct,7) = TM.(t).checks.trialtime;
    srateinfo(ct,8) = MC.(t).checks.trialtime;

    ct=ct+1;
end

% create table of data for easier displaying and saving
T = array2table(srateinfo,'VariableNames',...
    {'meanTMinterval','stdTMinterval','meanTMsrate',...
    'meanMCinterval','stdMCinterval','meanMCsrate',...
    'TMtrialtime','MCtrialtime'},'RowNames',trialTypes);

disp(T)

% save plots and table
if ~exist("savepath","var"), savepath=[]; end
if ~isempty(savepath)
    cd(savepath)
    print(hf,savepath+subj+"_check_srate",'-dpdf','-fillpage')
    saveas(hf, savepath+subj+"_check_srate.fig")
    writetable(T,savepath+subj+"_check_srate.xlsx")
end
