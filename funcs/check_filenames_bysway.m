function [TM, T] = check_filenames_bysway(TM,trialTypes,savepath)

% Reads treadmill sway data (SwayTheory, SwayActual) to check if filenames
% match condition
%
% INPUT VARIABLES:
%   TM              Treadmill data structure from import_dflow
%   trialsTypes     String array of trialTypes in TM
%   savepath        Folder of where to save table and figures. if it
%                   does not exist or is empty, nothing is saved.
%
% OUTPUT VARIABLES:
%   TM              TM data structure with sway additions to checks field
%   T               Table of the sway mean+/-std for TM
%
% written by Helen J. Huang (July 2023)
% updates tracked in git

initpath = pwd;

hf = 100;
figure(hf)
ct = 1;
for t = trialTypes
    % negative sign to make shifts positive to use findpeaks function
    [pks.the,locs.the] = findpeaks(-TM.(t).data.SwayTheory);
    [pks.act,locs.act] = findpeaks(-TM.(t).data.SwayActual);
    
    % remove "peaks" near zero
    i = find(pks.the < 0.005);
    if ~isempty(i)
        pks.the(i) = [];
        locs.the(i) = [];
    end

    % expect to find 200 peaks for perturbation shifts
    % 0.005 won't reject small shifts for chmag or bothch 
    % for shifts between 1-3, 3-5, 1-5
    % as a result, there will be more than 200 peaks
    i = find(pks.act < 0.005);
    if ~isempty(i)
        pks.act(i) = [];
        locs.act(i) = [];
    end

    subplot(5,1,ct)
    plot(TM.(t).data.stride,TM.(t).data.SwayTheory,'b');
    hold on
    plot(TM.(t).data.stride,TM.(t).data.SwayActual,'r');

    plot(TM.(t).data.stride(locs.the),-pks.the,'b+')
    plot(TM.(t).data.stride(locs.act),-pks.act,'ro')
    xlim([0 650])
    legend(["theory" "actual"]);

    if ct == 5, xlabel("strides"); end
    ylabel("sway (m)")

    strparts = split(TM.(t).file,"/");
    filename=strparts(end);
    title(filename)
    strparts = split(filename,"_");
    subj=strparts(1);

    intervals.the = diff(TM.(t).data.Time(locs.the));
    intervals.act = diff(TM.(t).data.Time(locs.act));

    TM.(t).checks.swaypks.the = pks.the;
    TM.(t).checks.swaypks.act = pks.act;

    TM.(t).checks.swaylocs.the = locs.the;
    TM.(t).checks.swaylocs.act = locs.act;

    TM.(t).checks.intervals.the = intervals.the;
    TM.(t).checks.intervals.act = intervals.act;

    % create table of data for easier displaying and saving
    if ~isempty(pks)
        swayinfo(ct,1) = mean(pks.the);
        swayinfo(ct,2) = std(pks.the);
        swayinfo(ct,3) = mean(pks.act);
        swayinfo(ct,4) = std(pks.act);
        swayinfo(ct,5) = mean(intervals.the);
        swayinfo(ct,6) = std(intervals.the);
        swayinfo(ct,7) = mean(intervals.act);
        swayinfo(ct,8) = std(intervals.act);
    end

    ct = ct + 1;
end
% creates table of info in code
T = array2table(swayinfo,'VariableNames',...
    {'PksTheory_mean','PksTheory_std','PksActual_mean','PksActual_std',...
    'IntervalsTheory_mean','IntervalsTheory_std',...
    'IntervalsActual_mean','IntervalsActual_std'},...
    'RowNames',trialTypes);
disp(T)

% save plots and table
if ~exist("savepath","var"), savepath=[]; end
if ~isempty(savepath)
    cd(savepath)
    saveas(hf, subj+"_check_filenames_bysway.pdf")
    saveas(hf, subj+"_check_filenames_bysway.fig")
    writetable(T,savepath+subj+"_check_filenames_bysway.xlsx")
    cd(initpath)
end