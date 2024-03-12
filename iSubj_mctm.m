function iSubj_mctm(subj)
% UCF BRaIN Lab starter code for treadmill projects
%
% This script, iSubj_mctm, processes *motion capture and treadmill data*
% for an individual subject. EMG and EEG data have separate sets of scripts
%
% Group analyses is performed in another script
%
% Data types (files)
%   Motion capture (.txt)
%   Treadmill, D-Flow (.txt)
%
% Steps for individual subject analyses
%
%  1 - Experiment info - project, subj codes, conditions, etc.
%  2 - Set up paths to raw data, folders for saving results, code
%  3 - Import data
%  4 - Sync and resample data
%  5 - Rotate lab axes to more traditional biomech axes
%  6 - Get gait events, perturbations events if applicable
%  7 - Check gait events
%  8 - Calculate any metrics where we want timeseries plots
%  9 - Split timeseries into strides
% 10 - Calculate and plot gait cycle mean+/-sd
% 11 - Calculate metrics and plot by steps or strides
%       a. step kinematics (step length, step width, step freq, step time,
%       double support time, single support, stance time, swing time, etc.)
%       b. stride kinematics
%       c. walking speed if self-paced walking
%       d. margin of stability, extrapolated center of mass, base of
%       support
% 11 - Save experiment info, stride timeseries data, and metrics for
%      individual subject
%
% Folder structure
%   Proj -> Rawdata -> MocapTreadmill -> Subjcode1 -> all the .txt files for all trials
%                   -> MocapTreadmill -> Subjcode2 -> all the .txt files for all trials
%                   -> MocapTreadmill -> Subjcode3 -> all the .txt files for all trials
%                           ...
%                   -> EMG -> SubjcodeX -> all the .txt files for all trials
%                   -> EEG -> SubjcodeX -> all the .txt files for all trials
%   Proj -> Results -> MocapTreadmill -> bySubj -> Subjcode1_results.mat (will include all trials)
%                   -> MocapTreadmill -> bySubj -> Subjcode2_results.mat
%                   -> MocapTreadmill -> bySubj -> Subjcode3_results.mat
%                           ...
%                       -> Group -> Group_results.mat
%   Proj -> Checks -> MocapTreadmill -> Subjcode1 -> any files, figures checking gait events or other
%   things
%                  -> MocapTreadmill -> Subjcode1 -> Subjcode1_checks.mat, Subjcode1_checks.fig
%                  -> MocapTreadmill -> Subjcode2 -> Subjcode2_checks.mat, Subjcode2_checks.fig
%                           ...
% Written by Helen J. Huang (2023.07)
%
% for version control - view git 
%
%% 0 - intialize

clearvars -except subj;
close all; clc;

%% 1 - Experiment info - project, subj codes, trials, etc.

projCode = "FE";
% if subj or var are not defined use SPP7 otherwise will used defined subj
% set to 1 to run checks
run_checks = 1;
% set to 1 to save checks
save_checks = 0;

if ~exist('subj','var') || isempty(subj), subj = "FE1"; else, subj = string(subj); end
disp("subj: "+subj);

% trials, aka conditions
trialTypes = ["bpm60"];

% adaptation studies
% blocks = ["pre" "pert" "post"];
% phases = ["early" "late"]; % not all will have this

% sampling rates of systems
% MC = mocap; TM = treadmill
MC.srate = 240; %hz, supposedly
TM.srate = 300; %hz, supposedly
des_srate = 240; %hz, for interpolation

%% 2 - Set up paths to raw data, folders for saving results, code.
% Then add to paths folders
% tutorial data is SPP

% on windows fs = \\ but on mac, fs = //
fs = string(filesep) + string(filesep);

fPath = string(pwd)+fs; % assumes script is ran from its folder

% path to project files and code in git
if ispc % for windows
    p2l.root = "D:\DATA\";
    p2l.git = "C:\Users\hjh45\Documents\Github\ucfbrainlab_mctm\";
elseif isunix % for mac and linux
    p2l.root = "/Users/lindseylee/Desktop/brainlab/"; %for spp data 
   % p2l.root = "/Volumes/Lexar/test_incline/"; %for incline testing
    p2l.git = "/Users/lindseylee/Documents/GitHub/FE";
end

p2l.proj = p2l.root + projCode + fs;
p2l.rawData = p2l.proj + "Rawdata" + fs + "MocapTreadmill" + fs + subj + fs;
p2l.results = p2l.proj + "Results" + fs + "MocapTreadmill" + fs + "bySubj" +fs + subj + fs;
p2l.checks = p2l.proj + "Checks" + fs + "MocapTreadmill" + fs + subj + fs;
p2l.events = p2l.proj + "Events" + fs + subj + fs; % Gait Events
% makes folders if they dont exist 
if ~isfolder(p2l.results), mkdir(p2l.results); end
if ~isfolder(p2l.checks), mkdir(p2l.checks); end
if ~isfolder(p2l.events), mkdir(p2l.events); end

addpath(genpath(fPath))

% create file paths to each trial's mocap & treadmill dflow files
for t = trialTypes
    f2l.mc.(t) = p2l.rawData + subj + "_" + t + "_0001.txt";
    f2l.tm.(t) = p2l.rawData + subj + "_" + t + "_treadmill0001.txt";
end

%%  3 - Import data
% loads mc and tm files & creates data structures that holds data and 
% varnames for each trial type
for t = trialTypes

    % so we'll know which file the data is from
    MC.(t).file = f2l.mc.(t);
    TM.(t).file = f2l.tm.(t);

    [MC.(t).data, MC.(t).varnames] = import_dflow(MC.(t).file);
    [TM.(t).data, TM.(t).varnames] = import_dflow(TM.(t).file);

end

%% Run treadmill file checks

if save_checks, savepath = p2l.checks; else savepath = []; end

if run_checks

    % Check treadmill filenames and conditions based on sway
    [TM, CHK.sway] = check_filenames_bysway(TM,trialTypes,savepath);

    % Check whether speed was clipped
    [TM, CHK.speed] = check_speed(TM,trialTypes,savepath);

    % Check sampling rates of marker and treadmill files
    [TM, MC, CHK.srate] = check_srate(TM,MC,trialTypes,savepath);

    % Check perturbation frequency onset
    [TM, CHK.freq] = check_freq(TM,trialTypes,savepath);
end

%%  4 - Sync and resample data
%
% Dflow sampling rate is 240 Hz but the treadmill file saves at an
% inconsistent rate, somewhat faster rate (265-299 hz).
% So we will interpolate both the MC and TM data to be at 240 Hz.

% use interp1 because the intervals between data points is not consistent,
% equidistant which is a requirement for using interpft. Here we really do
% need to interpolate points between the existing points, not just resample
% as interpft basically does.

for t = trialTypes
    tct = find(trialTypes == t);

    % take average of the 2 files, which should be very close
    trialtime = mean([CHK.srate.TMtrialtime(tct) CHK.srate.MCtrialtime(tct)]);

    timeq = linspace(0,trialtime,trialtime*des_srate);

    figure
    % time is in the first column of both files
    % sets first point of trial time to 0
    time0 = MC.(t).data.TimeStamp-MC.(t).data.TimeStamp(1);

    % interpolate mocap data
    MC.(t).idata = array2table(interp1(time0, MC.(t).data{:,:},timeq));
    MC.(t).idata.Properties.VariableNames = MC.(t).data.Properties.VariableNames;
    MC.(t).idata(end,:)=[]; % delete because last row is nan
    MC.(t).itime = timeq(1:end-1)'; % remove last value to keep same length as idata
    MC.(t).israte = 240;

    % column 16 is arbitrary here
    subplot(1,2,1)
    plot(time0,MC.(t).data.FP1_ForY, 'b', 'LineWidth',2)
    hold on
    plot(MC.(t).itime, MC.(t).idata.FP1_ForY,'r:','LineWidth',2);
    legend("original","interpolated")

    % treadmill file
    time0 = TM.(t).data.Time-TM.(t).data.Time(1);

    % interpolate treadmill data
    TM.(t).idata = array2table(interp1(time0,TM.(t).data{:,:},timeq));
    TM.(t).idata.Properties.VariableNames = TM.(t).data.Properties.VariableNames;
    TM.(t).idata(end,:)=[];  % delete because last row is nan
    TM.(t).itime = timeq(1:end-1)'; % remove last value to keep same length as idata
    TM.(t).israte = 240;

    % column 16 here is speed actual 1
    subplot(1,2,2)
    plot(time0,TM.(t).data.SpeedActual1, 'b', 'LineWidth',2)
    hold on
    plot(TM.(t).itime, TM.(t).idata.SpeedActual1,'r:','LineWidth',2);
    legend("original","interpolated")
    
end

%% 5 - Rotate lab axes to more traditional biomech axes
% rotate lab axes to more intuitive biomech axes
% lab axes, +x right, +y up, +z posterior
% biomech axes, +x right, +y anterior, +z up

% use interpolated data

for t = trialTypes

    % ** MARKERS
    % just want marker names.
    markers_x = MC.(t).varnames(find(contains(MC.(t).varnames,"PosX")));

    for m = markers_x
        idx = find(contains(MC.(t).varnames, m)); % find marker's PosX index

        temp = split(m,'_');
        marker = temp(1);

        % x,y,z is idx, idx+1, idx+2
        % use interpolated data
        % **** Rotate ****
        MC.(t).markers.(marker) = array2table(rotation_ucf2biomech(MC.(t).idata{:,idx:idx+2}));

        % handle dropped markers
        bad = find(abs(MC.(t).markers.(marker){:,2})>10)'; % use y
        disp(subj+ ", " +t+ ", "+marker+" dropped marker: "+num2str(bad));
        for i = 1:length(bad)
            while abs(MC.(t).markers.(marker){bad(i),:})>10
                MC.(t).markers.(marker){bad(i),:}=MC.(t).markers.(marker){bad(i)-1,:};
            end
        end
    end

    % ** FORCES
    for met = ["Cop" "For" "Mom"]
%         MC.(t).varnames(contains(MC.(t).varnames,met+"X"))
        % idx(1) for FP1; idx(2) for FP2
        idx = find(contains(MC.(t).varnames,met+"X"));

        % FP1 = left; FP2 = right;
        % use interpolated data
        % **** Rotate ****
        MC.(t).(met).left = array2table(rotation_ucf2biomech(MC.(t).idata{:,idx(1):idx(1)+2}));
        MC.(t).(met).left.Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];

        MC.(t).(met).right = array2table(rotation_ucf2biomech(MC.(t).idata{:,idx(2):idx(2)+2}));
        MC.(t).(met).right.Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];
    end
end

%% 6 - Get gait events, perturbations events if applicable

for t = trialTypes
%     idx = find(TM.(t).varnames == "stride");
%     window(1) = find(TM.(t).idata(:,idx)==1,1,'first');
% 
%     if t == "nopert"
%         maxstride = max(TM.(t).idata(:,idx));
%         window(2) = find(TM.(t).idata(:,idx)==maxstride,1,'last');
%     else
%         % 560 because there are 560 strides in the protocol
%         window(2) = find(TM.(t).idata(:,idx)==560,1,'last');
%     end

%chopped out first 10 sec (ramp up speed) and end at 420 seconds
    window(1) = find(TM.(t).idata.t < 10, 1, 'last') + 1;
    window(2) = find(TM.bpm60.idata.t < 420, 1, 'last') +1 ;

    % ****
    lpfc = 6; % lowpass cutoff frequency

    % ****
    disp("******* " + t + " *******")
    GE.(t) = get_gaitevents(MC.(t).markers,lpfc,des_srate,window);

end

%%  7 - Check gait events
for t = trialTypes

    % ~~ FILTER GRF DATA ~~
    lpfc = 6;
    [b, a] = butter(4,2.0*1/MC.srate*lpfc);
    MC.(t).For.filt.left = array2table(filtfilt(b,a,MC.(t).For.left{:,:}));
    MC.(t).For.filt.right = array2table(filtfilt(b,a,MC.(t).For.right{:,:}));
    MC.(t).For.filt.left.Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];
    MC.(t).For.filt.right.Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];

    % ~~
    % after rotation (above) col 3 = z, up
    % after interpolations
    vGRF_filt = [MC.(t).For.filt.left{:,3}, MC.(t).For.filt.right{:,3}];

    % remember to pass interpolated data. vGRF should be interpolated above
    % only plot goodstrides
    check_gaitevents(vGRF_filt,GE.(t).events(GE.(t).goodstrides,:),0,subj+", "+t)

end

%%  8 - Calculate any metrics where we want timeseries plots
% self-paced walking speed, COM, COMvel, xCOM
%
for t = trialTypes

    % ~~ FILTER DATA ~~
    lpfc = 6;
    [b, a] = butter(4,2.0*1/MC.srate*lpfc);
    TM.(t).fidata = array2table(filtfilt(b,a,TM.(t).idata{:,["SpeedActual1" "SpeedActual2" "SwayActual" "COMz"]}));
    TM.(t).fidata.Properties.VariableNames = ["SpeedActual1" "SpeedActual2" "SwayActual" "COMz"];

    for m = ["RASI" "LASI" "RPSI" "LPSI"]
        MC.(t).markers.filt.(m) = array2table(filtfilt(b,a,MC.(t).markers.(m){:,:}));
        MC.(t).markers.filt.(m).Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];
    end

    % ~~ END FILTER DATA ~~
    MC.(t).COM.pos = (MC.(t).markers.RASI{:,:} + MC.(t).markers.LASI{:,:} + ...
        MC.(t).markers.RPSI{:,:} + MC.(t).markers.LPSI{:,:})/4;

    % figure
    % plot(TM.(t).itime,-TM.(t).idata.COMz,'o');
    % hold on
    % plot(MC.(t).itime,MC.(t).COM(:,2),'+'); % prefer to use MC based one.
    MC.(t).COM.dot = diff(MC.(t).COM.pos)./diff(MC.(t).itime);
    MC.(t).COM.dot = [MC.(t).COM.dot;nan nan nan]; % add nan to keep same length as position
    MC.(t).wspd = mean([TM.(t).idata.SpeedActual1, TM.(t).idata.SpeedActual2],2) + MC.(t).COM.dot(:,2);

end
%%  9 - Split timeseries into strides
%
for t = trialTypes

    %     variables = TM.(t).data(:,["stride" "SpeedActual1" "SpeedActual2" "SwayActual"]);
    %     strides.TM.(t) = split_strides_normgc(variables, GE.(t).events, 200, 1);
    strides.MC.(t).COM.pos = split_strides_normgc(array2table(MC.(t).COM.pos), GE.(t).events, 200, 1);
    strides.MC.(t).COM.dot = split_strides_normgc(array2table(MC.(t).COM.dot), GE.(t).events, 200, 1);
   
    strides.MC.(t).wspd = split_strides_normgc(array2table(MC.(t).wspd), GE.(t).events, 200, 1);
   
    strides.MC.(t).For.right = split_strides_normgc(MC.(t).For.right, GE.(t).events, 200, 1);
    strides.MC.(t).For.left = split_strides_normgc(MC.(t).For.left, GE.(t).events, 200, 1);
   
   % variables = TM.(t).idata(:,["SwayActual" "SpeedActual1" "SpeedActual2" "stride" "COMz"]);
    variables = TM.(t).idata(:,["SwayActual" "SpeedActual1" "SpeedActual2" "COMz"]);
    strides.TM.(t) = split_strides_normgc(variables, GE.(t).events, 200, 1);
   
    pause
    close all
end

%% 10 - Calculate and plot gait cycle mean+/-sd
%
% clear plotinfo
% plotinfo.varname{1} = "X_ML";
% plotinfo.xlabel{1} = "gait cycle (%)";
% plotinfo.ylabel{1} = "left GRF (N)";
% plotinfo.color{1} = [1 0 0]; % red
%
% plotinfo.varname{2} = "Y_AP";
% plotinfo.xlabel{2} = "gait cycle (%)";
% plotinfo.ylabel{2} = "left GRF (N)";
% plotinfo.color{2} = [1 0 0]; % red
%
% plotinfo.varname{3} = "Z_UP";
% plotinfo.xlabel{3} = "gait cycle (%)";
% plotinfo.ylabel{3} = "left GRF (N)";
% plotinfo.color{3} = [1 0 0]; % red
% for t = trialTypes
%     GRP.(t) = plot_grpave_waveforms(strides.MC.(t).normgc, strides.MC.(t).For.left.normdata, plotinfo);
% end
figure
color2use = [1 0 0; 0 1 0; 0 0 1];

ct = 1; ntr = length(trialTypes);

for t = trialTypes
  
    % no perturbation stride
%     strides2plot = find(fix(min(strides.TM.(t).normdata.SwayActual)*100)~=0);
    %     if t == "nopert"
    %         varmean = mean(strides.MC.(t).For.left.normdata.(var2plot),2);
    %         varstd = std(strides.MC.(t).For.left.normdata.(var2plot),0,2);
    %     else
    %         varmean = mean(strides.MC.(t).For.left.normdata.(var2plot)(:,strides2plot),2);
    %         varstd = std(strides.MC.(t).For.left.normdata.(var2plot)(:,strides2plot),0,2);
    %     end
  
    wspddata = table2array(strides.MC.(t).wspd.normdata);
    strides.(t).pre = find(MC.(t).itime(GE.(t).events.RHS)>10 & MC.(t).itime(GE.(t).events.RHS)<70);
    strides.(t).pert = find(MC.(t).itime(GE.(t).events.RHS)>70 & MC.(t).itime(GE.(t).events.RHS)<370);
    strides.(t).post = find(MC.(t).itime(GE.(t).events.RHS)>370 & MC.(t).itime(GE.(t).events.RHS)<430);

    if t == "nopert"

        WSPD.nopert.mean = mean(wspddata,2);
        WSPD.nopert.std = std(wspddata,0,2);
       
        subplot(1,ntr,ct)
        hold on
      
        varmean = WSPD.(t).mean(:,1); varstd = WSPD.(t).std(:,1);

%         patch(patch_x, [varmean+varstd; flipud(varmean-varstd)], 'k', 'EdgeColor', 'k', 'facealpha', 0.5, 'edgealpha', 0.5);
        plot(strides.MC.(t).wspd.normgc,varmean, '-','Color','k','LineWidth',2)
        ylim([0 2])
        title(t)

    else
        WSPD.(t).mean(:,1) = mean(wspddata(:,1:strides.pre(end)),2);
        WSPD.(t).std(:,1) = std(wspddata(:,1:strides.pre(end)),0,2);
        WSPD.(t).mean(:,2) = mean(wspddata(:,81:480),2);
        WSPD.(t).std(:,2) = std(wspddata(:,81:480),0,2);
        WSPD.(t).mean(:,3) = mean(wspddata(:,481:end),2);
        WSPD.(t).std(:,3) = std(wspddata(:,481:end),0,2);
      
        patch_x = [strides.MC.(t).wspd.normgc, fliplr(strides.MC.(t).wspd.normgc)]';
        subplot(1,ntr,ct)
        hold on
       
        varmean = WSPD.(t).mean(:,1); varstd = WSPD.(t).std(:,1);
%         patch(patch_x, [varmean+varstd; flipud(varmean-varstd)], color2use(1,:), 'EdgeColor', color2use(1,:), 'facealpha', 0.5, 'edgealpha', 0.5);
        plot(strides.MC.(t).wspd.normgc,varmean, '-','Color', color2use(1,:),'LineWidth',2)
        varmean = WSPD.(t).mean(:,2); varstd = WSPD.(t).std(:,2);
%         patch(patch_x, [varmean+varstd; flipud(varmean-varstd)], color2use(2,:), 'EdgeColor', color2use(2,:), 'facealpha', 0.5, 'edgealpha', 0.5);
        plot(strides.MC.(t).wspd.normgc,varmean, '-','Color', color2use(2,:),'LineWidth',2)
        varmean = WSPD.(t).mean(:,3); varstd = WSPD.(t).std(:,3);
%         patch(patch_x, [varmean+varstd; flipud(varmean-varstd)], color2use(3,:), 'EdgeColor', color2use(3,:), 'facealpha', 0.5, 'edgealpha', 0.5);
        plot(strides.MC.(t).wspd.normgc,varmean, '-','Color',color2use(3,:),'LineWidth',2)
        legend(["pre" "pert" "post"])
        ylim([0 2])
        title(t)
    end
    ct=ct+1;
end
%% 11 - Calculate metrics and plot by steps or strides
%       a. step/stride kinematics (step length, step width, step freq, step time,
%       double support time, single support, stance time, swing time, etc.)
%       b. walking speed if self-paced walking
%       c. margin of stability, extrapolated center of mass, base of
%       support (will do later)

for t = trialTypes
    
%     idx(1) = find(TM.(t).varnames == "SpeedActual1");
%     idx(2) = find(TM.(t).varnames == "SpeedActual2");
    % GP paper says that we filtered treadmill data at a low pass of 6 hz
%     lpfc = 6;
%     [b, a] = butter(4,2.0*1/MC.srate*lpfc);
%     T.beltspd = filtfilt(b,a,TM.(t).idata{:,["SpeedActual1" "SpeedActual2"]});

    T.beltspd = TM.(t).fidata{:,["SpeedActual1" "SpeedActual2"]};
    T.pitch = 0; % in degrees
    T.itime = TM.(t).itime;

%     figure
%     plot(TM.(t).itime,TM.(t).idata{:,"SpeedActual1"},'LineWidth',1);
%     hold on
%     plot(TM.(t).itime,T.beltspd(:,1),'r:','LineWidth',2);
%     title(t + ", SpeedActual1"); legend(["unfilt" "filt"]);
%     xlabel("time (s)"); ylabel("belt speed")

    % ~~ FILTER DATA ~~
    lpfc = 6;
    [b, a] = butter(4,2.0*1/MC.srate*lpfc);

    for m = ["RHEE" "LHEE"]
        MC.(t).markers.filt.(m) = array2table(filtfilt(b,a,MC.(t).markers.(m){:,:}));
        MC.(t).markers.filt.(m).Properties.VariableNames = ["X_ML" "Y_AP" "Z_UP"];
    end

    % ~~ END FILTER DATA ~~
    SK.(t) = calc_stepkinematics(MC.(t).markers,GE.(t).events(GE.(t).goodstrides,:),T,des_srate);

end