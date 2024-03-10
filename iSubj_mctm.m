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

if ~exist('subj','var') || isempty(subj), subj = "FEp1"; else, subj = string(subj); end
disp("subj: "+subj);

% trials, aka conditions
trialTypes = ["bpm45"];

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
    time0 = MC.(t).data(:,1)-MC.(t).data(1,1);

    % interpolate mocap data
    MC.(t).idata = interp1(time0, MC.(t).data,timeq);
    MC.(t).idata(end,:)=[]; % delete because last row is nan
    MC.(t).itime = timeq(1:end-1)'; % remove last value to keep same length as idata
    MC.(t).israte = 240;

    % column 16 is arbitrary here
    subplot(1,2,1)
    plot(time0,MC.(t).data(:,16), 'b', 'LineWidth',2)
    hold on
    plot(MC.(t).itime, MC.(t).idata(:,16),'r:','LineWidth',2);
    legend("original","interpolated")

    % treadmill file
    time0 = TM.(t).data(:,1)-TM.(t).data(1,1);

    % interpolate treadmill data
    TM.(t).idata = interp1(time0,TM.(t).data,timeq);
    TM.(t).idata(end,:)=[];  % delete because last row is nan
    TM.(t).itime = timeq(1:end-1)'; % remove last value to keep same length as idata
    TM.(t).israte = 240;

    % column 16 here is speed actual 1
    subplot(1,2,2)
    plot(time0,TM.(t).data(:,16), 'b', 'LineWidth',2)
    hold on
    plot(TM.(t).itime, TM.(t).idata(:,16),'r:','LineWidth',2);
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
        MC.(t).markers.(marker) = rotation_ucf2biomech(MC.(t).idata(:,idx:idx+2));

        % handle dropped markers
        bad = find(abs(MC.(t).markers.(marker)(:,2))>10)'; % use y
        disp(subj+ ", " +t+ ", "+marker+" dropped marker: "+num2str(bad));
        for i = 1:length(bad)
            while abs(MC.(t).markers.(marker)(bad(i),:))>10
                MC.(t).markers.(marker)(bad(i),:)=MC.(t).markers.(marker)(bad(i)-1,:);
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
        MC.(t).(met).left = rotation_ucf2biomech(MC.(t).idata(:,idx(1):idx(1)+2));
        MC.(t).(met).right = rotation_ucf2biomech(MC.(t).idata(:,idx(2):idx(2)+2));
    end
end

%% 6 - Get gait events, perturbations events if applicable

for t = trialTypes
    idx = find(TM.(t).varnames == "stride");
    window(1) = find(TM.(t).idata(:,idx)==1,1,'first');

    if t == "nopert"
        maxstride = max(TM.(t).idata(:,idx));
        window(2) = find(TM.(t).idata(:,idx)==maxstride,1,'last');
    else
        % 560 because there are 560 strides in the protocol
        window(2) = find(TM.(t).idata(:,idx)==560,1,'last');
    end

    % ****
    lpfc = 6; % lowpass cutoff frequency
    % ****
    disp("******* " + t + " *******")
    GE.(t) = get_gaitevents(MC.(t).markers,lpfc,des_srate,window);

end

%%  7 - Check gait events

for t = trialTypes

    % after rotation (above) col 3 = z, up
    vGRF = [MC.(t).For.left(:,3), MC.(t).For.right(:,3)];

    % filter GRF
    lpfc = 6; 
    [b, a] = butter(4,2.0*1/MC.srate*lpfc);
    vGRF_filt = filtfilt(b,a,vGRF);

    % remember to pass interpolated data, idata
    % only plot goodstrides
    check_gaitevents(vGRF_filt,GE.(t).events(GE.(t).goodstrides,:),0,subj+", "+t)

end

%%  8 - Calculate any metrics where we want timeseries plots

% self-paced walking speed, COM, COMvel, xCOM
%

%%  9 - Split timeseries into strides
%

%% 10 - Calculate and plot gait cycle mean+/-sd
%

%% 11 - Calculate metrics and plot by steps or strides
%       a. step/stride kinematics (step length, step width, step freq, step time,
%       double support time, single support, stance time, swing time, etc.)
%       b. walking speed if self-paced walking
%       c. margin of stability, extrapolated center of mass, base of
%       support (will do later)

for t = trialTypes
    idx(1) = find(TM.(t).varnames == "SpeedActual1");
    idx(2) = find(TM.(t).varnames == "SpeedActual2");

    % GP paper says that we filtered treadmill data at a low pass of 6 hz
    lpfc = 6;
    [b, a] = butter(4,2.0*1/MC.srate*lpfc);
    T.beltspd = filtfilt(b,a,TM.(t).idata(:,[idx(1),idx(2)]));

    T.pitch = 0; % in degrees
    T.itime = TM.(t).itime;

    figure
    plot(TM.(t).itime,TM.(t).idata(:,idx(1)),'LineWidth',1);
    hold on
    plot(TM.(t).itime,T.beltspd(:,1),'r:','LineWidth',2);
    title(t + ", SpeedActual1"); legend(["unfilt" "filt"]); 
    xlabel("time (s)"); ylabel("belt speed")

    SK.(t) = calc_stepkinematics(MC.(t).markers,GE.(t).events(GE.(t).goodstrides,:),T,des_srate);
end