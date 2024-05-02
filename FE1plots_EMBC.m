% FE Plots for ASB

clearvars
clc
close all

% change to folder where .mat is
% /Users/hjhuang/Dropbox (Personal)/_DATA/FE/Results
cd('/Users/lindseylee/Desktop/brainlab/FE/Results/EMG/bySubj/test')
load test_isubj_EMGdatastruct.mat

cd('/Users/lindseylee/Desktop/brainlab/FE/Results/MocapTreadmill/bySubj/test')
load test_isubj_metrics.mat

projCode = "FE";
subj = "test";

trialTypes = ["nopert" "bpm15" "bpm30" "bpm45" "bpm60" "bpm75"];
% trialTypes = ["nopert" "bpm60" "bpm45" "bpm15" "bpm75" "bpm30"];
%order: 60, 45, 15, 75, 30

phases = ["pre" "pert" "post"];

saveplot = 1;

%%

clear GRP

muscles = string(fieldnames(strides.EMG.nopert))';
forces = ["X_ML" "Y_AP" "Z_UP"];

for t = trialTypes
    mct = 1; clear meanemg
    for m = muscles
        meanemg(:,mct) = mean(strides.EMG.(t).(m).normdata{:,:}(:,strides.(t).pert),2);
        mct=mct+1;
    end
    GRP.ave.(t).emg = array2table(meanemg);
    GRP.ave.(t).emg.Properties.VariableNames = muscles;

    fct = 1; clear meanforce
    for f = forces
        meanforce_R(:,fct) = mean(strides.MC.(t).For.right.normdata.(f)(:,strides.(t).pert),2);
        meanforce_L(:,fct) = mean(strides.MC.(t).For.left.normdata.(f)(:,strides.(t).pert),2);
        fct=fct+1;
    end
    GRP.ave.(t).grf.right = array2table(meanforce_R);
    GRP.ave.(t).grf.right.Properties.VariableNames = forces;

    GRP.ave.(t).grf.left = array2table(meanforce_L);
    GRP.ave.(t).grf.left.Properties.VariableNames = forces;
end

%%
figure

c = winter(5);
mycolors(1,:) = [1 0 0];
mycolors(2:6,:) = c;

% mycolors(1,:) = [255 0 0]/255;        % red
% mycolors(2,:) = [0 191 255]/255;    % deep sky blue
% mycolors(3,:) = [70 130 180]/255;   % steelblue
% mycolors(4,:) = [65 105 225]/255;   % royal blue
% mycolors(5,:) = [0 0 205]/255;      % medium blue
% mycolors(6,:) = [25 25 112]/255;    % midnight blue

emgtime = strides.EMG.(t).BICEPSFEMLTuV.normgc; % same for all EMG

mct=1;
for m = muscles([1 8 2 9 3 10 4 11 5 12 6 13 7 14])
    tct=1; clear data2plot
    for t = trialTypes
        data2plot(:,tct) = GRP.ave.(t).emg.(m);
        tct=tct+1;
    end
    bardata2plot(mct,:) = max(data2plot);
    subplot(7,2,mct)
    plot(emgtime,data2plot)
    ax = gca;
    ax.ColorOrder = mycolors;

    ylim([0 150])
    title(m)

    if mct==1, legend(trialTypes); end
    if mct >= 13, xlabel("gait cycle %"); end
    if mod(mct,2), ylabel("microvolts"); end
    mct=mct+1;
end

%%

grftime = strides.MC.(t).For.right.normgc; % same for all forces

figure

fct = 1; 
for f = forces
    tct=1; clear data2plot
    for t = trialTypes
        data2plot(:,tct) = GRP.ave.(t).grf.left.(f);
        tct=tct+1;
    end
    bardata2plot(fct,:) = max(data2plot);
    subplot(1,3,fct)
    plot(grftime,data2plot)
    ax = gca;
    ax.ColorOrder = mycolors;

%     ylim([0 150])
    title(f)

%     if mct==1, legend(trialTypes); end
    xlabel("gait cycle %");
    if fct == 1, ylabel("grf (N)"); end
    if fct == 3, legend(trialTypes); end
    fct=fct+1;
end

%%
figure
bar(bardata2plot)
xticklabels(muscles([1 8 2 9 3 10 4 11 5 12 6 13 7 14]))
ax = gca;
ax.ColorOrder = mycolors;

%%
tct = 1;

for t = trialTypes

    if t ~= "nopert"
        data2plot.wspd.mean(tct,1) = mean(mean(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).pre)));
        data2plot.wspd.mean(tct,2) = mean(mean(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).pert)));
        data2plot.wspd.mean(tct,3) = mean(mean(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).post)));

        data2plot.wspd.std(tct,1) = mean(std(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).pre)));
        data2plot.wspd.std(tct,2) = mean(std(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).pert)));
        data2plot.wspd.std(tct,3) = mean(std(strides.MC.(t).wspd.normdata{:,:}(:,strides.(t).post)));
    elseif t == "nopert"
        data2plot.wspd.mean(tct,1) = mean(mean(strides.MC.(t).wspd.normdata{:,:}));
        data2plot.wspd.mean(tct,2) = data2plot.wspd.mean(tct,1);
        data2plot.wspd.mean(tct,3) = data2plot.wspd.mean(tct,1);

        data2plot.wspd.std(tct,1) = mean(std(strides.MC.(t).wspd.normdata{:,:}));
        data2plot.wspd.std(tct,2) = data2plot.wspd.std(tct,1);
        data2plot.wspd.std(tct,3) = data2plot.wspd.std(tct,1);
    end

    tct = tct+1;
end

for v = ["width" "length" "freq" "time"]
    tct = 1;
    for t = trialTypes

        if v == "freq", c = 60; else c = 1; end
        if t ~= "nopert"
            data2plot.(v).mean(tct,1) = mean(abs(mean(SK.(t).step.(v)(strides.(t).pre,:)*c)));
            data2plot.(v).mean(tct,2) = mean(abs(mean(SK.(t).step.(v)(strides.(t).pert,:)*c)));
            data2plot.(v).mean(tct,3) = mean(abs(mean(SK.(t).step.(v)(strides.(t).post,:)*c)));

            data2plot.(v).std(tct,1) = mean(std(SK.(t).step.(v)(strides.(t).pre,:)*c));
            data2plot.(v).std(tct,2) = mean(std(SK.(t).step.(v)(strides.(t).pert,:)*c));
            data2plot.(v).std(tct,3) = mean(std(SK.(t).step.(v)(strides.(t).post,:)*c));
        else
            data2plot.(v).mean(tct,1) = mean(abs(mean(SK.(t).step.(v)*c)));
            data2plot.(v).mean(tct,2) = data2plot.(v).mean(tct,1);
            data2plot.(v).mean(tct,3) = data2plot.(v).mean(tct,1);

            data2plot.(v).std(tct,1) = mean(std(SK.(t).step.(v)));
            data2plot.(v).std(tct,2) = data2plot.(v).std(tct,1);
            data2plot.(v).std(tct,3) = data2plot.(v).std(tct,1);
        end

        tct = tct+1;
    end
end

%%
figure('Name','step kinematics averages', "Position",[100 100 500 700])
ct = 1;
var2plot = ["wspd"]; % "freq" "length" "time" "width"];
for v = var2plot
    for p = 1:3 %phases, pre/pert/post
        sh(ct) = subplot(length(var2plot),3,ct);
        hold on
        bar(data2plot.(v).mean(2:6,p))
        plot(xlim(gca),[1 1]*data2plot.(v).mean(1,1),'LineStyle',':','LineWidth',2)
        xticks(1:5);
        xticklabels([]);
        if ct <= 3, subtitle(phases(p)); end
        sgtitle(subj+": step kinematics ave")

        switch v
            case "wspd"
                ylabel('walking speed (m/s)')
                xlabel('condition')
                ylim([1,1.8])
            case "freq"
                ylabel(["step " v + " (bpm)"])
                ylim([100 140])
            case "length"
                ylabel(["step " v + " (m)"])
                ylim([0.5 0.8])
            case "time"
                ylabel(["step " v + " (s)"])
                ylim([0.3 0.6])
            case "width"
                ylabel(["step " v + " (m)"])
                ylim([0 0.15])
        end

        if ct >= length(var2plot)*3-2
            xticklabels(trialTypes(2:end));
            xlabel('condition');
        end

        ct = ct+1;
    end

end

if saveplot
    saveas(gcf,'test_step_kinematics_ave.fig')
    saveas(gcf,'test_step_kinematics_ave.pdf')
end
%%
figure('Name','step kinematics variability', "Position",[100 100 500 700])
ct = 1;
var2plot = ["wspd" "freq" "length" "time" "width"];
for v = var2plot
    for p = 1:3 %phases, pre/pert/post
        subplot(length(var2plot),3,ct)
        hold on
        bar(data2plot.(v).std(2:6,p))
        plot(xlim(gca),[1 1]*data2plot.(v).std(1,1),'LineStyle',':','LineWidth',2)
        xticks(1:5);
        xticklabels([]);
        if ct <= 3, subtitle(phases(p)); end
        sgtitle(subj+":step kinematics variability")

        switch v
            case "wspd"
                if p == 1, ylabel('sd walking speed (m/s)'); end
                ylim([0,0.2])
            case "freq"
                if p == 1, ylabel(["sd step " v + " (bpm)"]); end
                ylim([0 15])
            case "length"
                if p == 1, ylabel(["sd step " v + " (m)"]); end
                ylim([0 0.15])
            case "time"
                if p == 1, ylabel(["sd step " v + " (s)"]); end
                ylim([0 0.1])
            case "width"
                if p == 1, ylabel(["sd step " v + " (m)"]); end
                ylim([0 0.07])
        end

        if ct >= length(var2plot)*3-2
            xticklabels(trialTypes(2:end));
            xlabel('condition'); 
        end

        ct = ct+1;
    end

end

if saveplot
    saveas(gcf,'test_step_kinematics_variability.fig')
    saveas(gcf,'test_step_kinematics_variability.pdf')
end

%%
figure('Name','metabolic power', "Position",[100 100 600 300])
subplot(121)
hold on
bar(ME.ave.nmp(3:end))
plot(xlim(gca),[1 1]*ME.ave.nmp(2),'LineStyle',':','LineWidth',2)
xticks(1:5);
xticklabels(trialTypes(2:end));
xlabel('condition')
ylabel('net metabolic power (W/kg)')

subplot(122)
hold on
bar(ME.ave.cot(3:end))
plot(xlim(gca),[1 1]*ME.ave.cot(2),'LineStyle',':','LineWidth',2)
xticks(1:5);
xticklabels(trialTypes(2:end));
xlabel('condition')
ylabel('cost of transport (J/kg*m)')

sgtitle(subj+": metabolic energy")

% saveas(gcf,'FE1_metabolics.fig')
% saveas(gcf,'FE1_metabolics.pdf')

%% Entrainment

% figure

tct = 1;
for t = trialTypes(2:6)

    TM.(t).t_start_idx = find(TM.(t).idata.t < 0.0001, 1, 'last') + 1;
    TM.(t).t_start = TM.(t).itime(TM.(t).t_start_idx);

    ns = length(strides.(t).pert);

    itime = MC.(t).itime;
    potime = TM.(t).checks.freq.onset_time + TM.(t).t_start;

    GEtimes.RHS = itime(GE.(t).events.RHS);
    GEtimes.RHSn = itime(GE.(t).events.RHSn);

%     subplot(1,5,tct)
    figure("position",[300, 300, 300, 150])
    hold on
    plot(strides.MC.(t).For.right.normgc, strides.MC.(t).For.right.normdata.Z_UP(:,strides.(t).pert))
    xlabel('gait cycle (%)')
    ylabel('Right vGRF (N)')
    title(t + " pert onset times")
    ylim([0 1800])

    for p = 1:length(potime)
        strnum = find(GEtimes.RHS < potime(p),1,'last');
        onset_after_RHS_s = potime(p)-GEtimes.RHS(strnum);
        stridetime_s = GEtimes.RHSn(strnum)-GEtimes.RHS(strnum);
        onset_after_RHS_pgc = onset_after_RHS_s/stridetime_s*100;

        pertinfo.(t).strnum(p) = strnum;
        pertinfo.(t).onset_after_RHS_s(p) = onset_after_RHS_s;
        pertinfo.(t).stridetime(p) = stridetime_s;
        pertinfo.(t).onset_after_RHS_pgc(p) = onset_after_RHS_pgc;

        plot([1 1]*onset_after_RHS_pgc,[0 100],'k-','LineWidth',0.5)
    end
    tct=tct+1;
end
%%
figure
tct = 1;

for t = trialTypes(2:6)

    subplot(1,5,tct)
    histogram(pertinfo.(t).onset_after_RHS_pgc,10,'normalization','probability')

    xlabel("bins of percent gait cycle")
    if t == "bpm15", ylabel('probability'); end
    subtitle(t)
    if t == "bpm75", sgtitle("histogram of perturbation onsets within gait cycle"); end
    tct = tct+1;

    ylim([0 0.35]);
end

%%
figure
tct = 1;

edges = [0 10 30 50 60 75 85 100];
for t = trialTypes(2:6)

    subplot(1,5,tct)
    histogram(pertinfo.(t).onset_after_RHS_pgc,edges,'normalization','probability')

%     xlabel("bins of percent gait cycle")
    xlabel("gait phases")
    xticks([5 20 40 55 68 80 93]);
    xticklabels(["LR", "MS", "TS", "TO", "ESw", "MSw", "LSw"]);
    if t == "bpm15", ylabel('probability'); end
    subtitle(t)
    if t == "bpm75", sgtitle("histogram of perturbation onsets within gait cycle"); end
    tct = tct+1;

    ylim([0 0.35]);
    box off
end