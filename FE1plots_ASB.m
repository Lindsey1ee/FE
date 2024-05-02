% FE Plots for ASB

clearvars
clc
close all

% change to folder where .mat is
% /Users/hjhuang/Dropbox (Personal)/_DATA/FE/Results
cd('/Users/lindseylee/Desktop/brainlab/FE/Results/MocapTreadmill/bySubj/test')
load test_isubj_metrics.mat
load test_isubj_datastruct.mat

%cd('/Users/hjhuang/Dropbox (Personal)/_DATA/FE//Results//VO2//bySubj//FE1//')
%load FE1_isubj_ME.mat


projCode = "FE";
subj = "test";

BW = 89.34; %kg

trialTypes = ["nopert" "bpm15" "bpm30" "bpm45" "bpm60" "bpm75"];
% trialTypes = ["nopert" "bpm60" "bpm45" "bpm15" "bpm75" "bpm30"];
%order: 60, 45, 15, 75, 30

phases = ["pre" "pert" "post"];

saveplot = 1;
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
var2plot = ["wspd" "freq" "length" "time" "width"];
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
                ylim([0.8,1.8])
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
                ylim([0 0.25])
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
% figure('Name','metabolic power', "Position",[100 100 600 300])
% subplot(121)
% hold on
% bar(ME.ave.nmp(3:end))
% plot(xlim(gca),[1 1]*ME.ave.nmp(2),'LineStyle',':','LineWidth',2)
% xticks(1:5);
% xticklabels(trialTypes(2:end));
% xlabel('condition')
% ylabel('net metabolic power (W/kg)')
% 
% subplot(122)
% hold on
% bar(ME.ave.cot(3:end))
% plot(xlim(gca),[1 1]*ME.ave.cot(2),'LineStyle',':','LineWidth',2)
% xticks(1:5);
% xticklabels(trialTypes(2:end));
% xlabel('condition')
% ylabel('cost of transport (J/kg*m)')
% 
% sgtitle(subj+": metabolic energy")
% 
% % saveas(gcf,'FE1_metabolics.fig')
% % saveas(gcf,'FE1_metabolics.pdf')

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