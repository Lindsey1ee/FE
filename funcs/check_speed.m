function [TM, T] = check_speed(TM,trialTypes,savepath)

% checks if belt speeds were clipped in SPP. plots the speeds for each
% trialType and saves and xlsx of the belt speeds table.
%
% INPUT VARIABLES:
%   TM              Treadmill data structure from import_dflow
%   trialsTypes     String array of trialTypes in TM
%   savepath        Folder of where to save table and figures. if it
%                   doesn't exist or is empty, nothing is saved.
%
% OUTPUT VARIABLES:
%   TM              TM data structure with speed additions to checks field
%   T               Table of the belt speeds
%
% written by Helen J. Huang (July 2023)
% updates tracked in git

hf = 200;
figure(hf)
tct = 0;
for t = trialTypes

    strparts = split(TM.(t).file,"/");
    filename=strparts(end);

    strparts = split(filename,"_");
    subj=strparts(1);

    tct = tct+1;
    subplot(2,3,tct)
    plot(TM.(t).data.t,[TM.(t).data.SpeedActual1, TM.(t).data.SpeedActual2]);
    hold on
    plot([1 TM.(t).data.t(end)],[1.5 1.5], 'k', 'LineWidth',2);
    ylim([0 2.1])
    xlabel("time"); ylabel("belt speed m/s");
    title(subj+ ": "+ t)

    % Perturbation conditions include pre and post.
    % 80 strides, ~15% of all data. 80/(80+540+80)
    % Take the average from 20% (80/(80+400+80)+5% until the 85%=100-15%
    % is a quick estimate that avoids the initial increase in speed at
    % the beginning, and the pre and post
    %
    % For nopert, there was no pre, post. Take average for 15000 to
    % end-2400, ~10 seconds. 10 s* 240 samples/s.
    % to avoid initial increase in speed.
%     if contains(t,"nopert")
%         TM.(t).checks.beltspeed(1) = mean([TM.(t).data.SpeedActual1(15000:end-2400)]);
%         TM.(t).checks.beltspeed(2) = mean([TM.(t).data.SpeedActual2(15000:end-2400)]);
%     else
%         npts = length(TM.(t).data.t);
%         start = round(npts*0.20); % start of pert
%         stop = round(npts*0.85); % end of pert
%         TM.(t).checks.beltspeed(1) = mean([TM.(t).data.SpeedActual1(start:stop)]);
%         TM.(t).checks.beltspeed(2) = mean([TM.(t).data.SpeedActual2(start:stop)]);
%     end

%     speedinfo(tct,:) = [TM.(t).checks.beltspeed(1) TM.(t).checks.beltspeed(2)];

    % add lines of average speed
%     plot([1 TM.(t).data.stride(end)],[1 1]*TM.(t).checks.beltspeed(1), 'b:', 'LineWidth',2);
%     plot([1 TM.(t).data.stride(end)],[1 1]*TM.(t).checks.beltspeed(2), 'r:', 'LineWidth',2);
% 
%     % add text average speed
%     text(max(xlim(gca))*0.3, 0.8, "belt1, L:" + num2str(round(TM.(t).checks.beltspeed(1),2))+" m/s");
%     text(max(xlim(gca))*0.3, 0.6, "belt2, R:" + num2str(round(TM.(t).checks.beltspeed(2),2))+" m/s");
%     %
%     legend("belt1, L", "belt2, R",'location', 'southeast');
end

T = "Nothing for now";
% T = array2table(speedinfo,...
%     'VariableNames',{'belt1Speed','belt2Speed'},...
%     'RowNames',trialTypes);
% disp(T)

% save plots and table
if ~exist("savepath","var"), savepath=[]; end
if ~isempty(savepath)
    cd(savepath)
    saveas(hf, subj+"_check_speed.pdf")
    saveas(hf, subj+"_check_speed.fig")
   % writetable(T,savepath+subj+"_check_speed.xlsx")
end