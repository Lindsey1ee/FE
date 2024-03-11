function [TM, MC, T] = check_freq(TM,MC,trialTypes,savepath)

%% freq check

for t = trialTypes
%want to find the frequency the sway occurs at
%Find Onset
swaySign = sign(TM.(t).data.SwayActual);

%Find index where positive overshoots occur (where sign=1)
pos_idx = find(swaySign == 1);
swaySign(pos_idx) = 0 ;%set positive values to zero

time_idx = find(diff(0.03*swaySign) < 0); %finds where perturbation occurs and duration of perturbation

disp('time occured at: ')
disp(TM.(t).data.t(time_idx))

%find actual frequencies experienced 
pert_time = diff(TM.(t).data.t(time_idx));
disp('time between pert:')
disp(mean(pert_time))

avg_freq = 1/mean(pert_time);
disp('actual frequency:')
disp(avg_freq)

end

% create table of data for easier displaying and saving
T = array2table(avg_freq,'VariableNames',...
    {'avgFreq'},'RowNames',trialTypes);

disp(T)

% save plots and table
if ~exist("savepath","var"), savepath=[]; end
if ~isempty(savepath)
    cd(savepath)
    print(hf,savepath+subj+"_check_freq",'-dpdf','-fillpage')
    saveas(hf, savepath+subj+"_check_freq.fig")
    writetable(T,savepath+subj+"_check_freq.xlsx")
end




