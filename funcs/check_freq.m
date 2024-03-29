function [TM, T] = check_freq(TM,trialTypes,savepath)

%% freq check
ct=1;
if trialTypes == "nopert"
     pert_time = 0;
     avg_freq = 0;
     TM.(t).checks.freq.onset_time = 0;
else
    for t = trialTypes
    %want to find the frequency the sway occurs at
    %Find Onset
    %swaySign = sign(TM.(t).data.SwayActual);
    swaySign = zeros(size(TM.(t).data.SwayActual));
    swaySign(TM.(t).data.SwayActual < -0.002) = -1;
    
    %Find index where positive overshoots occur (where sign=1)
    % pos_idx = find(swaySign == 1);
    % swaySign(pos_idx) = 0 ;%set positive values to zero
    
    time_idx = find(diff(0.03*swaySign) < 0); %finds where perturbation occurs and duration of perturbation
    
    figure
    plot(TM.(t).data.t, TM.(t).data.SwayActual)
    hold on
    plot(TM.(t).data.t, 0.03*swaySign)
    
    % if swaySign == -1 & TM.(t).data.SwayActual >= -0.0001
    % z = find(swaySign == -1)
    % y = find(TM.(t).data.SwayActual >= -0.00001)
    % %c = find(swaySign == -1 & TM.(t).data.SwayActual >= -0.0001)
    % end 
    
    % disp('time occured at: ')
    % disp(TM.(t).data.t(time_idx))
    
    %find actual frequencies experienced 
    pert_time = diff(TM.(t).data.t(time_idx));
    disp('time between pert:')
    disp(mean(pert_time))
    
    avg_freq = 60/mean(pert_time);
    disp('actual frequency:')
    disp(avg_freq)
 
   TM.(t).checks.freq.pert_time = pert_time;
   TM.(t).checks.freq.avg_freq = avg_freq;
   TM.(t).checks.freq.onset_time = TM.(t).data.t(time_idx);

   freqinfo(ct,1) = TM.(t).checks.freq.avg_freq;

   ct=ct+1;
 end




end

% % create table of data for easier displaying and saving
T = array2table(freqinfo,'VariableNames',...
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




