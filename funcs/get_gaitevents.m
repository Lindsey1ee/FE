function GE = get_gaitevents(markers,lpfc,srate,window)

% We identified  heel strikes as the most anterior position of the
% calcaneus markers and toe-offs as the most posterior position of the
% second  metatarsal  head  markers  for  each  foot

% low pass filter the data, 4th order butterworth
[b, a] = butter(4,2.0*1/srate*lpfc);

for m = ["RHEE" "RTOE" "LHEE" "LTOE"]
    if istable(markers.(m))
            markers.(m) = table2array(markers.(m));
    end
    filteredmarkers.(m) = filtfilt(b,a,markers.(m));
end

% quick plot of a filtered vs unfiltered marker, whichever was last in the
% loop above

% works assuming lab axes have been rotated to be the more intuitive
% column 2 = y, a-p
figure
plot(markers.(m)(:,2),'LineWidth',2);
hold on
plot(filteredmarkers.(m)(:,2),'r:','LineWidth',2);
title(m); legend(["unfilt" "filt"]); xlabel("frames");
ylim([-2 2]);

% works assuming lab axes have been rotated to be the more intuitive
% column 2 = y, a-p
for m = ["RHEE" "LHEE"]
    [pks.(m),locs.(m)] = findpeaks(filteredmarkers.(m)(:,2), 'MinPeakDistance', srate*0.75);
end

% take negative so that most posterior is positive for find peaks function
% column 2 = y, a-p
for m = ["RTOE" "LTOE"]
    [pks.(m),locs.(m)] = findpeaks(-filteredmarkers.(m)(:,2), 'MinPeakDistance', srate*0.75);
end

% we don't need to have the gait events for the strides prior to or after
% the data to be analyzed. The window identifies the section of the data
% where we need to identify gait events and expect to have consistent
% strides and gait events. Plot marker trajectory and peaks
figure
ct=1;
for m = ["RHEE" "RTOE" "LHEE" "LTOE"]
    remove = find(locs.(m)<window(1));
    locs.(m)(remove)=[]; pks.(m)(remove)=[];

    remove = find(locs.(m)>window(2));
    locs.(m)(remove)=[]; pks.(m)(remove)=[];

%     remove = find(pks.(m)>mean(pks.(m))*4);

    subplot(2,2,ct)
    plot(filteredmarkers.(m)(:,2));
    hold on
    if m == "RHEE" || m == "LHEE"
        plot(locs.(m), pks.(m),'o','LineStyle','none')
    elseif m == "RTOE" || m == "LTOE"
        plot(locs.(m), -pks.(m),'o','LineStyle','none')
    end
    title(m); legend([m "peaks"]);
    ct=ct+1;
end

% associate locations/indices of peaks with gait events labels
RHS = locs.RHEE; LHS = locs.LHEE;
RTO = locs.RTOE; LTO = locs.LTOE;

% put all gait events into a column to sort. In real life, they must go in
% order with time so sort. 
temp = sort([RHS; LTO; LHS; RTO]);

% create a variable tags with the marker event associated with each value
% in temp
tags = strings(size(temp));
for m = ["RHEE" "LTOE" "LHEE" "RTOE"]
    [Lia,Locb] = ismember(temp,locs.(m));
    if m == "RHEE", event = "RHS";
    elseif m == "LTOE", event = "LTO";
    elseif m == "LHEE", event = "LHS";
    elseif m == "RTOE", event = "RTO";
    end
    tags(find(Locb~=0)) = event;
end

% If a gait event is not identified (maybe took a weird step) then there
% will not be a valid order. so a stride would need to have a sequence of
% RHS, LTO, LHS, RTO. Find that pattern to identify valid strides. 
% First find each instance of RHS
strides = find(ismember(tags, "RHS"));

% Save the valid strides into GE structure and the tags
GE.goodstrides = []; GE.badstrides = [];
for s = 1:length(strides)-1

    % a whole stride should have 5 events, RHS, LTO, LHS, RTO, RHSn
    % if there are not 5 elements, then not all events found
    if length(temp(strides(s):strides(s+1))) == 5 && ...
            all(tags(strides(s):strides(s+1))==["RHS"; "LTO"; "LHS"; "RTO"; "RHS"])
        GE.events(s,:) = temp(strides(s):strides(s+1));
        GE.tags(s,:) = tags(strides(s):strides(s+1));
        GE.goodstrides = [GE.goodstrides; s];
    else
        disp("full stride not identified, stride " + s)
        disp(tags(strides(s):strides(s+1))')
        GE.badstrides = [GE.badstrides; s];
    end 
end

GE.lpfc = lpfc;