function gaitevents_perts(GRF, GE, itime, ptime, threshold,label)
% plots gait events against the ground reaction forces. This is the gold
% standard check. Sometimes we can't use GRFs to identify the gait events
% because of weird steps (step on both belts. shuffle) and typically, the
% mocap way is more consistent
%
% INPUT VARIABLES:
%       GRF         Ground reaction force, expected to be vertical though
%                   technically doesn't have to be in a visual check
%       GE          Gait event structure
%       threshold   what GRF value to plot gait events at, often 20 N
%       label       adds title, perhaps of subj and condition being checked
%
% OUTPUT VARIABLES:
%
% written by Helen J. Huang (Aug 2023)
% updates tracked in git
nstrides = height(GE);
% plot gait events versus vertical GRF

figure
plot(itime, GRF);
hold on
plot(itime(GE.RHS),ones(nstrides,1)*threshold, 'rx','LineWidth',2);
plot(itime(GE.LTO),ones(nstrides,1)*threshold, 'bo','LineWidth',2);
plot(itime(GE.LHS),ones(nstrides,1)*threshold, 'bx','LineWidth',2);
plot(itime(GE.RTO),ones(nstrides,1)*threshold, 'ro','LineWidth',2);
legend(["GRF1" "GRF2" "RHS" "LTO" "LHS" "RTO"])
title(label)

% hold on 

plot(ptime,ones(length(ptime),1)*threshold, 'k^','LineWidth',2);


figure
plot(ptime,'o')
hold on 

x = itime(GE.RHS);
x = x(x>ptime(1));
x = x(x<ptime(end));

plot(x, 'x')

numstrides = height(GE);

temp = zeros(length(ptime),1);
GEtimes(:,1) = itime(GE.RHS);
GEtimes(:,2) = itime(GE.LTO);
GEtimes(:,3) = itime(GE.LHS);
GEtimes(:,4) = itime(GE.RTO);
GEtimes(:,5) = itime(GE.RHSn);

s = 1;
for p = 1:length(ptime)
    if ptime(p) >= itime(GE.RHS(s)) && ptime(p) <= itime(GE.RHSn(s))
        temp(s) = ptime(p);
        
    else
%         disp('here')
        temp(s) = -1;
    end
    s=s+1;
end

[correlation, lag] = xcorr(signal1, signal2);
