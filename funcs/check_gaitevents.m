function check_gaitevents(GRF, GE, threshold,label)
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
plot(GRF);
hold on
plot(GE.RHS,ones(nstrides,1)*threshold, 'rx','LineWidth',2);
plot(GE.LTO,ones(nstrides,1)*threshold, 'bo','LineWidth',2);
plot(GE.LHS,ones(nstrides,1)*threshold, 'bx','LineWidth',2);
plot(GE.RTO,ones(nstrides,1)*threshold, 'ro','LineWidth',2);
legend(["GRF1" "GRF2" "RHS" "LTO" "LHS" "RTO"])
title(label)