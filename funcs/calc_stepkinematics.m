function SK = calc_stepkinematics(M,GE,T,israte)
% data has been interpolated to the desired sampling rate
% biomech coords
% RHS = GE(:,1);
% LTO = GE(:,2);
% LHS = GE(:,3);
% RTO = GE(:,4);
% RHSn = GE(:,5);
SK.stride.time = (GE.RHSn-GE.RHS)*1/israte;
SK.step.time(:,1) = (GE.LHS-GE.RHS)*1/israte; % left
SK.step.time(:,2) = (GE.RHSn-GE.LHS)*1/israte; % right
% same as above
% SK.step.time0(:,1) = T.itime(LHS)-T.itime(RHS);
% SK.step.time0(:,2) = T.itime(RHSn)-T.itime(LHS);
SK.stride.freq = 1./SK.stride.time;
SK.step.freq = 1./SK.step.time;
SK.doublesupport.time(:,1) = (GE.LTO-GE.RHS)*1/israte; % first DS
SK.doublesupport.time(:,2) = (GE.RTO-GE.LHS)*1/israte; % second DS
% left step widths are negative
% column 1 = x
SK.step.width(:,1) = M.LHEE.X_ML(GE.LHS) - M.RHEE.X_ML(GE.LHS);
SK.step.width(:,2) = M.RHEE.X_ML(GE.RHSn) - M.LHEE.X_ML(GE.RHSn);
% average belt speed within a step
% belt 1 = left, belt 2 = right
for s = 1:length(GE.RHS)
    SK.step.beltspd(s,1) = mean(T.beltspd(GE.RHS(s):GE.LHS(s),1));
    SK.step.beltspd(s,2) = mean(T.beltspd(GE.LHS(s):GE.RHSn(s),2));
end
SK.step.dist = SK.step.beltspd.*SK.step.time;
% method Cesar used to calculate step length. takes the distance the belt
% traveled during the step time and the distance between the heel strikes
% column 2 = y, a-p
SK.step.length_cc(:,1) = ((M.LHEE.Y_AP(GE.LHS) - M.RHEE.Y_AP(GE.RHS))/cosd(T.pitch)) + SK.step.dist(:,1);
SK.step.length_cc(:,2) = ((M.RHEE.Y_AP(GE.RHSn) - M.LHEE.Y_AP(GE.LHS))/cosd(T.pitch)) + SK.step.dist(:,2);
% used Cesar's function here. same output as above. hooray.
SK.step.length_ccf(:,1) = steplength(M.LHEE.Y_AP, M.RHEE.Y_AP, GE.LHS, GE.RHS, SK.step.beltspd(:,1), SK.step.time(:,1), T.pitch);
SK.step.length_ccf(:,2) = steplength(M.RHEE.Y_AP, M.LHEE.Y_AP, GE.RHSn, GE.LHS, SK.step.beltspd(:,2), SK.step.time(:,2), T.pitch);
% common approach for calculating step length. takes the distance between
% the heel markers at a heel strike. which is what we said we did in the GP
% paper and will use in the correction.
% column 2 = y, a-p
SK.step.length(:,1) = (M.LHEE.Y_AP(GE.LHS) - M.RHEE.Y_AP(GE.LHS))/cosd(T.pitch);
SK.step.length(:,2) = (M.RHEE.Y_AP(GE.RHSn) - M.LHEE.Y_AP(GE.RHSn))/cosd(T.pitch);