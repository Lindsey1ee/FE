function SK = calc_stepkinematics(M,GE,T,israte)

% data has been interpolated to the desired sampling rate
% biomech coords

RHS = GE(:,1);
LTO = GE(:,2);
LHS = GE(:,3);
RTO = GE(:,4);
RHSn = GE(:,5);

SK.stride.time = (RHSn-RHS)*1/israte;
SK.step.time(:,1) = (LHS-RHS)*1/israte; % left
SK.step.time(:,2) = (RHSn-LHS)*1/israte; % right

% same as above
% SK.step.time0(:,1) = T.itime(LHS)-T.itime(RHS);
% SK.step.time0(:,2) = T.itime(RHSn)-T.itime(LHS);

SK.stride.freq = 1./SK.stride.time;
SK.step.freq = 1./SK.step.time;

SK.doublesupport.time(:,1) = (LTO-RHS)*1/israte; % first DS
SK.doublesupport.time(:,2) = (RTO-LHS)*1/israte; % second DS

% left step widths are negative
% column 1 = x
SK.step.width(:,1) = M.LHEE(LHS,1) - M.RHEE(LHS,1);
SK.step.width(:,2) = M.RHEE(RHSn,1) - M.LHEE(RHSn,1);

% average belt speed within a step
% belt 1 = left, belt 2 = right
for s = 1:length(RHS)
    SK.step.beltspd(s,1) = mean(T.beltspd(RHS(s):LHS(s),1));
    SK.step.beltspd(s,2) = mean(T.beltspd(LHS(s):RHSn(s),2));
end

SK.step.dist = SK.step.beltspd.*SK.step.time;

% method Cesar used to calculate step length. takes the distance the belt
% traveled during the step time and the distance between the heel strikes
% column 2 = y, a-p
SK.step.length_cc(:,1) = ((M.LHEE(LHS,2) - M.RHEE(RHS,2))/cosd(T.pitch)) + SK.step.dist(:,1);
SK.step.length_cc(:,2) = ((M.RHEE(RHSn,2) - M.LHEE(LHS,2))/cosd(T.pitch)) + SK.step.dist(:,2);

% used Cesar's function here. same output as above. hooray.
SK.step.length_ccf(:,1) = steplength(M.LHEE(:,2), M.RHEE(:,2), LHS, RHS, SK.step.beltspd(:,1), SK.step.time(:,1), T.pitch);
SK.step.length_ccf(:,2) = steplength(M.RHEE(:,2), M.LHEE(:,2), RHSn, LHS, SK.step.beltspd(:,2), SK.step.time(:,2), T.pitch);

% common approach for calculating step length. takes the distance between
% the heel markers at a heel strike. which is what we said we did in the GP
% paper and will use in the correction.
% column 2 = y, a-p 
SK.step.length(:,1) = (M.LHEE(LHS,2) - M.RHEE(LHS,2))/cosd(T.pitch);
SK.step.length(:,2) = (M.RHEE(RHSn,2) - M.LHEE(RHSn,2))/cosd(T.pitch);

% step variance



