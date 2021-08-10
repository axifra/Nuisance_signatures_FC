function [SLFOs, r_SLFOs_GS] = Generate_SLFOs(HR, resp, GS, Fs, ind_vol)
% Function to generate the SLFOs signal based on physiological recordings and the fMRI global signal

% Inputs:
% HR = heart rate signal extracted from the pulse oximeter
% resp = signal from the breathing belt
% GS = fMRI global signal (i.e. mean across brain voxels)
% Fs = sampling rate of the HR and resp signals
% ind_vol = instances where fMRI volumes were acquired expressed as indices at the timeline of the heart rate / respiration signals
% --> Example: assuming the heart rate and respiration are sampled at 10 Hz and the first two fMRI volumes occur at 5, 8, 11, ... seconds, the values of ind_vol would be ind_vol = [51, 81, 111, ...]

% Outputs: 
% SLFOs = systemic low frequency oscillations related to heart rate and breathing variations
% r_SLFOs_GS = correlation between SLFOs and the fMRI global signal

% Respiratory flow computation
resp = smooth(resp,Fs*1.5);
RF = diff(resp);
RF = [0;RF(:)]; 
RF = RF.^2;

% Computation of scan-specific physiological response functions to generate the SLFOs
options = optimoptions('fmincon','Display','off','Algorithm','interior-point','UseParallel',true,'MaxIterations',200,'MaxFunctionEvaluations',3000,'OptimalityTolerance',1e-8);    % 'PlotFcn','optimplotfval'
x0 = [  3.1136    2.4957   5.5743    0.9142    1.8940    2.9126   12.5214    0.4741 ];   % fmincon   -->  o.f.=47.20
ubMax = [10 3  10  3  15  3  20  3];
lbMin = [0 0.2  4  0  0  0.5  5  0.5];
bound_range = 2*ones(1,8);
lb = x0 - bound_range; lb = max(lb,lbMin);
ub = x0 + bound_range; ub = min(ub,ubMax);

h = @(P) Opt_HR_RF_GS(P,1/Fs,HR,RF,ind_vol,GS);
x_opt = fmincon(h,x0,[],[],[],[],lb,ub,[],options);        

[~, SLFOs] = h(x_opt);

r_SLFOs_GS = corr(SLFOs,GS);

end

function [output,yPred] = Opt_HR_RF_GS(P,Ts,HR,RF,ind_vol,GS)

    warning('off','all')

    t1=P(1);  s1=P(2);
    t2=P(3);  s2=P(4);
    t3=P(5);  s3=P(6);
    t4=P(7);  s4=P(8);

    NV = length(GS);
    t0 = 0; t_win= t0:Ts:60;

    a1= sqrt(t1)/s1; a2= sqrt(t1)*s1 ;
    IR = t_win.^a1.*exp(-t_win/a2); IR_cardF=IR/max(IR);   % plot(t_win,IR)
    a1= sqrt(t2)/s2; a2= sqrt(t2)*s2 ;
    IR = t_win.^a1.*exp(-t_win/a2); IR_cardS=IR/max(IR);   % plot(t_win,IR)

    a1= sqrt(t3)/s3; a2= sqrt(t3)*s3 ;
    IR = t_win.^a1.*exp(-t_win/a2); IR_respF=IR/max(IR);   % plot(t_win,IR)
    a1= sqrt(t4)/s4; a2= sqrt(t4)*s4 ;
    IR = t_win.^a1.*exp(-t_win/a2); IR_respS=IR/max(IR);   % plot(t_win,IR)

    HR_Fconv=conv(HR,IR_cardF); HR_Fconv=HR_Fconv(ind_vol);
    HR_Sconv=conv(HR,IR_cardS); HR_Sconv=HR_Sconv(ind_vol);
    RF_Fconv=conv(RF,IR_respF); RF_Fconv=RF_Fconv(ind_vol);
    RF_Sconv=conv(RF,IR_respS); RF_Sconv=RF_Sconv(ind_vol);
    regr=[ones(NV,1),HR_Fconv,HR_Sconv,RF_Fconv,RF_Sconv];

    B = regr\GS;     
    yPred = regr*B;
    output = 1 - corr(yPred,GS) ;
end