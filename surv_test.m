function [surv_significance_test] = surv_test(test_data, beta, cur_selected_idx, surv_days_test, censor_test, risk_threshold)
%Testing script for survival analysis
% indentify threshold and beta using surv_train
% Cheng Lu (10/27)
% Prateek Prasanna (10/29)

rs_val=test_data(:,cur_selected_idx)*beta;
%figure;subplot(2,1,1);bar(sort(rs_train)); xlabel('patient number/index in trainig cohort'); ylabel('image based risk score in trainig cohort');
%subplot(2,1,2);bar(sort(rs_val));xlabel('patient number/index in validation cohort');ylabel('image based risk score in validation cohort');
% set(gcf,'FontSize',16);

% opt_T=median(rs_val);
labels_pred=logical(rs_val>risk_threshold);
groupvar=[];
for i=1:length(logical(labels_pred))
    if labels_pred(i)==1
        groupvar{i}='high risk';
    else
        groupvar{i}='low risk';
    end
end
optionsKM.NoPlot=0;
% kick out the ones with survival time =0 
%idx_survival_usable=~((survival_time_LUSC==0)| isnan(survival_time_LUSC)|(survival_time_LUSC>60*30));
%[opt_p,~,opt_stats]=MatSurv(survival_time_LUSC(idx_survival_usable)/30, ~logical(survival_censor_LUSC(idx_survival_usable)), groupvar(idx_survival_usable),'Xstep',10,'Title','5-years survival in validation set',optionsKM);
group1=find(labels_pred);
group2=find(~labels_pred);
surv_significance_test=logrank([surv_days_test(group1) censor_test(group1)],[surv_days_test(group2) censor_test(group2)],0.05,0);

end

