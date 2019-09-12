function [surv_significance, risk_threshold,  beta, cur_selected_idx] = surv_train(train_data,surv_days,censor_train)

%Training script for survival analysis
%Creates radiomic risk score based on cox regression analysis, identifies
%optimal threshold to divide patients into low risk and high risk groups
% Cheng Lu (10/27)
% Prateek Prasanna (10/29)

%Inputs: train_data -> nxp feature matric
%surv_days -> nx1 vector of survival days
%censor_train -> nx1 vector, 0 for non-censored, 1 for censored

options.alpha = 1;
y = [surv_days, ~logical(censor_train)]; % load survival data 

feat_name=1:length(train_data);

% Fit Cox regression model in LASSO
cvfit = cvglmnet(train_data,y,'cox','nfolds',10);
cvglmnetPlot(cvfit)

% find top features corresponding to least errors
sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0)
sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_1se)~=0)

idx_min=find(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0);
idx_1se=find(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_1se)~=0);
%%% or get a fixed number of features
%num_top_feats=round(size(feat_data,1)/10);
num_top_feats=sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0);
tmp=find(cvfit.glmnet_fit.df<=num_top_feats);
idx_fix_top_feats=find(cvfit.glmnet_fit.beta(:,tmp(end)));

%Find beta values in regression model
beta_Lasso=cvfit.glmnet_fit.beta(:,tmp(end));
beta_Lasso=beta_Lasso(idx_fix_top_feats);
%% set the top features and show the feature name here
feat_name_idx_min=feat_name(idx_min)';
feat_name_idx_1se=feat_name(idx_1se)';
feat_name_fix_top_feats=feat_name(idx_fix_top_feats)';

% cur_selected_idx=idx_min;
cur_selected_idx=idx_fix_top_feats;

%% get the updated beta from the cox model again or use the beta from the elastic net
% cur_selected_idx
% [beta,logl,H,stats] = coxphfit(data_cLoCoM_train_truncated(:,cur_selected_idx),survival_time_train_CCF_YTMA,'Censoring',logical(survival_censor_train_CCF_YTMA));

beta=beta_Lasso;
%% get the risk score
rs_train=train_data(:,cur_selected_idx)*beta;
figure;bar(sort(rs_train));

% check the KM curve for training data
risk_threshold=median(rs_train);
labels_pred=logical(rs_train>risk_threshold);
groupvar=[];
groupvar_num=[];
for i=1:length(logical(labels_pred))
    if labels_pred(i)==1
        groupvar{i}='high risk';
    else
        groupvar{i}='low risk';
    end
end
optionsKM.NoPlot=0;
idx_survival_usable=~((surv_days==0)| isnan(surv_days)|(surv_days>0));
group1=find(labels_pred);
group2=find(~labels_pred);
surv_significance=logrank([surv_days(group1) censor_train(group1)],[surv_days(group2) censor_train(group2)],0.05,0);



%[risk_threshold,~,training_stats]=MatSurv(surv_days(idx_survival_usable), ~logical(censor_train(idx_survival_usable)), groupvar(idx_survival_usable),'Xstep',10,'Title','X years survival in training set',optionsKM);
end

