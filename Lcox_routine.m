%% %%%% D do the features selection using elastic cox model and get the model
addpath('F:\Nutstore\Nutstore\PathImAnalysis_Program\Program\Miscellaneous\glmnet_matlab');
% note that the right censor should be 0, not 1 as I normally did
% previously in the function

%data_cLoCoM_train_truncated is the data from training set 

options.alpha = 1;
y = [survival_time_train_CCF_YTMA, ~logical(survival_censor_train_CCF_YTMA)]; % load survival data 

feat_name=feature_list_truncated;

% cvfit = cvglmnet(data_cLoCoM_SC_CCF_YTMA,y,'cox',options);
cvfit = cvglmnet(data_cLoCoM_train_truncated,y,'cox','nfolds',10);
cvglmnetPlot(cvfit)

% find top features
sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0)
sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_1se)~=0)

idx_min=find(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0);
idx_1se=find(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_1se)~=0);
%%% or get a fixed number of features
num_top_feats=round(size(data_cLoCoM_train_truncated,1)/10);
% num_top_feats=sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0);
tmp=find(cvfit.glmnet_fit.df<=num_top_feats);
idx_fix_top_feats=find(cvfit.glmnet_fit.beta(:,tmp(end)));
beta_Lasso=cvfit.glmnet_fit.beta(:,tmp(end));
beta_Lasso=beta_Lasso(idx_fix_top_feats);
%% set the top features and show the feature name here
feat_name_idx_min=feat_name(idx_min)';
feat_name_idx_1se=feat_name(idx_1se)';
feat_name_fix_top_feats=feat_name(idx_fix_top_feats)';

% cur_selected_idx=idx_min;
cur_selected_idx=idx_fix_top_feats;
%% check feature distribution here
% for i=1:length(cur_selected_idx)
%     curidx=cur_selected_idx(i);
% %     mdl = fitlm(data_all(index_POI_indata,curidx),overallSurvival);
%     figure(i);
%     plot(survival_time_train_CCF_YTMA,data_cLoCoM_train_truncated(:,curidx),'o');
%     xlabel('survival days');
%     ylabel('feature value');
%     title(feat_name{curidx});
% end

%% get the updated beta from the cox model again or use the beta from the elastic net
% cur_selected_idx
% [beta,logl,H,stats] = coxphfit(data_cLoCoM_train_truncated(:,cur_selected_idx),survival_time_train_CCF_YTMA,'Censoring',logical(survival_censor_train_CCF_YTMA));

beta=beta_Lasso;
%% get the risk score
addpath('F:\Nutstore\Nutstore\PathImAnalysis_Program\Program\Miscellaneous\aebergl-MatSurv-d0f1877');
rs_train=data_cLoCoM_train_truncated(:,cur_selected_idx)*beta;
figure;bar(sort(rs_train));

% check the KM curve for training data
opt_T=median(rs_train);
labels_pred=logical(rs_train>opt_T);
groupvar=[];
for i=1:length(logical(labels_pred))
    if labels_pred(i)==1
        groupvar{i}='high risk';
    else
        groupvar{i}='low risk';
    end
end
optionsKM.NoPlot=0;
idx_survival_usable=~((survival_time_train_CCF_YTMA==0)| isnan(survival_time_train_CCF_YTMA)|(survival_time_train_CCF_YTMA>60));
[opt_p,~,opt_stats]=MatSurv(survival_time_train_CCF_YTMA(idx_survival_usable), ~logical(survival_censor_train_CCF_YTMA(idx_survival_usable)), groupvar(idx_survival_usable),'Xstep',10,'Title','5 years survival in training set',optionsKM);
%% %%%%  E find the optimal thresold from training set, I used the median or you can do the exhasutive search here
% %%% find best T here, based off the HR
% optionsKM.NoPlot=1;
% p=[];
% stats=[];
% HR=[];
% pred_train=rs_train;
% for t=round(length(pred_train)*0.1):round(length(pred_train)-length(pred_train)*0.1)
%     groupvar=[];
%     T=pred_train(t);
%     all_T(t)=T;
%     labels_pred=logical(pred_train<T);
%  
%     for i=1:length(logical(labels_pred))
%         if labels_pred(i)==1
%             groupvar{i}='long-term';
%         else
%             groupvar{i}='short-term';
%         end
%     end
%     [p(t),~,stats{t}]=MatSurv(survival_time_train_CCF_YTMA, ~logical(survival_censor_train_CCF_YTMA), groupvar,'Xstep',10,'Title','MatSurv KM-Plot',optionsKM);
%     HR(t)=stats{t}.HR_logrank;
% end
% p(p==0)=1;
% [val_min,idxmin]=min(p);
% 
% HR(HR==0)=1;
% [HR_min,idxmin]=min(HR);
% 
% 
% opt_T=all_T(idxmin);
% labels_pred=logical(pred_train<opt_T);
% groupvar=[];
% for i=1:length(logical(labels_pred))
%     if labels_pred(i)==1
%         groupvar{i}='long-term';
%     else
%         groupvar{i}='short-term';
%     end
% end
% optionsKM.NoPlot=0;
% [opt_p,~,opt_stats]=MatSurv(survival_time_train, ~logical(survival_censor_train), groupvar,'Xstep',10,'Title','model in training set',optionsKM);
%% %%%% F feed the validation data into the model and do analysis
%data_cLoCoM_LUAD_SC_truncated is the feature data from validation set

rs_val=data_cLoCoM_LUAD_SC_truncated(:,cur_selected_idx)*beta;
figure;subplot(2,1,1);bar(sort(rs_train)); xlabel('patient number/index in trainig cohort'); ylabel('image based risk score in trainig cohort');
subplot(2,1,2);bar(sort(rs_val));xlabel('patient number/index in validation cohort');ylabel('image based risk score in validation cohort');
% set(gcf,'FontSize',16);

% opt_T=median(rs_val);
labels_pred_TCGA=logical(rs_val>opt_T);
groupvar=[];
for i=1:length(logical(labels_pred_TCGA))
    if labels_pred_TCGA(i)==1
        groupvar{i}='high risk';
    else
        groupvar{i}='low risk';
    end
end
optionsKM.NoPlot=0;
% kick out the ones with survival time =0 
idx_survival_usable=~((survival_time_LUSC==0)| isnan(survival_time_LUSC)|(survival_time_LUSC>60*30));
[opt_p,~,opt_stats]=MatSurv(survival_time_LUSC(idx_survival_usable)/30, ~logical(survival_censor_LUSC(idx_survival_usable)), groupvar(idx_survival_usable),'Xstep',10,'Title','5-years survival in validation set',optionsKM);





