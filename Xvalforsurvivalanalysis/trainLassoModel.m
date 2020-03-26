function [opt_T,beta_Lasso,idx_top_feats,rs_train]=trainLassoModel(data,event_status,time_event)

cvfit = cvglmnet(data,[time_event,event_status],'cox','nfolds',10);
%cvglmnetPlot(cvfit);

num_idx_min=sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0);
num_idx_1se=sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_1se)~=0);
num_top_feats=max(num_idx_min,num_idx_1se);
if num_top_feats==0
    num_top_feats=round(size(data,1)/10);
end

tmp=find(cvfit.glmnet_fit.df<=num_top_feats);
idx_top_feats=find(cvfit.glmnet_fit.beta(:,tmp(end)));
beta_Lasso=cvfit.glmnet_fit.beta(:,tmp(end));
beta_Lasso=beta_Lasso(idx_top_feats);
%beta_Lasso = coxphfit(data(:,idx_top_feats),time_event,'Censoring',~event_status);

rs_train=data(:,idx_top_feats)*beta_Lasso;
opt_T=median(rs_train); % median risk score being used as threshold

end






