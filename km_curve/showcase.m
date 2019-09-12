% ground truth: grooups_total
% classifier label: predicted_label_SVM_total
% survival time in days: time_f

data_total=(feats_total(:,9));
for k=1:66
    if data_total(k,1)<1.01
        data_total(k,1)=0;
    else
        data_total(k,1)=1;
    end
end
    
%% 

group1=find(groups_total);
group2=find(~groups_total);
cens=zeros(66,1);
logrank([time_f(group1)' cens(group1)],[time_f(group2)' cens(group2)],0.05,0);

% figure;
% info = plotKaplanMeier(IDs_followup(:,1),labels_validation_65,ID_label_65(:,2));
