function [pval,hr,ci_up,ci_low,fig] = run_KM_XVal(features,event,eventTime,percentageTest,print,plot)
%RUN_KM_XVAL Runs a Kaplan-Meier cross-validation as suggested here: 
%Simon RM, Subramanian J, Li MC, Menezes S. Using cross-validation to 
%evaluate predictive accuracy of survival risk classifiers based on 
%high-dimensional data. Brief Bioinform. 2011;12(3):203–214. 
%doi:10.1093/bib/bbr001


if nargin<4
    percentageTest=30;
end
if nargin<5
    print=true;
end
if nargin<6
    plot=true;
end

numCases=length(event);

numTestCases=round(numCases*percentageTest/100);

fullList=1:numCases;
labels=nan(1,numCases);

i=1;
while i<=numCases
    
    lim=i+numTestCases-1;
    if lim>numCases
        lim=numCases;
    end
    idxTest=i:lim;
    idxTrain=fullList(~ismember(fullList,idxTest));
    
    [thrs,beta,idxTopFeats]=trainLassoModel(features(idxTrain,:),event(idxTrain),eventTime(idxTrain));
    
    riskTest = features(idxTest,idxTopFeats)*beta;
    labels(idxTest)=riskTest>thrs;
    
    i=i+numTestCases;
end
labels=logical(labels);

%% Plotting & Printing

groupvar=cell(length(labels),1);
groupvar(labels)={'high risk'};
groupvar(~labels)={'low risk'};

optionsKM.NoPlot = ~plot;
optionsKM.Print = false;

[pval,fig,stats]=MatSurv(eventTime,event,groupvar,optionsKM);

hr=stats.HR_logrank;
ci_low=stats.HR_95_CI_logrank(1);
ci_up=stats.HR_95_CI_logrank(2);

if print
    fprintf('Results: p-val = %.4f, HR = %.4f (%.4f - %.4f)\n',pval,hr,ci_low,ci_up);
end

end

