%cox elastic net test
% % % N=1000; p=30;
% % % nzc=p/3;
% % % %1000 observations with 30 predictors
% % % x=randn(N,p);
% % % beta=randn(nzc,1);
% % % fx=x(:,1:nzc)*beta/3;
% % % hx=exp(fx);
% % % %generates an N-by-1 array containing random numbers from the exponential
% % % %distribution with mean parameter 1./hx
% % % ty=exprnd(1./hx,N,1);
% % % %generates an N-by-1 array containing random numbers from the binomial
% % % %distribution with parameter, 1, nmber of trials, 0.3, probability of
% % % %success for each trial
% % % tcens=binornd(1,0.3,N,1);
% % % %contstructing time with censoring data
% % % y=cat(2,ty,1-tcens);
% % %
% % % fit=glmnet(x,y,'cox');
% % % glmnetPlot(fit);
clear, clc, close all
[~,~,raw] = xlsread('L:\project\PathRadioFuse\CCF_all.xlsx');
patientNamesInSheet = raw(:,1);
patientOSinMonth = raw(:,17);
patientCensoring = raw(:,18);

addpath('L:\project\predictiveOfChemo')
featureFolder = 'L:\project\predictiveOfChemo\featuresChemoReceived\';
matFiles = dir([featureFolder,'*.mat']);
matNames = {matFiles.name}';


[~,~,patientInfoETMA] = xlsread('L:\Data\TMA700\fromVVviaKaustav.xlsx');
excelNames = patientInfoETMA(1:815,2);
excelRecurrence = patientInfoETMA(1:815,3:5);
excelSuvival = patientInfoETMA(1:815,6:11);
overallSurvival = [];
recurrenceVector = [];
censoringVector = [];
counter2 = 0; counter3 = 0;
%retrieving the survival data
for i = 1:length(matNames)
    name = matNames{i};
    tmaID = regexp(name,'\d*', 'match', 'once');
    idx = contains(excelNames,tmaID);
    Index = find(idx==true);
    switch length(Index)
        case 0
            fprint('patients %s not in spread sheet', matNames{i})
        case 1
            
            if excelRecurrence{Index,3} ~= 1
                recurrenceVector = [recurrenceVector; excelRecurrence{Index,1}];
            else
                recurrenceVector = [recurrenceVector; NaN];
            end
            
            if excelRecurrence{Index,1} == 1 || ~isempty(excelSuvival{Index,3})
                censoringVector = [censoringVector; 0];
            else
                censoringVector = [censoringVector; 1];
            end
            
            
            osFromSheet = excelSuvival{Index,4};
            %if the os are not directly recorded from the spread sheet,
            %trying to calculate from individual date
            %all patients have the date of diagnosis
            if isnan(osFromSheet) || osFromSheet == 0
                %if the date of death is not NaN, used it with date of
                %diagnosis to calculate OS
                if ~isnan(excelSuvival{Index,3})
                    e = etime(datevec(excelSuvival{Index,3}),datevec(excelSuvival{Index,5}));
                    overallSurvival = [overallSurvival; e/365/24/3600];
                    %else if the last follow up time is not NaN, use it to
                    %calculate OS
                elseif ischar(excelSuvival{Index,6})
                    e = etime(datevec(excelSuvival{Index,6}),datevec(excelSuvival{Index,5}));
                    overallSurvival = [overallSurvival; e/365/24/3600];
                    %don't consider other situation since the current data
                    %already been handled by above condition
                else
                    overallSurvival = [overallSurvival; NaN];
                end
            else
                overallSurvival = [overallSurvival; osFromSheet];
            end
        case 2
            if excelRecurrence{Index(1),3} ~= 1
                recurrenceVector = [recurrenceVector; excelRecurrence{Index(1),1}];
            elseif excelRecurrence{Index(2),3} ~= 1
                recurrenceVector = [recurrenceVector; excelRecurrence{Index(2),1}];
            else
                recurrenceVector = [recurrenceVector; NaN];
            end
            
            if excelRecurrence{Index(1),1} == 1 || isempty(excelSuvival{Index(1),3})
                censoringVector = [censoringVector; 0];
            else
                censoringVector = [censoringVector; 1];
            end            
            
            osFromSheet1 = excelSuvival{Index(1),4};
            osFromSheet2 = excelSuvival{Index(2),4};
            %if the os are not directly recorded from the spread sheet,
            %trying to calculate from individual date
            %all patients have the date of diagnosis
            if isnan(osFromSheet1) || osFromSheet1 == 0
                %if the date of death is not NaN, used it with date of
                %diagnosis to calculate OS
                if ~isnan(excelSuvival{Index(1),3})
                    e = etime(datevec(excelSuvival{Index(1),3}),datevec(excelSuvival{Index(1),5}));
                    osFromSheet1 = e/365/24/3600;
                    %else if the last follow up time is not NaN, use it to
                    %calculate OS
                elseif ischar(excelSuvival{Index(1),6})
                    e = etime(datevec(excelSuvival{Index(1),6}),datevec(excelSuvival{Index(1),5}));
                    osFromSheet1 = e/365/24/3600;
                    %don't consider other situation since the current data
                    %already been handled by above condition
                else
                    osFromSheet1 = NaN;
                end
            else
                osFromSheet1 = osFromSheet1;
            end
            if isnan(osFromSheet2) || osFromSheet2 == 0
                %if the date of death is not NaN, used it with date of
                %diagnosis to calculate OS
                if ~isnan(excelSuvival{Index(2),3})
                    e = etime(datevec(excelSuvival{Index(2),3}),datevec(excelSuvival{Index(2),5}));
                    osFromSheet2 = e/365/24/3600;
                    %else if the last follow up time is not NaN, use it to
                    %calculate OS
                elseif ischar(excelSuvival{Index(2),6})
                    e = etime(datevec(excelSuvival{Index(2),6}),datevec(excelSuvival{Index(2),5}));
                    osFromSheet2 = e/365/24/3600;
                    %don't consider other situation since the current data
                    %already been handled by above condition
                else
                    osFromSheet2 = NaN;
                end
            else
                osFromSheet2 = osFromSheet2;
            end
            
            if ~isnan(osFromSheet1) && ~isnan(osFromSheet2)
                overallSurvival = [overallSurvival; mean([osFromSheet1 osFromSheet2])];
            elseif ~isnan(osFromSheet1)
                overallSurvival = [overallSurvival; osFromSheet1];
            elseif ~isnan(osFromSheet2)
                overallSurvival = [overallSurvival; osFromSheet2];
            else
                overallSurvival = [overallSurvival; NaN];
            end
    end
end

%ToDo: using the original excel sheet to generate the survival information,
%hoping to collect more patients with survival time
featurePatientMatrix = [];
a = 0;
for i = 1:length(matNames)
    nameIDNumber = regexp(matNames{i}, '\d*', 'match', 'once');
    idx = contains(patientNamesInSheet,nameIDNumber);
    Index = find(idx==true);
    if ~isempty(Index)
        a = a + 1;
    end
    
    load([featureFolder,matNames{i}])
    feature = averagingTileFeatures(featPerPatient);
    featurePatientMatrix = [featurePatientMatrix;feature];
end

% N=1000; p=30;
% nzc=p/3;
% %1000 observations with 30 predictors
% x=randn(N,p);^
% beta=randn(nzc,1);
% fx=x(:,1:nzc)*beta/3;
% hx=exp(fx);
% %generates an N-by-1 array containing random numbers from the exponential
% %distribution with mean parameter 1./hx
% ty=exprnd(1./hx,N,1);
% %generates an N-by-1 array containing random numbers from the binomial
% %distribution with parameter, 1, nmber of trials, 0.3, probability of
% %success for each trial
% tcens=binornd(1,0.3,N,1);
% %contstructing time with censoring data
% y=cat(2,ty,1-tcens);
featureMatrix = featurePatientMatrix;
iinf = isinf(featureMatrix);
featureMatrix(iinf) = NaN;
inan = isnan(featureMatrix);
featuresMean = nanmean(featureMatrix);
featuresMean = repmat(featuresMean,size(featureMatrix,1),1);
featureMatrix(inan) = featuresMean(inan);
featureMatrix(iinf) = featuresMean(iinf);
inan = isnan(featureMatrix);
featureMatrix(inan) = 0.0001;
tmp = repmat(mean(featureMatrix),size(featureMatrix,1),1);
featureMatrix(featureMatrix == 0) = tmp(featureMatrix == 0);
featureMean = mean(featureMatrix);
featureMin = min(featureMatrix);
featureMax = max(featureMatrix);
featureMatrix = (featureMatrix-featureMean)./(featureMax-featureMin);


% recurrenceVector(isnan(recurrenceVector)) = 0;

x = featureMatrix;
% y = [overallSurvival,censoringVector];
y = [overallSurvival, censoringVector];
cvfit = cvglmnet(x,y,'cox');
cvglmnetPlot(cvfit)

sum(cvfit.glmnet_fit.beta(:,cvfit.glmnet_fit.lambda==cvfit.lambda_min)~=0)


% cvglmnet(x,y,'cox');
% glmnetPlot(fit);
% glmnetPlot(fit, [],1);
% glmnetPrint(fit)

% pfit=glmnetPredict(fit,x(1,:),0.001713,'response')
% bk = 0;