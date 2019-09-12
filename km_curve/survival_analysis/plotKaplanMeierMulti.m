function info = plotKaplanMeierMulti(TTR,labelmask,predlabels)
% TTR - 1xn vector of months from treatment to last follow up
% labelmask - 1xn vector of outcomes (ie. -1, 0, 1)
% predlabels - 1xn vector of predicted labels (ie. -1, 0, 1)

X(:,1) = TTR;
X(:,2) = labelmask;

figure;
info = logrankMulti(X(predlabels == -1, :), X(predlabels == 0, :), X(predlabels == 1,:));
ylabel('Survival Rate','FontSize',14)
xlabel('Time (days)','FontSize',14)
title('Kaplan-Meier estimate of survival functions','FontSize',14)
set(gca,'FontSize',11)
%legend('Short Term','Long Term','Censored')
legend('Short Term','Mid Term','Long Term','Censored')


Y = X;
Z = predlabels';
Y(Z == -1, :) =8;
Z(Z == -1, :) =8;
Y(:, 2) = Z;
figure;
info = logrank_mod(Y(Z == 0,:),Y(Z == 1,:));
p1=info.p;
string = strcat(['p-value: ' num2str(info.p)]);
textbp(string, 'FontSize', 11)
legend('Mid Term','Long Term')
fig = gcf;
ax = get(fig,'children');
h = get(ax,'children');
set(h{2}(3),'Color','g')

Y = X;
Z = predlabels';
Y(Z == 1, :) = 8;
Z(Z == 1, :) = 8;
Z = Z + 1;
Y(:, 2) = Z;
figure;
info = logrank_mod(Y(Z == 0,:),Y(Z == 1,:));
p2=info.p;
string = strcat(['p-value: ' num2str(info.p)]);
textbp(string, 'FontSize', 11)
legend('Short Term','Mid Term')
fig = gcf;
ax = get(fig,'children');
h = get(ax,'children');
set(h{2}(2),'Color','g')

Y = X;
Z = predlabels';
Y(Z == 0, :) = 8;
Z(Z == 0, :) = 8;
Z = (Z + 1) / 2;
Y(:, 2) = Z;
figure;
info = logrank_mod(Y(Z == 0,:),Y(Z == 1,:));
p3=info.p;
string = strcat(['p-value: ' num2str(info.p)]);
textbp(string, 'FontSize', 11)
legend('Short Term','Long Term')
p=[p1,p2,p3];
disp(p);
end
