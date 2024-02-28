clc,clear,close all

% Q_estimation
TF_file = 'satellite.TF';
[period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult] = ...
    textread(TF_file,'%d %d %d %f %f %f %f %f %f','headerlines',1,'delimiter',' ');

dataQ = [period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult];

% KV2021
filename2 = 'Kuvshinov_2021_EPS_global_Q_responses.dat';
datakv=load(filename2);

% compare
fig = figure('Position',[400 200 600 500]);
subplot(2,1,1);set(gca,'position',[0.13 0.8 0.83 0.15]);
scatter(dataQ(:,4)/86400,dataQ(:,8),'ro','LineWidth',1.5);
hold on
scatter(datakv(:,1),datakv(:,5),'bd','LineWidth',1.5);
set(gca,'XScale','log');set(gca,'box','on');
ylabel('coh2')
xlim([1 160]);ylim([0.9 1]);set(gca, 'FontSize', 15, 'LineWidth', 1.0);

subplot(2,1,2);set(gca,'position',[0.13 0.15 0.83 0.55]);
Q_estimation_period_in_day=dataQ(:,4)/86400;
Q_estimation_C_Re=dataQ(:,5);
Q_estimation_C_error=dataQ(:,7);
errorbar(Q_estimation_period_in_day,Q_estimation_C_Re,Q_estimation_C_error,'ro','LineWidth',1.5)
hold on
KV_period_in_day=datakv(:,1);
KV_C_Re=datakv(:,2);
KV_C_error=datakv(:,4);
errorbar(KV_period_in_day,KV_C_Re,KV_C_error,'bd','LineWidth',1.5)
hold on
Q_estimation_C_Im=dataQ(:,6);
errorbar(Q_estimation_period_in_day,Q_estimation_C_Im,Q_estimation_C_error,'ro','LineWidth',1.5)
hold on
KV_C_Im=-datakv(:,3);
errorbar(KV_period_in_day,KV_C_Im,KV_C_error,'bd','LineWidth',1.5)
set(gca,'XScale','log')
xlabel('Period (days)')
ylabel('Q-response')
lgd = legend('MagTFs','Kuvshinov et al (2021)','NumColumns',2);
% set(lgd,'box','off','location','best')
set(gca, 'FontSize', 15, 'LineWidth', 1.0);
xlim([1 160]);
ylim([0 0.5]);
filename = strcat('Q_responses');
print(filename, '-dpdf', '-r300');