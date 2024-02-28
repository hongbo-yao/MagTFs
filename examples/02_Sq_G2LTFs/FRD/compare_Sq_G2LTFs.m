clc,clear,close all

% Sq_estimation
station_name='FRD';
filename1= [station_name '_Sq_G2LTFs.txt'];
fid1=fopen(filename1);
dataSq=textscan(fid1,' %f %f %f %f %f %f');
% MagTFs_BMT_Sq
coh2=load([station_name '_Sq_coh2.txt']);
Guzavina_Sq_G2LTFs=load(['Guzavina_' station_name '_Sq_G2LTF_data.txt']);

%BMT_Sq
Guzavina_Sq_G2LTFs_Re=Guzavina_Sq_G2LTFs(:,2);
Guzavina_Sq_G2LTFs_Im=Guzavina_Sq_G2LTFs(:,3);
Guzavina_Sq_G2LTFs_error=Guzavina_Sq_G2LTFs(:,4);

% compare
fig = figure('Position',[400 200 600 500]);
subplot(3,1,1);set(gca,'position',[0.12 0.83 0.82 0.12]);
plot([6,8,12,24],coh2(:,1),'ro','LineWidth',1.5);set(gca,'XScale','log');set(gca,'Xtick',[6 8 12 24]);
set(gca,'ylim',[0,1]);
set(gca,'xlim',[5.5,25]);
set(gca,'box','on');set(gca, 'FontSize', 15, 'LineWidth', 1.0);
ylabel('Coh2');title(station_name);

subplot(3,1,2);set(gca,'position',[0.12 0.49 0.82 0.25]);
Sq_estimation_period_in_hour=dataSq{1,3}/3600;
Sq_estimation_Sq_Re=dataSq{1,4};
Sq_estimation_Sq_error=dataSq{1,6};
errorbar(Sq_estimation_period_in_hour,Sq_estimation_Sq_Re,Sq_estimation_Sq_error,'ro','LineWidth',1.5)
hold on
Guzavina_period_in_hour=Sq_estimation_period_in_hour;
errorbar(Guzavina_period_in_hour,Guzavina_Sq_G2LTFs_Re,Guzavina_Sq_G2LTFs_error,'bd','LineWidth',1.5)
set(gca,'XScale','log')
ylabel('Real TFs');set(gca,'ylim',[-2,2]);
lgd = legend('MagTFs','Guzavina et al (2019)','NumColumns',2);
set(lgd,'location','best');set(gca, 'FontSize', 15, 'LineWidth', 1.0);
set(gca,'Xtick',[6 8 12 24]);set(gca,'xlim',[5.5,25]);

subplot(3,1,3);set(gca,'position',[0.12 0.15 0.82 0.25])
ylim=([-2 2]);
Sq_estimation_Sq_Im=dataSq{1,5};
errorbar(Sq_estimation_period_in_hour,Sq_estimation_Sq_Im,Sq_estimation_Sq_error,'ro','LineWidth',1.5)
hold on
errorbar(Guzavina_period_in_hour,Guzavina_Sq_G2LTFs_Im,Guzavina_Sq_G2LTFs_error,'bd','LineWidth',1.5)
set(gca,'XScale','log');set(gca, 'FontSize', 15, 'LineWidth', 1.0);
set(gca,'Xtick',[6 8 12 24]);set(gca,'xlim',[5.5,25]);
xlabel('Period (hours)')
ylabel('Imag TFs');set(gca,'ylim',[-2,2]);
filename = strcat(station_name, '_compare_Sq_G2LTFs');
print(filename, '-dpdf', '-r300');