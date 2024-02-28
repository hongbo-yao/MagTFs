clc,clear,close all

station_name = 'IRT';
TF_file = [station_name '.TF'];
[period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult] = ...
    textread(TF_file,'%d %d %d %f %f %f %f %f %f','headerlines',1,'delimiter',' ');

TF = [period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult];
Alexey_2011_PHD_C=load(['Semenov2011_' station_name '_C_responses.txt']);

theta_GM = load([station_name '_theta_GM.txt']);
a = 6371.2;
k = -a/2.0*tand(theta_GM);
C = k*(TF(:,5)+1j*TF(:,6));
C_error = abs(k*TF(:,7));
C_Re = real(C);
C_Im = imag(C);
coh2 = TF(:,8);
T = TF(:,4);
period_in_day = T/86400;

% Alexey_2011_PHD_C_FRD
Alexey_2011_PHD_C_Re = Alexey_2011_PHD_C(:,2);
Alexey_2011_PHD_C_Im = Alexey_2011_PHD_C(:,3);
Alexey_2011_PHD_error = Alexey_2011_PHD_C(:,4);
Alexey_2011_PHD_cho2=Alexey_2011_PHD_C(:,5);

fig = figure('Position',[400 200 600 500]);
subplot 211
semilogx(period_in_day,coh2,'ro','LineWidth',1.5)
hold on
semilogx(period_in_day,Alexey_2011_PHD_cho2,'bd','LineWidth',1.5)
ylabel('cho2')
ylim([0 1])
xlim([2.5 115]);title(station_name);
set(gca,'position',[0.14 0.8 0.84 0.15]);set(gca, 'FontSize', 15, 'LineWidth', 1.5);

subplot 212
errorbar(period_in_day,C_Re,C_error,'ro','LineWidth',1.5)
hold on
errorbar(period_in_day,Alexey_2011_PHD_C_Re,Alexey_2011_PHD_error,'bd','LineWidth',1.5)
hold on
errorbar(period_in_day,C_Im,C_error,'ro','LineWidth',1.5)
hold on
errorbar(period_in_day,Alexey_2011_PHD_C_Im,Alexey_2011_PHD_error,'bd','LineWidth',1.5)
set(gca,'XScale','log')
xlabel('Period (days)')
ylabel('C-responses (km)')
lgd = legend('MagTFs','Semenov (2011)','NumColumns',2);set(lgd,'location','best')
xlim([2.5 115])
set(gca,'position',[0.14 0.15 0.84 0.55]);set(gca, 'FontSize', 15, 'LineWidth', 1.5);

filename = [station_name '_observed_C_responses'];
print(filename, '-dpdf', '-r300');


fid = fopen([station_name '_observed_C_responses.txt'],'w');
for i=1:length(T)
    fprintf(fid,'%10.0f\t%10.2f\t%10.2f\t%10.2f\n',T(i),real(C(i)),imag(C(i)),C_error(i));
end
fclose(fid);

fid = fopen([station_name '_observed_C_responses_coh2.txt'],'w');
for i=1:length(period_in_day)
    fprintf(fid,'%8.2f\t%10.2f\t%10.2f\t%10.2f\t%10.2f\n',period_in_day(i),real(C(i)),imag(C(i)),C_error(i),coh2(i));
end
fclose(fid);