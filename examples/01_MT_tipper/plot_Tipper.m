clc,clear,close all

% MagTF
station_name = 'CKI';
TF_file = [station_name '.TF'];
[period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult] = ...
    textread(TF_file,'%d %d %d %f %f %f %f %f %f','headerlines',1,'delimiter',' ');

TF = [period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult];

% Tx
input_channel_id = 1;
loc = find(TF(:,3) == input_channel_id);
Tx = TF(loc,:);
Tx(:,6) = -Tx(:,6);
[n_data,~] = size(Tx);

% Ty
input_channel_id = 2;
loc = find(TF(:,3) == input_channel_id);
Ty = TF(loc,:);
Ty(:,6) = -Ty(:,6);

% coh_mult
loc = find(TF(:,3) == input_channel_id);
coh_mult = TF(loc,:);

% Rigaud et al (2021)
dat = load('Rigaud_etal_CKI2017Tipper.mat');

% Tx
Tx_re = reshape(real(dat.TF(:,1,:)),1,[]);
Tx_Im = reshape(imag(dat.TF(:,1,:)),1,[]);
Tx_err =reshape(real(dat.lim_TF(:,1,:)),1,[]);

% Ty
Ty_re = reshape(real(dat.TF(:,2,:)),1,[]);
Ty_Im = reshape(imag(dat.TF(:,2,:)),1,[]);
Ty_err =reshape(real(dat.lim_TF(:,2,:)),1,[]);

% coh_mult
Tx_coh=reshape(dat.coh2(:,1,:),1,[]);
Ty_coh=reshape(dat.coh2(:,2,:),1,[]);
% plot
fig = figure('Position',[400 200 600 500]);

subplot(3,1,1);
plot(Tx(:,4), Tx(:,8), 'ro','LineWidth', 1.5);
hold on;
plot(Ty(:,4), Ty(:,8), 'rd','LineWidth', 1.5);
hold on
plot(dat.T, Tx_coh, 'bo','LineWidth', 1.5);
hold on
plot(dat.T, Ty_coh, 'bd','LineWidth', 1.5);
ylabel('coh2');
xlim([0.8*Tx(1,4), 1.2 * Tx(n_data,4)]);
ylim([0,1]);
set(gca, 'XScale', 'log');set(gca,'position',[0.13 0.83 0.83 0.15]);
title(station_name);
set(gca, 'FontSize', 15, 'LineWidth', 1.5);

subplot(3,1,2);
h1=errorbar(Tx(:,4), Tx(:,5), Tx(:,7), 'ro', 'LineWidth', 1.5);
hold on;
h2=errorbar(Tx(:,4), Tx(:,6), Tx(:,7), 'rd', 'LineWidth', 1.5);
hold on
h3=errorbar(dat.T, Tx_re, Tx_err, 'bo', 'LineWidth', 1.5);
hold on
h4=errorbar(dat.T, Tx_Im, Tx_err, 'bd', 'LineWidth', 1.5);
ylabel('T_{zx}');
xlim([0.8*Tx(1,4), 1.2 * Tx(n_data,4)]);
ylim([-0.5,1.5]);
set(gca, 'XScale', 'log');set(gca,'position',[0.13 0.48 0.83 0.25]);
lgd=legend('MagTFs', 'MagTFs', 'Rigaud et al (2021)', 'Rigaud et al (2021)','NumColumns',2);
set(gca, 'FontSize', 15, 'LineWidth', 1.5);

subplot(3,1,3);
errorbar(Ty(:,4), Ty(:,5), Ty(:,7), 'ro', 'LineWidth', 1.5);
hold on;
errorbar(Ty(:,4), Ty(:,6), Ty(:,7), 'rd', 'LineWidth', 1.5);
errorbar(dat.T, Ty_re, Ty_err, 'bo', 'LineWidth', 1.5);
errorbar(dat.T, Ty_Im, Ty_err, 'bd', 'LineWidth', 1.5);
xlabel('Period (seconds)');
ylabel('T_{zy}');
xlim([0.8*Tx(1,4), 1.2 * Tx(n_data,4)]);
ylim([-0.25,0.5]);set(gca,'Ytick',[-0.25 0 0.25 0.5]);
set(gca, 'XScale', 'log');set(gca,'position',[0.13 0.15 0.83 0.25]);
set(gca, 'FontSize', 15, 'LineWidth', 1.5);

filename = strcat(station_name, '_tipper');
print(filename, '-dpdf', '-r300');
