clear all
close all
clc

load('CEAOptTmp')

figure('Units','inches','Position',[3,1,6.5,6.5])
font_size = 10;
titleStrings = {'Magnetization Evolution'};
y_labels = {'CEA-bb','VEA-const','VEA-max'};
subplot_labels = {'a)','b)','c)','d)','e)'};


subplot(2,3,2)
plot(excitation_angle_vect,abs(squeeze(std(fits,[],3)))/0.1*100,...
    'LineWidth',2)
xlabel('Flip Angle (Degrees)','FontSize',font_size,'FontWeight','bold')
ylabel('k_p_l Coeff. Varr. (%)','FontSize',font_size,'FontWeight','bold')
title('k_p_l Coeff. Varr.','FontSize',font_size,'FontWeight','bold')
grid on
set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches',...
    'Position',[3.7, 3.7949, 2.2176, 2.2176])
text(-5,12.75,subplot_labels{2},'FontSize',font_size,'FontWeight','bold')
subplot(2,3,1)
plot(excitation_angle_vect,abs(squeeze(mean(fits,3))-0.1)/0.1*100,...
    'LineWidth',2)
xlabel('Flip Angle (Degrees)','FontSize',font_size,'FontWeight','bold')
ylabel('k_p_l Error (%)','FontSize',font_size,'FontWeight','bold')
title('k_p_l Error','FontSize',font_size,'FontWeight','bold')
grid on
set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches',...
    'Position',[0.845, 3.7949, 2.2176, 2.2176])
text(-5,2.15,subplot_labels{1},'FontSize',font_size,'FontWeight','bold')
subplot(2,3,6)
plot(excitation_angle_vect,squeeze(mean(signals,3)),...
    'LineWidth',2)
xlabel('Flip Angle (Degrees)','FontSize',font_size,'FontWeight','bold')
title('Total Signal','FontSize',font_size,'FontWeight','bold')
grid on
set(gca,'FontSize',font_size,'FontWeight','bold')
text(-5,0.0735,subplot_labels{5},'FontSize',font_size,'FontWeight','bold')
subplot(2,3,4)
plot(excitation_angle_vect,squeeze(mean(pyrSig,3)),...
    'LineWidth',2)
xlabel('Flip Angle (Degrees)','FontSize',font_size,'FontWeight','bold')
ylabel('Signal (AU)','FontSize',font_size,'FontWeight','bold')
title('Pyruvate Signal','FontSize',font_size,'FontWeight','bold')
grid on
set(gca,'FontSize',font_size,'FontWeight','bold')
text(-5,0.053,subplot_labels{3},'FontSize',font_size,'FontWeight','bold')
subplot(2,3,5)
plot(excitation_angle_vect,squeeze(mean(lacSig,3)),...
    'LineWidth',2)
xlabel('Flip Angle (Degrees)','FontSize',font_size,'FontWeight','bold')
title('Lactate Signal','FontSize',font_size,'FontWeight','bold')
grid on
set(gca,'FontSize',font_size,'FontWeight','bold')
text(-5,0.0362,subplot_labels{4},'FontSize',font_size,'FontWeight','bold')