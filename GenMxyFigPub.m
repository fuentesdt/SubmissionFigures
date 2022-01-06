% GenMxyFigPub Generates the figure for the show Mxy figure
% Author: Chris Walker
% Date: 8/6/2018

clear all
close all
clc

load('tmpShowMxyPub.mat')

%% Show Flip angles
figure('Units','inches','Position',[3,1,6.5,6.5])
font_size = 8;
titleStrings = {'Magnetization Evolution'};
y_labels = {'CEA-bb','VEA-const','VEA-max'};
subplot_labels = {'a)','b)','c)','d)','e)','f)'};
x_axis = (1:N)*TR-TR;
%% Display Flip Angle Schedules
for i = 1:numel(FAType)
    subplot(numel(FAType),2,2*i-1)
    plot(x_axis,save_flip_angles{i}*180/pi,'LineWidth',2)
    ylabel(sprintf('%s\n Excitation Angle\n (Degrees)',y_labels{i}),...
        'FontSize',6,'FontWeight','bold')
    ylim([0,100]),xlim([0,44]),grid('on')
    if i == 1
        title('Excitation Angle Schedule','FontSize',font_size,'FontWeight','bold')
        tmpLegend = legend('Pyr','Lac','location','north');
        tmpLegend.FontSize = font_size;
    end
    xlabel(' Time (sec) ','FontSize',font_size,'FontWeight','bold')
    set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches')
    % set(gca,'Position',[1.8,3.5-(i-1)*1,0.9,0.8])
    text(-1,110,subplot_labels{i},'FontSize',font_size,'FontWeight','bold')
    subplot(numel(FAType),2,2*i)
    plot(x_axis,save_Mz{i}(1,:),'--b',x_axis,save_Mxy{i}(1,:),'-b',...
        x_axis,save_Mz{i}(2,:),'--r',x_axis,save_Mxy{i}(2,:),'-r',...
        x_axis,sum(save_Mxy{i}),'k','LineWidth',2)
    if (i == 2)
        tmpLegend = legend('Pyr Mz','Pyr Mxy','Lac Mz','Lac Mxy','Total Mxy','location','northeast');
        tmpLegend.FontSize = font_size;
        tmpLegend.FontWeight = 'bold';
        set(tmpLegend,'position',[0.757,0.49,0.1474,0.1298])
    end
    ylim([0,0.01]),xlim([0,44]),grid('on')
    set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches')
    ylabel('Magnetization (arb)',...
        'FontSize',font_size,'FontWeight','bold')
    if (i == 1)
        title(titleStrings,'FontSize',font_size,'FontWeight','bold')
    end
    xlabel(' Time (sec) ','FontSize',font_size,'FontWeight','bold')
    text(-1,.011,subplot_labels{i+3},'FontSize',font_size,'FontWeight','bold')
end
% %% Display Magnetizations
% for j = 1:numel(FAType)
%     tmp_ax(i) = subplot(numel(FAType),5,5*i-4+j+1);
%     plot(x_axis,save_Mz{i,j}(1,:),'--b',x_axis,save_Mxy{i,j}(1,:),'-b',...
%         x_axis,save_Mz{i,j}(2,:),'--r',x_axis,save_Mxy{i,j}(2,:),'-r',...
%         x_axis,sum(save_Mxy{i,j}),':k','LineWidth',2)
%     if (j == 2 && i == 1)
%         tmpLegend = legend('Pyr Mz','Pyr Mxy','Lac Mz','Lac Mxy','Total Mxy','location','northeast');
%         tmpLegend.FontSize = font_size;
%         tmpLegend.FontWeight = 'bold';
%         set(tmpLegend,'units','inches','position',[5.45,3.6275,0.6,0.65])
%     end
%     if (j == 3 && i == 1)
%         set(get(gca,'Children'),'Visible','off');
%         set(gca,'Visible','off')
%     end
%     ylim([0,0.01]),xlim([0,25]),grid('on')
%     set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches')
%     if (j==1)
%         ylabel('Magnetization (arb)',...
%             'FontSize',font_size,'FontWeight','bold')
%     else
%         set(gca,'YTickLabel',[])
%     end
%     if (i == 1 && any(j == [1,2])) || (i == 2 && j == 3)
%         title(titleStrings{j},'FontSize',font_size,'FontWeight','bold')
%     end
%     if i == numel(FAType)
%         xlabel(' Excitation  # ','FontSize',font_size,'FontWeight','bold')
%     else
%         set(gca,'XTickLabel',[])
%     end
%     if i == 2 && j == 3
%         text(-2,.01275,subplot_labels{i,j+2},'FontSize',font_size,'FontWeight','bold')
%     else
%         text(-1,.011,subplot_labels{i,j+2},'FontSize',font_size,'FontWeight','bold')
%     end
% end
% for i = 1:numel(FAType)
%     for j = 1:numel(perffused)
%         set(tmp_ax(i,j),'Position',[2.3+(1*j),3.5-(i-1)*1,0.9,0.8])
%     end
% end
%     
%     
%     
