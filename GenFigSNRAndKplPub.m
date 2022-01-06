clear all
close all
clc

% Global Variables
font_size = 8;
line_width = 2;
figure('units','inches','Position',[6,2,6.5,3])
files = 'tmpPerffusedPub2';
title_text = {'Pyruvate Signal (AU)','Lactate Signal (AU)',...
    'Total Signal (AU)'};
subplot_lables = {'a)','b)','c)'};
legend_labels = {'CEA-bb','VEA-const','VEA-max'};
load(files)
for row = 1:3
    subplot(1,3,row)
    switch row
        case 1
            plot(N_vect,mean(pyrSig(:,:,:),3),'LineWidth',line_width)
            ylim([0,0.07])
            title(title_text{row},...
                'FontSize',font_size,'FontWeight','bold')
        case 2
            plot(N_vect,mean(lacSig(:,:,:),3),'LineWidth',line_width)
            ylim([0,0.07])
            title(title_text{row},...
                'FontSize',font_size,'FontWeight','bold')
        case 3
            plot(N_vect,mean(signals(:,:,:),3),'LineWidth',line_width)
            ylim([0,0.07])
            title(title_text{row},...
                'FontSize',font_size,'FontWeight','bold')
    end
    xlabel(' Number of Excitations ',...
        'FontSize',font_size,'FontWeight','bold')
    grid on
    set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches')
    xlim([5,60])
    if row == 1
        tmpLegend = legend(legend_labels,'Location','northeast');
        tmpLegend.FontSize = font_size;
        tmpLegend.FontWeight = 'bold';
        ylabel('Signal (AU)','FontSize',font_size,'FontWeight','bold')
    else
        set(gca,'YTickLabel',[])
    end
    text(3,0.0735,subplot_lables{row},'FontSize',...
        font_size,'FontWeight','bold')
    set(gca,'Position',[0.6+(row-1)*2,0.4,1.7,2.3])
end
%% KPL Plot
% get line colors
 tmp = get(gca,'Children');
 tmp_inter = flip(2:numel(tmp));
for i = 1:numel(tmp_inter)
    line_color(i,:) = tmp(tmp_inter(i)).Color;
end
figure('units','inches','Position',[6,2,6.5,3])
title_text = legend_labels;
for row = 1:numel(FAType)
    tmp_ax(row) = subplot(1,3,row);
    plot(N_vect,zeros(size(N_vect))+Kpl,'g','LineWidth',line_width)
    hold on
    plot(N_vect,mean(fits(:,row,:),3),'color',line_color(row,:),'LineWidth',line_width)
    plot(N_vect,mean(fits(:,row,:),3)+std(fits(:,row,:)*1.96,[],3),...
        'lineStyle','--','color',line_color(row,:),'LineWidth',line_width)
        tmpLegend = legend({'True','Mean','95% C.I.'},'Location','southeast');
        tmpLegend.FontSize = font_size;
        tmpLegend.FontWeight = 'bold';
        tmpLegend.Units = 'inches';
        %tmpLegend.Position = [8.5,6.4-(row-1)*1.75,1.25,.75];
    plot(N_vect,mean(fits(:,row,:),3)-std(fits(:,row,:)*1.96,[],3),...
        'lineStyle','--','color',line_color(row,:),'LineWidth',line_width)
    ylim([0.075,0.125]),xlim([5,60])
    grid on
    set(gca,'FontSize',font_size,'FontWeight','bold','Units','inches')
    title(title_text{row},'FontSize',font_size,'FontWeight','bold')
    if row == 1
        ylabel('Fit k_P_L (sec^-^1)',...
            'FontSize',font_size,'FontWeight','bold')
    else
        set(gca,'YTickLabel',[])
    end
    xlabel(' No. Excitations ',...
        'FontSize',font_size,'FontWeight','bold')
    text(-0.25,.1275,subplot_lables{row},'FontSize',font_size,'FontWeight','bold')
    set(gca,'Position',[0.6+(row-1)*2,0.4,1.7,2.3])
end
for row = 1:numel(FAType)
   
end
