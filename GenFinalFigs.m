%% Generates all the final figures for the publication
clear all
close all
clc

resim_data = false; % this flag will re simulate all the data
if(resim_data)
    ShowMxyPub % Generates data to show signal evolution (Fig 1)
    save('tmpShowMxyPub')
    ShowMxyPubNoise % Generates data to show noisey signal evolution (S-Fig 1)
    save('tmpShowMxyPubNoise')
    NDependencePub % Generates data to show total signal and Kpl fit (Fig 2 & 3)
    save('tmpPerffusedPub2');
    NDependencePubPerf % Generates data to show Kpl fit with perfusion (Fig 4)
    save('tmpPerffusedPubPerf2');
    ConstantFlipAngleDetermination % Generates data to optimal CEA-bb (Fig S2)
    save('CEAOptTmp');
end
GenMxyFigPub
savefig(gcf,'Fig1')
GenMxyFigPubNoise
savefig(gcf,'FigS2')
GenFigSNRAndKplPub
savefig(gcf,'Fig3')
close
savefig(gcf,'Fig2')
GenFigSNRAndKplPerfPub
close, close, close
savefig(gcf,'Fig4')
GenFigS2
savefig(gcf,'FigS1')

StatisticalTesting