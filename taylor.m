% Micael Fernando Broggio
% Criacao: 5 de fevereiro de 2020
%
% FIGURA DIAGRAMA DE TAYLOR
%
% Para dados .nc.
%
% Scrip desenvolvido para analisar a variabilidade anual dos modelos CMIP6
% em relacao ao modelo BESM.
% 
% O intuito e a comparacao de variaveis estatisticas do modelo BESM em
% relacao aos demais modelos
% 
% Dados de input ja estao em medias mensais
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/SkillMetricsToolbox-master')
addpath('~/Documents/MATLAB/taylordiagram')

% var = {'tos','sos'}; %variaveis processadas
var = {'tos', 'sos', 'thetao', 'so'}; %variaveis processadas
labelvar = {'SST','SSS','\Theta','SO'};

dirin = '/home/micael/Documents/PROJETOS/ModCOSTA/DATA/CMIP6/'; %diretorio onde sao encontrados os dados anuais dos modelos
reg = 'AS'; %escala do processamento (AS ou GLOBAL)

%prepara os nomes dos modelos e referencia para que sejam carregados os
%dados---------------------------------------------------------------------
modelos = dir([dirin]); modelos = modelos(3:end); modelos = {modelos.name}; %encontra os modelos usados a partir dos nomes das pastas onde os dados estao armazenados
posbesm = find(strcmp(modelos,'BESM-OA2-5')); %encontra a posicao da nomenclatura do modelo BESM
posama = find(strcmp(modelos,'AMA')); %encontra a posicao da nomenclatura do modelo AMA
possma = find(strcmp(modelos,'SMA')); %encontra a posicao da nomenclatura do modelo SMA
modelos([posbesm posama possma])={'1' '2' '3'}; modelos = sort(modelos); %substitui o nome dos tres modelos acima para os numeros 1 2 3 para forçar a ordem de processamento desejada
modelos(1) = {'BESM-OA2-5'}; modelos(2) = {'AMA'}; modelos(3) = {'SMA'}; %retorna os nomes dos modelos aos seus respectivos valores
ref = dir([dirin 'REANALISES/']); ref = ref(3:end); ref = {ref.name}; %encontra os modelos usados a partir dos nomes das pastas onde os dados estao armazenados
%--------------------------------------------------------------------------

f1 = figure;
set(f1, 'Position', [1 1 1000 500]);

for i = 1:size(var,2)
    disp(var)
    nm = 1; dadomod = []; %prealoca o contador de modelos e a variavel dado
    for j = modelos %loop para todos os modelos
        dirmod = [dirin char(j) '/' reg '/1979-2014/']; %diretorio de cada modelo
        dadosdisp = dir([dirmod  '*' var{i} '*Omon*.nc']); dadosdisp = dadosdisp.name; %encontra o nome de cada arquivo .nc associado a variavel analisada
        levind = ncinfo([dirmod dadosdisp]); levind  = {levind.Dimensions.Name}; levind = cell2mat(strfind(levind,'lev'));
        
        %carrega os dados dos modelos--------------------------------------
        if isempty(levind)
            dadomod(:,:,:,nm) = ncread([dirmod dadosdisp], var{i} ); %carrega o netcdf
        else
            dadomod(:,:,:,:,nm) = ncread([dirmod dadosdisp], var{i} ); %carrega o netcdf
        end
        %------------------------------------------------------------------
        
        timemod = ncread([dirmod dadosdisp],'time'); %carrega o netcdf
        time2 = datetime(1979,01,01):calmonths:datetime(2014,12,01);
        nm = nm + 1; %contador de modelo 
    end
    
    %carrega os dados de reanalise-----------------------------------------
    dirref  = [dirin 'REANALISES/' ref{1} '/' reg '/']; %diretorio de cada modelo
    dadosdisp = dir([dirref  '*' var{i} '*.nc']); dadosdisp = dadosdisp.name;
    if isempty(levind)
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf
        dado = cat(4,dadoref,dadomod); ndado = cat(2,ref,modelos);
    else
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf
        dado = cat(5,dadoref,dadomod); ndado = cat(2,ref,modelos);
    end
    %----------------------------------------------------------------------
    
    dado = squeeze(nanmean(reshape(dado,1,1,[],size(timemod,1),size(modelos,2)+1),3)); %realiza a media territorial para cada timestep
    dado = movmean(dado,12); dado = dado(7:end-6,:); time2 = time2(:,7:end-6);
%     dado = squeeze(nanmean(reshape(dado,12,[],size(ndado,2)),1));
    clear dadoref dadomod
    
    [rho,pval] = corr(dado);
    rho = rho(1,:);
    erro = rmse2(dado(:,1),dado,1);
    desv = nanstd(dado,0,1);
    desvn = desv*2/desv(1);

    cor = {'k','b','r','c',[1 .5 0],[0 .5 0],'m',[.3 .3 .3]}; %cor para cada modelo
    label = {'Non-Dimensional Observation' modelos{:}};
    
    sp = subplot(2,2,i);
    sp.Position = sp.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot
    [A1, ht1, axl1] = taylor_diagram(desvn([1 2]),erro([1 2]),rho([1 2]), ...
        'numberPanels', 2, ...
        'markerSize',20, 'markerColor', cor{2}, ...
        'tickRMS',1:1:3,'styleRMS','--','titleRMS', 'off','colRMS','k', 'tickRMSangle', 140, ...
        'tickSTD',0:1:3,'limSTD',3,'colSTD','k', 'styleSTD', '-','widthSTD', 1.0, ...
        'tickCOR',[-1.0 -0.99 -0.95 -0.9 -0.8:0.2:-0.2 0:0.2:0.8 0.9 0.95 0.99 1.0],'titleCOR', 'off', 'colCOR','k', 'styleCOR', ':', 'widthCOR', 1.0, ...
        'styleOBS','-','colOBS','r','titleOBS',[labelvar{i} ' ORAS5'],'markerOBS','.','widthOBS',1.5);
    for j = 3:length(ndado)
        [B] = taylor_diagram(desvn([1 j]),erro([1 j]),rho([1 j]), ...
            'overlay','on', ...
            'markerSize',20, 'markerColor',cor{j});
        eval(['A' num2str(j-1) '= B']);
    end
    [B] = taylor_diagram(desvn([1 1]),erro([1 1]),rho([1 1]), ...
        'overlay','on', ...
        'markerSize',20, 'markerColor','k');
    if i==1
        legend([B,A1,A2,A3,A4,A5,A6,A7],ndado, 'Location', 'southwest')
    end
    set(gcf,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
    set(gcf,'color','w'); %figura emfundo branco
    %----------------------------------------------------------------------
%     close(f1)
end
% print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/taylor_historico_' reg],'-r1000')

