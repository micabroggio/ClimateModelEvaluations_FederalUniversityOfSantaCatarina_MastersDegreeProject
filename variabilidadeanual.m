% Micael Fernando Broggio
% Criacao: 21 de janeiro de 2020
%
% FIGURA VARIABILIDADE ANUAL
%
% Para dados .nc.
%
% Scrip desenvolvido para analisar a variabilidade anual dos modelos CMIP6
% em relacao ao modelo BESM
% 
% Dados de input ja estao em medias anuais
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')

% var = {'thetao'}; %variaveis processadas
var = {'tos', 'sos', 'thetao', 'so'}; %variaveis processadas
labelvar = {'SST - temperature (°C)','SSS - salinity','\Theta - temperature (°C)','SO - salinity'};

dirin = '/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/DATA/'; %diretorio onde sao encontrados os dados anuais dos modelos
reg = 'AS'; %escala do processamento (AS ou GLOBAL)

%prepara os nomes dos modelos e referencia para que sejam carregados os
%dados---------------------------------------------------------------------
modelos = dir([dirin 'MODELOS/']); modelos = modelos(3:end); modelos = {modelos.name}; %encontra os modelos usados a partir dos nomes das pastas onde os dados estao armazenados
posbesm = find(strcmp(modelos,'BESM-OA2-5')); %encontra a posicao da nomenclatura do modelo BESM
posama = find(strcmp(modelos,'AMA')); %encontra a posicao da nomenclatura do modelo AMA
possma = find(strcmp(modelos,'SMA')); %encontra a posicao da nomenclatura do modelo SMA
modelos([posbesm posama possma])={'1' '2' '3'}; modelos = sort(modelos); %substitui o nome dos tres modelos acima para os numeros 1 2 3 para forçar a ordem de processamento desejada
modelos(1) = {'BESM-OA2-5'}; modelos(2) = {'AMA'}; modelos(3) = {'SMA'}; %retorna os nomes dos modelos aos seus respectivos valores
ref = dir([dirin 'REANALISES/']); ref = ref(3:end); ref = {ref.name}; %encontra os modelos usados a partir dos nomes das pastas onde os dados estao armazenados
%--------------------------------------------------------------------------

f1 = figure
set(f1, 'Position', [1 1 700 900]); %set o tamanho da figura
set(gcf,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
set(gcf,'color','w'); %figura emfundo branco

for i = 1:size(var,2)
    disp(var{i}) %dispara a info de qual variavel e procesada
    nm = 1;, dadomod = []; %prealoca o contador de modelos e a variavel dado
    for j = modelos %loop para todos os modelos
        dirmod = [dirin 'MODELOS/' char(j) '/' reg '/1979-2014/']; %diretorio de cada modelo
        dadosdisp = dir([dirmod  '*' var{i} '*Omon*.nc']); dadosdisp = dadosdisp.name; %encontra o nome de cada arquivo .nc associado a variavel analisada
        levind = ncinfo([dirmod dadosdisp]); levind  = {levind.Dimensions.Name}; levind = cell2mat(strfind(levind,'lev'));
        
        timemod = ncread([dirmod dadosdisp],'time'); %carrega o netcdf
        time2 = datetime(1979,01,01):calmonths:datetime(2014,12,01); %define um vetor tempo feito a mao
        
        %carrega os dados dos modelos--------------------------------------
        for k = 1:12:length(timemod) %loop para variacao temporal ao passo 12
            if isempty(levind) %condicional de nivel supercial
                dadomod(k:k+11,nm) = double(squeeze(nanmean(reshape(ncread([dirmod dadosdisp], var{i}, [1 1 k], [inf inf 12], [1 1 1]),1,[],12),2))); %carrega o netcdf de 12 em 12 ja efetuando o calculo de media espacial
            else %condicional de varios niveis
                dadomod(k:k+11,nm) = double(squeeze(nanmean(reshape(ncread([dirmod dadosdisp], var{i}, [1 1 1 k], [inf inf inf 12], [1 1 1 1]),1,1,[],12),3))); %carrega o netcdf de 12 em 12 ja efetuando o calculo de media espacial
            end         
        end
        %------------------------------------------------------------------
        
        nm = nm + 1; %contador de modelo 
    end
    
    %carrega os dados de reanalise-----------------------------------------
    dirref  = [dirin 'REANALISES/' ref{1} '/' reg '/']; %diretorio de cada modelo
    dadosdisp = dir([dirref  '*' var{i} '*.nc']); dadosdisp = dadosdisp.name; %encontra o nome no arquivo nc que sera processado
    if isempty(levind) %condicional de nivel supercial
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf do dado de referencia
        dadoref = squeeze(nanmean(reshape(dadoref,1,[],size(dadoref,3)),2)); %media espacial para o dado ref
        dado = cat(2,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena o dado ref aos dados dos modelos
    else %condicional de varios niveis
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf do dado de referencia
        dadoref = squeeze(nanmean(reshape(dadoref,1,1,[],size(dadoref,4)),3)); %media espacial para o dado ref
        dado = cat(2,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena o dado ref aos dados dos modelos
    end
    %----------------------------------------------------------------------
    
    dado = movmean(dado,12); dado = dado(7:end-6,:); time2 = time2(:,7:end-6); %realiza a media movel com janela de 12 elementos
    dadodetrend = detrend(dado); %retira a tendencia dos dados
    dadonorm = dadodetrend - dadodetrend(1,:); %normaliza todos os modelos com o dado inicial em zero
    clear dadoref dadomod %elimina dos dados ref e mod para aliviar a memoria
    
    [rho1,pval1] = corr(dado); %calculo de correlacao aplicada aos dados das series temporais
    [rho2,pval2] = corr(dadonorm); %calculo de correlacao aplicada aos dados das series normalizadas

    cor = {'k','b','r','c',[1 .5 0],[0 .5 0],'m',[.3 .3 .3]}; %cor para cada modelo

    %PLOTA AS SERIES TEMPORAIS ANUAIS DOS MODELOS--------------------------
%     f1 = figure
%     set(f1, 'Position', [1 1 700 300]);
%     plot(time2,dado(:,3),'color',cor{3});
%     hold on
%     plot(time2,dado(:,4),'color',cor{4});
%     plot(time2,dado(:,5),'color',cor{5});
%     plot(time2,dado(:,6),'color',cor{6});
%     plot(time2,dado(:,7),'color',cor{7});
%     plot(time2,dado(:,8),'color',cor{8});
%     plot(time2,dado(:,2),'color',cor{2});
%     plot(time2,dado(:,1),'color',cor{1},'linewidth',1.5);
%     grid
%     xlabel('time (year)')
%     ylabel(labelvar{i})
%     print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/serietemporaldeyear' reg var{i}],'-r1000')
    %----------------------------------------------------------------------
    
    %PLOTA A VARIABILIDADE DOS DADOS SEM TENDENCIA-------------------------
    sp = subplot(4,1,i); %configura subplots em linhas
    sp.Position = sp.Position + [-.01 0 .04 .04]; %expande o tamanho do subplot
    p3=plot(time2,dadonorm(:,3),'color',cor{3}); %cada plot pn é um modelos
    hold on
    p4=plot(time2,dadonorm(:,4),'color',cor{4});
    p5=plot(time2,dadonorm(:,5),'color',cor{5});
    p6=plot(time2,dadonorm(:,6),'color',cor{6});
    p7=plot(time2,dadonorm(:,7),'color',cor{7});
    p8=plot(time2,dadonorm(:,8),'color',cor{8});
    p2=plot(time2,dadonorm(:,2),'color',cor{2});
    p1=plot(time2,dadonorm(:,1),'color',cor{1},'linewidth',1.5);
    grid % aplica a grade
    if i == 4 %condicional para o ultimo subplot
        xlabel('time (year)') %set de xlabel
        lx = legend([p1,p2,p3,p4,p5,p6,p7,p8],ndado,'Location','South','Orientation','horizontal','NumColumns',4); %aplica legenda e seta a legenda
        lx.Position = lx.Position + [0 -.09 0 0]; %reposiciona a legenda
    else %condicional para os demais plots
        xticklabels([]); %retira o xticklabel
    end
    ylabel(labelvar{i}) %set dos ylabels
    %----------------------------------------------------------------------
end
% print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/variabilidadeyear_historico_' reg],'-r1000') %salva a figura em png e 1000ppi de resolucao
