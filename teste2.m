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

var = {'tos'}; %variaveis processadas
% var = {'tos', 'sos', 'thetao', 'so'}; %variaveis processadas
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
end

%%

obs = dadodetrend(:,1);

models = dadodetrend(:,2:end);

[alpha,n,crps,AR,STD,it] = ODR(obs,models,100,40000);

% 
% 
% alpha = calculate_alpha(0.5,0.2,7);
% n = calculate_n(100,alpha);
% 
% STD = 0.25;
% contSTD = 0;
% for i = 1:40000
%     
%     if i<=5000
%         tolerancia = 0.005;
%     elseif i>5000 & i<=10000
%         tolerancia = 0.0025;
%     elseif i>10000 & i<=20000
%         tolerancia = 0.001;
%     elseif i>20000 & i<=30000
%         tolerancia = 0.0005;
%     elseif i>30000
%         tolerancia = 0.0001;
%     end
%     
%     if i==1
%         [a,b] = calculate_ab(n,obs,models);
%         crps = calculate_crps(a,b);
%     elseif i==2
%         alpha(:,i) = calculate_alpha(alpha(:,i-1),STD,7);
%         alphaused = alpha;
%         n(:,i) = calculate_n(100,alpha(:,i));
%         [a,b] = calculate_ab(n(:,i),obs,models);
%         crps(:,i) = calculate_crps(a,b);
%     else
%         if mean(alpha(:,i-1))==0
%             alpha(:,i) = calculate_alpha(alphaused(:,end),STD,7);
%             n(:,i) = calculate_n(100,alpha(:,i));
%             [a,b] = calculate_ab(n(:,i),obs,models);
%             crps(:,i) = calculate_crps(a,b);
%             dif=crps(:,i)-crps(:,i-1);
%             if dif<=tolerancia;
%                 alpha(:,i) = alpha(:,i);
%                 alphaused(:,end+1) = alpha(:,i);
%                 n(:,i) = n(:,1);
%             else
%                 alpha(:,i) = 0;
%                 n(:,i) = 0;
%                 crps(:,i) = 0;
%             end
%         else
%             alpha(:,i) = calculate_alpha(alpha(:,i-1),STD,7);
%             n(:,i) = calculate_n(100,alpha(:,i));
%             [a,b] = calculate_ab(n(:,i),obs,models);
%             crps(:,i) = calculate_crps(a,b);
%             dif=crps(:,i)-crps(:,i-1);
%             if dif<=tolerancia;
%                 alpha(:,i) = alpha(:,i);
%                 alphaused(:,end+1) = alpha(:,i);
%                 n(:,i) = n(:,1);
%             else
%                 alpha(:,i) = 0;
%                 n(:,i) = 0;
%                 crps(:,i) = 0;
%             end
%         end
%         AR = size(alphaused,2)/size(alpha,2);
%         if AR<0.1
%             STD = STD/2
%             contSTD = contSTD+1
%             if STD<0.03125
%                 STD = 0.25;
%             end
%         end             
%     end
% end
% 
