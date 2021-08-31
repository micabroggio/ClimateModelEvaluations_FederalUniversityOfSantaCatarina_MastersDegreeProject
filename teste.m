clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/nat/')

dirin = '/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/DATA/'; %diretorio onde sao encontrados os dados anuais dos modelos
reg = 'AS'; %escala do processamento (AS ou GLOBAL)

var = {'tos'}; %variaveis processadas
labelvar = {'SST'}; %nome dada a variavel no plot

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
    nm = 1; dadomod = []; %prealoca o contador de modelos e a variavel dado
    for j = modelos %loop para todos os modelos
        dirmod = [dirin 'MODELOS/' char(j) '/' reg '/1979-2014/']; %diretorio de cada modelo
        dadosdisp = dir([dirmod  '*' var{i} '*Omon*.nc']); dadosdisp = dadosdisp.name; %encontra o nome de cada arquivo .nc associado a variavel analisada
        levind = ncinfo([dirmod dadosdisp]); levind  = {levind.Dimensions.Name}; levind = cell2mat(strfind(levind,'lev')); %encontra de ha presenca de mais niveis ou so superficie
        
        %carrega os dados dos modelos--------------------------------------
        if isempty(levind) %dados superficiais (sem niveis)
            dadomod(:,:,:,nm) = ncread([dirmod dadosdisp], var{i} ); %carrega o netcdf
        else %dados subsuperficiais (com niveis)
            dadomod(:,:,:,:,nm) = ncread([dirmod dadosdisp], var{i} ); %carrega o netcdf
        end
        %------------------------------------------------------------------
        
        timemod = ncread([dirmod dadosdisp],'time'); %carrega o netcdf
        time2 = datetime(1979,01,01):calmonths:datetime(2014,12,01); %cria o vetor tempo na mao
        nm = nm + 1; %contador de modelo 
    end
    
    %carrega os dados de reanalise-----------------------------------------
    dirref  = [dirin 'REANALISES/' ref{1} '/' reg '/']; %diretorio de cada modelo
    dadosdisp = dir([dirref  '*' var{i} '*.nc']); dadosdisp = dadosdisp.name; %determina o nome do .nc que sera carregado
    if isempty(levind)
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf
        dado = cat(4,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena os dados sem nivel
    else
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf
        dado = cat(5,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena os dados com nivel
    end
    lat = double(ncread([dirref dadosdisp], 'lat' )); %carrega a latitude
    lon = double(ncread([dirref dadosdisp], 'lon' )); %carrega a longitude
    %----------------------------------------------------------------------
    
    dado = squeeze(nanmean(reshape(dado,size(lon,1),size(lat,1),[],size(timemod,1),size(modelos,2)+1),4)); %realiza a media territorial para cada timestep
    
    %filtro nans-----------------------------------------------------------
    NAN = nan(size(dado)); %cria uma matriz nan no tamanho e dimensoes da variavel dado
    NAN(~isnan(dado)) = 1; %a partir da variavel dado, se substitui tudo que nao é nan dentro de dados por 1 dentro da matriz nan
    NAN = mean(NAN,3); %cria uma matriz unica com a posicao de nans para todos os modelos
    dado = dado.*NAN; %aplica essa matriz nan dentro dentro de todos os modelos padronizando todas as areas que nao contem dados
    %----------------------------------------------------------------------
    
end

%%

obs = dado(:,:,1);
obs = permute(reshape(obs,1,[],size(obs,3)),[2 1]);

models = dado(:,:,2:end);
models = permute(reshape(models,1,[],size(models,3)),[2 3 1]);

[alpha,n,crps,AR,STD,it] = ODR(obs,models,40000);

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
