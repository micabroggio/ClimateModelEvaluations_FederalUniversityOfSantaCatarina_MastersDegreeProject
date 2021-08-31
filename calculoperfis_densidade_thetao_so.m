% Micael Fernando Broggio
% Criacao: 27 de janeiro de 2020
%
% Calculo para DENSIDADE CMIP6
%
% Para dados .nc.
%
% Scrip desenvolvido para calcular e densidade nos modelos CMIP6 e BESM
% 
% Dados de input estao em medias mensais
%
% As variaveis de nivel do ORAS foram alteradas para lev
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/m_map')

dirin = '/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/DATA/'; %diretorio onde sao encontrados os dados anuais dos modelos
reg = 'AS'; %escala do processamento (AS ou GLOBAL)

var = {'thetao','so'}; %variaveis processadas

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

crz = {'A','B','C'};
latA = -25;
latB = -30;
latC = -35;

for i = 1:size(crz,2)
    for j = 1:size(var,2)
        disp(var{j})
        nm = 1; dadomod = []; %prealoca o contador de modelos e a variavel dado
        for k = modelos %loop para todos os modelos
            dirmod = [dirin 'MODELOS/' char(k) '/' reg '/1979-2014/']; %diretorio de cada modelo
            dadosdisp = dir([dirmod  '*' var{j} '*Omon*.nc']); dadosdisp = dadosdisp.name; %encontra o nome de cada arquivo .nc associado a variavel analisada
            levind = ncinfo([dirmod dadosdisp]); levind  = {levind.Dimensions.Name}; levind = cell2mat(strfind(levind,'lev')); %encontra de ha presenca de mais niveis ou so superficie

            %carrega os dados dos modelos--------------------------------------
            if isempty(levind) %dados superficiais (sem niveis)
                dadomod(:,:,:,nm) = ncread([dirmod dadosdisp], var{j} ); %carrega o netcdf
            else %dados subsuperficiais (com niveis)
                dadomod(:,:,:,:,nm) = nanmean(ncread([dirmod dadosdisp], var{j} ),4); %carrega o netcdf
            end
            %------------------------------------------------------------------

            timemod = ncread([dirmod dadosdisp],'time'); %carrega o netcdf
            time2 = datetime(1979,01,01):calmonths:datetime(2014,12,01); %cria o vetor tempo na mao
            nm = nm + 1; %contador de modelo 
        end

        %carrega os dados de reanalise-----------------------------------------
        dirref  = [dirin 'REANALISES/' ref{1} '/' reg '/']; %diretorio de cada modelo
        dadosdisp = dir([dirref  '*' var{j} '*.nc']); dadosdisp = dadosdisp.name; %determina o nome do .nc que sera carregado
        if isempty(levind)
            dadoref = ncread([dirref dadosdisp], var{j} ); %carrega o netcdf
            dado = cat(4,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena os dados sem nivel
        else
            dadoref = nanmean(ncread([dirref dadosdisp], var{j} ),4); %carrega o netcdf
            dado = cat(5,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena os dados com nivel
        end
        
        lat = double(ncread([dirref dadosdisp], 'lat' )); %carrega a latitude
        lon = double(ncread([dirref dadosdisp], 'lon' )); %carrega a longitude
        lev = double(ncread([dirref dadosdisp], 'lev' )); %carrega os niveis
        [Y,X,Z] = meshgrid(lat,lon,lev); %cria uma grade dos modelos para interpolacao
        %----------------------------------------------------------------------
        
        latperf = eval(['lat' crz{i}]);
        [YY,XX,ZZ] = meshgrid(latperf,lon,lev); %cria uma grade dos modelos para interpolacao

        dado = squeeze(dado); %organiza os dados

        for k = 1:size(dado,4)
            dadoz(:,:,k,j) = squeeze(griddata(X,Y,Z,dado(:,:,:,k),XX,YY,ZZ)); %interpolacao dos dados dos modelos para as coordenadas do GOSHIP
        end
    end

    for j = 1:size(dadoz,3)
        d(:,:,j) = sw_dens0(dadoz(:,:,j,2),dadoz(:,:,j,1))-1000; %Calcula a densidade apartir dos dados de sal e temp dos modelos e do ORAS5
        dens(j).modelos = ndado{j}; %Cria a estrutura para ser exportada
        eval(['dens(j).' crz{i} '= d(:,:,j);'])
    end
    clear dadoz dadomod d
end
save(['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/OUTPUTS/densidade_perfis_hist'], 'dens', '-v7.3'); %salva a struct gerada do diretorio de saida)
