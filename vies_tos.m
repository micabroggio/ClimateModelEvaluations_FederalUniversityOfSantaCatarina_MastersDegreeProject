% Micael Fernando Broggio
% Criacao: 27 de janeiro de 2020
%
% FIGURA VIES SUPERFICIAL
%
% Para dados .nc.
%
% Scrip desenvolvido para plotar o VIES dos modelos CMIP6 e BESM em relação
% ao dado de referencia ORAS5
% 
% Dados de input estao em medias mensais
%
% As variaveis de salidade e temperatura do ORAS foram alteradas para tos,
% sos, thetao e so
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/m_map')

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
    disp(var)
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
    clear dadoref dadomod
    
    bias = dado - dado(:,:,1); %calculo do VIES
    
    maxdado = max(max(max(dado(:,:,1)))); %encontra o maior valor no dado ORAS5
    mindado = min(min(min(dado(:,:,1)))); %encontra o menor valor no dado ORAS5
    cbarlim1 = floor(mindado):2:ceil(maxdado)+1; %cria o vetor para delimitar a colorbar
    l1 = length(cbarlim1); %comprimento do vetor da colorbar
    
    maxbias = max(max(max(bias))); %encontra o maior VIES em todos os dados
    minbias = min(min(min(bias))); %encontra o menor VIES em todos os dados
    cbarlim2 = ceil(max(abs(maxbias),abs(minbias))); %cria o vetor para delimitar a colorbar VIES
    l2 = length(-cbarlim2:cbarlim2); %comprimento do vetor da colorbar VIES
    
    %FIGURA VIES
    f1 = figure;
    set(f1, 'Position', [1 1 1000 500]); %set de tamanho da figura
    
    %plot do dado do ORAS (referencia)-------------------------------------
    sp = subplot(2,4,1);
    sp.Position = sp.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot
    m_proj('miller','lon',[round(min(lon)) round(max(lon))],'lat',[ceil(min(lat)) floor(max(lat))]); %set para a projecao do mapa  
    [CS,CH]=m_contourf(lon,lat,dado(:,:,1)',cbarlim1,'edgecolor','none'); %plot utilizando linhas de contorno
    hold on
    [CS2,CH2]=m_contour(lon,lat,dado(:,:,1)',[-10 18],'color','k','linestyle','--'); %linha de 18 de temperatura
    colormap(sp,(m_colmap('jet',l1-1))); %set na coloracao do mapa
    m_coast('patch',[.8 .8 .8],'edgecolor','k'); %linha de costa
    m_grid('tickdir','in','linewi',1,'xtick',[-180:30:180],'ytick',[-90:15:90],'xticklabels',[],'yticklabels',[]); %grade de latitude e longitude
    title(['a) ' labelvar{1} ' ORAS5']) %titulo do mapa a)
    
    ax = colorbar; %insere colobar
    ax.Location = 'west'; %set a colorbar do lado esquerdo do subplot
    ax.Position = ax.Position + [-.035 .0 .0 .0]; %ajuste de posicao da colorbar
    ax.AxisLocation = 'out'; %set dos numeros da colobar para o lado de fora do subplot 
    ax.TickDirection = 'out'; %set da orientacao dos marcadores para o lado de fora do subplot
    ax.Ticks = [cbarlim1(1):4:cbarlim1(end)]%set da orientacao dos marcadores para o lado de fora do subplot
    caxis([cbarlim1(1) cbarlim1(end)]) %limites da colorbar
    %----------------------------------------------------------------------
    
    %plot dos VIESES de todos os modelos-----------------------------------
    letra = {'b','c','d','e','f','g','h'}; %vetor para a variacao de letras nos titulos de cada mapa
    for j=2:size(ndado,2) %loop para todos os VIESES
        sp = subplot(2,4,j);
        sp.Position = sp.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot
        m_proj('miller','lon',[round(min(lon)) round(max(lon))],'lat',[ceil(min(lat)) floor(max(lat))]); %set para a projecao do mapa
        [CS,CH]=m_contourf(lon,lat,bias(:,:,j)',-cbarlim2:cbarlim2,'edgecolor','none'); %plot utilizando linhas de contorno
        hold on
        [CS2,CH2]=m_contour(lon,lat,dado(:,:,1)',[-10 18],'color','k','linestyle','--'); %linha de 18 de temperatura
        [CS2,CH2]=m_contour(lon,lat,dado(:,:,j)',[-10 18],'color','k','linestyle','-'); %linha de 18 de temperatura
        colormap(sp,(m_colmap('div',l2-1))); %set na coloracao do mapa
        m_coast('patch',[.8 .8 .8],'edgecolor','k'); %linha de costa
        if j == size(ndado,2)/2+1
            m_grid('tickdir','in','linewi',1,'xtick',[-180:30:180],'ytick',[-90:15:90]); %grade de latitude e longitude. Neste caso os labels serao visiveis
            xlabel('longitude') %label do eixo x          
            yl = ylabel('latitude'); %label do eixo y
            yl.Position = yl.Position + [-.08 0 0]; %ajuste de posicao do label y
        else
            m_grid('tickdir','in','linewi',1,'xtick',[-180:30:180],'ytick',[-90:15:90],'xticklabels',[],'yticklabels',[]); %grade de latitude e longitude.
        end
        title([letra{j-1} ') ' ndado{j} ' x ORAS5']) %titulos dos demais mapas
                
        if j == size(ndado,2)/2
            cbarpos1 = sp.Position; %posicao do mapa superior direito
        elseif j == size(ndado,2) 
            cbarpos2 = sp.Position; %posicao do mapa inferior direito
            ax = colorbar; %colorbar somente no mapa inferior direito
            ax.Location = 'east'; %set de colobar para o lado direito do subplot
            ax.Position = ax.Position + [.035 -(ax.Position(2)-cbarpos2(2)) .0 ((cbarpos1(2)+cbarpos1(4))-cbarpos2(2))-ax.Position(4)]; %ajuste de posicao da colorbar com referencia aos mapas superior e inferior direito
            ax.AxisLocation = 'out'; %set dos numeros da colobar para o lado de fora do subplot
            ax.TickDirection = 'out'; %set da orientacao dos marcadores para o lado de fora do subplot
            ax.Ticks = [-cbarlim2:3:cbarlim2]; %set da orientacao dos marcadores para o lado de fora do subplot
        end
        caxis([-cbarlim2 cbarlim2]) %limites da colorbar    
    end
    set(gcf,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
    set(gcf,'color','w'); %figura emfundo branco
    print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/vies_superficie_' reg '_' var{i}],'-r1000') %save da figura em png e 1000ppi de resolucao
    %----------------------------------------------------------------------    
end
