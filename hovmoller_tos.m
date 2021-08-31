% Micael Fernando Broggio
% Criacao: 27 de janeiro de 2020
%
% FIGURA DIAGRAMA HOVMOLLER
%
% Para dados .nc.
%
% Scrip desenvolvido para plotar o VIES dos diagramas Hovmoller de cada
% modelo em comparacao com o dado de referencia
% 
% Dados de input estao em medias mensais
%
% As variaveis de salidade e temperatura do ORAS foram alteradas para tos,
% sos, thetao e so
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/m_map') %adiciona o path m_map para o processamento

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
        
        lat = double(ncread([dirmod dadosdisp],'lat')); %carrega lat dos modelos
        lon = double(ncread([dirmod dadosdisp],'lon')); %carrega lon dos modelos
        timemod = ncread([dirmod dadosdisp],'time'); %carrega o netcdf
        time2 = datetime(1979,01,01):calmonths:datetime(2014,12,01); %define um vetor tempo feito a mao
        
        %carrega os dados dos modelos--------------------------------------
        ano = 1;
        for k = 1:12:length(timemod) %loop para variacao temporal ao passo 12
            if isempty(levind) %condicional de nivel supercial                
                dadomod(:,1:12,ano,nm) = double(squeeze(nanmean(ncread([dirmod dadosdisp], ...
                    var{i}, [1 1 k], [inf inf 12], [1 1 1]),1))); %carrega o netcdf de 12 em 12 ja efetuando o calculo de mediazonal
            else %condicional de varios niveis
                dadomod(:,1:12,ano,nm) = double(squeeze(nanmean(reshape(permute(ncread([dirmod dadosdisp], ...
                    var{i}, [1 1 1 k], [inf inf inf 12], [1 1 1 1]),[2 1 3 4]),length(lat),1,[],12),3))); %carrega o netcdf de 12 em 12 ja efetuando o calculo de media zonal
            end
            ano = ano + 1; %contador de anos
        end
        %------------------------------------------------------------------
        nm = nm + 1; %contador de modelo 
    end
    
    %carrega os dados de reanalise-----------------------------------------
    dirref  = [dirin 'REANALISES/' ref{1} '/' reg '/']; %diretorio de cada modelo
    dadosdisp = dir([dirref  '*' var{i} '*.nc']); dadosdisp = dadosdisp.name; %encontra o nome no arquivo nc que sera processado
    if isempty(levind) %condicional de nivel supercial
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf do dado de referencia
        dadoref = reshape(squeeze(nanmean(dadoref,1)),length(lat),12,[]); %calcula a media zonal e rearranja os dados de referencia para se encaixarem na matrix dos modelos
        dado = cat(4,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena o dado ref aos dados dos modelos
    else %condicional de varios niveis
        dadoref = ncread([dirref dadosdisp], var{i} ); %carrega o netcdf do dado de referencia
        dadoref = squeeze(nanmean(reshape(permute(dadoref,[2 1 3 4]),length(lat),1,[],size(dadoref,4)),3)); %organiza os dados e realiza as media zonal de profundidade para dados de subsurficie
        dadoref = reshape(dadoref,length(lat),12,[]); %calcula a media zonal e rearranja os dados de referencia para se encaixarem na matrix dos modelos
        dado = cat(4,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena o dado ref aos dados dos modelos
    end
    %----------------------------------------------------------------------
    
    dadoSTD = squeeze(std(dado,0,3)); %realiza o desvio padrao do periodo para todos os modelos
    dado = squeeze(nanmean(dado,3)); %realiza a media anual dos dados
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
    [CS,CH]=contourf(1:12,lat,dado(:,:,1),cbarlim1,'edgecolor','none'); %plot utilizando linhas de contorno
    hold on
    colormap(sp,(m_colmap('jet',l1-1))); %set na coloracao do mapa
    xticks([1:12]); xticklabels([]); %set do eixo x
    yticks([-90:15:90]); yticklabels([]); %set do eixo y
    grid %adiciona grade
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
    for j=2:size(ndado,2) %loop para todos os VIESES (modelos)
        sp = subplot(2,4,j); %configura os subplots segundo os modelos
        sp.Position = sp.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot
        [CS,CH]=contourf(1:12,lat,bias(:,:,j),-cbarlim2:cbarlim2,'edgecolor','none'); %plot utilizando linhas de contorno
        hold on
        colormap(sp,(m_colmap('div',l2-1))); %set na coloracao do mapa
        if j == size(ndado,2)/2+1
            xticks([1:12]); xticklabels ({'J','F','M','A','M','J','J','A','S','O','N','D'}); %set do eixo x
            yticks([-90:15:90]); %set do eixo y
            xlabel('months') %label do eixo x          
            yl = ylabel('latitude'); %label do eixo y
            yl.Position = yl.Position + [-.08 0 0]; %ajuste de posicao do label y
        else
            xticks([1:12]); xticklabels([]); %set do eixo x
            yticks([-90:15:90]); yticklabels([]); %set do eixo y
        end
        title([letra{j-1} ') ' ndado{j} ' x ORAS5']); %titulos dos demais mapas
                
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
        caxis([-cbarlim2 cbarlim2]); %limites da colorbar 
        grid %adiciona grade
    end
    set(gcf,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
    set(gcf,'color','w'); %figura emfundo branco
    print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/hovmoller_' reg '_' var{i}],'-r1000') %save da figura em png e 1000ppi de resolucao
    %---------------------------------------------------------------------- 
    
    %FIGURA STD
    
    maxSTD = max(max(max(dadoSTD(:,:,1)))); %encontra o maior valor no dado ORAS5
    minSTD = min(min(min(dadoSTD(:,:,1)))); %encontra o menor valor no dado ORAS5
    cbarlim1 = floor(minSTD*10):.5:ceil(maxSTD*10); %cria o vetor para delimitar a colorbar mudando a casa decimal para o correto funcionamento das funcoes floor e ceil
    cbarlim1 = cbarlim1/10; %retorna os dados para a casa decimal original
    l1 = length(cbarlim1); %comprimento do vetor da colorbar
    
    f2 = figure;
    set(f2, 'Position', [1 1 1000 500]); %set de tamanho da figura
    
    %plot dos VIESES de todos os modelos-----------------------------------
    letra = {'a','b','c','d','e','f','g','h'}; %vetor para a variacao de letras nos titulos de cada mapa
    for j=1:size(ndado,2) %loop para todos os VIESES (modelos)
        sp = subplot(2,4,j); %configura os subplots segundo os modelos
        sp.Position = sp.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot
        [CS,CH]=contourf(1:12,lat,dadoSTD(:,:,j),cbarlim1,'edgecolor','none'); %plot utilizando linhas de contorno
        hold on
        colormap(sp,(m_colmap('jet',l1-1))); %set na coloracao do mapa
        if j == size(ndado,2)/2+1
            xticks([1:12]); xticklabels ({'J','F','M','A','M','J','J','A','S','O','N','D'}); %set do eixo x
            yticks([-90:15:90]); %set do eixo y
            xlabel('months') %label do eixo x          
            yl = ylabel('latitude'); %label do eixo y
            yl.Position = yl.Position + [-.08 0 0]; %ajuste de posicao do label y
        else
            xticks([1:12]); xticklabels([]); %set do eixo x
            yticks([-90:15:90]); yticklabels([]); %set do eixo y
        end
        title([letra{j} ') ' ndado{j} ' STD']); %titulos dos demais mapas
                
        if j == size(ndado,2)/2
            cbarpos1 = sp.Position; %posicao do mapa superior direito
        elseif j == size(ndado,2) 
            cbarpos2 = sp.Position; %posicao do mapa inferior direito
            ax = colorbar; %colorbar somente no mapa inferior direito
            ax.Location = 'east'; %set de colobar para o lado direito do subplot
            ax.Position = ax.Position + [.035 -(ax.Position(2)-cbarpos2(2)) .0 ((cbarpos1(2)+cbarpos1(4))-cbarpos2(2))-ax.Position(4)]; %ajuste de posicao da colorbar com referencia aos mapas superior e inferior direito
            ax.AxisLocation = 'out'; %set dos numeros da colobar para o lado de fora do subplot
            ax.TickDirection = 'out'; %set da orientacao dos marcadores para o lado de fora do subplot
            ax.Ticks = [cbarlim1(1):3:cbarlim1(end)]; %set da orientacao dos marcadores para o lado de fora do subplot
        end
        caxis([cbarlim1(1) cbarlim1(end)]); %limites da colorbar 
        grid %adiciona grade
    end
    set(gcf,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
    set(gcf,'color','w'); %figura emfundo branco
%     print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/hovmollerSTD_' reg '_' var{i}],'-r1000') %save da figura em png e 1000ppi de resolucao
    %---------------------------------------------------------------------- 
end
