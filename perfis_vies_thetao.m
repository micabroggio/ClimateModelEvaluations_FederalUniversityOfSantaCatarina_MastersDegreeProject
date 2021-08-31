% Micael Fernando Broggio
% Criacao: 27 de janeiro de 2020
%
% FIGURA VIES SUBSUPERFICIAL
%
% Para dados .nc.
%
% Scrip desenvolvido para plotar o VIES dos modelos CMIP6 e BESM em relação
% ao dado de referencia ORAS5
% 
% Dados de input estao em medias mensais
%
% As variaveis de nivel do ORAS foram alteradas para lev
%
% Ha a necessidade de se ter as matrizes de densidade ja prontas para o
% plot das figuras.
%
%==========================================================================

clc, clear all, close all %limpa workspace, command window e fecha todo popup matlab

addpath('~/Documents/MATLAB/')
addpath('~/Documents/MATLAB/m_map')

dirin = '/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/DATA/'; %diretorio onde sao encontrados os dados anuais dos modelos
reg = 'AS'; %escala do processamento (AS ou GLOBAL)

var = {'thetao'}; %variaveis processadas
labelvar = {'\Theta'}; %nome dada a variavel no plot

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

% Perfil A_________________________________________________________________

load('/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/OUTPUTS/densidade_perfis_hist'); %carrega as matrizes preprocessadas de densidade
denss = {dens(:).A}; dens = reshape(cell2mat(denss),size(denss{1},1),size(denss{1},2),[]); clear denss %cria a variavel dens importando de celulas e transformando em matriz, ja dividindo entre todos os modelos
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
            dadomod(:,:,:,:,nm) = nanmean(ncread([dirmod dadosdisp], var{i} ),4); %carrega o netcdf
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
        dadoref = nanmean(ncread([dirref dadosdisp], var{i} ),4); %carrega o netcdf
        dado = cat(5,dadoref,dadomod); ndado = cat(2,ref,modelos); %concatena os dados com nivel
    end
    lat = double(ncread([dirref dadosdisp], 'lat' )); %carrega a latitude
    lon = double(ncread([dirref dadosdisp], 'lon' )); %carrega a longitude
    lev = double(ncread([dirref dadosdisp], 'lev' )); %carrega os niveis
    [Y,X,Z] = meshgrid(lat,lon,lev); %cria grade das variaveis lat lon lev
    %----------------------------------------------------------------------
    
    latperf = -25.5;
    [YY,XX,ZZ] = meshgrid(latperf,lon,lev); %cria uma grade dos modelos para interpolacao
    
    dado = squeeze(dado); %organiza os dados
    
    for j = 1:size(dado,4)
        dadoz(:,:,j) = squeeze(griddata(X,Y,Z,dado(:,:,:,j),YY,XX,ZZ)); %interpolacao dos dados dos modelos segundo as coordenadas do GOSHIP
    end
    YY = squeeze(YY); XX = squeeze(XX); ZZ = squeeze(ZZ);
    ZZZ = mean(ZZ,1);
    
    clear dadoref dadomod
    
    bias = dadoz - dadoz(:,:,1); %calculo do VIES
    
    maxdado = max(max(max(dadoz(:,:,1)))); %encontra o maior valor no dado ORAS5
    mindado = min(min(min(dadoz(:,:,1)))); %encontra o menor valor no dado ORAS5
    cbarlim1 = floor(mindado):1:ceil(maxdado); %cria o vetor para delimitar a colorbar
    l1 = length(cbarlim1); %comprimento do vetor da colorbar
    
    maxbias = max(max(max(bias))); %encontra o maior VIES em todos os dados
    minbias = min(min(min(bias))); %encontra o menor VIES em todos os dados
    cbarlim2 = ceil(max(abs(maxbias),abs(minbias))); %encontra o delimitador da colorbar VIES
    cbarlim2 = -cbarlim2+.5:.5:cbarlim2-.5; %cria o vetor da colorbar VIES
    l2 = length(cbarlim2); %comprimento do vetor da colorbar VIES
    
    %FIGURA VIES
    f1 = figure;
    set(f1, 'Position', [1 1 1000 500]); %set de tamanho da figura
    set(f1,'color','w'); %figura emfundo branco
    set(f1,'renderer','painters'); %determina o tipo de renderizador. Paiters tem melhor aplicacao para graficos em 2D
    
    %plot do dado do ORAS (referencia)-------------------------------------
    sp1 = subplot(2,4,1);
    sp1.Position = sp1.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot 
    [CS,CH]=contourf(XX(:,1)',-ZZ(1,ZZZ>=1000),dadoz(:,ZZZ>=1000,1),cbarlim1,'edgecolor','none'); %plot utilizando linhas de contorno para a profundidade maior que 1000 metros
    colormap(sp1,(m_colmap('jet',l1-1))); %configura o map de cores
    hold on
    [C,H] = contour(XX(:,1),-ZZ(1,ZZZ>1000),dens(:,ZZZ>1000,1),[25.7,26.8,27.53],'linestyle','-','facecolor','none','edgecolor','w'); %plot utilizando linhas de contorno para as densidades das massas de agua destacadas
    clabel(C,H,'FontSize',10,'Color','w') %label das linhas de densidade em braco
    sp1b = sp1.Position; %backup das proporcoes do subplot
    sp1.Position = sp1.Position + [0 (sp1.Position(4)*0.5) 0 -(sp1.Position(4)*0.5)]; %edita as dimencoes de cada subplot
    xticks(-50:10:20); xticklabels([]); %set de ticks x
    yticks(-5000:2000:-100); yticklabels([]); %set de ticks y
    set(gca,'Color',[.7 .7 .7]); %fundo do plot em cinza
    set(gca,'TickDir','in'); %colorca os ticks na area interna dos plots
    caxis([cbarlim1(1) cbarlim1(end)]); %aplica limits no colorbar
    grid %aplica o grid
    
    ax(2) = copyobj(sp1, gcf); %copia o subplot sp1     
    set(ax(2), 'position', sp1b + [0 0 0 -(sp1b(4)*0.5)]); %coloca a copia em outra posicao, logo abaixo do novo plot que esta por vir
    set(gca,'Yticklabel',[]); %apaga o ylabel
    ax(1) = axes('position', sp1.Position); %set novos eixos na posicao anterior de sp1, ou seja, a posicao antes da copia

    [CS,CH]=contourf(XX(:,1),-ZZ(1,ZZZ<1000),dadoz(:,ZZZ<1000,1),cbarlim1,'edgecolor','none'); %plot utilizando linhas de contorno para a profundidade maior que 1000 metros na nova regiao
    colormap(gca,(m_colmap('jet',l1-1))); %configura o map de cores
    hold on
    [C,H] = contour(XX(:,1),-ZZ(1,ZZZ<1000),dens(:,ZZZ<1000,1),[25.7,26.8,27.53],'linestyle','-','facecolor','none','edgecolor','w'); %plot utilizando linhas de contorno para as densidades das massas de agua destacadas
    clabel(C,H,'FontSize',10,'Color','w') %labels dos contornos em branco
    xticks(-100:10:10); xticklabels([]); %set de xticks   
    ylim([-1000 0]); yticks(-1000:200:0); yticklabels([]); %set de yticks
    set(gca,'Xticklabel',[]); %retira o label de x
    set(gca,'Color',[.7 .7 .7]); %fundo do plot em cinza
    set(gca,'TickDir','in'); %coloca os ticks na area interna dos plots
    linkaxes(ax, 'x'); %vincula os eixos da copia com a nova figura
    grid %aplica grade
    title(['a) ' labelvar{1} ' ORAS5 - A10']) %titulo do mapa a)
    
    ax = colorbar; %insere colobar
    ax.Location = 'west'; %set a colorbar do lado esquerdo do subplot
    ax.Position = ax.Position + [-.035 -sp1b(4)*0.5 .0 sp1b(4)*0.5]; %ajuste de posicao da colorbar
    ax.AxisLocation = 'out'; %set dos numeros da colobar para o lado de fora do subplot 
    ax.TickDirection = 'out'; %set da orientacao dos marcadores para o lado de fora do subplot
    ax.Ticks = [round(cbarlim1(1)):4:round(cbarlim1(end))];%set da orientacao dos marcadores para o lado de fora do subplot
    caxis([cbarlim1(1) cbarlim1(end)]) %limites da colorbar
    
    letra = {'b','c','d','e','f','g','h'}; %vetor para a variacao de letras nos titulos de cada mapa 
    for j=2:size(ndado,2) %loop para todos os VIESES
        sp1 = subplot(2,4,j);
        sp1.Position = sp1.Position + [-.04 -.04 .04 .06]; %expande o tamanho do subplot 
        [CS,CH]=contourf(XX(:,1),-ZZ(1,ZZZ>=1000),dadoz(:,ZZZ>=1000,j),cbarlim2,'edgecolor','none'); %plot utilizando linhas de contorno para a profundidade maior que 1000 metros
        colormap(sp1,(m_colmap('div',l2-1))); %configura o mapa de cores
        hold on
        [C,H] = contour(XX(:,1),-ZZ(1,ZZZ>1000),dens(:,ZZZ>1000,j),[25.7,26.8,27.53],'linestyle','-','facecolor','none','edgecolor','k'); %plot utilizando linhas de contorno para as densidades das massas de agua destacadas
        clabel(C,H,'FontSize',10,'Color','k') %labels dos contornos em preto
        if j == size(ndado,2)/2+1 %condicional para o subplot 5 pois e o unico que exibe os eixos x e y nomeados
            xticks(-50:10:20); %set de xticks
            yticks(-5000:2000:-100); yticklabels({'5000', '3000','1000'});  %set de yticks
        else %os demais subplots sem eixos nomeados
            xticks(-50:10:20); xticklabels([]); %set de xticks
            yticks(-5000:2000:-100); yticklabels([]); %set de yticks
        end
        sp1b = sp1.Position; %copia de posicao do subplot
        sp1.Position = sp1.Position + [0 (sp1.Position(4)*0.5) 0 -(sp1.Position(4)*0.5)]; %set de novo posicionamento para o subplot atual
        set(gca,'Color',[.7 .7 .7]); %cor do fundo de plot em cinza
        set(gca,'TickDir','in'); %ticks para o interior do plot
        grid %aplica grade
        caxis([cbarlim2(1) cbarlim2(end)]) %limites da colorbar 

        ax(2) = copyobj(sp1, gcf); %copia o subplot sp1        
        set(ax(2), 'position', sp1b + [0 0 0 -(sp1b(4)*0.5)]); %coloca a copia em outra posicao, logo abaixo do novo plot que esta por vir
        set(gca,'Yticklabel',[]);  %apaga o ylabel
        ax(1) = axes('position', sp1.Position); %set novos eixos na posicao anterior de sp1, ou seja, a posicao antes da copia

        [CS,CH]=contourf(XX(:,1),-ZZ(1,ZZZ<1000),dadoz(:,ZZZ<1000,j),cbarlim2,'edgecolor','none'); %plot utilizando linhas de contorno para a profundidade maior que 1000 metros
        colormap(gca,(m_colmap('div',l2-1))); %configura o map de cores
        hold on
        [C,H] = contour(XX(:,1),-ZZ(1,ZZZ<1000),dens(:,ZZZ<1000,j),[25.7,26.8,27.53],'linestyle','-','facecolor','none','edgecolor','k'); %plot utilizando linhas de contorno para as densidades das massas de agua destacadas
        clabel(C,H,'FontSize',10,'Color','k'); %label dos contornos em preto
        if j == size(ndado,2)/2+1 %condicional para o subplot 5 pois e o unico que exibe os eixos x e y nomeados
            xticks(-100:10:10); xticklabels([]); %set de xticks
            ylim([-1000 0]); yticks(-1000:200:0); yticklabels({' ','800','600','400','200'}); %set de yticks
            xl = xlabel('longitude'); %xlabel
            set(xl,'Units','normalized'); %normaliza as unidades do eixo x
            phl = get(xl,'position'); %copia da posicao do xlabel
            set(xl, 'position', [phl(1) -1.15 phl(3)]); %correcao de posicao do xlabel
            yl = ylabel('depth'); %ylabel
            set(yl,'Units','normalized'); %normaliza as unidade do eixo y
            phl = get(yl,'position'); %copia da posicao do eixo y
            set(yl, 'position', [phl(1)-0.03 0 phl(3)]); %correcao de posicao do ylabel
        else %os demais subplots sem eixos nomeados
            xticks(-100:10:10); xticklabels([]); %set de xtick
            ylim([-1000 0]); yticks(-1000:200:0); yticklabels([]); %set de ytick
        end
        sp1b = sp1.Position; %copia da posicao do subplot
        set(gca,'Color',[.7 .7 .7]); %set fundo cinza no plot
        set(gca,'TickDir','in'); %set de tick interno ao plot
        linkaxes(ax, 'x'); %vincula os eixos da copia com a nova figura
        grid %aplica grade
        title([letra{j-1} ') ' ndado{j} ' x ORAS5']) %titulos dos demais mapas
        
        if j == size(ndado,2)/2 %condicional para subplot 4 para setar os tamanho da colorbar VIES
            cbarpos1 = sp1b; %posicao do mapa superior direito
        elseif j == size(ndado,2) % demais subplots
            cbarpos2 = sp1b; %posicao do mapa inferior direito
            ax = colorbar; %colorbar somente no mapa inferior direito
            ax.Location = 'east'; %set de colobar para o lado direito do subplot
            ax.Position = ax.Position + [.035 -(ax.Position(2)-cbarpos2(2)+cbarpos2(4)) .0 ((cbarpos1(2)+cbarpos1(4))-cbarpos2(2))-ax.Position(4)+cbarpos2(4)]; %ajuste de posicao da colorbar com referencia aos mapas superior e inferior direito
            ax.AxisLocation = 'out'; %set dos numeros da colobar para o lado de fora do subplot
            ax.TickDirection = 'out'; %set da orientacao dos marcadores para o lado de fora do subplot
            ax.Ticks = [round(cbarlim2(1)):1:round(cbarlim2(end))]; %set da orientacao dos marcadores para o lado de fora do subplot
        end
        caxis([cbarlim2(1) cbarlim2(end)]) %limites da colorbar 
    end
    set(gcf,'InvertHardcopy','off') %altera a copia impressa e permite salvamento com fundos coloridos
%     print('-dpng',['/home/micael/Documents/PROJETOS/ModCOSTA/CMIP6/FIGURAS/vies_subsuperficieA10_' reg '_' var{i}],'-r1000') %save da figura em png e 1000ppi de resolucao
end
clear dadoz ax cbarlim1 cbarlim2 cbarpos1 cbarpos2 H CH CS xl yl
close(f1)
%__________________________________________________________________________