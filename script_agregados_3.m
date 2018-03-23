% Este script hace representaciones espacio temporales de velocidad, flujo
% y densidad a partir de medidas en R(t) y R(t,s):

smax        = param.smax;
tmax        = param.tmax;
smin        = 400;

% Submuestreo del campo para 
tiempos = 0:1:param.tmax;
I       = ismember(stv_gt(:,2),tiempos);
stv_gt_n = stv_gt(I,:);

%Las regiones R(t,s) están separadas G y tienen unas dimensiones de Tmp x K

max_pac = 2;                % Número de coches máximo por segmento tamaño K
K = max_pac*(param.l+param.s0);  % K en función del número de coches máximo
D = 100;                         % Distancia entre zonas de medida
G = D - K + 1;                   % Separación en función de K
Tmp = 60;                        % Promedio temporal   
tmin = Tmp;

%Calcula parámetros de la rejilla
franjas = floor((smax-smin)/K);
smax    = franjas*K+smin;
NT      = floor((tmax-tmin)/Tmp);    
tmax    = NT*Tmp+tmin;
T       = (tmin:tmax)';
S       = (smin:smax)';
Tc      = (tmin:Tmp:tmax)';
Sc      = (smin:K+G:smax)';

disp('------------------------------');
disp(['Número de regiones R(t,:): ' num2str(numel(Sc))]);
disp(['Número de regiones R(:,s): ' num2str(NT)])
disp(['K:                         ' num2str(K)])
disp(['Tmp:                       ' num2str(Tmp)]);
disp(['Empaquetamiento máximo:    ' num2str(max_pac)]);
disp(['Tmax:                      ' num2str(max(T(:)))]);
disp(['Smax:                      ' num2str(max(S(:)))]);

Rts = zeros((length(Tc)-1)*length(Sc),4);
Rmat = zeros(length(Sc),length(Tc)-1,3);

kk = 1;

for i = 1:length(Sc)
    
    I = stv_gt_n(:,1) >= Sc(i) & stv_gt_n(:,1) <= Sc(i)+K;
    temp1 = stv_gt_n(I,:);
    
    for j = 1:length(Tc)-1
        I = temp1(:,2) >= Tc(j) & temp1(:,2) < Tc(j+1);
        temp2   = temp1(I,:);
        
        %%%%% CALCULOS %%%%%
        tiempos   = Tc(j):Tc(j+1)-1;
        rho       = 1000/(Tmp*K)*sum(histc(temp2(:,2),tiempos));
        q1        = zeros(Tmp,1);
        
        for k = 1:Tmp
            q1(k) = 1/K*sum(temp2(temp2(:,2) == tiempos(k),3));  
        end
        
        q         = (3600/Tmp)*sum(q1);
        
        Rts(kk,:) = [i j rho q];
        Rmat(i,j,:) = [rho q q/rho]; 
        kk        = kk + 1;
    end
end
clear i j k kk q rho I temp1 temp2 tiempos q1

%% Medidas de espiras puestas en línea

figure;
axes('YTickLabel',{'0','1','2','3','4','5','6','7','8'},...
     'YTick',[0 1000 2000 3000 4000 5000 6000 7000 8000],...
     'XTickLabel',{'0','22','44','66','88','110','132'},...
     'XTick',[0 1330 2660 3990 5320 6650 7980],...
     'ZTickLabel',{'0','20','40','60','80','100','120'},...
     'ZTick',[0 20 40 60 80 100 120],...
     'FontWeight','bold',...
     'ZDir','reverse',...
     'FontSize',14);

grid on         
axis([0 param.tmax 0 param.smax 0 120]);
hold on

for i = 1:length(Sc)
    I = Rts(:,1) == i;
    plot3(Tc(Rts(I,2)),repmat(Sc(i)+K/2,1,sum(I)),Rts(I,4)./Rts(I,3));
end 
hold off

%% Medidas de espiras puestas en grid
h = figure;
axes1 = axes('Parent',h,...
             'FontWeight','bold',...
             'ZDir','reverse',...
             'FontSize',14);

axis([0 param.tmax 0 param.smax 0 120]);
hold on

[Tmesh,Smesh] = meshgrid(Tc(1:end-1)-Tmp/2,Sc);

surf(Tmesh,Smesh,Rmat(:,:,3), 'edgecolor', 'none');
%contour3(Tmesh,Smesh,Rmat(:,:,3));
%surface(Tmesh,Smesh,Rmat(:,:,3), 'facecolor', 'white')
view(-36,36);
grid on