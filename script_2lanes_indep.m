clear all
close all
%dt, tiempo máximo, longitud carretera
dt      = 0.25;         %(s) 
tmax    = 8000;         %(s)
%smax   = 4500;         %(m) 
smax    = 8000;
%Parámetros IDM (malos conductores)

%T       = 1.2;          %(s)       %timehead demasiado pequeño
%a       = 0.8;          %(m/s^2)   %aceleración demasiado pequeña
%b       = 2.5;          %(m/s^2)   %decelación brusca
%ab      = 2*sqrt(a*b);  %(m/s^2)
%s0      = 2;            %(m)
%l       = 5;            %(m)
%v0      = 120/3.6;      %(m/s)

%Parámetros IDM (buenos conductores)

s0      = 2;            %(m)        
T       = 1.6;          %(s)       %timehead recomendado internacional
a       = 1;            %(m/s^2)   %aceleración suave pero suficiente
b       = 1.65;         %(m/s^2)   %deceleración suave
ab      = 2*sqrt(a*b);  %(m/s^2)
l       = 5;            %(m)
v0      = 120/3.6;      %(m/s)

%Aleatoriedad velocidad deseada
desv    = 0/3.6;         %(m/s)

%Parámetros entrada vehículos lanes:
lambda1 = 5; %1 coche por lambda s lane 1 
lambda2 = 2; %1 coche por lambda s lane 2

%Introducción de vehículos para cumplir condición de contorno:
llegadas1 = zeros(1,numel(0:dt:tmax));
llegadas2 = zeros(1,numel(0:dt:tmax));

vin = @(va,s) min(v0,real((sqrt(ab.^2.*T.^2-2.*ab.*T.*va+va.^2-4.*ab.*(s0-s))-ab.*T+va)./2));
da  = @(va) (s0+va*T)/sqrt(1-(va/v0)^4);

%Generación velocidades deseadas aleatorias
V0    = v0+desv.*randn(1,ceil(2/lambda1*tmax)+ceil(2/lambda2*tmax)); %esto asegura que la vmax=v0
I     = V0 <= v0-0.7*desv | V0 > v0+0.7*desv;
V0(I) = v0; %Y que la velocidad mínima no va a ser menor que mu-n*sigma

%Variables de salida y control CARRIL 1:
%Ax = [S T V a Ds Dv id s0 T a ab v0 l lane_x];

%carril 1
id = 1;
A1 = [10000 0 0 0 0 0 -1 0 0 0 0 0 0 1];
A1 = [0 0 V0(id) 0 A1(1,1)-l 0 1 s0 T a ab v0 l 1; A1];

%carril 2
id = 2;
A2 = [10000 0 0 0 0 0 -2 0 0 0 0 0 0 2];
A2 = [0 0 V0(id) 0 A2(1,1)-l 0 1 s0 T a ab v0 l 2; A2];
id = id + 1;

numt = numel(0:dt:tmax);
nums = ceil(smax/(s0+l));
stv_gt = zeros(2*numt*nums,8);

j = 1;
k1 = 1;
k2 = 1;

tic
for t = 0:dt:tmax
   
   %CARRIL 1 
   %Llegadas Poissonianas, parámetro lambda1
   t0 = 0;
   arrivaltimes1=[];
   while t0<=dt %
       t0 = t0 - log(rand)*lambda1; %tiempo a la proxima llegada
       arrivaltimes1 = [arrivaltimes1 t0];
   end
   umbral1 = sum(arrivaltimes1<=dt);
   %Introduzco el coche, si cabe
   if umbral1 >= 1
       if A1(1,1) > s0+l
            %A = [S T V a Ds Dv id s0 T a ab v0 l lane_x];
            v_temp = vin(A1(1,3),A1(1,1)-A1(1,13));
            %v_temp = vin(A(1,3),0);
            %v_temp = A(1,3);
            A1 = [0 t v_temp 0 A1(1,1)-A1(1,13) v_temp-A1(1,3) id s0 T a ab V0(id) l 1; A1];
            %A = [0 t v_temp 0 A(1,1)-A(1,13) 0 id s0 T a ab V0(id) l; A];
            id = id + 1;
            llegadas1(k1) = 1;
       end
   end
   k1 = k1+1;
   
   %Paso de simulación
   A1 = IDM_engine(A1,dt);
   
   %Elimina coches de la carretera cuando alcanzan el final
   I = A1(:,1) > smax & A1(:,7) > 0;
   A1(I,:) = [];
   
   %CARRIL 2:
   %Llegadas Poissonianas, parámetro lambda2
   t0 = 0;
   arrivaltimes2=[];
   while t0<=dt
       t0 = t0 - log(rand)*lambda2; %tiempo a la proxima llegada
       arrivaltimes2 = [arrivaltimes2 t0];
   end
   umbral2 = sum(arrivaltimes2<=dt);
   %Introduzco el coche, si cabe
   if umbral2 >= 1
       if A2(1,1) > s0+l
            %A = [S T V a Ds Dv id s0 T a ab v0 l lane_x];
            v_temp = vin(A2(1,3),A2(1,1)-A2(1,13));
            %v_temp = vin(A(1,3),0);
            %v_temp = A(1,3);
            A2 = [0 t v_temp 0 A2(1,1)-A2(1,13) v_temp-A2(1,3) id s0 T a ab V0(id) l 2; A2];
            %A = [0 t v_temp 0 A(1,1)-A(1,13) 0 id s0 T a ab V0(id) l; A];
            id = id + 1;
            llegadas2(k2) = 1;
       end
   end
   k2 = k2+1;
   
   %Paso de simulación
   A2 = IDM_engine(A2,dt);
   
   %Elimina coches de la carretera cuando alcanzan el final
   I = A2(:,1) > smax & A2(:,7) > 0;
   A2(I,:) = [];
   
   %Guarda en el campo
   stv_gt(j:(j-1)+size(A1,1),:) = A1(:,[1:7 14]);
   j = j+size(A1,1);
   
   stv_gt(j:(j-1)+size(A2,1),:) = A2(:,[1:7 14]);
   j = j+size(A2,1);
   
end

%Limpieza campo STV
stv_gt = stv_gt(stv_gt(:,7) > 0,:);

param.s0    = s0;
param.T     = T;
param.a     = a;
param.b     = b;
param.ab    = ab;
param.v0    = v0;
param.desv  = desv;
param.l     = l;
param.V0    = V0;

param.dt = dt;
param.tmax = tmax; %en s.
param.t0 = t0;
param.smax = smax;
param.lambda1 = lambda1;
param.llegadas1 = llegadas1;
param.lambda2 = lambda2;
param.llegadas2 = llegadas2;

clear s0 Tt Ss a b ab v0 desv l V0 dt tmax smax T I1 I2 I3 I4 In dIn t0 vin n v_temp
clear umbral lambda A I id k t j numt nums sih plot_rt m Tin  Tf lih arrivaltimes llegadas

toc
%% DIBUJA DISTINTOS CARRILES

% fracción de puntos a representar
prob = 0.05;  %fracción 
tipo = 0;     %puntos al azar u ordenados por vehículo
deltax = 100;
deltat = 120;

if prob == 1; 
    stv_dibuja = stv_gt;
elseif prob < 1
    [stv_dibuja, ~, ~] = sensordata10(stv_gt,prob,tipo,7);
end

wr = [0 + deltax param.smax - deltax 0 + deltat param.tmax-deltat];
%wr = [0 400 0 1000];
I = stv_dibuja(:,1) > wr(1) & stv_dibuja(:,1) <= wr(2) &...
    stv_dibuja(:,2) > wr(3) & stv_dibuja(:,2) <= wr(4);

stv_dibuja = stv_dibuja(I,:);

%ejes de la figura
h = figure;
%axes1 = axes('Parent',h,...
%    'YTickLabel',{'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5'},...
%    'YTick',[0 500 1000 1500 2000 2500 3000 3500 4000 4500],...
%    'XTickLabel',{'0','15','30','45','60','75','90','105','120','135','150'},...
%    'XTick',[0 900 1800 2700 3600 4500 5400 6300 7200 8100 9000],...
%    'FontWeight','bold',...
%    'FontSize',14);
hold on
colores = ['-k';'-r'];
%pinta:
%A = [S T V a Ds Dv id s0 T a ab v0 l];
ids = unique(stv_dibuja(:,7));
for i = 1:numel(ids)
    I = stv_dibuja(:,7) == ids(i);
    c = stv_dibuja(I,8);
    plot(stv_dibuja(I,2),stv_dibuja(I,1),colores(c(1),:));
end

axis([0 param.tmax 0 param.smax]);
xlabel('time (min)','fontsize',14);
ylabel('space (km)','fontsize',14);
box on
%caxis([0 3.6*param.v0]);
%colorbar('NorthOutside');
title('speed (km/h)','fontsize',14);
