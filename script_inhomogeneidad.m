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

%Parámetros inhomogeneidad
%sih = 2400;             %Posicion (m)
sih  = 6000;
lih = 600;               %Longitud (m)
Tf  = 2;                 %T final  (s)

%Parámetros entrada vehículos
lambda = 1.5;            %1 coche por lambda s 

%Entrada de vehículos: (esto en realidad se llama beta, pero bueno)
lambda = 1/lambda;       %(s)

%Introducción de vehículos para cumplir condición de contorno:
llegadas = zeros(1,numel(0:dt:tmax));
vin = @(va,s) min(v0,real((sqrt(ab.^2.*T.^2-2.*ab.*T.*va+va.^2-4.*ab.*(s0-s))-ab.*T+va)./2));
da = @(va) (s0+va*T)/sqrt(1-(va/v0)^4);
%Generación velocidades deseadas aleatorias
V0    = v0+desv.*randn(1,ceil(2*lambda*tmax)); %esto asegura que la vmax=v0
I     = V0 <= v0-0.7*desv | V0 > v0+0.7*desv;
V0(I) = v0; %Y que la velocidad mínima no va a ser menor que mu-n*sigma
%figure, hist(V0);

%Inhomogeneidad:
m     = (Tf-T)/lih;
n     = T-sih/lih*(Tf-T);
Tin   = @(x) m*x+n;

%Variables de salida y control
%A = [S T V a Ds Dv id s0 T a ab v0 l];

id = 1;
A = [10000 0 0 0 0 0 -1 0 0 0 0 0 0];
A = [0 0 V0(id) 0 A(1,1)-l 0 1 s0 T a ab v0 l; A];
id = id +1;

numt = numel(0:dt:tmax);
nums = ceil(smax/(s0+l));
stv_gt = zeros(numt*nums,7);

j = 1;
k = 1;

tic
for t = 0:dt:tmax
   
   %Llegadas Poissonianas, parámetro lambda
   t0 = 0;
   arrivaltimes=[];
   while t0<=dt
       t0 = t0 - log(rand)/lambda; %tiempo a la proxima llegada
       arrivaltimes = [arrivaltimes t0];
   end
   umbral = sum(arrivaltimes<=dt);
   %deltaxe = da(A(1,3));
   %Introduzco el coche, si cabe
   if umbral >= 1
       if A(1,1) > s0+l
            %A = [S T V a Ds Dv id s0 T a ab v0 l];
            v_temp = vin(A(1,3),A(1,1)-A(1,13));
            %v_temp = vin(A(1,3),0);
            %v_temp = A(1,3);
            A = [0 t v_temp 0 A(1,1)-A(1,13) v_temp-A(1,3) id s0 T a ab V0(id) l; A];
            %A = [0 t v_temp 0 A(1,1)-A(1,13) 0 id s0 T a ab V0(id) l; A];
            id = id + 1;
            llegadas(k) = 1;
       end
   end
   k = k+1;
   
   %Inhomogeneidad timehead
   I1 = A(:,1) > sih;
   A(I1,9) = min(Tin(A(I1,1)),Tf);
   
   %Paso de simulación
   A = IDM_engine(A,dt);
   
   %Elimina coches de la carretera cuando alcanzan el final
   I = A(:,1) > smax & A(:,7) > 0;
   A(I,:) = [];
   
   %Guarda en el campo
   stv_gt(j:(j-1)+size(A,1),:) = A(:,1:7);
   j = j+size(A,1);
end

%Limpieza campo STV
stv_gt = stv_gt(stv_gt(:,7) > 0,:);

param.s0   = s0;
param.T    = [T Tf];
param.sih  = sih;
param.lih  = lih;
param.a    = a;
param.b    = b;
param.ab   = ab;
param.v0   = v0;
param.desv = desv;
param.l    = l;
param.V0   = V0;

param.dt = dt;
param.tmax = tmax; %en s.
param.t0 = t0;
param.smax = smax;
param.lambda = 1/lambda;
param.llegadas = llegadas;

clear s0 Tt Ss a b ab v0 desv l V0 dt tmax smax T I1 I2 I3 I4 In dIn t0 vin n v_temp
clear umbral lambda A I id k t j numt nums sih plot_rt m Tin  Tf lih arrivaltimes llegadas

toc
%% dibuja el campo, haciendo diezmado de puntos (para que no tarde mil años):

% fracción de puntos a representar
prob = 0.01;  %fracción 
tipo = 1;     %puntos al azar u ordenados por vehículo
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

%pinta:
%A = [S T V a Ds Dv id s0 T a ab v0 l];
scatter(stv_dibuja(:,2),stv_dibuja(:,1),1,3.6.*stv_dibuja(:,3),'filled');

axis([0 param.tmax 0 param.smax]);
xlabel('time (min)','fontsize',14);
ylabel('space (km)','fontsize',14);
box on
caxis([0 3.6*param.v0]);
colorbar('NorthOutside');
title('speed (km/h)','fontsize',14);

%% graba si nos gusta
clear I h prob stv_dibuja n w
n = 1;
nombre = ['lambda_' num2str(10*param.lambda) '_tf_' num2str(10*param.T(2)) '_' num2str(n) '.mat'];
clear n;
save(nombre) 
%% ST-rho
rho = 1000./(stv_dibuja(:,5)+param.l);

figure
%A = [S T V a Ds Dv id s0 T a ab v0 l];

scatter(stv_dibuja(:,2),stv_dibuja(:,1),5,rho(:));
axis([0 param.tmax 0 param.smax]);
xlabel('Time (s)');
ylabel('Space (m)');
%zlabel('Speed (m/s)');
grid on
caxis([0 150]);
colorbar

%ST-Q
figure,
scatter(stv_dibuja(:,2),stv_dibuja(:,1),5,3.6.*stv_dibuja(:,3).*rho);
xlabel('Time (s)');
ylabel('Space (m)');
grid on
caxis([0 2000]);
colorbar;

clear i Ts Tss h ids;

% Scatter RHO-V
figure,
%A = [S T V a Ds Dv id s0 T a ab v0 l];
scatter(rho,3.6.*stv_dibuja(:,3),1,stv_dibuja(:,4),'filled')
axis([0 150 0 120]);
colorbar

% Scatter RHO-Q
figure,
%A = [S T V a Ds Dv id s0 T a ab v0 l];
scatter(rho,3.6.*stv_dibuja(:,3).*rho,1,stv_dibuja(:,4),'filled')
axis([0 150 0 2500]);
colorbar
%caxis([-5 5]);

figure,
scatter(stv_dibuja(:,6),stv_dibuja(:,5),1,stv_dibuja(:,3),'filled')
axis([-10 33 0 200]);
%% Diezmado a 1s 
tiempos = 0:param.tmax;
I = ismember(stv_gt(:,2),tiempos);
stv_gt_n = stv_gt(I,:);
%%
h = figure;
%A = [S T V a Ds Dv id s0 T a ab v0 l];

subplot(1,2,1), scatter(stv_gt_n(:,2),stv_gt_n(:,1),5,3.6.*stv_gt_n(:,3),'filled');
axis([1000 2000 2000 3000]);
xlabel('Time (s)');
ylabel('Space (m)');
%zlabel('Speed (m/s)');
grid on
caxis([0 3.6*param.v0]);
colorbar

subplot(1,2,2), scatter(stv_gt(:,2),stv_gt(:,1),5,3.6.*stv_gt(:,3),'filled');
axis([1000 2000 2000 3000]);
xlabel('Time (s)');
ylabel('Space (m)');
%zlabel('Speed (m/s)');
grid on
caxis([0 3.6*param.v0]);
colorbar
%% Llegadas
tiempos = 0:param.dt:param.tmax;
llegadas = tiempos(logical(param.llegadas));
figure
hold on
%t = 0:param.lambda:param.tmax;
%distrib = [];
%for j = 1:length(t)-1
%    I = llegadas > t(j) & llegadas <= t(j+1); 
%    distrib = [distrib sum(I)]; 
%end
%fdx = histc(distrib,0:10);
%figure
%plot(0:10,poisspdf(0:10,1/param.lambda),'.-g');
%hold on
%plot(0:10,fdx./sum(fdx),'.-r');
%title(['Lambda= ' num2str(param.lambda)]);

% Comparación con distribución exponencial
ll = 0:0.25:20;
ll2 = 0:0.05:20;

fdemp = histc(diff(llegadas),ll);
fdte  = histc(exprnd(param.lambda,100000,1),ll);
fdte2 = histc(exprnd(mean(diff(llegadas)),100000,1),ll2);
plot(ll,cumsum(fdemp./sum(fdemp)),'-g','LineWidth',2);
%plot(ll,cumsum(fdte./sum(fdte)),[symb(k,:) 'r'],'LineWidth',2);
plot(ll2,cumsum(fdte2./sum(fdte2)),'--r','LineWidth',2);
axis([0 20 0 1]);