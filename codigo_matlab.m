clear;      
close all;  
clc;        

%% 1. TREN DE PULSOS PERIODICO
% -------------------------------------------------------------------------
% En esta sección se define una señal rectangular periódica.

A  = 1;                  % Amplitud del pulso.
f0 = 1000;               % Frecuencia fundamental [Hz].
T  = 1/f0;               % Periodo [s].
D  = 0.50;               % Ciclo de trabajo: 50%.
K  = 20;                 % Se analizarán armónicos desde -20 hasta 20.
w0 = 2*pi*f0;            % Frecuencia angular fundamental [rad/s].

% Se crean 5000 puntos dentro de un periodo.
% Mientras más puntos tengamos, más precisa será la representación
Nt = 5000;

% Genera un vector de tiempo desde 0 hasta T.
t = linspace(0,T,Nt);

% -------------------------------------------------------------------------
% DEFINICION DEL TREN DE PULSOS
% Como D = 0.5, el pulso ocupa la mitad del periodo.
% -------------------------------------------------------------------------

x = A*(t < D*T);

% Se crea una figura para observar la señal en el dominio del tiempo.
figure;

% Se grafica el tiempo en milisegundos para que sea más fácil de leer.
plot(t*1e3,x,'b','LineWidth',1.5);
grid on;

xlabel('Tiempo [ms]');
ylabel('Amplitud');

title(sprintf('Tren de pulsos: D = %.2f',D));

% Se limita el eje vertical para visualizar mejor el pulso.
ylim([-0.1 1.1]);


%% 2. COEFICIENTES DE FOURIER ANALITICOS Y NUMERICOS
% -------------------------------------------------------------------------

% Vector de índices de armónicos.
% Como K = 20:
% k = [-20 -19 ... -1 0 1 ... 19 20]
k = -K:K;

% -------------------------------------------------------------------------
% CALCULO ANALITICO DE LOS COEFICIENTES
% -------------------------------------------------------------------------


Ck_an = A*D*(k==0) + ...
        A*sinc(k*D).*exp(-1j*pi*k*D).*(k~=0);

Ck_an(k==0) = A*D;


% -------------------------------------------------------------------------
% CALCULO NUMERICO DE LOS COEFICIENTES
% -------------------------------------------------------------------------

% Se crea un vector inicialmente lleno de ceros para guardar los
% coeficientes numéricos.
Ck_num = zeros(size(k));

% Se recorre cada armónico.
for ii = 1:length(k)

    Ck_num(ii) = (1/T)*trapz(t,...
        x.*exp(-1j*k(ii)*w0*t));

end


% -------------------------------------------------------------------------
% COMPARACION ENTRE RESULTADO ANALITICO Y NUMERICO
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

fprintf('Error maximo analitico-numerico = %.3e\n',...
    max(abs(Ck_an-Ck_num)));


%% 3. ESPECTRO BILATERAL
% -------------------------------------------------------------------------

% angle() obtiene la fase de cada coeficiente complejo en radianes.
fase = angle(Ck_an);

% Cuando la magnitud es prácticamente cero, la fase no tiene significado.
% En esos puntos colocamos NaN para que MATLAB no los grafique.
fase(abs(Ck_an)<1e-10) = NaN;


% -------------------------------------------------------------------------
% GRAFICA DEL ESPECTRO
% -------------------------------------------------------------------------

figure;

% ---------------- MAGNITUD ----------------
subplot(2,1,1);

% stem() genera un gráfico discreto mediante "palitos".
stem(k,abs(Ck_an),'filled');

grid on;
xlabel('k');
ylabel('|C_k|');

title('Espectro bilateral: magnitud');


% ---------------- FASE ----------------
subplot(2,1,2);

stem(k,fase,'filled');

grid on;
xlabel('k');
ylabel('Fase [rad]');

title('Espectro bilateral: fase');

% La fase se limita entre -pi y pi.
ylim([-pi pi]);


%% 4. ENVOLVENTE SINC, EFECTO DE D Y CANCELACION DE ARMONICOS

% Se crea un vector con muchos puntos para dibujar suavemente
% la envolvente.
k_env = linspace(-K,K,2000);

% Se calcula la envolvente:
%
%       A*D*|sinc(kD)|
%
env = A*D*abs(sinc(k_env*D));


% -------------------------------------------------------------------------
% GRAFICA DE LA ENVOLVENTE Y LOS COEFICIENTES
% -------------------------------------------------------------------------

figure;

subplot(2,1,1);

% Se grafica la envolvente de manera continua.
plot(k_env,env,'r','LineWidth',1.5);
hold on;

% Se superponen los valores reales de los coeficientes de Fourier.
stem(k,abs(Ck_an),'b','filled');

grid on;
xlabel('k');
ylabel('Magnitud');

title('Envolvente |AD sinc(kD)|');

legend('Envolvente','|C_k|');


% -------------------------------------------------------------------------
% EFECTO DEL CICLO DE TRABAJO D
% -------------------------------------------------------------------------

D_values = [0.25 0.50 0.75];

subplot(2,1,2);
hold on;

% Se analiza cada valor de D.
for Dv = D_values

    % Se calculan los coeficientes de Fourier correspondientes
    % al nuevo ciclo de trabajo.
    Cv = A*Dv*(k==0) + ...
         A*sinc(k*Dv).*exp(-1j*pi*k*Dv).*(k~=0);

    % Se asegura el valor correcto de C0.
    Cv(k==0) = A*Dv;

    % Se grafica la magnitud de los coeficientes.
    stem(k,abs(Cv),'filled');

end

grid on;
xlabel('k');
ylabel('|C_k|');

title('Efecto del ciclo de trabajo D');

legend('D=0.25','D=0.50','D=0.75');


% -------------------------------------------------------------------------
% CANCELACION DE ARMONICOS
% -------------------------------------------------------------------------

fprintf('\nCancelacion para D = %.2f:\n',D);

% Se buscan los armónicos entre 1 y K donde:
%
% sin(pi*k*D) aproximadamente = 0

canc = find(abs(sin(pi*(1:K)*D))<1e-10);

fprintf('Armonicos cancelados: ');

disp(canc);


%% 5. RECONSTRUCCION, GIBBS Y ERROR RMS
% -------------------------------------------------------------------------

K_values = [3 5 20 50];

% Se analiza cada cantidad de armónicos.
for ii = 1:length(K_values)

    Kr = K_values(ii);

    % -------------------------------------------------------------
    % VECTOR DE ARMÓNICOS
    % -------------------------------------------------------------
    kk = -Kr:Kr;

    % -------------------------------------------------------------
    % CALCULO DE LOS COEFICIENTES
    % -------------------------------------------------------------
    Ck = A*D*(kk==0) + ...
         A*sinc(kk*D).*exp(-1j*pi*kk*D).*(kk~=0);

    % Se asegura el valor correcto de C0.
    Ck(kk==0) = A*D;

    % -------------------------------------------------------------
    % RECONSTRUCCION DE LA SEÑAL
    % -------------------------------------------------------------
    xr = real(Ck*exp(1j*w0*kk'*t));
    % -------------------------------------------------------------
    % GRAFICA INDIVIDUAL
    % -------------------------------------------------------------
    figure;
    
    plot(t*1e3,x,'b','LineWidth',1.5);
    hold on;
    
    plot(t*1e3,xr,'r','LineWidth',1.3);
    
    grid on;
    
    xlabel('Tiempo [ms]');
    ylabel('Amplitud');
    
    title(sprintf('Reconstrucción con K = %d - Fenómeno de Gibbs',Kr));
    
    legend('Señal original','Señal reconstruida');
    
    ylim([-0.8 1.8]);


end


% -------------------------------------------------------------------------
% ERROR RMS DE LA RECONSTRUCCION
% -------------------------------------------------------------------------

% Se prueban diferentes cantidades impares de armónicos:
% 1, 3, 5, ..., 99

K_err = 1:2:100;

% Vector donde se almacenará el RMSE de cada caso.
RMSE = zeros(size(K_err));


% Se repite el procedimiento para cada K.
for ii = 1:length(K_err)

    Kr = K_err(ii);

    % Vector de armónicos.
    kk = -Kr:Kr;


    % Coeficientes de Fourier.
    Ck = A*D*(kk==0) + ...
         A*sinc(kk*D).*exp(-1j*pi*kk*D).*(kk~=0);

    Ck(kk==0) = A*D;


    % Reconstrucción.
    xr = real(Ck*exp(1j*w0*kk'*t));


    % Cálculo del error RMS.
    RMSE(ii) = sqrt(mean((x-xr).^2));

end


% -------------------------------------------------------------------------
% GRAFICA DEL ERROR
% -------------------------------------------------------------------------

figure;

plot(K_err,RMSE,'o-','LineWidth',1.3);

grid on;

xlabel('Numero de armonicos K');
ylabel('RMSE');

title('Error de reconstruccion vs K');


%% 6. MUESTREO Y DTFS
% -------------------------------------------------------------------------

N = 64;

% Índices de las muestras:
% n = 0,1,2,...,63

n = 0:N-1;


% -------------------------------------------------------------------------
% NUMERO DE MUESTRAS QUE PERTENECEN AL PULSO
% -------------------------------------------------------------------------

M = round(D*N);


% Se crea inicialmente una secuencia de N ceros.
xn = zeros(1,N);

% Las primeras M muestras toman el valor A.
xn(1:M) = A;


% -------------------------------------------------------------------------
% CALCULO DE LA DTFS
% -------------------------------------------------------------------------

Xk = (1/N)*xn*exp(...
    -1j*2*pi*(0:N-1)'*(0:N-1)/N);


% -------------------------------------------------------------------------
% GRAFICAS DE LA SEÑAL DISCRETA Y SU ESPECTRO
% -------------------------------------------------------------------------

figure;

% ---------------- SECUENCIA DISCRETA ----------------
subplot(2,1,1);

stem(n,xn,'filled');

grid on;

xlabel('n');
ylabel('x[n]');

title('Secuencia periodica muestreada');


% ---------------- DTFS ----------------
subplot(2,1,2);

stem(0:N-1,abs(Xk),'filled');

grid on;

xlabel('k');
ylabel('|X[k]|');

title('Magnitud de la DTFS');


%% 7. RELACION ENTRE FOURIER CONTINUO Y DTFS
% -------------------------------------------------------------------------

Kc = 10;

kc = -Kc:Kc;


% -------------------------------------------------------------------------
% COEFICIENTES DE FOURIER CONTINUO
% -------------------------------------------------------------------------

Cc = A*D*(kc==0) + ...
     A*sinc(kc*D).*exp(-1j*pi*kc*D).*(kc~=0);

% Se asegura el valor de C0.
Cc(kc==0) = A*D;


% -------------------------------------------------------------------------
% EXTRACCION DE LOS COEFICIENTES CORRESPONDIENTES DE LA DTFS
% -------------------------------------------------------------------------

Xc = Xk(mod(kc,N)+1);


% -------------------------------------------------------------------------
% COMPARACION GRAFICA
% -------------------------------------------------------------------------

figure;

% Fourier continuo.
stem(kc,abs(Cc),'b','filled');

hold on;

% DTFS.
stem(kc,abs(Xc),'r--');

grid on;

xlabel('k');
ylabel('Magnitud');

title('Fourier continuo vs DTFS');

legend('|C_k| continuo','|X[k]| DTFS');


% -------------------------------------------------------------------------
% CALCULO DEL ERROR ENTRE AMBOS
% -------------------------------------------------------------------------

fprintf('\nError maximo Fourier continuo-DTFS = %.3e\n',...
    max(abs(Cc-Xc)));
