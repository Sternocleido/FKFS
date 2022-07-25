%% Auswertung EMG-Daten 
clear all; clc;
%% Daten laden
Tabelle = xlsread('EMG_Trials.xlsx.');
EMG_Nr = Tabelle(2:6,1:21);
%% 
s = 9;          % Anzahl der Probanden 
for p=1:s;      
    for x=1:5;
        if EMG_Nr(x,p)<10
        data(x,p)=load(['F:\FKFS\Simulationsdaten\Proband' num2str(p) '\Versuch000' num2str(EMG_Nr(x,p)) '.mat']); % Proband 1-> Versuch 10 (1.)
        else data(x,p)=load(['F:\FKFS\Simulationsdaten\Proband' num2str(p) '\Versuch00' num2str(EMG_Nr(x,p)) '.mat']); 
        end
    end 
end
%% Daten auswählen  
for p=1:s;
    for x=1:5;
        STR_L{x,p} = data(x,p).data(:,1);       % [mV],O-EMG M. Sternocleidomastoideus links
        STR_R{x,p} = data(x,p).data(:,2);       % [mV],O-EMG M. Sternocleidomastoideus rechts
        TRO_L{x,p} = data(x,p).data(:,3);       % [mV],O-EMG M.Trapezius pars descendens links
        TRO_R{x,p} = data(x,p).data(:,4);       % [mV],O-EMG M.Trapezius pars descendens rechts
        TRU_L{x,p} = data(x,p).data(:,5);       % [mV],O-EMG M.Trapezius pars transversa links
        TRU_R{x,p} = data(x,p).data(:,6);       % [mV],O-EMG M.Trapezius pars transversa rechts
        EKG{x,p} = data(x,p).data(:,7);         % [mV],EKG
        Ax_car{x,p} = data(x,p).data(:,8);      % [m/s²], x-Beschleunigung Kopf
        Ay_car{x,p} = data(x,p).data(:,9);      % [m/s²], y-Beschleunigung Kopf
        Az_car{x,p} = data(x,p).data(:,10);      % [m/s²], z-Beschleunigung Kopf
        Ax_head{x,p} = data(x,p).data(:,12);     % [m/s²], x-Beschleunigung Auto
        Ay_head{x,p} = data(x,p).data(:,13);     % [m/s²], y-Beschleunigung Auto
        Az_head{x,p} = data(x,p).data(:,14);     % [m/s²], z-Beschleunigung Auto
        t{x,p} = (0:0.0005:(length(data(x,p).data)/2000)-0.0005);     % [s], Zeit
    end 
end 

%%
figure(1)
sgtitle('Kopf- und Autobeschleunigung')
subplot(3,2,1);
plot(t{1,7},Ax_head{1,7})
title('x-head')
subplot(3,2,2);
plot(t{1,7},Ax_car{1,7})
title('x-car')
subplot(3,2,3);
plot(t{1,7},Ay_head{1,7})
title('y-head')
subplot(3,2,4);
plot(t{1,7},Ay_car{1,7})
title('y-car')
subplot(3,2,5);
plot(t{1,7},Az_head{1,7})
title('z-head')
subplot(3,2,6);
plot(t{1,7},Az_car{1,7})
title('z-car')

%% Beschleunigung 
% Beschleunigung car 
for p=1:s;
    for x=1:5;
        [Ay_car_max{x,p},t_car_max_index{x,p}] = max(abs(Ay_car{x,p}(:,1)));        % maximale Beschleunigung finden und dazugehöriger X-Index
        t_car_max{x,p} = t{x,p}(t_car_max_index{x,p});                                        % Zeit bei maximaler Beschleunigung 
        Ay_car_max_half{x,p} = Ay_car_max{x,p}/2;                                        % Hälfte der maximalen Beschleunigung 
        t_car_half{x,p} = t_car_max_index{x,p}/2;                                        % Index bei der Hälfte der Zeit zum Erreichen des Maxima 
        Ay_car_half_mean{x,p} = mean(abs(Ay_car{x,p}(1:t_car_half{x,p})));                         % Mittelwert der Beschleunigung von Beginn bis zur Hälfte der Zeit 
        Ay_car_half_std{x,p} = std(abs(Ay_car{x,p}(1:t_car_half{x,p})));                           % Standardabweichung des Mittelwerts 

% Berechnung Beginn Beschleunigung car
        n=1;
        while abs(Ay_car{x,p}(n))<= Ay_car_half_mean{x,p}+5*Ay_car_half_std{x,p}   % erster Wert > Mittelwert + 5*SD 
        n=n+1;
        end 
        Beginn_car_Beschleunigung{x,p}=n;
    end 
end 

% Beschleunigung head 
for p=1:s;
    for x=1:5;
        [Ay_head_max{x,p},t_head_max_index{x,p}] = max(abs(Ay_head{x,p}(:,1))); % maximale Beschleunigung finden und dazugehöriger X-Index           
        t_head_max{x,p} = t{x,p}(t_head_max_index{x,p});                                   % Zeit bei maximaler Beschleunigung 
        Ay_head_max_half{x,p} = Ay_head_max{x,p}/2;                                   % Hälfte der maximalen Beschleunigung 
        t_head_half{x,p} = t_head_max_index{x,p}/2;                                   % Index bei der Hälfte der Zeit zum Erreichen des Maxima 
        Ay_head_half_mean{x,p} = mean(abs(Ay_head{x,p}(1:t_head_half{x,p})));              % Mittelwert der Beschleunigung von Beginn bis zur Hälfte der Zeit 
        Ay_head_half_std{x,p} = std(abs(Ay_head{x,p}(1:t_head_half{x,p})));                % Standardabweichung des Mittelwerts 

% Berechnung Beginn Beschleunigung head
        z=1;
        while abs(Ay_head{x,p}(z))<= Ay_head_half_mean{x,p}+5*Ay_head_half_std{x,p}        % erster Wert > Mittelwert -6*SD 
        z=z+1;
        end 
        Beginn_head_Beschleunigung{x,p}=z;
    end
end 
%% 
% figure(2)
plot(abs(Ay_car{3,9}),'b'),title('Beschleunigung Kopf/Auto'), xlabel('Samples'),ylabel('Beschleunigung[m/s²]')
hold on 
plot(abs(Ay_head{3,9}),'g')
hold on
xline(Beginn_car_Beschleunigung{3,9},'r');
hold on
xline(24303,'m');
xline(Beginn_head_Beschleunigung{3,9},'s');
legend('Beschleunigung Auto', 'Beschleunigung Kopf')
%% Signalverarbeitung 
% EMG-Signal nullen 
for p=1:s;
    for x=1:5;
        STR_L_null{x,p} = STR_L{x,p}-mean(STR_L{x,p}(1:5000,1));
        STR_R_null{x,p} = STR_R{x,p}-mean(STR_R{x,p}(1:5000,1));
        TRO_L_null{x,p} = STR_L{x,p}-mean(TRO_L{x,p}(1:5000,1));
        TRO_R_null{x,p} = STR_L{x,p}-mean(TRO_R{x,p}(1:5000,1));
        TRU_L_null{x,p} = STR_L{x,p}-mean(TRU_L{x,p}(1:5000,1));
        TRU_R_null{x,p} = STR_L{x,p}-mean(TRU_R{x,p}(1:5000,1));
    end 
end 

%% 
% R Zacken des EKG Signals extrahieren
for p=1:s;
    for x=1:5;
    peaks{x,p}=[]; locs{x,p}=[];vek{x,p}=[];
    [peaks{x,p}(:,1) locs{x,p}(:,1)]=findpeaks(EKG{x,p},'MinPeakDistance',1000, 'MinPeakHeight',1); % Diese Zeile müsste auf Basisi des EKG Signals bei den einzelnen Probanden angepasst werden

    % Bau des Vektors zur Detektion des EKG Signals
    le{x,p}=length(EKG{x,p});
    vek{x,p}(1:le{x,p},1)=0;

    for j=1:length(locs{x,p})
        k=locs{x,p}(j,1);
        if k-70<0
           k=locs{x,p}(j+1,1);
        end
        vek{x,p}(k-80:k+80,1)=1;    % Fenstebreite des EKG Signals %140, 70
                                      % Diese Zeile müsste auf Basis des EKG Signals bei den einzelnen Probanden angepasst werden
    end
    end 
end 
%%
% Filtern des EMG Signals mit einem Butterworth Tiefpassfilter. 
fc=50; n=4; fs=2000; form='low';
for p=1:s;
    for x=1:5;
        [STR_L_FiltSignal{x,p}]=butterw(STR_L_null{x,p}, fc, fs, n, form);
        [STR_R_FiltSignal{x,p}]=butterw(STR_R_null{x,p}, fc, fs, n, form);
        [TRO_L_FiltSignal{x,p}]=butterw(TRO_L_null{x,p}, fc, fs, n, form);
        [TRO_R_FiltSignal{x,p}]=butterw(TRO_R_null{x,p}, fc, fs, n, form);
        [TRU_L_FiltSignal{x,p}]=butterw(TRU_L_null{x,p}, fc, fs, n, form);
        [TRU_R_FiltSignal{x,p}]=butterw(TRU_R_null{x,p}, fc, fs, n, form);
    end 
end 
%%
% Eliminierung der EKG Artefakte an den spezifischen Artefaktstellen aus dem EMG-Signal
for p=1:s;
    for x=1:5;
        for j=1:length(STR_L{x,p});
            STR_L_cor{x,p}(j,1)=STR_L_null{x,p}(j,1)-STR_L_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
            STR_R_cor{x,p}(j,1)=STR_R_null{x,p}(j,1)-STR_R_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
            TRO_L_cor{x,p}(j,1)=TRO_L_null{x,p}(j,1)-TRO_L_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
            TRO_R_cor{x,p}(j,1)=TRO_R_null{x,p}(j,1)-TRO_R_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
            TRU_L_cor{x,p}(j,1)=TRU_L_null{x,p}(j,1)-TRU_L_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
            TRU_R_cor{x,p}(j,1)=TRU_R_null{x,p}(j,1)-TRU_R_FiltSignal{x,p}(j,1)*vek{x,p}(j,1);
        end 
    end 
end 
%%
% Tiefpassfilter 500Hz
% F_T =500;                                      % [Hz], Frequenz 
% Fs_T =2000;                                    % Abtrastfrequnez
% [b,a]= butter(4,F_T/(Fs_T/2),'low');                  % 4nd Ordnung butterworth filter 
% for p=1:s;
%     for x=1:5;
%         STR_L_filter_T{x,p}= filtfilt(b,a,STR_L_cor{x,p});
%         STR_R_filter_T{x,p}= filtfilt(b,a,STR_R_cor{x,p});
%         TRO_L_filter_T{x,p}= filtfilt(b,a,TRO_L_cor{x,p});
%         TRO_R_filter_T{x,p}= filtfilt(b,a,TRO_R_cor{x,p});
%         TRU_L_filter_T{x,p}= filtfilt(b,a,TRU_L_cor{x,p});
%         TRU_R_filter_T{x,p}= filtfilt(b,a,TRU_R_cor{x,p});
%     end 
% end 

%%
% % Hochpassfilter 10Hz
% F_H =10;                                      % [Hz], Frequenz 
% Fs_H =2000;                                   % Abtrastfrequnez
% [b,a]= butter(4,F_T/(Fs_H/2),'high');                  % 4nd Ordnung butterworth filter 
% for p=1:s;
%     for x=1:5;
%         STR_L_filter_T_H{x,p}= filtfilt(b,a,STR_L_filter_T{x,p});
%         STR_R_filter_T_H{x,p}= filtfilt(b,a,STR_R_filter_T{x,p});
%         TRO_L_filter_T_H{x,p}= filtfilt(b,a,TRO_L_filter_T{x,p});
%         TRO_R_filter_T_H{x,p}= filtfilt(b,a,TRO_R_filter_T{x,p});
%         TRU_L_filter_T_H{x,p}= filtfilt(b,a,TRU_L_filter_T{x,p});
%         TRU_R_filter_T_H{x,p}= filtfilt(b,a,TRU_R_filter_T{x,p});
%     end 
% end 
%% 
% Gleichrichten 
for p=1:s;
    for x=1:5;
        STR_L_abs{x,p} = abs(STR_L_cor{x,p});
        STR_R_abs{x,p} = abs(STR_R_cor{x,p});
        TRO_L_abs{x,p} = abs(TRO_L_cor{x,p});
        TRO_R_abs{x,p} = abs(TRO_R_cor{x,p});
        TRU_L_abs{x,p} = abs(TRU_L_cor{x,p});
        TRU_R_abs{x,p} = abs(TRU_R_cor{x,p});
    end 
end 
%%
% Glätten
for p=1:s;
    for x=1:5;
        STR_L_rms{x,p}=rms(STR_L_abs{x,p},250,249,0);
        STR_R_rms{x,p}=rms(STR_R_abs{x,p},250,249,0);
        TRO_L_rms{x,p}=rms(TRO_L_abs{x,p},250,249,0);
        TRO_R_rms{x,p}=rms(TRO_R_abs{x,p},250,249,0);
        TRU_L_rms{x,p}=rms(TRU_L_abs{x,p},250,249,0);
        TRU_R_rms{x,p}=rms(TRU_R_abs{x,p},250,249,0);
    end 
end 
%%
figure(1)
sgtitle('EMG')
subplot(3,2,1);
plot(t{2,3},STR_L{2,3}); title('STR_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{2,3},STR_R{2,3}); title('STR_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{2,3},TRO_L{2,3}); title('TRO_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{2,3},TRO_R{2,3}); title('TRO_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{2,3},TRU_L{2,3}); title('TRU_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{2,3},TRU_R{2,3}); title('TRU_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')

%%
figure(2)
sgtitle('Roh-EMG')
subplot(3,2,1);
plot(t{2,3},STR_L{2,3}); title('STR_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{2,3},STR_R{2,3}); title('STR_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{2,3},TRO_L{2,3}); title('TRO_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{2,3},TRO_R{2,3}); title('TRO_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{2,3},TRU_L{2,3}); title('TRU_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{2,3},TRU_R{2,3}); title('TRU_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
%%
figure(3)
sgtitle('Gleichgerichtetes EMG')
subplot(3,2,1);
plot(t{2,3},STR_L_abs{2,3}); title('STR_Labs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{2,3},STR_R_abs{2,3}); title('STR_Rabs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{2,3},TRO_L_abs{2,3}); title('TRO_Labs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{2,3},TRO_R_abs{2,3}); title('TRO_Rabs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{2,3},TRU_L_abs{2,3}); title('TRU_Labs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{2,3},TRU_R_abs{2,3}); title('TRU_Rabs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')



%%


% F =10;                                      % [Hz], Frequenz 
% Fs =2000;                                   % Abtrastfrequnez
% [b,a]= butter(4,F/(Fs/2));                  % 4nd Ordnung butterworth filter 
% for p=1:s;
%     for x=1:5;
%     STR_L_filter{x,p}= filtfilt(b,a,STR_L_abs{x,p});
%     end 
% end 
% STR_R_filter = filtfilt(b,a,STR_R_abs);
% TRO_L_filter = filtfilt(b,a,TRO_L_abs);
% TRO_R_filter = filtfilt(b,a,TRO_R_abs);
% TRU_L_filter = filtfilt(b,a,TRU_L_abs);
% TRU_R_filter = filtfilt(b,a,TRU_R_abs);

 
%%

%%
% 
% F =30;                                      % [Hz], Frequenz 
% Fs =2000;                                   % Abtrastfrequnez
% 
% STR_L_EKG_filter= highpass(STR_L_abs,30,2000);
% plot(STR_L_EKG_filter)
% 
% 
% 
% 
% 
% 
% %%
% % Filtern 
% F =10;                                      % [Hz], Frequenz 
% Fs =2000;                                   % Abtrastfrequnez
% [b,a]= butter(4,F/(Fs/2));                  % 4nd Ordnung butterworth filter 
% STR_L_filter = filtfilt(b,a,STR_L_abs);
% STR_R_filter = filtfilt(b,a,STR_R_abs);
% TRO_L_filter = filtfilt(b,a,TRO_L_abs);
% TRO_R_filter = filtfilt(b,a,TRO_R_abs);
% TRU_L_filter = filtfilt(b,a,TRU_L_abs);
% TRU_R_filter = filtfilt(b,a,TRU_R_abs);
% 
%% Bandpassfilter;
fpass = [10 500];
fs = 2000; 
for p=1:s;
    for x=1:5;
    STR_L_bandpass{x,p} = bandpass(STR_L_abs{x,p},fpass,fs);
    end 
end 
