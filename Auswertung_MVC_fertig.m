%% Auswertung MVC-Daten 
clear all; clc;
%% Daten laden
Tabelle = xlsread('MVC_Trials.xlsx.');
MVC_Nr = Tabelle(2:9,1:21);
%%
s = 1;          % Anzahl der Probanden 
for p=1:s;      % Anzahl der Probanden 
    for x=1:8;
    if MVC_Nr(x,p)<10
    data(x,p)=load(['C:\Users\Lisa Broß\Documents\6. Semester\FKFS\MVC_Daten\Proband' num2str(p) '\Versuch0000000' num2str(MVC_Nr(x,p)) '.mat']); % Proband 1-> Versuch 10 (1.)
    else data(x,p)=load(['C:\Users\Lisa Broß\Documents\6. Semester\FKFS\MVC_Daten\Proband' num2str(p) '\Versuch000000' num2str(MVC_Nr(x,p)) '.mat']); 
    end
    end 
end

%% Daten auswählen  
for p=1:s;
    for x=1:8;
        STR_L{x,p}=data(x,p).data(:,1);    % [mV],O-EMG M. Sternocleidomastoideus links
        STR_R{x,p}=data(x,p).data(:,2);    % [mV],O-EMG M. Sternocleidomastoideus rechts
        TRO_L{x,p}=data(x,p).data(:,3);    % [mV],O-EMG M.Trapezius pars descendens links
        TRO_R{x,p}=data(x,p).data(:,4);    % [mV],O-EMG M.Trapezius pars descendens rechts
        TRU_L{x,p}=data(x,p).data(:,5);    % [mV],O-EMG M.Trapezius pars transversa links
        TRU_R{x,p}=data(x,p).data(:,6);    % [mV],O-EMG M.Trapezius pars transversa rechts
        EKG{x,p}=data(x,p).data(:,7);           % [mV],EKG
        t{x,p} = (0:0.0005:(length(data(x,p).data)/2000)-0.0005);     % [s], Zeit
    end
end 
%%
x=1;
p=1;
figure(1)
sgtitle('Roh-EMG')
subplot(3,2,1);
plot(t{x,p},STR_L{x,p}); title('STR_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{x,p},STR_R{x,p}); title('STR_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{x,p},TRO_L{x,p}); title('TRO_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{x,p},TRO_R{x,p}); title('TRO_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{x,p},TRU_L{x,p}); title('TRU_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{x,p},TRU_R{x,p}); title('TRU_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')

%% Signalverarbeitung
% EMG-Signal nullen 
for p=1:s;
    for x=1:8;
        STR_L_null{x,p} = STR_L{x,p}-mean(STR_L{x,p}(1:1000,1));
        STR_R_null{x,p} = STR_R{x,p}-mean(STR_R{x,p}(1:1000,1));
        TRO_L_null{x,p} = TRO_L{x,p}-mean(TRO_L{x,p}(1:1000,1));
        TRO_R_null{x,p} = TRO_R{x,p}-mean(TRO_R{x,p}(1:1000,1));
        TRU_L_null{x,p} = TRU_L{x,p}-mean(TRU_L{x,p}(1:1000,1));
        TRU_R_null{x,p} = TRU_R{x,p}-mean(TRU_R{x,p}(1:1000,1));
    end 
end
%%
% Gleichrichten roh EMG 
for p=1:s;
    for x=1:8;
        STR_L_abs_roh{x,p} = abs(STR_L_null{x,p});
        STR_R_abs_roh{x,p} = abs(STR_R_null{x,p});
        TRO_L_abs_roh{x,p} = abs(TRO_L_null{x,p});
        TRO_R_abs_roh{x,p} = abs(TRO_R_null{x,p});
        TRU_L_abs_roh{x,p} = abs(TRU_L_null{x,p});
        TRU_R_abs_roh{x,p} = abs(TRU_R_null{x,p});
    end 
end
%%
x=1;
p=1;
figure(2)
sgtitle('genulltes EMG')
subplot(3,2,1);
plot(t{x,p},STR_L_null{x,p}); title('STR_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{x,p},STR_R_null{x,p}); title('STR_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{x,p},TRO_L_null{x,p}); title('TRO_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{x,p},TRO_R_null{x,p}); title('TRO_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{x,p},TRU_L_null{x,p}); title('TRU_L'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{x,p},TRU_R_null{x,p}); title('TRU_R'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
%% 
% R Zacken des EKG Signals extrahieren
for p=1:s;
    for x=1:8;
    peaks{x,p}=[]; locs{x,p}=[];vek{x,p}=[];
    [peaks{x,p}(:,1) locs{x,p}(:,1)]=findpeaks(EKG{x,p},'MinPeakDistance',1000, 'MinPeakHeight',1); % Diese Zeile müsste auf Basisi des EKG Signals bei den einzelnen Probanden angepasst werden

    % Bau des Vektors zur Detektion des EKG Signals
    le{x,p}=length(EKG{x,p});
    vek{x,p}(1:le{x,p})=0;

    for j=1:length(locs{x,p})
        k=locs{x,p}(j);
        if k-80<0
           k=locs{x,p}(j+1);
        end
        vek{x,p}(k-80:k+80)=1;        % Fenstebreite des EKG Signals %140, 70
                                      % Diese Zeile müsste auf Basis des EKG Signals bei den einzelnen Probanden angepasst werden
    end
    end 
end
%%
x=1;
p=1;
figure(3)
title('EKG mit Peak')
plot(EKG{x,p},'b')
hold on
plot(locs{x,p},peaks{x,p},'r*')
xlabel('Samples'), ylabel('Aktivität [mV]')
%% 
x=1;
p=1;
figure(4)
plot(t{x,p},TRO_L_abs_roh{x,p})
hold on 
plot(t{x,p},vek{x,p})

%%
% Filtern des EMG Signals mit einem Butterworth Tiefpassfilter. 
fc=30; n=4; fs=2000; form='low';
for p=1:s;
    for x=1:8;
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
    for x=1:8;
        for j=1:length(STR_L{x,p});
            STR_L_cor{x,p}(j,:)=STR_L_null{x,p}(j)-STR_L_FiltSignal{x,p}(j)*(vek{x,p}(j))';
            STR_R_cor{x,p}(j,:)=STR_R_null{x,p}(j)-STR_R_FiltSignal{x,p}(j)*(vek{x,p}(j))';
            TRO_L_cor{x,p}(j,:)=TRO_L_null{x,p}(j)-TRO_L_FiltSignal{x,p}(j)*(vek{x,p}(j))';
            TRO_R_cor{x,p}(j,:)=TRO_R_null{x,p}(j)-TRO_R_FiltSignal{x,p}(j)*(vek{x,p}(j))';
            TRU_L_cor{x,p}(j,:)=TRU_L_null{x,p}(j)-TRU_L_FiltSignal{x,p}(j)*(vek{x,p}(j))';
            TRU_R_cor{x,p}(j,:)=TRU_R_null{x,p}(j)-TRU_R_FiltSignal{x,p}(j)*(vek{x,p}(j))';
        end 
    end 
end
%%
% Gleichrichten 
for p=1:s;
    for x=1:8;
        STR_L_abs{x,p} = abs(STR_L_cor{x,p});
        STR_R_abs{x,p} = abs(STR_R_cor{x,p});
        TRO_L_abs{x,p} = abs(TRO_L_cor{x,p});
        TRO_R_abs{x,p} = abs(TRO_R_cor{x,p});
        TRU_L_abs{x,p} = abs(TRU_L_cor{x,p});
        TRU_R_abs{x,p} = abs(TRU_R_cor{x,p});
    end 
end 
%%
% Tiefpassfilter 20Hz
F_T =20;                                      % [Hz], Frequenz 
Fs_T =2000;                                    % Abtrastfrequnez
[b,a]= butter(4,F_T/(Fs_T/2),'low');                  % 4nd Ordnung butterworth filter 
for p=1:s;
    for x=1:8;
        STR_L_filter_T{x,p}= filtfilt(b,a,STR_L_abs{x,p});
        STR_R_filter_T{x,p}= filtfilt(b,a,STR_R_abs{x,p});
        TRO_L_filter_T{x,p}= filtfilt(b,a,TRO_L_abs{x,p});
        TRO_R_filter_T{x,p}= filtfilt(b,a,TRO_R_abs{x,p});
        TRU_L_filter_T{x,p}= filtfilt(b,a,TRU_L_abs{x,p});
        TRU_R_filter_T{x,p}= filtfilt(b,a,TRU_R_abs{x,p});
    end 
end 

%%
x=1;
p=1;
figure(5)
sgtitle('Gefiltertes EMG und roh-EMG')
subplot(3,2,1);
plot(t{x,p},STR_L_abs_roh{x,p},'b',t{x,p},STR_L_filter_T{x,p},'r'); title('STR_L abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,2);
plot(t{x,p},STR_R_abs_roh{x,p},'b',t{x,p},STR_R_filter_T{x,p},'r'); title('STR_R abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,3);
plot(t{x,p},TRO_L_abs_roh{x,p},'b',t{x,p},TRO_L_filter_T{x,p},'r'); title('TRO_L abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,4);
plot(t{x,p},TRO_R_abs_roh{x,p},'b',t{x,p},TRO_R_filter_T{x,p},'r'); title('TRO_R abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,5);
plot(t{x,p},TRU_L_abs_roh{x,p},'b',t{x,p},TRU_L_filter_T{x,p},'r'); title('TRU_L abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
subplot(3,2,6);
plot(t{x,p},TRU_R_abs_roh{x,p},'b',t{x,p},TRU_R_filter_T{x,p},'r'); title('TRU_R abs'),xlabel('Zeit[s]'), ylabel('EMG[mV]')
legend('roh-EMG','gefiltertes ENG')
%% maximale EMG-Amplitude jedes Muskels 
for p=1:s;
    for x=1:8
    STR_L_max{x,p}=max(STR_L_filter_T{x,p}(:,1));
    STR_R_max{x,p}=max(STR_R_filter_T{x,p}(:,1));
    TRO_L_max{x,p}=max(TRO_L_filter_T{x,p}(:,1));
    TRO_R_max{x,p}=max(TRO_R_filter_T{x,p}(:,1));
    TRU_L_max{x,p}=max(TRU_L_filter_T{x,p}(:,1));
    TRU_R_max{x,p}=max(TRU_R_filter_T{x,p}(:,1));
    end 
end

%% Cell zu Double 
STR_L_max_d = cell2mat(STR_L_max);
STR_R_max_d = cell2mat(STR_R_max);
TRO_L_max_d = cell2mat(TRO_L_max);
TRO_R_max_d = cell2mat(TRO_R_max);
TRU_L_max_d = cell2mat(TRU_L_max);
TRU_R_max_d = cell2mat(TRU_R_max);

%% Mittelwert von den zwei Werten je Richtung
for p=1:s;
        STR_L_max_mean(:,p) = arrayfun(@(i) mean(STR_L_max_d(i:i+1,p)),1:2:(length(STR_L_max_d)-1))';
        STR_R_max_mean(:,p) = arrayfun(@(i) mean(STR_R_max_d(i:i+1,p)),1:2:(length(STR_R_max_d)-1))';
        TRO_L_max_mean(:,p) = arrayfun(@(i) mean(TRO_L_max_d(i:i+1,p)),1:2:(length(TRO_L_max_d)-1))';
        TRO_R_max_mean(:,p) = arrayfun(@(i) mean(TRO_R_max_d(i:i+1,p)),1:2:(length(TRO_R_max_d)-1))';
        TRU_L_max_mean(:,p) = arrayfun(@(i) mean(TRU_L_max_d(i:i+1,p)),1:2:(length(TRU_L_max_d)-1))';
        TRU_R_max_mean(:,p) = arrayfun(@(i) mean(TRU_R_max_d(i:i+1,p)),1:2:(length(TRU_R_max_d)-1))';
end 
%% Maximale MVC jedes Muskels aus allen vier Richtungen 
for p=1:s;
        STR_L_max_max(p) = max(STR_L_max_mean(:,p));
        STR_R_max_max(p) = max(STR_R_max_mean(:,p));
        TRO_L_max_max(p) = max(TRO_L_max_mean(:,p));
        TRO_R_max_max(p) = max(TRO_R_max_mean(:,p));
        TRU_L_max_max(p) = max(TRU_L_max_mean(:,p));
        TRU_R_max_max(p) = max(TRU_R_max_mean(:,p));
end 

%% Tabelle erstellen 
for p=1:s;
    MVC_Proband(:,p) = [STR_L_max_max(:,p); STR_R_max_max(:,p);TRO_L_max_max(:,p);TRO_R_max_max(:,p);TRU_L_max_max(:,p);TRU_R_max_max(:,p)];
end  
%%
Muskeln = ["STR_L";"STR_R";"TRO_L";"TRO_R";"TRU_L";"TRU_R"];
Probanden = ['1','2'];
MVC = table(['STR_L';'STR_R';'TRO_L';'TRO_R';'TRU_L';'TRU_R'],MVC_Proband);
%
save('C:\Users\Lisa Broß\Documents\6. Semester\FKFS\MVC_Daten\data_Output_MVC\data_Output.mat','MVC_Proband');
