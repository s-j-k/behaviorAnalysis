%%To run this script, you need the following code:
% bpod_data_preprocessing(rawdata_folder
% 

% datapath='T:\su\DATA\behaviorData\';
datapath='S:\su\';
% cohortname='BPOD4';
cohortname='BPOD4';
dpath = [datapath cohortname '\'];
folders=dir(fullfile(dpath,'SK*'));
Animals={{folders.name}};Animals=Animals{1,1};Animals=[Animals];
Animals=string(Animals);
rawdata_folder=dpath;
% GNG_protocol='LickingGNG_SJK'; %'LickingGNG_Reverse_cam';
GNG_protocol='GNG_LaserSineWavePhasicTone_Celine_ZZAW';
data=bpod_data_preprocessing(rawdata_folder,GNG_protocol,Animals);
GNG_protocol='LGNG_SJK';
data1=bpod_data_preprocessing(rawdata_folder,GNG_protocol,Animals);
data2=data;
data=[data data1];
%%
% plot the d' over time for each animal
figure;
n=2;
for ii=1:(length(data)/(n-1))
    dprime1=[data(ii).dprime data(ii+n).dprime];
    dprime1=cell2mat(dprime1);
    subplot(1,(length(data)/(n-1)),ii);
            hold on;
    plot(1:length(dprime1),dprime1(1,:));
    plot(1:length(dprime1),dprime1(2,:))
    legend('Reinforced','Probe')
    xlabel('days of training');ylabel('d`');
    title(data(1,ii).mouse)
    hold off
end
%% if you ran one protocol for everyone
        for ii=1:(length(data))
        dprime1=[data(ii).dprime];
        dprime1=[data(ii).dprime];
        dprime1=cell2mat(dprime1);
        subplot(1,(length(data)),ii);
                hold on;
        plot(1:length(dprime1),dprime1(1,:));
        plot(1:length(dprime1),dprime1(2,:))
        legend('Reinforced','Probe')
        xlabel('days of training');ylabel('d`');
        title(data(1,ii).mouse)
        hold off
        end
%%

% now look at the first lick time
figure(10);
hold on;
for jj=1:length(data)
% for jj=1 % 1 is SK09
    for kk=1:length(data(jj).firstlickonset)
        subplot(2,length(data)/2,kk);
        title(['Day ' num2str(kk) ]);
        hold on
%         numlicks=length(data(jj).firstlickonset{1,kk}.RT);
%         scatter(1:numlicks,data(jj).firstlickonset{1,kk}.RT); % for target cues
        numlicks=length(data(jj).firstlickonset{1,kk}.RF);
        scatter(1:numlicks,data(jj).firstlickonset{1,kk}.RF);
        hold off
    end
end


