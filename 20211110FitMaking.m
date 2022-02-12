for ii=1:length(data)
    clear ReinforcedTrials ProbeTrials TtrialHistR TtrialHistP TrialResponseR TrialResponseP
    ReinforcedTrials=find(BigMat{ii}(1,:)<=2);
    ProbeTrials=find(BigMat{ii}(1,:)>2);
    TtrialHistR=BigMat{ii}(1,ReinforcedTrials);
    TrialResponseR=BigMat{ii}(2,ReinforcedTrials);
    TtrialHistP=BigMat{ii}(1,ProbeTrials);
    TrialResponseP=BigMat{ii}(2,ProbeTrials);
    NumTones(ii,:)=[sum(TtrialHistR==1)/length(TtrialHistR) sum(TtrialHistR==2)/length(TtrialHistR)];
    sliwi=0; %for reinforced

    for sli=1:((max(ReinforcedTrials)-winincrease3)/winincrease3)
    til=find(ReinforcedTrials>sliwi & ReinforcedTrials<=sliwi+winincrease3);
    xtrialsR(ii,sli)=round(nanmean(ReinforcedTrials(til)));
    end
end


%% FIGURE 3: Fits for individual animals and finds plateau performance

colors=[[187, 186, 183]; [255, 195, 185]]./255;

winlength3=100;
winincrease3=100;
winlength4=20;
winincrease4=20;
maxsizemat=500;

clear dp_sli400R dp_sli400P dp_slitouseP
xtrialsP=NaN(length(data),maxsizemat);
xtrialsR=NaN(length(data),maxsizemat);
dp_slitouseP=NaN(length(data),maxsizemat);
dp_sli100P=NaN(length(data),maxsizemat);
dp_sli100R=NaN(length(data),maxsizemat);
c_sli100R=NaN(length(data),maxsizemat);
pHsli100R=NaN(length(data),maxsizemat);
pFAsli100R=NaN(length(data),maxsizemat);
pHsli100P=NaN(length(data),maxsizemat);
pFAsli100P=NaN(length(data),maxsizemat);
maxdpsli100R=NaN(length(data),maxsizemat);
maxcsli100R=NaN(length(data),maxsizemat);
pH_sli100R=NaN(length(data),maxsizemat);
pFA_sli100R=NaN(length(data),maxsizemat);

for ii=1:length(data)
    clear ReinforcedTrials ProbeTrials TtrialHistR TtrialHistP TrialResponseR TrialResponseP
    ReinforcedTrials=find(BigMat{ii}(1,:)<=2);
    ProbeTrials=find(BigMat{ii}(1,:)>2);
    TtrialHistR=BigMat{ii}(1,ReinforcedTrials);
    TrialResponseR=BigMat{ii}(2,ReinforcedTrials);
    TtrialHistP=BigMat{ii}(1,ProbeTrials);
    TrialResponseP=BigMat{ii}(2,ProbeTrials);
    NumTones(ii,:)=[sum(TtrialHistR==1)/length(TtrialHistR) sum(TtrialHistR==2)/length(TtrialHistR)];
    sliwi=0; %for reinforced
    for sli=1:((max(ReinforcedTrials)-winincrease3)/winincrease3)
        clear til
        til=find(ReinforcedTrials>sliwi & ReinforcedTrials<=sliwi+winincrease3);
        tillen{ii}(sli,:)=length(til);
        xtrialsR(ii,sli)=round(nanmean(ReinforcedTrials(til)));
        [maxdpsli100R(ii,sli) maxcsli100R(ii,sli)]=dprime_criterion_calc(sum(TtrialHistR(til)==1), 0, sum(TtrialHistR(til)==2), 0);
        Rhitssli=sum(TrialResponseR(til)==1);
        Rmisssli=sum(TrialResponseR(til)==2);
        Rcrsli=sum(TrialResponseR(til)==3);
        Rfasli=sum(TrialResponseR(til)==4);
        [dp_sli100R(ii,sli), c_sli100R(ii,sli)]=dprime_criterion_calc(Rhitssli, Rmisssli, Rcrsli, Rfasli);
        [bla1,bla2,pH_sli100R(ii,sli), pFA_sli100R(ii,sli)]=dprime_criterion_pH_pFA_calc(Rhitssli, Rmisssli, Rcrsli, Rfasli);
        pHsli100R(ii,sli)=Rhitssli/(Rhitssli+Rmisssli);
        pFAsli100R(ii,sli)=Rfasli/(Rfasli+Rcrsli);
        sliwi=sliwi+winincrease3;
    end
    rrr(ii,:)=min(tillen{ii});
     
    sliwi=0;
    slidwindow4=[1:winlength4]; %for probe
    for sli=1:length(TtrialHistP)/winincrease4
        clear til        
        xtrialsP(ii,sli)=round(nanmean(ProbeTrials(slidwindow4)));
        Phitssli=sum(TrialResponseP(slidwindow4)==1);
        Pmisssli=sum(TrialResponseP(slidwindow4)==2);
        Pcrsli=sum(TrialResponseP(slidwindow4)==3);
        Pfasli=sum(TrialResponseP(slidwindow4)==4);
        [dp_sli100P(ii,sli), c_sli100P(ii,sli)]=dprime_criterion_calc(Phitssli, Pmisssli, Pcrsli, Pfasli);
        pHsli100P(ii,sli)=Phitssli/(Phitssli+Pmisssli);
        pFAsli100P(ii,sli)=Pfasli/(Pfasli+Pcrsli);
        slidwindow4=slidwindow4+winincrease4;
        if sli>1
            if dp_sli100P(ii,sli)>=max(dp_sli100P(ii,:))
                dp_slitouseP(ii,sli)=dp_sli100P(ii,sli);
            else
                dp_slitouseP(ii,sli)=max(dp_sli100P(ii,:));
            end
        else
            dp_slitouseP(ii,sli)=dp_sli100P(ii,sli);
        end
    end
end
%% plot individual fitted d' in one figure
colors=[[187, 186, 183]; [255, 195, 185]]./255;
plotyes=1; %1 if you want to get figure with individual fits, 0 if not
exampleanimalsR=0; %to plot example Reinforced animals
clear xf yf xfplat yfplat yfplathalf xfplathalf mousename
i=0;
exi=0;
% micetoeval=[WTtouseyoung MUTtouseyoung];
micetoeval=[MUTtouseyoung];
fignum=2;
% micetoeval=[WTtouseyoung];
for ii=micetoeval
    i=i+1;
    
    %Reinforced
    xf=[1 xtrialsR(ii,1:maxblock)]'; %x data for fit
    yf=[0 dp_sli100R(ii,1:maxblock)]'; %y data for fit

    %Probe
    xpf=[1 xtrialsP(ii,1:maxblock)]'; %x data for fit
    ypf=[0 dp_slitouseP(ii,1:maxblock)]'; %y data for fit
    
    if sum(~isnan(xf))>=4
        [fitresultR{ii}, gofR{ii}] = mysigmfit(xf, yf);
        %find plateau in reinforced
        clear yplat2 ord closestIndex
        yfit_dp_sli100R{ii}=feval(fitresultR{ii}, xf);
        ord=max(yfit_dp_sli100R{ii})-max(yfit_dp_sli100R{ii})*0.05;
        [minValue,closestIndex]=min(abs(yfit_dp_sli100R{ii}-ord)); %sort(abs(yplat-ordmean))
        if xf(closestIndex)>1000
            xfplat(ii,:)=xf(closestIndex);
            yfplat(ii,:)=yfit_dp_sli100R{ii}(closestIndex);
        else
            vct=xf(~isnan(xf)); vct2=yfit_dp_sli100R{ii}(~isnan(yfit_dp_sli100R{ii}));
            xfplat(ii,:)=vct(end);
            yfplat(ii,:)=vct2(end);
        end
        mousename(ii,:)=char(mouse(ii));

        %find half to plateau
        ordhalf=nanmean(yfit_dp_sli100R{ii});
        [minValuehalf,closestIndexhalf]=min(abs(yfit_dp_sli100R{ii}-ordhalf)); %sort(abs(yplat-ordmean))
        yfplathalf(ii,:)=yfit_dp_sli100R{ii}(closestIndexhalf);
        xfplathalf(ii,:)=xf(closestIndexhalf);        
    end
    if sum(~isnan(xpf))>=4
        %probe
        [fitresultP{ii}, gofP{ii}] = mysigmfit(xpf, ypf);
        yfit_dp_sli100P{ii}=feval(fitresultP{ii}, xpf); 
        
        ordp=max(yfit_dp_sli100P{ii})-max(yfit_dp_sli100P{ii})*0.05;
        [minValuep,closestIndexp]=min(abs(yfit_dp_sli100P{ii}-ordp)); %sort(abs(yplat-ordmean))
        if xpf(closestIndexp)>1000
            xfplatP(ii,:)=xpf(closestIndexp);
            yfplatP(ii,:)=yfit_dp_sli100P{ii}(closestIndexp);
        else
            vctp=xpf(~isnan(xpf)); vct2p=yfit_dp_sli100P{ii}(~isnan(yfit_dp_sli100P{ii}));
            xfplatP(ii,:)=vctp(end);
            yfplatP(ii,:)=vct2p(end);
        end
    else
        yfit_dp_sli100P{ii}=NaN;
    end
    
    if exampleanimalsR==1 && strcmp(mouse(ii),'KF023') || strcmp(mouse(ii),'FZ010') || strcmp(mouse(ii),'AW031') || strcmp(mouse(ii),'AW023')
        gi=figure(i+1);
        exi=exi+1;
        subplot(4,length(micetoeval)/4,exi)
        hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 250, 500]);
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
        set(gca, 'Xtick', [1:4000:8000], 'Xticklabel', [1:4000:8000]);
        xlim([0 10000]); ylim([-1 5]); xlabel('Trials'); ylabel('d'''); title(['(' char(mouse(ii)) ')']); legend off
%         text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
%         text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        mmi=char(mouse(ii));
        dpexamps.(mmi)=[xf yf yfit_dp_sli100R{ii}];
    end

    if plotyes==1
        gi=figure(fignum);
        subplot(4,length(micetoeval)/4,i)
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 553, 736]);
        hold all
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
            plot(xf(2:end),yf(2:end),'k.');
        h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.', 'MarkerEdgeColor', colors(1,:));
%         plot(xfplatP(ii,:), yfplatP(ii,:), 'ro') %for plateau
%         plot(xfplathalf(ii,:), yfplathalf(ii,:), 'go') %for halfplateau       
        set(h1, 'linewidth', 2);
            set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
        set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
       
        if ismember(ii,WTtouseyoung)
            title(['\color{Blue}' char(mouse(ii))])
        else
            title(char(mouse(ii)))
        end
        
        set(gca, 'Xtick', [0:4000:maxblock*100], 'Xticklabel', [0:4000:maxblock*100]);
        xlim([1 maxblock*100])
        ylim([-1 5])
        xlabel('Trials')
        ylabel('d''')
        legend off
        if i==1
            text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
            text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        end     
        
%         figure(1012)
%         subplot(1,2,1)
%         hold all
%         h1=plot(fitresultR{ii},'-',xf(2:end),yf(2:end),'.');
%         subplot(1,2,2)
%         hold all
%         h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.');
    end
end



%%

%% plot individual fitted d' in one figure
colors=[[187, 186, 183]; [255, 195, 185]]./255;
plotyes=1; %1 if you want to get figure with individual fits, 0 if not
exampleanimalsR=0; %to plot example Reinforced animals
clear xf yf xfplat yfplat yfplathalf xfplathalf mousename
i=0;
exi=0;
% micetoeval=[WTtouseyoung MUTtouseyoung];
micetoeval=[MUTtouseyoung];
fignum=2;
% micetoeval=[WTtouseyoung];
for ii=micetoeval
    i=i+1;
    
    %Reinforced
    xf=[1 xtrialsR(ii,1:maxblock)]'; %x data for fit
    yf=[0 dp_sli100R(ii,1:maxblock)]'; %y data for fit

    %Probe
    xpf=[1 xtrialsP(ii,1:maxblock)]'; %x data for fit
    ypf=[0 dp_slitouseP(ii,1:maxblock)]'; %y data for fit
    
    if sum(~isnan(xf))>=4
        [fitresultR{ii}, gofR{ii}] = mysigmfit(xf, yf);
        %find plateau in reinforced
        clear yplat2 ord closestIndex
        yfit_dp_sli100R{ii}=feval(fitresultR{ii}, xf);
        ord=max(yfit_dp_sli100R{ii})-max(yfit_dp_sli100R{ii})*0.05;
        [minValue,closestIndex]=min(abs(yfit_dp_sli100R{ii}-ord)); %sort(abs(yplat-ordmean))
        if xf(closestIndex)>1000
            xfplat(ii,:)=xf(closestIndex);
            yfplat(ii,:)=yfit_dp_sli100R{ii}(closestIndex);
        else
            vct=xf(~isnan(xf)); vct2=yfit_dp_sli100R{ii}(~isnan(yfit_dp_sli100R{ii}));
            xfplat(ii,:)=vct(end);
            yfplat(ii,:)=vct2(end);
        end
        mousename(ii,:)=char(mouse(ii));

        %find half to plateau
        ordhalf=nanmean(yfit_dp_sli100R{ii});
        [minValuehalf,closestIndexhalf]=min(abs(yfit_dp_sli100R{ii}-ordhalf)); %sort(abs(yplat-ordmean))
        yfplathalf(ii,:)=yfit_dp_sli100R{ii}(closestIndexhalf);
        xfplathalf(ii,:)=xf(closestIndexhalf);        
    end
    if sum(~isnan(xpf))>=4
        %probe
        [fitresultP{ii}, gofP{ii}] = mysigmfit(xpf, ypf);
        yfit_dp_sli100P{ii}=feval(fitresultP{ii}, xpf); 
        
        ordp=max(yfit_dp_sli100P{ii})-max(yfit_dp_sli100P{ii})*0.05;
        [minValuep,closestIndexp]=min(abs(yfit_dp_sli100P{ii}-ordp)); %sort(abs(yplat-ordmean))
        if xpf(closestIndexp)>1000
            xfplatP(ii,:)=xpf(closestIndexp);
            yfplatP(ii,:)=yfit_dp_sli100P{ii}(closestIndexp);
        else
            vctp=xpf(~isnan(xpf)); vct2p=yfit_dp_sli100P{ii}(~isnan(yfit_dp_sli100P{ii}));
            xfplatP(ii,:)=vctp(end);
            yfplatP(ii,:)=vct2p(end);
        end
    else
        yfit_dp_sli100P{ii}=NaN;
    end
    
    if exampleanimalsR==1 && strcmp(mouse(ii),'KF023') || strcmp(mouse(ii),'FZ010') || strcmp(mouse(ii),'AW031') || strcmp(mouse(ii),'AW023')
        gi=figure(i+1);
        exi=exi+1;
%         subplot(4,length(micetoeval)/4,exi)
        hold all
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 250, 500]);
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
        set(gca, 'Xtick', [1:4000:8000], 'Xticklabel', [1:4000:8000]);
        xlim([0 10000]); ylim([-1 5]); xlabel('Trials'); ylabel('d'''); title(['(' char(mouse(ii)) ')']); legend off
%         text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
%         text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        mmi=char(mouse(ii));
        dpexamps.(mmi)=[xf yf yfit_dp_sli100R{ii}];
    end

    if plotyes==1
        gi=figure();
%         subplot(4,length(micetoeval)/4,i)
        set(gca,'TickDir','out','FontSize',12,'fontname','arial');
        set(gi, 'position', [1547, 7, 553, 736]);
        hold all
        h1=plot(fitresultR{ii},'k-',xf(2:end),yf(2:end),'k.');
            plot(xf(2:end),yf(2:end),'k.');
        h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.', 'MarkerEdgeColor', colors(1,:));
%         plot(xfplatP(ii,:), yfplatP(ii,:), 'ro') %for plateau
%         plot(xfplathalf(ii,:), yfplathalf(ii,:), 'go') %for halfplateau       
        set(h1, 'linewidth', 2);
            set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
        set(h2, {'color'}, num2cell(colors(1,:), 2), 'linewidth', 2);
       
        if ismember(ii,WTtouseyoung)
            title(['\color{Blue}' char(mouse(ii))])
        else
            title(char(mouse(ii)))
        end
        
        set(gca, 'Xtick', [0:4000:maxblock*100], 'Xticklabel', [0:4000:maxblock*100]);
        xlim([1 maxblock*100])
        ylim([-1 5])
        xlabel('Trials')
        ylabel('d''')
        legend off
        if i==1
            text([300 300],[4.5 4.5], 'Reinforced', 'Color', 'k','FontSize', 12, 'fontname','arial')
            text([300 300],[3.8 3.8], 'Probe', 'Color', colors(1,:),'FontSize', 12, 'fontname','arial')
        end     
        
%         figure(1012)
%         subplot(1,2,1)
%         hold all
%         h1=plot(fitresultR{ii},'-',xf(2:end),yf(2:end),'.');
%         subplot(1,2,2)
%         hold all
%         h2=plot(fitresultP{ii},'-'); plot(xpf(2:end),ypf(2:end),'.');
    end
end
