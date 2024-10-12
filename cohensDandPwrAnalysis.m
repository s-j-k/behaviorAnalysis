%% cohens D test    
% compute cohen's d
% medium effect of 0.5 is visible to the naked eye
cohensDTest={}; % this will give you effect sizes we see in the pilot data
% using the mean of each session
cohensDTest{1,1}={'Animal Name'};
cohensDTest{1,2}={'MGB Full Opto Hit'};
cohensDTest{1,3}={'MGB Full Opto FA'};
cohensDTest{1,4}={'MGB Tone Opto Hit'};
cohensDTest{1,5}={'MGB Tone Opto FA'};
cohensDTest{1,6}={'MGB Choice Opto Hit'};
cohensDTest{1,7}={'MGB Choice Opto FA'};
pp=2;

% correct for high and low values
% newBinRates is uncorrected, binRates is corrected

% newBinRates=binRates;
uu=2;yy=3;
for uu=2:size(binRates,1)
    for yy=3:size(binRates,2)
        dimSub=length(binRates{uu,yy});
        if find(cell2mat(binRates(uu,yy))==1)
            idx=find(cell2mat(binRates(uu,yy))==1);
            binRates{uu,yy}(idx)=(dimSub-1)/dimSub;
        elseif find(cell2mat(binRates(uu,yy))==0)
            idx=find(cell2mat(binRates(uu,yy))==0);
            binRates{uu,yy}(idx)=1/dimSub;
        else
        end
    end
end

for tt=2:size(binRates,1)
    cohensDTest{pp,1}=binRates{tt,1};
    cohensDTest{pp,2}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,5}))),...
        binRates{tt,5},'paired'); % Full Opto H
    mgbRHitStd=std(binRates{tt,3}(1:(length(binRates{tt,5}))));
    cohensDTest{pp,3}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,6}))),...
        binRates{tt,6},'paired'); % Full Opto FA
    mgbRFaStd=std(binRates{tt,4}(1:(length(binRates{tt,6}))));
    cohensDTest{pp,4}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,7}))),...
        binRates{tt,7},'paired'); % Tone Opto H
    cohensDTest{pp,5}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,8}))),...
        binRates{tt,8},'paired'); % Tone Opto FA
    cohensDTest{pp,6}=computeCohen_d(binRates{tt,3}(1:(length(binRates{tt,9}))),...
        binRates{tt,9},'paired'); % Choice Opto H
    cohensDTest{pp,7}=computeCohen_d(binRates{tt,4}(1:(length(binRates{tt,10}))),...
        binRates{tt,10},'paired'); % Choice Opto FA
    pp=pp+1;
end
%%    % how many samples do i need for each value??? % Inputs for sample size est.
nn=1:20; count=-1;
tt=2;
figure(tt);
for tt=2:size(binRates,1)
     %hits first
    subplot(3,3,tt+count);
    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
        mean(binRates{tt,5}),0.80); % % (type of test, [mean1 stdev], mean2, power)
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
        mean(binRates{tt,5}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),[],nn);

    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,5})))) mgbRHitStd], ...
    mean(binRates{tt,5}),[],nn);
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');ylim([0.65 1]);
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    if tt==4
        xlim([0 10]);
    else
    end
    title([char(binRates{tt,1}) 'MGB Hit Full Power vs. Sample size']);
    xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,7})))) mgbRHitStd], ...
        mean(binRates{tt,7}),[],nn);
    count=count+1; subplot(3,3,tt+count);plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');xlim([0 10]);
    title('Hit Tone Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');

    nout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,3}(1:(length(binRates{tt,9})))) mgbRHitStd], ...
        mean(binRates{tt,9}),[],nn); 
    count=count+1;subplot(3,3,tt+count); plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');xlim([0 10]);
    title('Hit Choice Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
end

nn=1:55; count=-1;
figure(tt+1);
for tt=2:size(binRates,1)
    subplot(3,3,tt+count);
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.80);
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
    mean(binRates{tt,6}),[],nn); 
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');ylim([0.65 1]);
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    title([char(binRates{tt,1}) 'MGB FA Full Power vs. Sample size']);
    xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,8}),[],nn);
    count=count+1; subplot(3,3,tt+count);plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');
    title('MGB FA Tone Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
    
    nout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.80); 
    pwrout8=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn); 
    nout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.90); 
    pwrout9=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn);
    nout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),0.95); 
    pwrout95=sampsizepwr('t',[mean(binRates{tt,4}(10:(length(binRates{tt,6})+10))) mgbRFaStd], ...
        mean(binRates{tt,10}),[],nn); 
    count=count+1;subplot(3,3,tt+count); 
    plot(nn,pwrout8,'b-',nout8,0.8,'ro');hold on;
    plot(nn,pwrout9,'b-',nout9,0.9,'ro');
    plot(nn,pwrout95,'b-',nout95,0.95,'ro');
    title('MGB FA Choice Power vs. Sample size');xlabel('sample size (bins of 5 trials)');ylabel('power');
end


%%

% for a given animal, take effect sizes from cohensDTest variable and then
% calculate the sample sizes needed. 
for uu=2:size(cohensDTest,1)-2
    effect_sizes = [cohensDTest{uu,2}, cohensDTest{uu,3},cohensDTest{uu,4},cohensDTest{uu,5},...
        cohensDTest{uu,6},cohensDTest{uu,7}];  
    effect_str   = {'MGB Full Hit', 'MGB Full FA','MGB Tone Hit','MGB Tone FA',...
        'MGB Choice Hit','MGB Choice FA'};
    alpha        = 0.05; % 0.05, 0.01, 0.005, 0.001                   
    stat_power   = 0.85; % 0.80, 0.85, 0.90                    

    % Sample size estimates
    sample_sizes = arrayfun(@(effect) calculateSampleSize(effect, alpha, stat_power), effect_sizes);
%     sample_sizes = arrayfun(@(effect) sampsizepwr('t2', [effect 1], alpha, stat_power), effect_sizes);

    for i = 1:length(effect_sizes)
        fprintf('For shift between %s conditions an effect size of %.4f requires a sample size of %d\n', effect_str{i}, effect_sizes(i), ceil(sample_sizes(i)));
    end

end
