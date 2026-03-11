function lickLatOptoPerAnimal2(mgbTempTestsOnly,allDataTestsOnly,icDataTestsOnly,allDataCtlOnly,reinfcolor,optocolor)

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;
allDataTestsOnly(:,42:50)=[];
allDataTestsOnly2=vertcat(allDataTestsOnly,icDataTestsOnly);
allDataTestsOnly2(8,:)=[];
% tests for MGB for each animal
nbsubj=2;
for nbsubj=2:size(allDataTestsOnly,1)
    matrix=allDataTestsOnly{nbsubj,26};
    for i=1:max(matrix(:,SESS))
        mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
        sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
        mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
        mo_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        so_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
        mto_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        sto_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
        mco_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        sco_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
        
        mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
        mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
        mo_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        so_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
        mto_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        sto_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
        mco_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        sco_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
        
        mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
                        
        mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
        mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
        mlo_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        slo_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
        mlto_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        slto_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        mlco_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        slco_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        
        mcco_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==4 & matrix(:,OUTCOME)==3,LICKL));
        scco_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==4 & matrix(:,OUTCOME)==3,LICKL));
        mlcco_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==4 & matrix(:,SESS)==i,LICKR));
        slcco_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==4 & matrix(:,SESS)==i,LICKR));
        
    end
    switch nbsubj
        case 2 % sk163
            mgbDays = [0 0 0 0 1 1 0 0 1 0 1 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 1 1 0 0 1 1 0 0 0 1 1];icDays=logical(icDays);
            expRange=1:11; 
            mgbDays=logical(mgbDays);icDays=logical(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        case 3 %sk175
        expRange=22:32; % cohort 2
            mgbDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                1 1 0 0 1 1 0 0 0 0 0];mgbDays=logical(mgbDays);
            icDays = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 1 1 0 0 1 1 ...
                0 0 0];
            icDays=logical(icDays); 
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            MGBmr_hit(MGBmr_hit==0)=NaN;
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICmr_hit(ICmr_hit==0)=NaN;
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays);            
            
        case 4 % sk176
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                0 0 0 1 0 1 1 0 0 0 0 1];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ... 
                0 1 0 0 1 0 0 1 0 1 1 0];
            icDays=logical(icDays);
            expRange=21:length(mgbDays); % cohort 3 sk177
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            MGBmr_hit(MGBmr_hit==0)=NaN;
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICmr_hit(ICmr_hit==0)=NaN;
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays);            

        case 5 % sk198
            mgbDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                1 1 0 1 0 1 0 ...
                0 0 0]; %... % catch mgb, catch ic, lob
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                1 0 1 0 1 0 1 ...
                0 0 0]; %...
            icDays=logical(icDays);
            expRange=38:length(mgbDays);    
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        
        case 6 % sk203
            mgbDays=[0 0 0 0 0 0 0 0 ...
                1 1 0 0 0 1 1 0 0 ...
                0 0 0]; % catch MGB, catch IC, LOB
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 ...
                0 0 0 1 1 0 0 1 1 ...
                 0 0 0]; % catch MGB, catch IC, LOB
            icDays=logical(icDays);
            expRange=9:length(mgbDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
    
        case 7 % sk204
            mgbDays=[0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ... %15
                0 0 0 0 0 0 0 0 0 ...
                0 0 ... %26
                0 0 1 1 0 1 0 1 ...
                0 0 ...  % no opto, then ic tone
                0 0 0 0]; % catch ic, catch mgb, lob
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 ...
                0 0 ...
                1 0 0 0 1 0 1 0 ...
                0 1 ...
                0 0 1 0]; % catch ic, catch mgb, ic tone, lob
            icDays=logical(icDays);
            expRange=27:length(mgbDays); 
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        case 8 % sk258
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ... 
                0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0]; % catch, fc
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ... 
                0 0 0 0 0 0 0 0 1 1 1 1 0 0 1]; 
            icDays=logical(icDays);
            expRange=24:length(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        case 9 % sk259
            mgbDays=[0 0 0 0 0 0 0 0 0 0 ... 
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0]; % catch, fc
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 ...
                0 0 0 0 0 0 0 0 0 0 ...
                0 0 0 0 0 0 0 1 1 0 ...
                1 0 1 1]; % catch
            icDays=logical(icDays);
            expRange=23:length(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        case 10 % sk260
            mgbDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 0 ... 
                0 0 0 0 0 ... % forgot task?
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 0]; % catch, fc
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 0 ... 
                0 0 0 0 0 ... % forgot task?
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...% 45
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                1 1 1 0 1 ...
                1 0];
            icDays=logical(icDays);
            expRange=50:length(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
        case 11 % sk262
            mgbDays=[0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0];
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 1 1 1 ...
                1];
            icDays=logical(icDays);
            expRange=24:length(icDays); % Start with 60dB 
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
    case 12 % sk276
            mgbDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 0 ...
                0 0 0 0 0 ...
                0 0]; % catch, fc
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 0 ...
                0 1 1 0 1 ...
                1 0]; % catch, fc
            icDays=logical(icDays);
            expRange=17:length(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
    case 13 % sk278
           mgbDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 0 ...
                0 0 0 0 0]; % catch, fc
            mgbDays=logical(mgbDays);
            icDays=[0 0 0 0 0 ...
                0 0 0 0 0 ... 
                0 0 0 0 1 ...
                1 1 0 1 0 ...
                0]; % last one is catch
            icDays=logical(icDays);
            expRange=15:length(icDays);
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
            MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
            MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
            MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
            MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
            MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
            MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
            MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
            MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
            MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
            MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
            
            MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
            MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
            MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
            MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
            MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
            MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
            MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
            MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
            MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
            MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);
            
            MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
            MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
            MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
            MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
            MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
            MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
            MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
            MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
            MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
            MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

            MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
            MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
            MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
            MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
            MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
            MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
            MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
            MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
            MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
            MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);
            
            % do the IC
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); ICsr_hit(ICsr_hit==0)=NaN;
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays); ICmp_hit(ICmp_hit==0)=NaN;
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays); ICsp_hit(ICsp_hit==0)=NaN;
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays); ICmo_hit(ICmo_hit==0)=NaN;
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays); ICso_hit(ICso_hit==0)=NaN;
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);ICmto_hit(ICmto_hit==0)=NaN;
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays); ICsto_hit(ICsto_hit==0)=NaN;
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays); ICmco_hit(ICmco_hit==0)=NaN;
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays); ICsco_hit(ICsco_hit==0)=NaN;
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            ICmr_hit(nbsubj,icDays) = mr_hit(nbsubj,icDays); % reinf hit
            ICsr_hit(nbsubj,icDays) = sr_hit(nbsubj,icDays); 
            ICmp_hit(nbsubj,icDays) = mp_hit(nbsubj,icDays);
            ICsp_hit(nbsubj,icDays) = sp_hit(nbsubj,icDays);
            ICmo_hit(nbsubj,icDays) = mo_hit(nbsubj,icDays);
            ICso_hit(nbsubj,icDays) = so_hit(nbsubj,icDays);
            ICmto_hit(nbsubj,icDays) = mto_hit(nbsubj,icDays);
            ICsto_hit(nbsubj,icDays) = sto_hit(nbsubj,icDays);
            ICmco_hit(nbsubj,icDays) = mco_hit(nbsubj,icDays);
            ICsco_hit(nbsubj,icDays) = sco_hit(nbsubj,icDays);
            
            ICmr_fa(nbsubj,icDays) = mr_fa(nbsubj,icDays);
            ICsr_fa(nbsubj,icDays) = sr_fa(nbsubj,icDays);
            ICmp_fa(nbsubj,icDays) = mp_fa(nbsubj,icDays);
            ICsp_fa(nbsubj,icDays) = sp_fa(nbsubj,icDays);
            ICmo_fa(nbsubj,icDays) = mo_fa(nbsubj,icDays);
            ICso_fa(nbsubj,icDays) = so_fa(nbsubj,icDays);
            ICmto_fa(nbsubj,icDays) = mto_fa(nbsubj,icDays);
            ICsto_fa(nbsubj,icDays) = sto_fa(nbsubj,icDays);
            ICmco_fa(nbsubj,icDays) = mco_fa(nbsubj,icDays);
            ICsco_fa(nbsubj,icDays) = sco_fa(nbsubj,icDays);
            
            ICmlr_hit(nbsubj,icDays) = mlr_hit(nbsubj,icDays);
            ICslr_hit(nbsubj,icDays) = slr_hit(nbsubj,icDays);
            ICmlp_hit(nbsubj,icDays) = mlp_hit(nbsubj,icDays);
            ICslp_hit(nbsubj,icDays) = slp_hit(nbsubj,icDays);
            ICmlo_hit(nbsubj,icDays) = mlo_hit(nbsubj,icDays);
            ICslo_hit(nbsubj,icDays) = slo_hit(nbsubj,icDays);
            ICmlto_hit(nbsubj,icDays) = mlto_hit(nbsubj,icDays);
            ICslto_hit(nbsubj,icDays) = slto_hit(nbsubj,icDays);
            ICmlco_hit(nbsubj,icDays) = mlco_hit(nbsubj,icDays);
            ICslco_hit(nbsubj,icDays) = slco_hit(nbsubj,icDays);

            ICmlr_fa(nbsubj,icDays) = mlr_fa(nbsubj,icDays);
            ICslr_fa(nbsubj,icDays) = slr_fa(nbsubj,icDays);
            ICmlp_fa(nbsubj,icDays) = mlp_fa(nbsubj,icDays);
            ICslp_fa(nbsubj,icDays) = slp_fa(nbsubj,icDays);
            ICmlo_fa(nbsubj,icDays) = mlo_fa(nbsubj,icDays);
            ICslo_fa(nbsubj,icDays) = slo_fa(nbsubj,icDays);
            ICmlto_fa(nbsubj,icDays) = mlto_fa(nbsubj,icDays);
            ICslto_fa(nbsubj,icDays) = slto_fa(nbsubj,icDays);
            ICmlco_fa(nbsubj,icDays) = mlco_fa(nbsubj,icDays);
            ICslco_fa(nbsubj,icDays) = slco_fa(nbsubj,icDays); 
    end
end

MGBmr_hit(MGBmr_hit==0)=NaN;
MGBsr_hit(MGBsr_hit==0)=NaN;
MGBmp_hit(MGBmp_hit==0)=NaN;
MGBsp_hit(MGBsp_hit==0)=NaN;
MGBmo_hit(MGBmo_hit==0)=NaN;
MGBso_hit(MGBso_hit==0)=NaN;
MGBmto_hit(MGBmto_hit==0)=NaN;
MGBsto_hit(MGBsto_hit==0)=NaN;
MGBmco_hit(MGBmco_hit==0)=NaN;
MGBsco_hit(MGBsco_hit==0)=NaN;

MGBmr_fa(MGBmr_fa==0) = NaN;
MGBsr_fa(MGBsr_fa==0) = NaN;
MGBmp_fa(MGBmp_fa==0) = NaN;
MGBsp_fa(MGBsp_fa==0) = NaN;
MGBmo_fa(MGBmo_fa==0) = NaN;
MGBso_fa(MGBso_fa==0) = NaN;
MGBmto_fa(MGBmto_fa==0) = NaN;
MGBsto_fa(MGBsto_fa==0) = NaN;
MGBmco_fa(MGBmco_fa==0) = NaN;
MGBsco_fa(MGBsco_fa==0) = NaN;

MGBmlr_hit(MGBmlr_hit==0) = NaN;
MGBslr_hit(MGBslr_hit==0) = NaN;
MGBmlp_hit(MGBmlp_hit==0) = NaN;
MGBslp_hit(MGBslp_hit==0) = NaN;
MGBmlo_hit(MGBmlo_hit==0) = NaN;
MGBslo_hit(MGBslo_hit==0) = NaN;
MGBmlto_hit(MGBmlto_hit==0) = NaN;
MGBslto_hit(MGBslto_hit==0) = NaN;
MGBmlco_hit(MGBmlco_hit==0) = NaN;
MGBslco_hit(MGBslco_hit==0) = NaN;

MGBmlr_fa(MGBmlr_fa==0) = NaN;
MGBslr_fa(MGBslr_fa==0) = NaN;
MGBmlp_fa(MGBmlp_fa==0) = NaN;
MGBslp_fa(MGBslp_fa==0) = NaN;
MGBmlo_fa(MGBmlo_fa==0) = NaN;
MGBslo_fa(MGBslo_fa==0) = NaN;
MGBmlto_fa(MGBmlto_fa==0) = NaN;
MGBslto_fa(MGBslto_fa==0) = NaN;
MGBmlco_fa(MGBmlco_fa==0) = NaN;
MGBslco_fa(MGBslco_fa==0) = NaN;

ICmr_hit(ICmr_hit==0)=NaN;
ICsr_hit(ICsr_hit==0)=NaN;
ICmp_hit(ICmp_hit==0)=NaN;
ICsp_hit(ICsp_hit==0)=NaN;
ICmo_hit(ICmo_hit==0)=NaN; 
ICso_hit(ICso_hit==0)=NaN; 
ICmto_hit(ICmto_hit==0)=NaN;
ICsto_hit(ICsto_hit==0)=NaN; 
ICmco_hit(ICmco_hit==0)=NaN; 
ICsco_hit(ICsco_hit==0)=NaN;

ICmr_fa(ICmr_fa==0) = NaN;
ICsr_fa(ICsr_fa==0) = NaN;
ICmp_fa(ICmp_fa==0) = NaN;
ICsp_fa(ICsp_fa==0) = NaN;
ICmo_fa(ICmo_fa==0) = NaN;
ICso_fa(ICso_fa==0) = NaN;
ICmto_fa(ICmto_fa==0) = NaN;
ICsto_fa(ICsto_fa==0) = NaN;
ICmco_fa(ICmco_fa==0) = NaN;
ICsco_fa(ICsco_fa==0) = NaN;

ICmlr_hit(ICmlr_hit==0) = NaN;
ICslr_hit(ICslr_hit==0) = NaN;
ICmlp_hit(ICmlp_hit==0) = NaN;
ICslp_hit(ICslp_hit==0) = NaN;
ICmlo_hit(ICmlo_hit==0) = NaN;
ICslo_hit(ICslo_hit==0) = NaN;
ICmlto_hit(ICmlto_hit==0) = NaN;
ICslto_hit(ICslto_hit==0) = NaN;
ICmlco_hit(ICmlco_hit==0) = NaN;
ICslco_hit(ICslco_hit==0) = NaN;

ICmlr_fa(ICmlr_fa==0) = NaN;
ICslr_fa(ICslr_fa==0) = NaN;
ICmlp_fa(ICmlp_fa==0) = NaN;
ICslp_fa(ICslp_fa==0) = NaN;
ICmlo_fa(ICmlo_fa==0) = NaN;
ICslo_fa(ICslo_fa==0) = NaN;
ICmlto_fa(ICmlto_fa==0) = NaN;
ICslto_fa(ICslto_fa==0) = NaN;
ICmlco_fa(ICmlco_fa==0) = NaN;
ICslco_fa(ICslco_fa==0) = NaN; 

close all
% % plots by session
bySess=1;
if bySess==1
    jj=2;     
    for jj=2:size(allDataTestsOnly,1)
        lickFig=figure(jj+103); 
        subplot(2,3,1); hold on; %MGB FUll
        lickbar=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmo_hit(jj,:)) nanmean(MGBmo_fa(jj,:))]); hold on;
        error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBso_hit(jj,:)) nanmean(MGBso_fa(jj,:))];
        lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        for kk=1:numel(lickbar)
            xtips=lickbar(kk).XEndPoints;
            ytips=lickbar(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        scatter(repmat(lickbar(1).XEndPoints(1),size(MGBmr_hit,1),1), ...
            MGBmr_hit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(2),size(MGBmr_fa,1),1), ...
            MGBmr_fa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(3),size(MGBmo_hit,1),1), ...
            MGBmo_hit,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbar(1).XEndPoints(4),size(MGBmo_fa,1),1), ...
            MGBmo_fa,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmo_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmo_fa(jj,:));
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
        xticklabels({'hit', 'fa','hit','fa'});

        subplot(2,3,2) % tone MGB
        lickbarT=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmto_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmto_fa(jj,:))]); hold on;
        error=[nanmean(MGBsr_hit(jj,:)) mean(MGBsto_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsto_fa(jj,:))];
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmto_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmto_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,3); % choice MGB
        lickbarC=bar([nanmean(MGBmr_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(MGBmco_hit(jj,:)) nanmean(MGBmco_fa(jj,:))]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        error=[nanmean(MGBsr_hit(jj,:)) nanmean(MGBsr_fa(jj,:)) nanmean(MGBsco_hit(jj,:)) nanmean(MGBsco_fa(jj,:))];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        scatter(repmat(lickbarC(1).XEndPoints(1),size(MGBmr_hit(jj,:),1),1), ...
            MGBmr_hit(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(2),size(MGBmr_fa(jj,:),1),1), ...
            MGBmr_fa(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(3),size(MGBmco_hit(jj,:),1),1), ...
            MGBmco_hit(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarT(1).XEndPoints(4),size(MGBmco_fa(jj,:),1),1), ...
            MGBmco_fa(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [MGBmr_hit(1) MGBmr_fa(1) MGBmr_hit(2) MGBmr_fa(2)...
            MGBmr_hit(3) MGBmr_fa(3) MGBmr_hit(4) MGBmr_fa(4) ...
            MGBmr_hit(5) MGBmr_fa(5) MGBmr_hit(6) MGBmr_fa(6)];
        oLickRate=[MGBmco_hit(1) MGBmco_fa(1) MGBmco_hit(2) MGBmco_fa(1)...
            MGBmco_hit(3) MGBmco_fa(3) MGBmco_hit(4) MGBmco_fa(4)...
            MGBmco_hit(5) MGBmco_fa(5) MGBmco_hit(6) MGBmco_fa(6)];
        line(xVector,lickRate, 'LineWidth', 0.5, 'Color', [0 0 0]);
        line(xoVector,oLickRate,'LineWidth', 0.5, 'Color', [0 0 0]);
        [h,pHit,ci,stats] = ttest2(MGBmr_hit(jj,:),MGBmco_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmr_fa(jj,:),MGBmco_fa(jj,:));
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});


        subplot(2,3,4); hold on; %IC Full
        lickbar=bar([nanmean(ICmr_hit(jj,:)) nanmean(ICmr_fa(jj,:)) nanmean(ICmo_hit(jj,:)) nanmean(ICmo_fa(jj,:))]); hold on;
        error=[nanmean(ICsr_hit(jj,:)) nanmean(ICsr_fa(jj,:)) nanmean(ICso_hit(jj,:)) nanmean(ICso_fa(jj,:))];
        lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        for kk=1:numel(lickbar)
            xtips=lickbar(kk).XEndPoints;
            ytips=lickbar(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmo_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmo_fa(jj,:));
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,5) % tone IC
        lickbarT=bar([nanmean(ICmr_hit(jj,:)) nanmean(MGBmr_fa(jj,:)) nanmean(ICmto_hit(jj,:)) nanmean(ICmto_fa(jj,:))]); hold on;
        error=[nanmean(ICsr_hit(jj,:)) nanmean(ICsr_fa(jj,:)) mean(ICsto_hit(jj,:)) nanmean(ICsto_fa(jj,:))];
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmto_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmto_fa(jj,:));
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,6); % choice IC
        lickbarC=bar([nanmean(ICmr_hit(jj,:)) nanmean(ICmr_fa(jj,:)) nanmean(ICmco_hit(jj,:)) nanmean(ICmco_fa(jj,:))]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        error=[nanmean(ICsr_hit(jj,:)) nanmean(ICsr_fa(jj,:)) nanmean(ICsco_hit(jj,:)) nanmean(ICsco_fa(jj,:))];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(ICmr_hit(jj,:),ICmco_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmr_fa(jj,:),ICmco_fa(jj,:));
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        lickFig.Position(3:4)=[725 475];
        saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickLat_Opto']);
        saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickLat_Opto.png']);

        % RATE
        lickFigl=figure(jj+13); 
        subplot(2,3,1); % full trial MGB
        lickbarL=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlo_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlo_fa(jj,:))]); hold on;
        lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslo_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslo_fa(jj,:))];
        for kk=1:numel(lickbarL)
            xtips=lickbarL(kk).XEndPoints;
            ytips=lickbarL(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlo_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlo_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' MGB Full Trial Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,2) % tone MGB
        lickbarT=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlto_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlto_fa(jj,:))]); hold on;
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslto_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslto_fa(jj,:))];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlto_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlto_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,3); % choice MGB
        lickbarC=bar([nanmean(MGBmlr_hit(jj,:)) nanmean(MGBmlco_hit(jj,:)) nanmean(MGBmlr_fa(jj,:)) nanmean(MGBmlco_fa(jj,:))]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(MGBslr_hit(jj,:)) nanmean(MGBslco_hit(jj,:)) nanmean(MGBslr_fa(jj,:)) nanmean(MGBslco_fa(jj,:))];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(MGBmlr_hit(jj,:),MGBmlco_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(MGBmlr_fa(jj,:),MGBmlco_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' MGB Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,4); % full trial IC
        lickbar=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlo_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlo_fa(jj,:))]); hold on;
        error=[nanmean(ICslr_hit(jj,:)) nanmean(ICslo_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslo_fa(jj,:))];
        lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        for kk=1:numel(lickbar)
            xtips=lickbar(kk).XEndPoints;
            ytips=lickbar(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end

        [~,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlo_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlo_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' IC Full Trial Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,5) % tone IC
        lickbarT=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlto_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlto_fa(jj,:))]); hold on;
        error=[nanmean(ICslr_hit(jj,:)) mean(ICslto_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslto_fa(jj,:))];
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlto_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlto_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,6); % choice IC
        lickbarC=bar([nanmean(ICmlr_hit(jj,:)) nanmean(ICmlr_fa(jj,:)) nanmean(ICmlco_hit(jj,:)) nanmean(ICmlco_fa(jj,:))]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
        error=[nanmean(ICslr_hit(jj,:)) nanmean(ICslr_fa(jj,:)) nanmean(ICslco_hit(jj,:)) nanmean(ICslco_fa(jj,:))];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        scatter(repmat(lickbar(1).XEndPoints(1),size(ICmlr_hit(jj,:),1),1), ...
            ICmlr_hit(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(2),size(ICslr_fa(jj,:),1),1), ...
            ICslr_fa(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(3),size(ICslco_hit(jj,:),1),1), ...
            ICslco_hit(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbar(1).XEndPoints(4),size(ICslco_fa(jj,:),1),1), ...
            ICslco_fa(jj,:),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[ohitLickRate(1) ofaLickRate(1) ohitLickRate(2) ofaLickRate(1)...
            ohitLickRate(3) ofaLickRate(3) ohitLickRate(4) ofaLickRate(4)...
            ohitLickRate(5) ofaLickRate(5) ohitLickRate(6) ICslco_hit(6)];
        line(xVector,lickRate, 'LineWidth', 0.5, 'Color', [0 0 0]);
        line(xoVector,oLickRate,'LineWidth', 0.5, 'Color', [0 0 0]);
        [h,pHit,ci,stats] = ttest2(ICmlr_hit(jj,:),ICmlco_hit(jj,:));
        [h,pFA,ci,stats] = ttest2(ICmlr_fa(jj,:),ICmlco_fa(jj,:));
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title([allDataTestsOnly{jj,1} ' IC Tone Inactivation']);
        xticklabels({'hit', 'hit','fa','fa'});

        lickFig.Position(3:4)=[725 475];
        saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickRate_Opto']);
        saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_IC_LickRate_Opto.png']);  
    end
    else
end

%%
makeLickLat=1;
if makeLickLat==0
else
    jj=2;
    close all
    % LICK LATENCY, AVERAGED ACROSS ANIMALS
    % change these variable names, theyre swapped
    lickFig=figure(jj+103); 
    xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
    xoVector = xVector+2;
    MGBmr_hit1=MGBmr_hit;
    MGBmr_fa1=MGBmr_fa;
    for jj=2:size(mgbTempTestsOnly,1)
        MGBmr_hit1(isnan(MGBmo_hit))=nan;
        MGBmr_fa1(isnan(MGBmo_fa))=nan;
        hitLickRate(jj-1)=nanmean(MGBmr_hit1(jj,:));
        ohitLickRate(jj-1)=nanmean(MGBmo_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmr_fa1(jj,:));
        ofaLickRate(jj-1)=nanmean(MGBmo_fa(jj,:));

        MGBsr_hit(isnan(MGBso_hit))=nan;
        MGBsr_fa(isnan(MGBso_fa))=nan;
        hitLickRateE(jj-1)=nanmean(MGBsr_hit(jj,:));
        ohitLickRateE(jj-1)=nanmean(MGBso_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBsr_fa(jj,:));
        ofaLickRateE(jj-1)=nanmean(MGBso_fa(jj,:));
    end

        subplot(2,3,1); % full trial MGB
    %     [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBarPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate,hitLickRateE,faLickRateE,ohitLickRateE,ofaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate);
        mgbTempTestsOnly{2,42} = [hitLickRate; ohitLickRate; faLickRate; ofaLickRate];
        ylim([0 2.5]);
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Animal MGB Full Trial']);
        xticklabels({'hit', 'fa','hit','fa'});

    clear hitLickRate faLickRate pHit pFA
    MGBmr_hit1=MGBmr_hit;
    MGBmr_fa1=MGBmr_fa;
    for jj=2:size(mgbTempTestsOnly,1)
        MGBmr_hit1(isnan(MGBmto_hit))=nan;
        MGBmr_fa1(isnan(MGBmto_fa))=nan;
        hitLickRate(jj-1)=nanmean(MGBmr_hit1(jj,:));
        othitLickRate(jj-1)=nanmean(MGBmto_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmr_fa1(jj,:));
        otfaLickRate(jj-1)=nanmean(MGBmto_fa(jj,:));

        MGBsr_hit(isnan(MGBsto_hit))=nan;
        MGBsr_fa(isnan(MGBsto_fa))=nan;
        hitLickRateE(jj-1)=nanmean(MGBsr_hit(jj,:));
        othitLickRateE(jj-1)=nanmean(MGBsto_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBsr_fa(jj,:));
        otfaLickRateE(jj-1)=nanmean(MGBsto_fa(jj,:));
    end
    subplot(2,3,2) % tone MGB
    %     [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBarPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate,hitLickRateE,faLickRateE,othitLickRateE,otfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate);
        mgbTempTestsOnly{2,43} = [hitLickRate; othitLickRate; faLickRate; otfaLickRate];
        ylim([0 2.5]);
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Animal MGB Tone']);
        xticklabels({'hit', 'fa','hit','fa'});

    clear hitLickRate faLickRate pHit pFA
    MGBmr_hit1=MGBmr_hit;
    MGBmr_fa1=MGBmr_fa;
    for jj=2:size(mgbTempTestsOnly,1)
        MGBmr_hit1(isnan(MGBmco_hit))=nan;
        MGBmr_fa1(isnan(MGBmco_fa))=nan;
        hitLickRate(jj-1)=nanmean(MGBmr_hit1(jj,:));
        ochitLickRate(jj-1)=nanmean(MGBmco_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmr_fa1(jj,:));
        ocfaLickRate(jj-1)=nanmean(MGBmco_fa(jj,:));

        MGBsr_hit(isnan(MGBsco_hit))=nan;
        MGBsr_fa(isnan(MGBsco_fa))=nan;
        hitLickRateE(jj-1)=nanmean(MGBsr_hit(jj,:));
        ochitLickRateE(jj-1)=nanmean(MGBsco_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBsr_fa(jj,:));
        ocfaLickRateE(jj-1)=nanmean(MGBsco_fa(jj,:));
    end
    subplot(2,3,3); % choice MGB
    %     [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBarPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate,hitLickRateE,faLickRateE,ochitLickRateE,ocfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate);
        mgbTempTestsOnly{2,44} = [hitLickRate; ochitLickRate; faLickRate; ocfaLickRate];
        ylim([0 2.5]);
        sigstar({[1,3],[2,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Animal MGB Choice']);
        xticklabels({'hit', 'fa','hit','fa'});


    icIdx=[8:13];
    clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ochitLickRate ocfaLickRate pHit pFA
    MGBmr_hit1=MGBmr_hit;
    MGBmr_fa1=MGBmr_fa;
    for jj=1:length(icIdx)
        ICmr_hit2=ICmr_hit;
        ICmr_fa2=ICmr_fa;
        ICsr_hit2=ICsr_hit;
        ICsr_fa2=ICsr_fa;
        ICmr_hit2(isnan(ICmo_hit))=nan;
        ICmr_fa2(isnan(ICmo_fa))=nan;
        hitLickRate(jj)=nanmean(ICmr_hit2(icIdx(jj),:));
        ohitLickRate(jj)=nanmean(ICmo_hit(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa2(icIdx(jj),:));
        ofaLickRate(jj)=nanmean(ICmo_fa(icIdx(jj),:));

        ICsr_hit2(isnan(ICso_hit))=nan;
        ICsr_fa2(isnan(ICso_fa))=nan;
        hitLickRateE(jj)=nanmean(ICsr_hit2(icIdx(jj),:));
        ohitLickRateE(jj)=nanmean(ICso_hit(icIdx(jj),:));
        faLickRateE(jj)=nanmean(ICsr_fa2(icIdx(jj),:));
        ofaLickRateE(jj)=nanmean(ICso_fa(icIdx(jj),:));
    end

    subplot(2,3,4); % full trial IC
    %     [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBarPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate,hitLickRateE,faLickRateE,ohitLickRateE,ofaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate);

    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2.5]);
    title(['By Animal IC Full Trial']);
    xticklabels({'hit', 'hit','fa','fa'});


    clear hitLickRate faLickRate pHit pFA
    for jj=1:length(icIdx)
        ICmr_hit3=ICmr_hit;
        ICmr_fa3=ICmr_fa;
        ICsr_hit3=ICsr_hit;
        ICsr_fa3=ICsr_fa;

        ICmr_hit3(isnan(ICmto_hit))=nan;
        ICmr_fa3(isnan(ICmto_fa))=nan;
        ICsr_hit3(isnan(ICsto_hit))=nan;
        ICsr_fa3(isnan(ICsto_fa))=nan;

        hitLickRate(jj)=nanmean(ICmr_hit3(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa3(icIdx(jj),:));
        othitLickRate(jj)=nanmean(ICmto_hit(icIdx(jj),:));
        otfaLickRate(jj)=nanmean(ICmto_fa(icIdx(jj),:));

        hitLickRateE(jj)=nanmean(ICsr_hit3(icIdx(jj),:));    
        faLickRateE(jj)=nanmean(ICsr_fa3(icIdx(jj),:));
        othitLickRateE(jj)=nanmean(ICsto_hit(icIdx(jj),:));
        otfaLickRateE(jj)=nanmean(ICsto_fa(icIdx(jj),:));
    end

    subplot(2,3,5) % tone IC
    %     [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBarPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate,hitLickRateE,faLickRateE,othitLickRateE,otfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate);

    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2.5]);
    title(['By Animal IC Tone']);
    xticklabels({'hit', 'hit','fa','fa'});

    clear hitLickRate faLickRate pHit pFA
    for jj=1:length(icIdx)
        ICmr_hit4=ICmr_hit;
        ICmr_fa4=ICmr_fa;
        ICsr_hit4=ICsr_hit;
        ICsr_fa4=ICsr_fa;    
        ICmr_hit4(isnan(ICmco_hit))=nan;
        ICmr_fa4(isnan(ICmco_fa))=nan;
        ICsr_hit4(isnan(ICsco_hit))=nan;
        ICsr_fa4(isnan(ICsco_fa))=nan;

        hitLickRate(jj)=nanmean(ICmr_hit4(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa4(icIdx(jj),:));
        ochitLickRate(jj)=nanmean(ICmco_hit(icIdx(jj),:));
        ocfaLickRate(jj)=nanmean(ICmco_fa(icIdx(jj),:));

        hitLickRateE(jj)=nanmean(ICsr_hit4(icIdx(jj),:));    
        faLickRateE(jj)=nanmean(ICsr_fa4(icIdx(jj),:));    
        ochitLickRateE(jj)=nanmean(ICsco_hit(icIdx(jj),:));
        ocfaLickRateE(jj)=nanmean(ICsco_fa(icIdx(jj),:));
    end

    subplot(2,3,6); % choice IC
    %     [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBarPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate,hitLickRateE,faLickRateE,ochitLickRateE,ocfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate);

    sigstar({[1,3],[2,4]}, [pHit pFA]); ylim([0 2.5]);
    ylabel('mean lick latency');
    title(['By Animal IC Choice']);
    xticklabels({'hit', 'fa','hit','fa'});

    lickFig.Position(3:4)=[725 475];
    saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto']);
    saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto.pdf']);
    close all


    % RATE
    lickFigl=figure(jj+13); 
    clear hitLickRate faLickRate pHit pFA
    for jj=2:size(mgbTempTestsOnly,1)
        hitLickRate(jj-1)=nanmean(MGBmlr_hit(jj,:));
        ohitLickRate(jj-1)=nanmean(MGBmlo_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmlr_fa(jj,:));
        ofaLickRate(jj-1)=nanmean(MGBmlo_fa(jj,:));
        othitLickRate(jj-1)=nanmean(MGBmlto_hit(jj,:));
        otfaLickRate(jj-1)=nanmean(MGBmlto_fa(jj,:));
        ochitLickRate(jj-1)=nanmean(MGBmlco_hit(jj,:));
        ocfaLickRate(jj-1)=nanmean(MGBmlco_fa(jj,:));

        hitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
        ohitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBslr_fa(jj,:));
        ofaLickRateE(jj-1)=nanmean(MGBslo_fa(jj,:));
        othitLickRateE(jj-1)=nanmean(MGBslto_hit(jj,:));
        otfaLickRateE(jj-1)=nanmean(MGBslto_fa(jj,:));
        ochitLickRateE(jj-1)=nanmean(MGBslco_hit(jj,:));
        ocfaLickRateE(jj-1)=nanmean(MGBslco_fa(jj,:));    

    end

    subplot(2,3,1); % full trial MGB
    %     [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBarPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate,hitLickRateE,faLickRateE,ohitLickRateE,ofaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate);
    mgbTempTestsOnly{2,46} = [hitLickRate; ohitLickRate; faLickRate; ofaLickRate];
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Full Trial']);
    xticklabels({'hit', 'fa','hit','fa'});

    clear pHit pFA
    subplot(2,3,2) % tone MGB
    %     [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBarPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate,hitLickRateE,faLickRateE,othitLickRateE,otfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate);    
    [h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
    mgbTempTestsOnly{2,46} = [hitLickRate; othitLickRate; faLickRate; otfaLickRate];
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Tone']);
    xticklabels({'hit', 'fa','hit','fa'});

    clear pHit pFA
    subplot(2,3,3); % choice MGB
    %     [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBarPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate,hitLickRateE,faLickRateE,ochitLickRateE,ocfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate);
    mgbTempTestsOnly{2,47} = [hitLickRate; ochitLickRate; faLickRate; ocfaLickRate];
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Choice']);
    xticklabels({'hit', 'fa','hit','fa'});

    % do the IC
    clear pFA pHit hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ochitLickRate ocfaLickRate lickbar lickbarL lickbarC lickbarT
    for jj=1:length(icIdx)
        ICmr_hit2=ICmlr_hit;
        ICmr_fa2=ICmlr_fa;
        ICsr_hit2=ICslr_hit;
        ICsr_fa2=ICslr_fa;
        ICmr_hit2(isnan(ICmlo_hit))=nan;
        ICmr_fa2(isnan(ICmlo_fa))=nan;
        hitLickRate(jj)=nanmean(ICmr_hit2(icIdx(jj),:));
        ohitLickRate(jj)=nanmean(ICmlo_hit(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa2(icIdx(jj),:));
        ofaLickRate(jj)=nanmean(ICmlo_fa(icIdx(jj),:));

        ICsr_hit2(isnan(ICslo_hit))=nan;
        ICsr_fa2(isnan(ICslo_fa))=nan;
        hitLickRateE(jj)=nanmean(ICsr_hit2(icIdx(jj),:));
        ohitLickRateE(jj)=nanmean(ICslo_hit(icIdx(jj),:));
        faLickRateE(jj)=nanmean(ICsr_fa2(icIdx(jj),:));
        ofaLickRateE(jj)=nanmean(ICslo_fa(icIdx(jj),:));
    end
    subplot(2,3,4); % full trial IC
    %     [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBarPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate,hitLickRateE,faLickRateE,ohitLickRateE,ofaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ohitLickRate,ofaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ohitLickRate,ofaLickRate);

    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal IC Full Trial']);
    xticklabels({'hit', 'fa','hit','fa'});

    clear hitLickRate faLickRate pHit pFA
    for jj=1:length(icIdx)
        ICmr_hit2=ICmlr_hit;
        ICmr_fa2=ICmlr_fa;
        ICsr_hit2=ICslr_hit;
        ICsr_fa2=ICslr_fa;
        ICmr_hit2(isnan(ICmlto_hit))=nan;
        ICmr_fa2(isnan(ICmlto_fa))=nan;
        hitLickRate(jj)=nanmean(ICmr_hit2(icIdx(jj),:));
        othitLickRate(jj)=nanmean(ICmlto_hit(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa2(icIdx(jj),:));
        otfaLickRate(jj)=nanmean(ICmlto_fa(icIdx(jj),:));

        ICsr_hit2(isnan(ICslto_hit))=nan;
        ICsr_fa2(isnan(ICslto_fa))=nan;
        hitLickRateE(jj)=nanmean(ICsr_hit2(icIdx(jj),:));
        othitLickRateE(jj)=nanmean(ICslto_hit(icIdx(jj),:));
        faLickRateE(jj)=nanmean(ICsr_fa2(icIdx(jj),:));
        otfaLickRateE(jj)=nanmean(ICslto_fa(icIdx(jj),:));
    end
    subplot(2,3,5) % tone IC
    %     [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBarPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate,hitLickRateE,faLickRateE,othitLickRateE,otfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,othitLickRate,otfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,othitLickRate,otfaLickRate);    
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');ylim([0 8]);
    title(['By Animal IC Tone']);
    xticklabels({'hit', 'fa','hit','fa'});

    clear hitLickRate faLickRate pHit pFA
    for jj=1:length(icIdx)
        ICmr_hit2=ICmlr_hit;
        ICmr_fa2=ICmlr_fa;
        ICsr_hit2=ICslr_hit;
        ICsr_fa2=ICslr_fa;
        ICmr_hit2(isnan(ICmlco_hit))=nan;
        ICmr_fa2(isnan(ICmlco_fa))=nan;
        hitLickRate(jj)=nanmean(ICmr_hit2(icIdx(jj),:));
        ochitLickRate(jj)=nanmean(ICmlco_hit(icIdx(jj),:));
        faLickRate(jj)=nanmean(ICmr_fa2(icIdx(jj),:));
        ocfaLickRate(jj)=nanmean(ICmlco_fa(icIdx(jj),:));

        ICsr_hit2(isnan(ICslto_hit))=nan;
        ICsr_fa2(isnan(ICslto_fa))=nan;
        hitLickRateE(jj)=nanmean(ICsr_hit2(icIdx(jj),:));
        ochitLickRateE(jj)=nanmean(ICslco_hit(icIdx(jj),:));
        faLickRateE(jj)=nanmean(ICsr_fa2(icIdx(jj),:));
        ocfaLickRateE(jj)=nanmean(ICslco_fa(icIdx(jj),:));  
    end
    subplot(2,3,6); % choice IC
    %     [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBarPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate,hitLickRateE,faLickRateE,ochitLickRateE,ocfaLickRateE);
        [pHit,pFA,hitLickRate,faLickRate,ochitLickRate,ocfaLickRate]=lickBoxPlot(hitLickRate,faLickRate,ochitLickRate,ocfaLickRate);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');ylim([0 8]);
    title(['By Animal IC Choice']);
    xticklabels({'hit','fa', 'hit','fa'});

    lickFigl.Position(3:4)=[725 475];
    saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto']);
    saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto.pdf']);  
    close all
end

%% by session
bySession=0;
    if bySession==1
        lickFig1=figure(jj+108); 
        % latency
        clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ...
            ochitLickRate ocfaLickRate hitLickRateE ohitLickRateE faLickRateE ...
            ofaLickRateE othitLickRateE otfaLickRateE ochitLickRateE ocfaLickRateE 
        hitLickRate=NaN; ohitLickRate=NaN; faLickRate=NaN; ofaLickRate=NaN; othitLickRate=NaN;
        otfaLickRate=NaN; ochitLickRate=NaN; ocfaLickRate=NaN; hitLickRateE=NaN; ohitLickRateE=NaN;
        faLickRateE=NaN; ofaLickRateE=NaN; othitLickRateE=NaN; otfaLickRateE=NaN; ochitLickRateE=NaN;
        ocfaLickRateE=NaN; 
        for jj=2:size(mgbTempTestsOnly,1)
            hitLickRate=cat(2,hitLickRate,MGBmr_hit(jj,:));
            ohitLickRate=cat(2,ohitLickRate,MGBmo_hit(jj,:));
            faLickRate=cat(2,faLickRate,MGBmr_fa(jj,:));
            ofaLickRate=cat(2,ofaLickRate,MGBmo_fa(jj,:));
            othitLickRate=cat(2,othitLickRate,MGBmto_hit(jj,:));
            otfaLickRate=cat(2,otfaLickRate,MGBmto_fa(jj,:));
            ochitLickRate=cat(2,ochitLickRate,MGBmco_hit(jj,:));
            ocfaLickRate=cat(2,ocfaLickRate,MGBmco_fa(jj,:));

            hitLickRateE=cat(2,hitLickRateE,MGBsr_hit(jj,:));
            ohitLickRateE=cat(2,ohitLickRateE,MGBso_hit(jj,:));
            faLickRateE=cat(2,faLickRateE,MGBsr_fa(jj,:));
            ofaLickRateE=cat(2,ofaLickRateE,MGBso_fa(jj,:));
            othitLickRateE=cat(2,othitLickRateE,MGBsto_hit(jj,:));
            otfaLickRateE=cat(2,otfaLickRateE,MGBsto_fa(jj,:));
            ochitLickRateE=cat(2,ochitLickRateE,MGBsco_hit(jj,:));
            ocfaLickRateE=cat(2,ocfaLickRateE,MGBsco_fa(jj,:));

        end

        subplot(2,3,1); % full trial MGB
        lickbar=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
        error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
        lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        for kk=1:numel(lickbar)
            xtips=lickbar(kk).XEndPoints;
            ytips=lickbar(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Sess MGB Full Trial']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,2) % tone MGB
        lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
        error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Sess MGB Tone']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,3); % choice MGB
        lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick latency');
        title(['By Sess MGB Choice']);
        xticklabels({'hit', 'hit','fa','fa'});

        % for rate
        clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ...
            ochitLickRate ocfaLickRate hitLickRateE ohitLickRateE faLickRateE ...
            ofaLickRateE othitLickRateE otfaLickRateE ochitLickRateE ocfaLickRateE 
        hitLickRate=NaN; ohitLickRate=NaN; faLickRate=NaN; ofaLickRate=NaN; othitLickRate=NaN;
        otfaLickRate=NaN; ochitLickRate=NaN; ocfaLickRate=NaN; hitLickRateE=NaN; ohitLickRateE=NaN;
        faLickRateE=NaN; ofaLickRateE=NaN; othitLickRateE=NaN; otfaLickRateE=NaN; ochitLickRateE=NaN;
        ocfaLickRateE=NaN;
        for jj=2:size(mgbTempTestsOnly,1)
            hitLickRate=cat(2,hitLickRate,MGBmlr_hit(jj,:));
            ohitLickRate=cat(2,ohitLickRate,MGBmlo_hit(jj,:));
            faLickRate=cat(2,faLickRate,MGBmlr_fa(jj,:));
            ofaLickRate=cat(2,ofaLickRate,MGBmlo_fa(jj,:));
            othitLickRate=cat(2,othitLickRate,MGBmlto_hit(jj,:));
            otfaLickRate=cat(2,otfaLickRate,MGBmlto_fa(jj,:));
            ochitLickRate=cat(2,ochitLickRate,MGBmlco_hit(jj,:));
            ocfaLickRate=cat(2,ocfaLickRate,MGBmlco_fa(jj,:));

            hitLickRateE=cat(2,hitLickRateE,MGBslr_hit(jj,:));
            ohitLickRateE=cat(2,ohitLickRateE,MGBslo_hit(jj,:));
            faLickRateE=cat(2,faLickRateE,MGBslr_fa(jj,:));
            ofaLickRateE=cat(2,ofaLickRateE,MGBslo_fa(jj,:));
            othitLickRateE=cat(2,othitLickRateE,MGBslto_hit(jj,:));
            otfaLickRateE=cat(2,otfaLickRateE,MGBslto_fa(jj,:));
            ochitLickRateE=cat(2,ochitLickRateE,MGBslco_hit(jj,:));
            ocfaLickRateE=cat(2,ocfaLickRateE,MGBslco_fa(jj,:));
        end

        subplot(2,3,4); % full trial MGB
        lickbarL=bar([nanmean(hitLickRate) nanmean(ohitLickRate) nanmean(faLickRate) nanmean(ofaLickRate)]); hold on;
        lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(hitLickRateE) nanmean(ohitLickRateE) nanmean(faLickRateE) nanmean(ofaLickRateE)];
        for kk=1:numel(lickbarL)
            xtips=lickbarL(kk).XEndPoints;
            ytips=lickbarL(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title(['By Sess MGB Full Trial']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,5) % tone MGB
        lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
        lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
        for kk=1:numel(lickbarT)
            xtips=lickbarT(kk).XEndPoints;
            ytips=lickbarT(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title(['By Sess MGB Tone']);
        xticklabels({'hit', 'hit','fa','fa'});

        subplot(2,3,6); % choice MGB
        lickbarC=bar([nanmean(hitLickRate) nanmean(ochitLickRate) nanmean(faLickRate) nanmean(ocfaLickRate)]); hold on;
        lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
        error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
        for kk=1:numel(lickbarC)
            xtips=lickbarC(kk).XEndPoints;
            ytips=lickbarC(kk).YEndPoints;
            errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
        end
        [h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
        [h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
        sigstar({[1,2],[3,4]}, [pHit pFA])
        ylabel('mean lick rate');
        title(['By Sess MGB Choice']);
        xticklabels({'hit', 'hit','fa','fa'});

        saveas(gcf,['BySess_T_MGB_LickRateLat_Opto']);
        saveas(gcf,['BySess_T_MGB_LickRateLat_Opto.png']);  
    else
    end

%% control animals
ctlFlag=0;
if ctlFlag==0
else    
    clearvars -except allDataCtlOnly reinfcolor optocolor

    SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
    START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

    % ctl for MGB for each animal
    nbsubj=2;
    for nbsubj=2:size(allDataCtlOnly,1)
        matrix=allDataCtlOnly{nbsubj,26};
        for i=1:max(matrix(:,SESS))
            mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
            sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
            mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            mo_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            so_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            mto_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            sto_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            mco_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
            sco_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));

            mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            mo_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            so_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            mto_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            sto_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            mco_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
            sco_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));

            mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mlo_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            slo_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mlto_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            slto_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mlco_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            slco_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));

            mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mlo_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            slo_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mlto_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            slto_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mlco_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            slco_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));        

        end
        switch nbsubj
            case 2 % sk177
                mgbDays=[0 0 0 0 0 0 0 0 1 1 0 0 0 1 1 0 0];
                mgbDays=logical(mgbDays);
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

            case 3 %sk179
                mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 0 1 1 1 0 0 1 0 0];
                mgbDays=logical(mgbDays);
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
                MGBmr_hit(MGBmr_hit==0)=NaN;
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

            case 4 % sk183
                mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                    0 0 1 1 0 0 1 1 0 0];
                mgbDays=logical(mgbDays);
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);
                MGBmr_hit(MGBmr_hit==0)=NaN;
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

        case 5 % sk189
                mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                    0 0 1 1 0 0 1 1 1 0];
                mgbDays=logical(mgbDays);  
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

        case 6 % sk290
                 mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 0 0 ...
                    0 0 0 1 1];
                mgbDays=logical(mgbDays);
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

        case 7 % sk191
                mgbDays=[0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 0 0 0 0 0 0 0 ...
                    0 0 0 1 1 1 1 0 0 0 ...
                    0];
                mgbDays=logical(mgbDays);
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); MGBsr_hit(MGBsr_hit==0)=NaN;
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays); MGBmp_hit(MGBmp_hit==0)=NaN;
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays); MGBsp_hit(MGBsp_hit==0)=NaN;
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays); MGBmo_hit(MGBmo_hit==0)=NaN;
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays); MGBso_hit(MGBso_hit==0)=NaN;
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);MGBmto_hit(MGBmto_hit==0)=NaN;
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays); MGBsto_hit(MGBsto_hit==0)=NaN;
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays); MGBmco_hit(MGBmco_hit==0)=NaN;
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays); MGBsco_hit(MGBsco_hit==0)=NaN;
                MGBmr_hit(nbsubj,mgbDays) = mr_hit(nbsubj,mgbDays); % reinf hit
                MGBsr_hit(nbsubj,mgbDays) = sr_hit(nbsubj,mgbDays); 
                MGBmp_hit(nbsubj,mgbDays) = mp_hit(nbsubj,mgbDays);
                MGBsp_hit(nbsubj,mgbDays) = sp_hit(nbsubj,mgbDays);
                MGBmo_hit(nbsubj,mgbDays) = mo_hit(nbsubj,mgbDays);
                MGBso_hit(nbsubj,mgbDays) = so_hit(nbsubj,mgbDays);
                MGBmto_hit(nbsubj,mgbDays) = mto_hit(nbsubj,mgbDays);
                MGBsto_hit(nbsubj,mgbDays) = sto_hit(nbsubj,mgbDays);
                MGBmco_hit(nbsubj,mgbDays) = mco_hit(nbsubj,mgbDays);
                MGBsco_hit(nbsubj,mgbDays) = sco_hit(nbsubj,mgbDays);

                MGBmr_fa(nbsubj,mgbDays) = mr_fa(nbsubj,mgbDays);
                MGBsr_fa(nbsubj,mgbDays) = sr_fa(nbsubj,mgbDays);
                MGBmp_fa(nbsubj,mgbDays) = mp_fa(nbsubj,mgbDays);
                MGBsp_fa(nbsubj,mgbDays) = sp_fa(nbsubj,mgbDays);
                MGBmo_fa(nbsubj,mgbDays) = mo_fa(nbsubj,mgbDays);
                MGBso_fa(nbsubj,mgbDays) = so_fa(nbsubj,mgbDays);
                MGBmto_fa(nbsubj,mgbDays) = mto_fa(nbsubj,mgbDays);
                MGBsto_fa(nbsubj,mgbDays) = sto_fa(nbsubj,mgbDays);
                MGBmco_fa(nbsubj,mgbDays) = mco_fa(nbsubj,mgbDays);
                MGBsco_fa(nbsubj,mgbDays) = sco_fa(nbsubj,mgbDays);

                MGBmlr_hit(nbsubj,mgbDays) = mlr_hit(nbsubj,mgbDays);
                MGBslr_hit(nbsubj,mgbDays) = slr_hit(nbsubj,mgbDays);
                MGBmlp_hit(nbsubj,mgbDays) = mlp_hit(nbsubj,mgbDays);
                MGBslp_hit(nbsubj,mgbDays) = slp_hit(nbsubj,mgbDays);
                MGBmlo_hit(nbsubj,mgbDays) = mlo_hit(nbsubj,mgbDays);
                MGBslo_hit(nbsubj,mgbDays) = slo_hit(nbsubj,mgbDays);
                MGBmlto_hit(nbsubj,mgbDays) = mlto_hit(nbsubj,mgbDays);
                MGBslto_hit(nbsubj,mgbDays) = slto_hit(nbsubj,mgbDays);
                MGBmlco_hit(nbsubj,mgbDays) = mlco_hit(nbsubj,mgbDays);
                MGBslco_hit(nbsubj,mgbDays) = slco_hit(nbsubj,mgbDays);

                MGBmlr_fa(nbsubj,mgbDays) = mlr_fa(nbsubj,mgbDays);
                MGBslr_fa(nbsubj,mgbDays) = slr_fa(nbsubj,mgbDays);
                MGBmlp_fa(nbsubj,mgbDays) = mlp_fa(nbsubj,mgbDays);
                MGBslp_fa(nbsubj,mgbDays) = slp_fa(nbsubj,mgbDays);
                MGBmlo_fa(nbsubj,mgbDays) = mlo_fa(nbsubj,mgbDays);
                MGBslo_fa(nbsubj,mgbDays) = slo_fa(nbsubj,mgbDays);
                MGBmlto_fa(nbsubj,mgbDays) = mlto_fa(nbsubj,mgbDays);
                MGBslto_fa(nbsubj,mgbDays) = slto_fa(nbsubj,mgbDays);
                MGBmlco_fa(nbsubj,mgbDays) = mlco_fa(nbsubj,mgbDays);
                MGBslco_fa(nbsubj,mgbDays) = slco_fa(nbsubj,mgbDays);

        end
    end

    MGBmr_hit(MGBmr_hit==0)=NaN;
    MGBsr_hit(MGBsr_hit==0)=NaN;
    MGBmp_hit(MGBmp_hit==0)=NaN;
    MGBsp_hit(MGBsp_hit==0)=NaN;
    MGBmo_hit(MGBmo_hit==0)=NaN;
    MGBso_hit(MGBso_hit==0)=NaN;
    MGBmto_hit(MGBmto_hit==0)=NaN;
    MGBsto_hit(MGBsto_hit==0)=NaN;
    MGBmco_hit(MGBmco_hit==0)=NaN;
    MGBsco_hit(MGBsco_hit==0)=NaN;

    MGBmr_fa(MGBmr_fa==0) = NaN;
    MGBsr_fa(MGBsr_fa==0) = NaN;
    MGBmp_fa(MGBmp_fa==0) = NaN;
    MGBsp_fa(MGBsp_fa==0) = NaN;
    MGBmo_fa(MGBmo_fa==0) = NaN;
    MGBso_fa(MGBso_fa==0) = NaN;
    MGBmto_fa(MGBmto_fa==0) = NaN;
    MGBsto_fa(MGBsto_fa==0) = NaN;
    MGBmco_fa(MGBmco_fa==0) = NaN;
    MGBsco_fa(MGBsco_fa==0) = NaN;

    MGBmlr_hit(MGBmlr_hit==0) = NaN;
    MGBslr_hit(MGBslr_hit==0) = NaN;
    MGBmlp_hit(MGBmlp_hit==0) = NaN;
    MGBslp_hit(MGBslp_hit==0) = NaN;
    MGBmlo_hit(MGBmlo_hit==0) = NaN;
    MGBslo_hit(MGBslo_hit==0) = NaN;
    MGBmlto_hit(MGBmlto_hit==0) = NaN;
    MGBslto_hit(MGBslto_hit==0) = NaN;
    MGBmlco_hit(MGBmlco_hit==0) = NaN;
    MGBslco_hit(MGBslco_hit==0) = NaN;

    MGBmlr_fa(MGBmlr_fa==0) = NaN;
    MGBslr_fa(MGBslr_fa==0) = NaN;
    MGBmlp_fa(MGBmlp_fa==0) = NaN;
    MGBslp_fa(MGBslp_fa==0) = NaN;
    MGBmlo_fa(MGBmlo_fa==0) = NaN;
    MGBslo_fa(MGBslo_fa==0) = NaN;
    MGBmlto_fa(MGBmlto_fa==0) = NaN;
    MGBslto_fa(MGBslto_fa==0) = NaN;
    MGBmlco_fa(MGBmlco_fa==0) = NaN;
    MGBslco_fa(MGBslco_fa==0) = NaN;


    close all
    % LICK LATENCY, AVERAGED ACROSS ANIMALS
    % change these variable names, theyre swapped
    jj=2;
    lickFig=figure(jj+103); 
    xVector = [1 2 1 2 1 2 1 2 1 2 1 2];
    xoVector = xVector+2;
    for jj=2:size(allDataCtlOnly,1)
        hitLickRate(jj-1)=nanmean(MGBmr_hit(jj,:));
        ohitLickRate(jj-1)=nanmean(MGBmo_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmr_fa(jj,:));
        ofaLickRate(jj-1)=nanmean(MGBmo_fa(jj,:));
        othitLickRate(jj-1)=nanmean(MGBmto_hit(jj,:));
        otfaLickRate(jj-1)=nanmean(MGBmto_fa(jj,:));
        ochitLickRate(jj-1)=nanmean(MGBmco_hit(jj,:));
        ocfaLickRate(jj-1)=nanmean(MGBmco_fa(jj,:));

        hitLickRateE(jj-1)=nanmean(MGBsr_hit(jj,:));
        ohitLickRateE(jj-1)=nanmean(MGBso_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBsr_fa(jj,:));
        ofaLickRateE(jj-1)=nanmean(MGBso_fa(jj,:));
        othitLickRateE(jj-1)=nanmean(MGBsto_hit(jj,:));
        otfaLickRateE(jj-1)=nanmean(MGBsto_fa(jj,:));
        ochitLickRateE(jj-1)=nanmean(MGBsco_hit(jj,:));
        ocfaLickRateE(jj-1)=nanmean(MGBsco_fa(jj,:));

    end

    subplot(2,3,1); % full trial MGB
    lickbar=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(ohitLickRate) nanmean(ofaLickRate)]); hold on;
    error=[nanmean(hitLickRateE) nanmean(faLickRateE) nanmean(ohitLickRateE) nanmean(ofaLickRateE)];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    % for kk=1:numel(lickbar)
    %     xtips=lickbar(kk).XEndPoints;
    %     ytips=lickbar(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
        scatter(repmat(lickbar(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbar(1).XEndPoints(3),size(ohitLickRate,1),1), ...
            ohitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbar(1).XEndPoints(4),size(ofaLickRate,1),1), ...
            ofaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[ohitLickRate(1) ofaLickRate(1) ohitLickRate(2) ofaLickRate(2)...
            ohitLickRate(3) ofaLickRate(3) ohitLickRate(4) ofaLickRate(4)...
            ohitLickRate(5) ofaLickRate(5) ohitLickRate(6) ofaLickRate(6)];
            line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
    ylim([0 2.5]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title(['By Animal MGB Full Trial']);
    xticklabels({'hit', 'fa','hit','fa'});

    subplot(2,3,2) % tone MGB
    lickbarT=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(othitLickRate) nanmean(otfaLickRate)]); hold on;
    error=[nanmean(hitLickRateE) nanmean(faLickRateE) nanmean(othitLickRateE) nanmean(otfaLickRateE)];
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    % for kk=1:numel(lickbarT)
    %     xtips=lickbarT(kk).XEndPoints;
    %     ytips=lickbarT(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
    scatter(repmat(lickbarT(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(3),size(othitLickRate,1),1), ...
            othitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarT(1).XEndPoints(4),size(otfaLickRate,1),1), ...
            otfaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[othitLickRate(1) otfaLickRate(1) othitLickRate(2) otfaLickRate(2)...
            othitLickRate(3) otfaLickRate(3) othitLickRate(4) otfaLickRate(4)...
            othitLickRate(5) otfaLickRate(5) othitLickRate(6) otfaLickRate(6)];
        line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
    ylim([0 2.5]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title(['By Animal MGB Tone']);
    xticklabels({'hit', 'fa','hit','fa'});

    subplot(2,3,3); % choice MGB
    lickbarC=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(ochitLickRate) nanmean(ocfaLickRate)]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    error=[nanmean(hitLickRateE)  nanmean(faLickRateE) nanmean(ochitLickRateE) nanmean(ocfaLickRateE)];
    % for kk=1:numel(lickbarC)
    %     xtips=lickbarC(kk).XEndPoints;
    %     ytips=lickbarC(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
    scatter(repmat(lickbarC(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarC(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarC(1).XEndPoints(3),size(ochitLickRate,1),1), ...
            ochitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarC(1).XEndPoints(4),size(ocfaLickRate,1),1), ...
            ocfaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[ochitLickRate(1) ocfaLickRate(1) ochitLickRate(2) ocfaLickRate(2)...
            ochitLickRate(3) ocfaLickRate(3) ochitLickRate(4) ocfaLickRate(4)...
            ochitLickRate(5) ocfaLickRate(5) ochitLickRate(6) ocfaLickRate(6)];
        line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
    ylim([0 2.5]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick latency');
    title(['By Animal MGB Choice']);
    xticklabels({'hit', 'fa','hit','fa'});

    lickFig.Position(3:4)=[725 475];
    saveas(gcf,['ByAnimal_C_MGB_LickLat_Opto']);
    saveas(gcf,['ByAnimal_C_MGB_LickLat_Opto.pdf']);
    close all


    % control animals, RATE
    lickFigl=figure(jj+13); 
    for jj=2:size(allDataCtlOnly,1)
        hitLickRate(jj-1)=nanmean(MGBmlr_hit(jj,:));
        ohitLickRate(jj-1)=nanmean(MGBmlo_hit(jj,:));
        faLickRate(jj-1)=nanmean(MGBmlr_fa(jj,:));
        ofaLickRate(jj-1)=nanmean(MGBmlo_fa(jj,:));
        othitLickRate(jj-1)=nanmean(MGBmlto_hit(jj,:));
        otfaLickRate(jj-1)=nanmean(MGBmlto_fa(jj,:));
        ochitLickRate(jj-1)=nanmean(MGBmlco_hit(jj,:));
        ocfaLickRate(jj-1)=nanmean(MGBmlco_fa(jj,:));

        hitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
        ohitLickRateE(jj-1)=nanmean(MGBslo_hit(jj,:));
        faLickRateE(jj-1)=nanmean(MGBslr_fa(jj,:));
        ofaLickRateE(jj-1)=nanmean(MGBslo_fa(jj,:));
        othitLickRateE(jj-1)=nanmean(MGBslto_hit(jj,:));
        otfaLickRateE(jj-1)=nanmean(MGBslto_fa(jj,:));
        ochitLickRateE(jj-1)=nanmean(MGBslco_hit(jj,:));
        ocfaLickRateE(jj-1)=nanmean(MGBslco_fa(jj,:));    

    end

    subplot(2,3,1); % full trial MGB
    lickbarL=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(ohitLickRate) nanmean(ofaLickRate)]); hold on;
    lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    error=[nanmean(hitLickRateE) nanmean(faLickRateE) nanmean(ohitLickRateE) nanmean(ofaLickRateE)];
    % for kk=1:numel(lickbarL)
    %     xtips=lickbarL(kk).XEndPoints;
    %     ytips=lickbarL(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
        scatter(repmat(lickbarL(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarL(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarL(1).XEndPoints(3),size(ohitLickRate,1),1), ...
            ohitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarL(1).XEndPoints(4),size(ofaLickRate,1),1), ...
            ofaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[ohitLickRate(1) ofaLickRate(1) ohitLickRate(2) ofaLickRate(2)...
            ohitLickRate(3) ofaLickRate(3) ohitLickRate(4) ofaLickRate(4)...
            ohitLickRate(5) ofaLickRate(5) ohitLickRate(6) ofaLickRate(6)];
        line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,ohitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,ofaLickRate);
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Full Trial']);
    xticklabels({'hit', 'fa','hit','fa'});

    subplot(2,3,2) % tone MGB
    lickbarT=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(othitLickRate) nanmean(otfaLickRate)]); hold on;
    lickbarT(1).FaceColor='flat'; lickbarT(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    error=[nanmean(hitLickRateE) nanmean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
    % for kk=1:numel(lickbarT)
    %     xtips=lickbarT(kk).XEndPoints;
    %     ytips=lickbarT(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
    scatter(repmat(lickbarT(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarT(1).XEndPoints(3),size(othitLickRate,1),1), ...
            othitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarT(1).XEndPoints(4),size(otfaLickRate,1),1), ...
            otfaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[othitLickRate(1) otfaLickRate(1) othitLickRate(2) otfaLickRate(2)...
            othitLickRate(3) otfaLickRate(3) othitLickRate(4) otfaLickRate(4)...
            othitLickRate(5) otfaLickRate(5) othitLickRate(6) otfaLickRate(6)];
        line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,othitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Tone']);
    xticklabels({'hit', 'fa','hit','fa'});

    subplot(2,3,3); % choice MGB
    lickbarC=bar([nanmean(hitLickRate) nanmean(faLickRate) nanmean(ochitLickRate) nanmean(ocfaLickRate)]); hold on;
    lickbarC(1).FaceColor='flat'; lickbarC(1).CData=[reinfcolor;reinfcolor;optocolor;optocolor];
    error=[nanmean(hitLickRateE) nanmean(ochitLickRateE) nanmean(faLickRateE) nanmean(ocfaLickRateE)];
    % for kk=1:numel(lickbarC)
    %     xtips=lickbarC(kk).XEndPoints;
    %     ytips=lickbarC(kk).YEndPoints;
    %     errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    % end
    scatter(repmat(lickbarC(1).XEndPoints(1),size(hitLickRate,1),1), ...
            hitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarC(1).XEndPoints(2),size(faLickRate,1),1), ...
            faLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',reinfcolor);
        scatter(repmat(lickbarC(1).XEndPoints(3),size(ochitLickRate,1),1), ...
            ochitLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        scatter(repmat(lickbarC(1).XEndPoints(4),size(ocfaLickRate,1),1), ...
            ocfaLickRate,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',optocolor);
        lickRate = [hitLickRate(1) faLickRate(1) hitLickRate(2) faLickRate(2)...
            hitLickRate(3) faLickRate(3) hitLickRate(4) faLickRate(4) ...
            hitLickRate(5) faLickRate(5) hitLickRate(6) faLickRate(6)];
        oLickRate=[ochitLickRate(1) ocfaLickRate(1) ochitLickRate(2) ocfaLickRate(2)...
            ochitLickRate(3) ocfaLickRate(3) ochitLickRate(4) ocfaLickRate(4)...
            ochitLickRate(5) ocfaLickRate(5) ochitLickRate(6) ocfaLickRate(6)];
        line([1 2],lickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([1 2],lickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(1:2), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(3:4), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(5:6), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(7:8), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(9:10), 'LineWidth', 0.5, 'Color', [0 0 0]);
        line([3 4],oLickRate(11:12), 'LineWidth', 0.5, 'Color', [0 0 0]);
    [h,pHit,ci,stats] = ttest2(hitLickRate,ochitLickRate);
    [h,pFA,ci,stats] = ttest2(faLickRate,ocfaLickRate);
    ylim([0 8]);
    sigstar({[1,3],[2,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title(['By Animal MGB Choice']);
    xticklabels({'hit', 'fa','hit','fa'});

    lickFigl.Position(3:4)=[725 475];
    saveas(gcf,['ByAnimal_C_MGB_LickRate_Opto']);
    saveas(gcf,['ByAnimal_C_MGB_LickRate_Opto.pdf']);  
    close all
end
%% lick latency by percent correct, by session

close all
% LICK LATENCY, AVERAGED ACROSS ANIMALS
% change these variable names, theyre swapped
lickFig=figure(103); 
clear rpc opc
counter=1;
cCounter=1;
hitLickRate=MGBmr_hit(~isnan(MGBmr_hit));
faLickRate=MGBmr_fa(~isnan(MGBmr_fa));
ochitLickRate=MGBmco_hit(~isnan(MGBmco_hit));
ocfaLickRate=MGBmco_fa(~isnan(MGBmco_fa));

ohitLickRate=MGBmo_hit(~isnan(MGBmo_hit));
ofaLickRate=MGBmo_fa(~isnan(MGBmo_fa));
hitFLickRate=MGBmr_hit(~isnan(MGBmo_hit));
faFLickRate=MGBmr_fa(~isnan(MGBmo_fa));

hitTLickRate=MGBmr_hit(~isnan(MGBmto_hit));
faTLickRate=MGBmr_fa(~isnan(MGBmto_fa));
othitLickRate=MGBmto_hit(~isnan(MGBmto_hit));
otfaLickRate=MGBmto_fa(~isnan(MGBmto_fa));

for jj=2:size(allDataTestsOnly,1)
    rpc(counter:counter+1)=allDataTestsOnly{jj,27}(~isnan(allDataTestsOnly{jj,27}));
    opc(counter:counter+1)=allDataTestsOnly{jj,28}(~isnan(allDataTestsOnly{jj,28}));
    rpcTone(counter:counter+1)=allDataTestsOnly{jj,29}(~isnan(allDataTestsOnly{jj,29}));
    opcTone(counter:counter+1)=allDataTestsOnly{jj,30}(~isnan(allDataTestsOnly{jj,30}));
    rpcChoice(cCounter:cCounter+3)=allDataTestsOnly{jj,31}(~isnan(allDataTestsOnly{jj,31}));
    if jj==5
        opcChoice(cCounter:cCounter+3)=allDataTestsOnly{jj,32}(~isnan(allDataTestsOnly{jj,32}));
        opcChoice(cCounter:cCounter+2)=[opcChoice(cCounter) opcChoice(cCounter+2) opcChoice(cCounter+3)];
        opcChoice(cCounter+3)=[];
        counter=counter+2;cCounter=cCounter+4;
        %sk198 has 0 hit rate for one choice session, so can't compute lick latency or lick rate; remove
    else
        opcChoice(cCounter:cCounter+3)=allDataTestsOnly{jj,32}(~isnan(allDataTestsOnly{jj,32}));
        counter=counter+2;cCounter=cCounter+4;
    end
    
end

opcChoice(:,16)=[];
subplot(1,3,1); % full trial, hit
scatter(hitFLickRate, rpc,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ohitLickRate, opc,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitFLickRate, rpc,1); 
Bfit = polyval(Fit, hitFLickRate);
plot(hitFLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ohitLickRate, opc,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ohitLickRate);
plot(ohitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency');
title('Full Trial Hit'); legend('off', 'on','Location', 'Best');


subplot(1,3,2); % tone, hit
scatter(hitTLickRate, rpcTone,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(othitLickRate, opcTone,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitTLickRate, rpcTone,1); 
Bfit = polyval(Fit, hitTLickRate);
plot(hitTLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(othitLickRate, opcTone,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, othitLickRate);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
plot(othitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
ylabel('percent correct'); xlabel('lick latency');
title('Tone Hit'); 


subplot(1,3,3); % choice, hit
scatter(hitLickRate, rpcChoice,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ochitLickRate, opcChoice,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitLickRate, rpcChoice,1); 
Bfit = polyval(Fit, hitLickRate);
plot(hitLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ochitLickRate, opcChoice,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ochitLickRate);
plot(ochitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency');
title('Choice Hit'); lickFig.Position(3:4)=[725 225];
saveas(gcf,['BySess_T_MGB_Hit_LickLatPC_Opto']);
saveas(gcf,['BySess_T_MGB_Hit_LickLatPC_Opto.png']);


close all
lickFig=figure(9);
subplot(1,3,1) % full trial, fa
scatter(faFLickRate, rpc,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ofaLickRate, opc,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faFLickRate, rpc,1); 
Bfit = polyval(Fit, faFLickRate);
plot(faFLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ofaLickRate, opc,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ofaLickRate);
plot(ofaLickRate,Ofit,'Color',optocolor)
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
xlim([0 1.75]);
ylabel('percent correct'); xlabel('lick latency');
title('Full Trial FA'); legend('off', 'on','Location', 'Best');
lickFig.Position(3:4)=[725 225];


subplot(1,3,2) % tone, fa
scatter(faTLickRate, rpcTone,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(otfaLickRate, opcTone,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faTLickRate, rpcTone,1); 
Bfit = polyval(Fit, faTLickRate);
plot(faTLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(otfaLickRate, opcTone,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, otfaLickRate);
plot(otfaLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency');
title('Tone'); 
lickFig.Position(3:4)=[425 225];

subplot(1,3,3) % choice, fa
scatter(faLickRate, rpcChoice,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ocfaLickRate, opcChoice,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faLickRate, rpcChoice,1); 
Bfit = polyval(Fit, faLickRate);
plot(faLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ocfaLickRate, opcChoice,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ocfaLickRate);
plot(ocfaLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency');
title('Choice FA'); 
lickFig.Position(3:4)=[725 225];
saveas(gcf,['BySess_T_MGB_FA_LickLatPC_Opto']);
saveas(gcf,['BySess_T_MGB_FA_LickLatPC_Opto.png']);
close all


%%
% by animal


close all
% LICK LATENCY, AVERAGED ACROSS ANIMALS
% change these variable names, theyre swapped
lickFig=figure(103); 
clear rpc opc rpcTone opcTone rpcChoice opcChoice hitLickRate ohitLickRate faLickRate othitLickRate
clear ofaLickRate otfaLickRate ochitLickRate ocfaLickRate

for jj=2:size(allDataTestsOnly,1)
    hitLickRate(jj-1)=nanmean(MGBmr_hit(jj,:));
    ohitLickRate(jj-1)=nanmean(MGBmo_hit(jj,:));
    faLickRate(jj-1)=nanmean(MGBmr_fa(jj,:));
    ofaLickRate(jj-1)=nanmean(MGBmo_fa(jj,:));
    othitLickRate(jj-1)=nanmean(MGBmto_hit(jj,:));
    otfaLickRate(jj-1)=nanmean(MGBmto_fa(jj,:));
    ochitLickRate(jj-1)=nanmean(MGBmco_hit(jj,:));
    ocfaLickRate(jj-1)=nanmean(MGBmco_fa(jj,:));

    rpc(jj-1,:)=nanmean(allDataTestsOnly{jj,27}(1:4));
    opc(jj-1,:)=nanmean(allDataTestsOnly{jj,28}(1:4));
    rpcTone(jj-1,:)=nanmean(allDataTestsOnly{jj,29}(1:4));
    opcTone(jj-1,:)=nanmean(allDataTestsOnly{jj,30}(1:4));
    rpcChoice(jj-1,:)=nanmean(allDataTestsOnly{jj,31}(1:4));
    opcChoice(jj-1,:)=nanmean(allDataTestsOnly{jj,32}(1:4));

end

subplot(1,3,1); % full trial, hit
scatter(hitLickRate, rpc,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ohitLickRate, opc,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitLickRate, rpc,1); 
Bfit = polyval(Fit, hitLickRate);
plot(hitLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ohitLickRate, opc,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ohitLickRate);
plot(ohitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency');
title('Full Trial Hit'); legend('off', 'on','Location', 'Best');


subplot(1,3,2); % tone, hit
scatter(hitLickRate, rpcTone,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(othitLickRate, opcTone,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitLickRate, rpcTone,1); 
Bfit = polyval(Fit, hitLickRate);
plot(hitLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(othitLickRate, opcTone,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, othitLickRate);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
plot(othitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
ylabel('percent correct'); xlabel('lick latency');
title('Tone Hit'); 


subplot(1,3,3); % choice, hit
scatter(hitLickRate, rpcChoice,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ochitLickRate, opcChoice,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(hitLickRate, rpcChoice,1); 
Bfit = polyval(Fit, hitLickRate);
plot(hitLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ochitLickRate, opcChoice,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ochitLickRate);
plot(ochitLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency (s)');
title('Choice Hit'); lickFig.Position(3:4)=[725 225];
saveas(gcf,['ByA_T_MGB_Hit_LickLatPC_Opto']);
saveas(gcf,['ByA_T_MGB_Hit_LickLatPC_Opto.png']);
saveas(gcf,['ByA_T_MGB_Hit_LickLatPC_Opto.pdf']);

close all
lickFig=figure(9);
subplot(1,3,1) % full trial, fa
scatter(faLickRate, rpc,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ofaLickRate, opc,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faLickRate, rpc,1); 
Bfit = polyval(Fit, faLickRate);
plot(faLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ofaLickRate, opc,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ofaLickRate);
plot(ofaLickRate,Ofit,'Color',optocolor)
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
xlim([0 1.75]);
ylabel('percent correct'); xlabel('lick latency (s)');
title('Full Trial FA'); legend('off', 'on','Location', 'Best');
lickFig.Position(3:4)=[725 225];


subplot(1,3,2) % tone, fa
scatter(faLickRate, rpcTone,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(otfaLickRate, opcTone,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faLickRate, rpcTone,1); 
Bfit = polyval(Fit, faLickRate);
plot(faLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(otfaLickRate, opcTone,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, otfaLickRate);
plot(otfaLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency (s)');
title('Tone'); 
lickFig.Position(3:4)=[425 225];

subplot(1,3,3) % choice, fa
scatter(faLickRate, rpcChoice,'MarkerEdgeColor',reinfcolor); hold on; 
scatter(ocfaLickRate, opcChoice,'MarkerEdgeColor',optocolor); hold on; 
Fit = polyfit(faLickRate, rpcChoice,1); 
Bfit = polyval(Fit, faLickRate);
plot(faLickRate,Bfit,'Color',reinfcolor)
Fit = polyfit(ocfaLickRate, opcChoice,1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Ofit = polyval(Fit, ocfaLickRate);
plot(ocfaLickRate,Ofit,'Color',optocolor);xlim([0 1.75]);
[h,p,ci,stats] = ttest2(Bfit,Ofit);
sigstar({[0.5,1]}, [p])
ylabel('percent correct'); xlabel('lick latency (s)');
title('Choice FA'); 
lickFig.Position(3:4)=[725 225];
saveas(gcf,['ByA_T_MGB_FA_LickLatPC_Opto']);
saveas(gcf,['ByA_T_MGB_FA_LickLatPC_Opto.png']);
saveas(gcf,['ByA_T_MGB_FA_LickLatPC_Opto.pdf']);
close all
%% within subject anova1
% model on light on vs. light off, check animal contribution
% t=table(
% rm = fitrm(t,'year1-year6 ~ meanEmploy','WithinDesign',Year);
% anovatbl=anova(rm,'WithinModel',WM);
