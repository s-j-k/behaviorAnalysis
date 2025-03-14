function lickLatOptoPerAnimalDead(allDataTestsOnly,reinfcolor,optocolor)

SESS = 1; CTXT = 2; TONE = 3; OUTCOME = 4; 
START = 5; STOP = 6; TONE_T = 7; LICKL = 8; LICKR = 9;

% tests for MGB for each animal
nbsubj=2;
for nbsubj=2:size(allDataTestsOnly,1)
    matrix=allDataTestsOnly{nbsubj,49};
    for i=1:max(matrix(:,SESS))
        if i==1 % protocol 1, delay 1 = ctx 6; delay 2 = ctx 5; ctx 2 reinf; ctx 0 probe
            mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
            sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
            mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            md2o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            sd2o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            md1o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
            sd1o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
            mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            md2o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            sd2o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            md1o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
            sd1o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));

            mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld2o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld2o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mld1o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            sld1o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            
            mld2o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld2o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mld1o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            sld1o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
                   
        elseif i==2 || i==3
            mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
            sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
            mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            md3o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            sd3o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            md4o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            sd4o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));      
            mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            md3o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            sd3o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            md4o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            sd4o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            
            mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld3o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld3o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld4o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld4o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld3o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld3o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld4o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld4o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        elseif i==4
            mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
            sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
            mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            md5o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            sd5o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            md1o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
            sd1o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==1,LICKL));
            mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));      
            mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            md5o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            sd5o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            md1o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
            sd1o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==6 & matrix(:,OUTCOME)==3,LICKL));
            
            mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld5o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld5o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld1o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            sld1o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld5o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld5o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld1o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
            sld1o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==6 & matrix(:,SESS)==i,LICKR));
        elseif i==5
            mr_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); % reinf hit
            sr_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==1,LICKL)); 
            mp_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            sp_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==1,LICKL));
            md5o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            sd5o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==1,LICKL));
            md2o_hit(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            sd2o_hit(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==1,LICKL));
            mr_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));
            sr_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==2 & matrix(:,OUTCOME)==3,LICKL));      
            mp_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            sp_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==0 & matrix(:,OUTCOME)==3,LICKL));
            md5o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            sd5o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==1 & matrix(:,OUTCOME)==3,LICKL));
            md2o_fa(nbsubj,i) = nanmean(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            sd2o_fa(nbsubj,i) = nansem(matrix(matrix(:,SESS)==i & matrix(:,CTXT)==5 & matrix(:,OUTCOME)==3,LICKL));
            
            mlr_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld5o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld5o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld2o_hit(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld2o_hit(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==1 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            mlr_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            slr_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==2 & matrix(:,SESS)==i,LICKR));
            mlp_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            slp_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==0 & matrix(:,SESS)==i,LICKR));
            mld5o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            sld5o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==1 & matrix(:,SESS)==i,LICKR));
            mld2o_fa(nbsubj,i) = nanmean(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
            sld2o_fa(nbsubj,i) = nansem(matrix(matrix(:,OUTCOME)==3 & matrix(:,CTXT)==5 & matrix(:,SESS)==i,LICKR));
        end
%             
%         MGBsr_hit(nbsubj,i) = sr_hit(nbsubj,i); MGBsr_hit(MGBsr_hit==0)=NaN;
%         MGBmp_hit(nbsubj,i) = mp_hit(nbsubj,i); MGBmp_hit(MGBmp_hit==0)=NaN;
%         MGBsp_hit(nbsubj,i) = sp_hit(nbsubj,i); MGBsp_hit(MGBsp_hit==0)=NaN;
%         MGBmo_hit(nbsubj,i) = mo_hit(nbsubj,i); MGBmo_hit(MGBmo_hit==0)=NaN;
%         MGBso_hit(nbsubj,i) = so_hit(nbsubj,i); MGBso_hit(MGBso_hit==0)=NaN;
%         MGBmto_hit(nbsubj,i) = mto_hit(nbsubj,i);MGBmto_hit(MGBmto_hit==0)=NaN;
%         MGBsto_hit(nbsubj,i) = sto_hit(nbsubj,i); MGBsto_hit(MGBsto_hit==0)=NaN;
%         MGBmco_hit(nbsubj,i) = mco_hit(nbsubj,i); MGBmco_hit(MGBmco_hit==0)=NaN;
%         MGBsco_hit(nbsubj,i) = sco_hit(nbsubj,i); MGBsco_hit(MGBsco_hit==0)=NaN;
%         MGBmr_hit(nbsubj,i) = mr_hit(nbsubj,i); % reinf hit
%         MGBsr_hit(nbsubj,i) = sr_hit(nbsubj,i); 
%         MGBmp_hit(nbsubj,i) = mp_hit(nbsubj,i);
%         MGBsp_hit(nbsubj,i) = sp_hit(nbsubj,i);
%         MGBmo_hit(nbsubj,i) = mo_hit(nbsubj,i);
%         MGBso_hit(nbsubj,i) = so_hit(nbsubj,i);
%         MGBmto_hit(nbsubj,i) = mto_hit(nbsubj,i);
%         MGBsto_hit(nbsubj,i) = sto_hit(nbsubj,i);
%         MGBmco_hit(nbsubj,i) = mco_hit(nbsubj,i);
%         MGBsco_hit(nbsubj,i) = sco_hit(nbsubj,i);
% 
%         MGBmr_fa(nbsubj,i) = mr_fa(nbsubj,i);
%         MGBsr_fa(nbsubj,i) = sr_fa(nbsubj,i);
%         MGBmp_fa(nbsubj,i) = mp_fa(nbsubj,i);
%         MGBsp_fa(nbsubj,i) = sp_fa(nbsubj,i);
%         MGBmo_fa(nbsubj,i) = mo_fa(nbsubj,i);
%         MGBso_fa(nbsubj,i) = so_fa(nbsubj,i);
%         MGBmto_fa(nbsubj,i) = mto_fa(nbsubj,i);
%         MGBsto_fa(nbsubj,i) = sto_fa(nbsubj,i);
%         MGBmco_fa(nbsubj,i) = mco_fa(nbsubj,i);
%         MGBsco_fa(nbsubj,i) = sco_fa(nbsubj,i);
% 
%         MGBmlr_hit(nbsubj,i) = mlr_hit(nbsubj,i);
%         MGBslr_hit(nbsubj,i) = slr_hit(nbsubj,i);
%         MGBmlp_hit(nbsubj,i) = mlp_hit(nbsubj,i);
%         MGBslp_hit(nbsubj,i) = slp_hit(nbsubj,i);
%         MGBmlo_hit(nbsubj,i) = mlo_hit(nbsubj,i);
%         MGBslo_hit(nbsubj,i) = slo_hit(nbsubj,i);
%         MGBmlto_hit(nbsubj,i) = mlto_hit(nbsubj,i);
%         MGBslto_hit(nbsubj,i) = slto_hit(nbsubj,i);
%         MGBmlco_hit(nbsubj,i) = mlco_hit(nbsubj,i);
%         MGBslco_hit(nbsubj,i) = slco_hit(nbsubj,i);
% 
%         MGBmlr_fa(nbsubj,i) = mlr_fa(nbsubj,i);
%         MGBslr_fa(nbsubj,i) = slr_fa(nbsubj,i);
%         MGBmlp_fa(nbsubj,i) = mlp_fa(nbsubj,i);
%         MGBslp_fa(nbsubj,i) = slp_fa(nbsubj,i);
%         MGBmlo_fa(nbsubj,i) = mlo_fa(nbsubj,i);
%         MGBslo_fa(nbsubj,i) = slo_fa(nbsubj,i);
%         MGBmlto_fa(nbsubj,i) = mlto_fa(nbsubj,i);
%         MGBslto_fa(nbsubj,i) = slto_fa(nbsubj,i);
%         MGBmlco_fa(nbsubj,i) = mlco_fa(nbsubj,i);
%         MGBslco_fa(nbsubj,i) = slco_fa(nbsubj,i);
    end      
end

mr_hit(mr_hit==0) = NaN;
sr_hit(sr_hit==0) = NaN;
mp_hit(mp_hit==0) = NaN;
sp_hit(sp_hit==0) =NaN;
md5o_hit(md5o_hit==0) = NaN;
sd5o_hit(sd5o_hit==0) =NaN;
md2o_hit(md2o_hit==0) = NaN;
sd2o_hit(sd2o_hit==0) =NaN;
md3o_hit(md3o_hit==0) = NaN;
sd3o_hit(sd3o_hit==0) =NaN;
md4o_hit(md4o_hit==0) = NaN;
sd4o_hit(sd4o_hit==0) =NaN;
md1o_hit(md1o_hit==0) = NaN;
sd1o_hit(sd1o_hit==0) =NaN;
mr_fa(mr_fa==0) =NaN;
sr_fa(sr_fa==0) = NaN;
mp_fa(mp_fa==0) = NaN;
sp_fa(sp_fa==0) = NaN;
md5o_fa(md5o_fa==0) = NaN;
sd5o_fa(sd5o_fa==0) = NaN;
md2o_fa(md2o_fa==0) = NaN;
sd2o_fa(sd2o_fa==0) = NaN;
md3o_fa(md3o_fa==0) = NaN;
sd3o_fa(sd3o_fa==0) = NaN;
md4o_fa(md4o_fa==0) = NaN;
sd4o_fa(sd4o_fa==0) = NaN;
md1o_fa(md1o_fa==0) = NaN;
sd1o_fa(sd1o_fa==0) = NaN;

mlr_hit(mlr_hit==0) = NaN;
slr_hit(slr_hit==0) = NaN;
mlp_hit(mlp_hit==0) = NaN;
slp_hit(slp_hit==0) =NaN;
mld5o_hit(mld5o_hit==0) = NaN;
sld5o_hit(sld5o_hit==0) =NaN;
mld2o_hit(mld2o_hit==0) = NaN;
sld2o_hit(sld2o_hit==0) =NaN;
mld3o_hit(mld3o_hit==0) = NaN;
sld3o_hit(sld3o_hit==0) =NaN;
mld4o_hit(mld4o_hit==0) = NaN;
sld4o_hit(sld4o_hit==0) =NaN;
mld1o_hit(mld4o_hit==0) = NaN;
sld1o_hit(sld4o_hit==0) =NaN;
mlr_fa(mlr_fa==0) =NaN;
slr_fa(slr_fa==0) = NaN;
mlp_fa(mlp_fa==0) = NaN;
slp_fa(slp_fa==0) = NaN;
mld5o_fa(mld5o_fa==0) = NaN;
sld5o_fa(sld5o_fa==0) = NaN;
mld2o_fa(mld2o_fa==0) = NaN;
sld2o_fa(sld2o_fa==0) = NaN;
mld3o_fa(mld3o_fa==0) = NaN;
sld3o_fa(sld3o_fa==0) = NaN;
mld4o_fa(mld4o_fa==0) = NaN;
sld4o_fa(sld4o_fa==0) = NaN;
mld1o_fa(mld1o_fa==0) = NaN;
sld1o_fa(sld1o_fa==0) = NaN;
% 
% ICmr_hit(ICmr_hit==0)=NaN;
% ICsr_hit(ICsr_hit==0)=NaN;
% ICmp_hit(ICmp_hit==0)=NaN;
% ICsp_hit(ICsp_hit==0)=NaN;
% ICmo_hit(ICmo_hit==0)=NaN; 
% ICso_hit(ICso_hit==0)=NaN; 
% ICmto_hit(ICmto_hit==0)=NaN;
% ICsto_hit(ICsto_hit==0)=NaN; 
% ICmco_hit(ICmco_hit==0)=NaN; 
% ICsco_hit(ICsco_hit==0)=NaN;
% 
% ICmr_fa(ICmr_fa==0) = NaN;
% ICsr_fa(ICsr_fa==0) = NaN;
% ICmp_fa(ICmp_fa==0) = NaN;
% ICsp_fa(ICsp_fa==0) = NaN;
% ICmo_fa(ICmo_fa==0) = NaN;
% ICso_fa(ICso_fa==0) = NaN;
% ICmto_fa(ICmto_fa==0) = NaN;
% ICsto_fa(ICsto_fa==0) = NaN;
% ICmco_fa(ICmco_fa==0) = NaN;
% ICsco_fa(ICsco_fa==0) = NaN;
% 
% ICmlr_hit(ICmlr_hit==0) = NaN;
% ICslr_hit(ICslr_hit==0) = NaN;
% ICmlp_hit(ICmlp_hit==0) = NaN;
% ICslp_hit(ICslp_hit==0) = NaN;
% ICmlo_hit(ICmlo_hit==0) = NaN;
% ICslo_hit(ICslo_hit==0) = NaN;
% ICmlto_hit(ICmlto_hit==0) = NaN;
% ICslto_hit(ICslto_hit==0) = NaN;
% ICmlco_hit(ICmlco_hit==0) = NaN;
% ICslco_hit(ICslco_hit==0) = NaN;
% 
% ICmlr_fa(ICmlr_fa==0) = NaN;
% ICslr_fa(ICslr_fa==0) = NaN;
% ICmlp_fa(ICmlp_fa==0) = NaN;
% ICslp_fa(ICslp_fa==0) = NaN;
% ICmlo_fa(ICmlo_fa==0) = NaN;
% ICslo_fa(ICslo_fa==0) = NaN;
% ICmlto_fa(ICmlto_fa==0) = NaN;
% ICslto_fa(ICslto_fa==0) = NaN;
% ICmlco_fa(ICmlco_fa==0) = NaN;
% ICslco_fa(ICslco_fa==0) = NaN; 

%% plots for individual animals
jj=2;     
for jj=2:size(allDataTestsOnly,1)
    lickFig=figure(jj+103); 
    subplot(2,3,1); hold on; %MGB delay 1
    mr_hit_temp = mr_hit(~isnan(md1o_hit));
    mr_fa_temp = mr_fa(~isnan(md1o_fa));
    sr_hit_temp = sr_hit(~isnan(sd1o_hit));
    sr_fa_temp = sr_fa(~isnan(sd1o_fa));
    lickbar=bar([nanmean(mr_hit_temp(jj,:)) nanmean(md1o_hit(jj,:)) nanmean(mr_fa_temp(jj,:)) nanmean(md1o_fa(jj,:))]); hold on;
    error=[nanmean(sr_hit_temp(jj,:)) nanmean(sd1o_hit(jj,:)) nanmean(sr_fa_temp(jj,:)) nanmean(sd1o_fa(jj,:))];
    lickbar(1).FaceColor='flat'; lickbar(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     for kk=1:numel(lickbar)
%         xtips=lickbar(kk).XEndPoints;
%         ytips=lickbar(kk).YEndPoints;
%         errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
%     end
    [h,pHit,ci,stats] = ttest2(mr_hit_temp(jj,:),md1o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mr_fa_temp(jj,:),md1o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2]);
    box off;
    title([allDataTestsOnly{jj,1} ' MGB Delay 1 Inactivation']);
    xticklabels({'','hit', 'fa','hit','fa'});

    subplot(2,3,2) % delay 2 MGB
    mr_hit_temp = mr_hit(~isnan(md2o_hit));
    mr_fa_temp = mr_fa(~isnan(md2o_fa));
    sr_hit_temp = sr_hit(~isnan(sd2o_hit));
    sr_fa_temp = sr_fa(~isnan(sd2o_fa));
    lickbar2=bar([nanmean(mr_hit_temp(jj,:)) nanmean(md2o_hit(jj,:)) nanmean(mr_fa_temp(jj,:)) nanmean(md2o_fa(jj,:))]); hold on;
    error=[nanmean(sr_hit_temp(jj,:)) nanmean(sd2o_hit(jj,:)) nanmean(sr_fa_temp(jj,:)) nanmean(sd2o_fa(jj,:))];
    lickbar2(1).FaceColor='flat'; lickbar2(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     for kk=1:numel(lickbar2)
%         xtips=lickbar2(kk).XEndPoints;
%         ytips=lickbar2(kk).YEndPoints;
%         errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
%     end
    [h,pHit,ci,stats] = ttest2(mr_hit_temp(jj,:),md2o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mr_fa_temp(jj,:),md2o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2]);box off;
    title([allDataTestsOnly{jj,1} ' Delay 2']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,3); % delay 3
    mr_hit_temp = mr_hit(~isnan(md3o_hit));
    mr_fa_temp = mr_fa(~isnan(md3o_fa));
    sr_hit_temp = sr_hit(~isnan(sd3o_hit));
    sr_fa_temp = sr_fa(~isnan(sd3o_fa));
    lickbar3=bar([nanmean(mr_hit_temp(jj,:)) nanmean(md3o_hit(jj,:)) nanmean(mr_fa_temp(jj,:)) nanmean(md3o_fa(jj,:))]); hold on;
    error=[nanmean(sr_hit_temp(jj,:)) nanmean(sd3o_hit(jj,:)) nanmean(sr_fa_temp(jj,:)) nanmean(sd3o_fa(jj,:))];
    lickbar3(1).FaceColor='flat'; lickbar3(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     for kk=1:numel(lickbar3)
%         xtips=lickbar3(kk).XEndPoints;
%         ytips=lickbar3(kk).YEndPoints;
%         errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
%     end
    [h,pHit,ci,stats] = ttest2(mr_hit_temp(jj,:),md3o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mr_fa_temp(jj,:),md3o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2]);box off;
    title([allDataTestsOnly{jj,1} ' Delay 3']);
    xticklabels({'hit', 'hit','fa','fa'});


    subplot(2,3,4); hold on; %delay 4
    mr_hit_temp = mr_hit(~isnan(md4o_hit));
    mr_fa_temp = mr_fa(~isnan(md4o_fa));
    sr_hit_temp = sr_hit(~isnan(sd4o_hit));
    sr_fa_temp = sr_fa(~isnan(sd4o_fa));
    lickbar4=bar([nanmean(mr_hit_temp(jj,:)) nanmean(md4o_hit(jj,:)) nanmean(mr_fa_temp(jj,:)) nanmean(md4o_fa(jj,:))]); hold on;
    error=[nanmean(sr_hit_temp(jj,:)) nanmean(sd4o_hit(jj,:)) nanmean(sr_fa_temp(jj,:)) nanmean(sd4o_fa(jj,:))];
    lickbar4(1).FaceColor='flat'; lickbar4(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     for kk=1:numel(lickbar4)
%         xtips=lickbar4(kk).XEndPoints;
%         ytips=lickbar4(kk).YEndPoints;
%         errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
%     end
    [h,pHit,ci,stats] = ttest2(mr_hit_temp(jj,:),md4o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mr_fa_temp(jj,:),md4o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylim([0 2]);
    ylabel('mean lick latency');box off;
    title([allDataTestsOnly{jj,1} ' Delay 4']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,5) % delay 5
    mr_hit_temp = mr_hit(~isnan(md5o_hit));
    mr_fa_temp = mr_fa(~isnan(md5o_fa));
    sr_hit_temp = sr_hit(~isnan(sd5o_hit));
    sr_fa_temp = sr_fa(~isnan(sd5o_fa));
    lickbar5=bar([nanmean(mr_hit_temp(jj,:)) nanmean(md5o_hit(jj,:)) nanmean(mr_fa_temp(jj,:)) nanmean(md5o_fa(jj,:))]); hold on;
    error=[nanmean(sr_hit_temp(jj,:)) nanmean(sd5o_hit(jj,:)) nanmean(sr_fa_temp(jj,:)) nanmean(sd5o_fa(jj,:))];
    lickbar5(1).FaceColor='flat'; lickbar5(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
%     for kk=1:numel(lickbar5)
%         xtips=lickbar5(kk).XEndPoints;
%         ytips=lickbar5(kk).YEndPoints;
%         errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
%     end
    [h,pHit,ci,stats] = ttest2(mr_hit_temp(jj,:),md5o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mr_fa_temp(jj,:),md5o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick latency');ylim([0 2]);
    box off;
    title([allDataTestsOnly{jj,1} ' Delay 5']);
    xticklabels({'hit', 'hit','fa','fa'});


    lickFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Delay_LickLat_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_MGB_Delay_LickLat_Opto.png']);

    % RATE
    lickFigl=figure(jj+13); subplot(2,3,1); % full trial MGB
    mlr_hit_temp = mlr_hit(~isnan(mld1o_hit));
    mlr_fa_temp = mlr_fa(~isnan(mld1o_fa));
    slr_hit_temp = slr_hit(~isnan(sld1o_hit));
    slr_fa_temp = slr_fa(~isnan(sld1o_fa));
    lickbarL=bar([nanmean(mlr_hit_temp(jj,:)) nanmean(mld1o_hit(jj,:)) nanmean(mlr_fa_temp(jj,:)) nanmean(mld1o_fa(jj,:))]); hold on;
    lickbarL(1).FaceColor='flat'; lickbarL(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(slr_hit_temp(jj,:)) nanmean(sld1o_hit(jj,:)) nanmean(slr_fa_temp(jj,:)) nanmean(sld1o_fa(jj,:))];
    for kk=1:numel(lickbarL)
        xtips=lickbarL(kk).XEndPoints;
        ytips=lickbarL(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(mlr_hit_temp(jj,:),mld1o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mlr_fa_temp(jj,:),mld1o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    box off;ylim([0 10]);
    title([allDataTestsOnly{jj,1} ' Delay 1 Inactivation']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,2) % delay 2
    mlr_hit_temp = mlr_hit(~isnan(mld2o_hit));
    mlr_fa_temp = mlr_fa(~isnan(mld2o_fa));
    slr_hit_temp = slr_hit(~isnan(sld2o_hit));
    slr_fa_temp = slr_fa(~isnan(sld2o_fa));
    lickbar2L=bar([nanmean(mlr_hit_temp(jj,:)) nanmean(mld2o_hit(jj,:)) nanmean(mlr_fa_temp(jj,:)) nanmean(mld2o_fa(jj,:))]); hold on;
    lickbar2L(1).FaceColor='flat'; lickbar2L(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(slr_hit_temp(jj,:)) nanmean(sld2o_hit(jj,:)) nanmean(slr_fa_temp(jj,:)) nanmean(sld2o_fa(jj,:))];
    for kk=1:numel(lickbar2L)
        xtips=lickbar2L(kk).XEndPoints;
        ytips=lickbar2L(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(mlr_hit_temp(jj,:),mld2o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mlr_fa_temp(jj,:),mld2o_fa(jj,:));
    box off;ylim([0 10]);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' Delay 2']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,3); % delay 3
    mlr_hit_temp = mlr_hit(~isnan(mld3o_hit));
    mlr_fa_temp = mlr_fa(~isnan(mld3o_fa));
    slr_hit_temp = slr_hit(~isnan(sld3o_hit));
    slr_fa_temp = slr_fa(~isnan(sld3o_fa));
    lickbar3L=bar([nanmean(mlr_hit_temp(jj,:)) nanmean(mld3o_hit(jj,:)) nanmean(mlr_fa_temp(jj,:)) nanmean(mld3o_fa(jj,:))]); hold on;
    lickbar3L(1).FaceColor='flat'; lickbar3L(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(slr_hit_temp(jj,:)) nanmean(sld3o_hit(jj,:)) nanmean(slr_fa_temp(jj,:)) nanmean(sld3o_fa(jj,:))];
    for kk=1:numel(lickbar3L)
        xtips=lickbar3L(kk).XEndPoints;
        ytips=lickbar3L(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(mlr_hit_temp(jj,:),mld2o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mlr_fa_temp(jj,:),mld2o_fa(jj,:));
    sigstar({[1,2],[3,4]}, [pHit pFA])
    box off;ylim([0 10]);
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' Delay 3']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,4); % Delay 4
    mlr_hit_temp = mlr_hit(~isnan(mld4o_hit));
    mlr_fa_temp = mlr_fa(~isnan(mld4o_fa));
    slr_hit_temp = slr_hit(~isnan(sld4o_hit));
    slr_fa_temp = slr_fa(~isnan(sld4o_fa));
    lickbar4L=bar([nanmean(mlr_hit_temp(jj,:)) nanmean(mld4o_hit(jj,:)) nanmean(mlr_fa_temp(jj,:)) nanmean(mld4o_fa(jj,:))]); hold on;
    lickbar4L(1).FaceColor='flat'; lickbar4L(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(slr_hit_temp(jj,:)) nanmean(sld4o_hit(jj,:)) nanmean(slr_fa_temp(jj,:)) nanmean(sld4o_fa(jj,:))];
    for kk=1:numel(lickbar4L)
        xtips=lickbar4L(kk).XEndPoints;
        ytips=lickbar4L(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(mlr_hit_temp(jj,:),mld4o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mlr_fa_temp(jj,:),mld4o_fa(jj,:));
    box off;ylim([0 10]);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' Delay 4']);
    xticklabels({'hit', 'hit','fa','fa'});

    subplot(2,3,5) % Delay 5
    mlr_hit_temp = mlr_hit(~isnan(mld5o_hit));
    mlr_fa_temp = mlr_fa(~isnan(mld5o_fa));
    slr_hit_temp = slr_hit(~isnan(sld5o_hit));
    slr_fa_temp = slr_fa(~isnan(sld5o_fa));
    lickbar5L=bar([nanmean(mlr_hit_temp(jj,:)) nanmean(mld5o_hit(jj,:)) nanmean(mlr_fa_temp(jj,:)) nanmean(mld5o_fa(jj,:))]); hold on;
    lickbar5L(1).FaceColor='flat'; lickbar5L(1).CData=[reinfcolor;optocolor;reinfcolor;optocolor];
    error=[nanmean(slr_hit_temp(jj,:)) nanmean(sld5o_hit(jj,:)) nanmean(slr_fa_temp(jj,:)) nanmean(sld5o_fa(jj,:))];
    for kk=1:numel(lickbar5L)
        xtips=lickbar5L(kk).XEndPoints;
        ytips=lickbar5L(kk).YEndPoints;
        errorbar(xtips,ytips,error,'.k','MarkerSize',0.1);
    end
    [h,pHit,ci,stats] = ttest2(mlr_hit_temp(jj,:),mld5o_hit(jj,:));
    [h,pFA,ci,stats] = ttest2(mlr_fa_temp(jj,:),mld5o_fa(jj,:));
    box off;ylim([0 10]);
    sigstar({[1,2],[3,4]}, [pHit pFA])
    ylabel('mean lick rate');
    title([allDataTestsOnly{jj,1} ' Delay 5']);
    xticklabels({'hit', 'hit','fa','fa'});

    lickFig.Position(3:4)=[725 475];
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_Delay_LickRate_Opto']);
    saveas(gcf,[char(allDataTestsOnly{jj,1}) '_T_Delay_LickRate_Opto.png']);  
end


close all
% LICK LATENCY, AVERAGED ACROSS ANIMALS
% change these variable names, theyre swapped
lickFig=figure(jj+103); 
for jj=2:size(mgbTempTestsOnly,1)
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
title(['By Animal MGB Full Trial']);
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
[h,pHit,~,stats] = ttest2(hitLickRate,othitLickRate);
[h,pFA,ci,stats] = ttest2(faLickRate,otfaLickRate);
sigstar({[1,2],[3,4]}, [pHit pFA])
ylabel('mean lick latency');
title(['By Animal MGB Tone']);
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
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});

icIdx=[3,5,6,7];icIdx=icIdx-1;
clear hitLickRate ohitLickRate faLickRate ofaLickRate othitLickRate otfaLickRate ochitLickRate ocfaLickRate
for jj=2:length(icIdx)
    hitLickRate(jj)=nanmean(ICmr_hit(icIdx(jj),:));
    ohitLickRate(jj)=nanmean(ICmo_hit(icIdx(jj),:));
    faLickRate(jj)=nanmean(ICmr_fa(icIdx(jj),:));
    ofaLickRate(jj)=nanmean(ICmo_fa(icIdx(jj),:));
    othitLickRate(jj)=nanmean(ICmto_hit(icIdx(jj),:));
    otfaLickRate(jj)=nanmean(ICmto_fa(icIdx(jj),:));
    ochitLickRate(jj)=nanmean(ICmco_hit(icIdx(jj),:));
    ocfaLickRate(jj)=nanmean(ICmco_fa(icIdx(jj),:));
    
    hitLickRateE(jj)=nanmean(ICsr_hit(icIdx(jj),:));
    ohitLickRateE(jj)=nanmean(ICso_hit(icIdx(jj),:));
    faLickRateE(jj)=nanmean(ICsr_fa(icIdx(jj),:));
    ofaLickRateE(jj)=nanmean(ICso_fa(icIdx(jj),:));
    othitLickRateE(jj)=nanmean(ICsto_hit(icIdx(jj),:));
    otfaLickRateE(jj)=nanmean(ICsto_fa(icIdx(jj),:));
    ochitLickRateE(jj)=nanmean(ICsco_hit(icIdx(jj),:));
    ocfaLickRateE(jj)=nanmean(ICsco_fa(icIdx(jj),:));
    
end

subplot(2,3,4); % full trial IC
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
title(['By Animal IC Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,5) % tone IC
lickbarT=bar([nanmean(hitLickRate) nanmean(othitLickRate) nanmean(faLickRate) nanmean(otfaLickRate)]); hold on;
error=[nanmean(hitLickRateE) mean(othitLickRateE) nanmean(faLickRateE) nanmean(otfaLickRateE)];
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
title(['By Animal IC Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,6); % choice IC
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
title(['By Animal IC Choice']);
xticklabels({'hit', 'hit','fa','fa'});

lickFig.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto']);
saveas(gcf,['ByAnimal_T_MGB_IC_LickLat_Opto.png']);
close all
% RATE
% need to edit this part! 10/18/24
lickFigl=figure(jj+13); 
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
title(['By Animal MGB Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,2) % tone MGB
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
title(['By Animal MGB Tone']);
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
ylabel('mean lick rate');
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});

% do the IC

icIdx=[3,5,6,7];
for jj=1:length(icIdx)
    hitLickRate(jj)=nanmean(ICmlo_hit(icIdx(jj),:));
    ohitLickRate(jj)=nanmean(ICmlo_hit(icIdx(jj),:));
    faLickRate(jj)=nanmean(ICmlr_fa(icIdx(jj),:));
    ofaLickRate(jj)=nanmean(ICmlo_fa(icIdx(jj),:));
    othitLickRate(jj)=nanmean(ICmlto_hit(icIdx(jj),:));
    otfaLickRate(jj)=nanmean(ICmlto_fa(icIdx(jj),:));
    ochitLickRate(jj)=nanmean(ICmlco_hit(icIdx(jj),:));
    ocfaLickRate(jj)=nanmean(ICmlco_fa(icIdx(jj),:));
    %error
    hitLickRateE(jj)=nanmean(ICslo_hit(icIdx(jj),:));
    ohitLickRateE(jj)=nanmean(ICslo_hit(icIdx(jj),:));
    faLickRateE(jj)=nanmean(ICslr_fa(icIdx(jj),:));
    ofaLickRateE(jj)=nanmean(ICslo_fa(icIdx(jj),:));
    othitLickRateE(jj)=nanmean(ICslto_hit(icIdx(jj),:));
    otfaLickRateE(jj)=nanmean(ICslto_fa(icIdx(jj),:));
    ochitLickRateE(jj)=nanmean(ICslco_hit(icIdx(jj),:));
    ocfaLickRateE(jj)=nanmean(ICslco_fa(icIdx(jj),:));    
    
end

subplot(2,3,4); % full trial IC
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
title(['By Animal IC Full Trial']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,5) % tone IC
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
title(['By Animal IC Tone']);
xticklabels({'hit', 'hit','fa','fa'});

subplot(2,3,6); % choice IC
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
title(['By Animal IC Choice']);
xticklabels({'hit', 'hit','fa','fa'});

lickFigl.Position(3:4)=[725 475];
saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto']);
saveas(gcf,['ByAnimal_T_MGB_IC_LickRate_Opto.png']);  
close all

% by session
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
title(['By Animal MGB Full Trial']);
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
title(['By Animal MGB Tone']);
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
title(['By Animal MGB Choice']);
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
title(['By Animal MGB Full Trial']);
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
title(['By Animal MGB Tone']);
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
title(['By Animal MGB Choice']);
xticklabels({'hit', 'hit','fa','fa'});


saveas(gcf,['BySess_T_MGB_LickRateLat_Opto']);
saveas(gcf,['BySess_T_MGB_LickRateLat_Opto.png']);  


end