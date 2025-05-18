function [subjlist, explist, pathsave]=getCohort(cohort)

switch cohort
    case 1
    subjlist={'sk158','sk160','sk161', 'sk162','sk163','sk164'}; %PV Cre
    explist=[2 2 2 2 1 1]; % 2 is ctl
    pathsave='O:\sjk\Behavior\optoMGBIC_1\';
    
    case 2
    subjlist={'sk170','sk175'}; %GTACR
    explist=[2 1];
    pathsave='O:\sjk\Behavior\OptoMGBIC_2\';

    case 3
    pathsave='O:\sjk\Behavior\optoMGBIC_3\';
    subjlist={'sk176','sk177','sk178','sk179'}; %GTACR
    explist=[1 2 1 2];
    
    case 4
    pathsave='O:\sjk\Behavior\MGBIC_4\';
    subjlist={'sk182','sk183','sk184','sk185','sk186','sk187','sk188','sk189'}; %GTACR
    explist=[1 2 1 1 1 1 1 2];
    
    case 5
    pathsave='O:\sjk\Behavior\MGBIC_5\';
    subjlist={'sk190','sk191','sk192','sk193','sk194','sk195','sk196','sk197'}; %GTACR
%     explist=[2 2 1 1 1 1 1 1];
    explist=[1 1 1 1 1 2 1 1];
    
    case 6
    pathsave='O:\sjk\Behavior\MGBIC_6\';
    subjlist={'sk198','sk199','sk200','sk201','sk202','sk203','sk204','sk205'}; %GTACR
    explist=[1 1 1 1 2 1 1 1];    
    
   case 7
    pathsave='O:\sjk\Behavior\IC_cohort_1\';
    subjlist={'sk218','sk219','sk220'}; %GTACR
    explist=[1 1 1];    

    case 8
    pathsave='O:\sjk\Behavior\IC_cohort_2\';
    subjlist={'sk223','sk224','sk225','sk226'}; %GTACR
    explist=[1 1 1 1];
   
    case 9
    pathsave='O:\sjk\Behavior\cohort_7\';
    subjlist={'sk230','sk231','sk232','sk233','sk234','sk235','sk236','sk237'}; %GTACR
    explist=[1 1 1 1 1 1 1 1];
    
    case 10
    pathsave='O:\sjk\Behavior\cohort_8\';
    subjlist={'sk246','sk247','sk248','sk249','sk250','sk251','sk252','sk253'}; %dreads - 5, window - 1, fibers - 3
    explist=[1 1 1 1 1 1 1 1];
    
    
    otherwise
        disp('Cohort not found');    
end