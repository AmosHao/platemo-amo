function p=testmops(problemName,pd,od)
% This function is coded by Zhang, Hu. It includes many multiobjective
% optimization problems. In the function, all the input variables x must be
% row vectors.
%
% INPUT: problemName - the function name.
% OUTPUT: problem - the problem to be invoked.
%
% Author: Zhang Hu, Harbin Institute of Technology Last modified: September
% 25, 2013 Copyright (C) 2013-2016 by Zhang Hu (e-mail: jxzhanghu@126.com)

%%
% Define a structure 'p'
p = struct('name',[],'od',[],'pd',[],'pf',[],'domain',[],'func',[]);
k=2*(od-1);
l=pd-k;
switch lower(problemName)
    % LZ09 test suite
    case 'kno1'
        p=kno1(p,pd);
    case 'planner'
        p=planner(p,pd);
    case 'tec09_f1'
        p=tec09_f1(p,pd);
    case 'tec09_f2'
        p=tec09_f2(p,pd);
    case 'tec09_f3'
        p=tec09_f3(p,pd);
    case 'tec09_f4'
        p=tec09_f4(p,pd);
    case 'tec09_f5'
        p=tec09_f5(p,pd);
    case 'tec09_f6'
        p=tec09_f6(p,pd);
    case 'tec09_f7'
        p=tec09_f7(p,pd);
    case 'tec09_f8'
        p=tec09_f8(p,pd);
    case 'tec09_f9'
        p=tec09_f9(p,pd);
        % DTLZ test suite
    case 'dtlz1'
        p=dtlz1(p,pd,od);
    case 'dtlz2'
        p=dtlz2(p,pd,od);
    case 'dtlz3'
        p=dtlz3(p,pd,od);
    case 'dtlz4'
        p=dtlz4(p,pd,od);
    case 'dtlz5'
        p=dtlz5(p,pd,od);
    case 'dtlz5i'
        p=dtlz5i(p,pd,od);
    case 'dtlz6'
        p=dtlz6(p,pd,od);
    case 'dtlz7'
        p=dtlz7(p,pd,od);
%     case 'viennet2'
%         p=viennet2(p,pdension);
%     case 'viennet3'
%         p=viennet3(p,pdension);
%     case 'fonseca'
%         p=fonseca(p,pd);
%     case 'schaffer'
%         p=schaffer(p,pdension);
%     case 'kursawe'
%         p=kursawe(p,pd);
    % ZDT test suite
    case 'zdt1'
        p=zdt1(p,pd);
    case 'zdt2'
        p=zdt2(p,pd);
    case 'zdt3'
        p=zdt3(p,pd);
    case 'zdt4'
        p=zdt4(p,pd);
    case 'zdt6'
        p=zdt6(p,pd);
    % DMOEA/D test suite
    case 'smeaglt1'
        p=smeaglt1(p,pd);
    case 'smeaglt2'
        p=smeaglt2(p,pd);
    case 'smeaglt3'
        p=smeaglt3(p,pd);
    case 'smeaglt4'
        p=smeaglt4(p,pd);
    case 'smeaglt5'
        p=smeaglt5(p,pd);
    case 'smeaglt6'
        p=smeaglt6(p,pd);
    case 'wfg1'
        p=wfg1(p,od,k,l);
    case 'wfg2'
        p=wfg2(p,od,k,l);        
    case 'wfg3'
        p=wfg3(p,od,k,l);        
    case 'wfg4'
        p=wfg4(p,od,k,l);        
    case 'wfg5'
        p=wfg5(p,od,k,l);        
    case 'wfg6'
        p=wfg6(p,od,k,l);        
    case 'wfg7'
        p=wfg7(p,od,k,l);        
    case 'wfg8'
        p=wfg8(p,od,k,l);        
    case 'wfg9'
        p=wfg9(p,od,k,l);  
    %
    % The GL benchmark is used in Liu Hailin's Paper "A novel weight
    % design in multi-objective evolutionary algorithm"
    case 'gl1'
        p=gl1(p,pd);
    case 'gl2'
        p=gl2(p,pd);
    case 'gl3'
        p=gl3(p,pd);
    case 'gl4'
        p=gl4(p,pd);
    % M2M test suite
    % The test plem is from Hai-Lin Liu's paper "Decomposition of a
    % Multiobjective Optimization plem into a Number of Simple
    % Multiobjective Subplems"
    case 'm2m_f1'
        p=m2m_f1(p,pd);
    case 'm2m_f2'
        p=m2m_f2(p,pd);
    case 'm2m_f3'
        p=m2m_f3(p,pd);
    case 'm2m_f4'
        p=m2m_f4(p,pd);
    case 'm2m_f5'
        p=m2m_f5(p,pd);
    case 'm2m_f6'
        p=m2m_f6(p,pd);
    case 'm2m_f7'
        p=m2m_f7(p,pd);
    % RM-MEDA test suite
    case 'zzj_f1'
        p=zzj_f1(p,pd);
    case 'zzj_f2'
        p=zzj_f2(p,pd);
    case 'zzj_f3'
        p=zzj_f3(p,pd);
    case 'zzj_f4'
        p=zzj_f4(p,pd);
    case 'zzj_f5'
        p=zzj_f5(p,pd);
    case 'zzj_f6'
        p=zzj_f6(p,pd);
    case 'zzj_f7'
        p=zzj_f7(p,pd);
    case 'zzj_f8'
        p=zzj_f8(p,pd);
    case 'zzj_f9'
        p=zzj_f9(p,pd);
    case 'zzj_f10'
        p=zzj_f10(p,pd);
    % IRM-MEDA test suite
    case 'wxc_f1'
        p=wxc_f1(p,pd);
    case 'wxc_f2'
        p=wxc_f2(p,pd);
    case 'wxc_f3'
        p=wxc_f3(p,pd);
    case 'wxc_f4'
        p=wxc_f4(p,pd);
    case 'wxc_f5'
        p=wxc_f5(p,pd);
    case 'wxc_f6'
        p=wxc_f6(p,pd);
    case 'wxc_f7'
        p=wxc_f7(p,pd);
    case 'wxc_f8'
        p=wxc_f8(p,pd);
    case 'wxc_f9'
        p=wxc_f9(p,pd);
    case 'wxc_f10'
        p=wxc_f10(p,pd);
    % CEC09 unconstraint instances
    case 'uf1'
        p = uf1(p,pd);
    case 'uf2'
        p = uf2(p,pd);
    case 'uf3'
        p = uf3(p,pd);
    case 'uf4'
        p = uf4(p,pd);
    case 'uf5'
        p = uf5(p,pd);
    case 'uf6'
        p = uf6(p,pd);
    case 'uf7'
        p = uf7(p,pd);
    case 'uf8'
        p = uf8(p,pd);
    case 'uf9'
        p = uf9(p,pd);
    case 'uf10'
        p = uf10(p,pd);
    case 'eff'
        p=eff(p,pd);
    case 'ate'
        p=ate(p,pd);
    case 'jy_f1'
        p=jy_f1(p,pd);
    case 'jy_f2'
        p=jy_f2(p,pd);
    case 'jy_f3'
        p=jy_f3(p,pd);
    case 'jy_f4'
        p=jy_f4(p,pd);
    case 'jy_mf4'
        p=jy_mf4(p,pd);
    case 'jy_f5'
        p=jy_f5(p,pd);
    case 'jy_f6'
        p=jy_f6(p,pd);
    case 'jy_convexdtlz2'
        p=jy_convexdtlz2(p,pd);
    case 'mjy1'
        p=mjy1(p,pd);
    case 'mjy2'
        p=mjy2(p,pd);
    case 'mjy3'
        p=mjy3(p,pd);
    case 'mjy4'
        p=mjy4(p,pd);
    case 'mjy5'
        p=mjy5(p,pd);
    case 'mjy6'
        p=mjy6(p,pd);
    case 'mzdt3'
        p=mzdt3(p,pd);
    case 'mdtlz5'
        p=mdtlz5(p,pd);
    case 'mdtlz7'
        p=mdtlz7(p,pd);
    otherwise
        error('Undefined test problem name');
end
end