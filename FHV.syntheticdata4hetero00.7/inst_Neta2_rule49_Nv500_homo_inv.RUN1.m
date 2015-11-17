function inst_Neta2_rule49_Nv500_homo_inv()
%function inst_Neta2_rule49_Nv500_homo_inv()

%%%%%%%%%%%%%%%%%%%%
outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.homo.inv';
fprintf(1,'inst_Neta2_rule49_Nv500_homo_inv: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%
Neta=2;
deltachi=[4.7 4.7]; %image sampling intervals in Angstroms
R2=197.4; %outer radius R2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.  To do nothing, specify cmd=[];.  To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
%%%%%%%%%%%%%%%%%%%%
cmd=cell(100,1); %Preallocate for 100 operators.
ii=0;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_setpseudorandomnumberseed';
cmd{ii}.pseudorandomnumberseed=29831;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_readimagestack';
cmd{ii}.imagestackformat='mrc';
cmd{ii}.fn_imagestack='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.fw.imagestack.mrc';
cmd{ii}.startSlice=1;
cmd{ii}.numSlices=500;
%%%%%%%%%%%%%%%%%%%%NOT IN pre4hetero_inv
ii=ii+1;
cmd{ii}.operator='basic_setsizeof2Drealspaceimagesfromimagestack';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_set2Drealspacesamplingintervals';
cmd{ii}.samplingintervals=deltachi;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_annulusstatistics';
cmd{ii}.radius01=[R2+deltachi(1) R2+10*deltachi(1)];
%%%%%%%%%%%%%%%%%%%%NOT IN pre4hetero_inv
ii=ii+1;
cmd{ii}.operator='basic_compute2Dreciprocalspaceproperties';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='realrecip_2DFFT';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}='FHV.ic.lmax0pmax1.clnp.c001iszero.txt';
cmd{ii}.fn_clnp{2}='FHV.ic.lmax0pmax1.clnp.perturbed.c001iszero.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}='FHV.ic.lmax0pmax1.nu.homogeneous.txt';
cmd{ii}.fn_nu{2}='FHV.ic.lmax0pmax1.nu.homogeneous.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='FHV.ic.lmax0pmax1.q.equalclassprobs.txt';
cmd{ii}.fn_q{2}='FHV.ic.lmax0pmax1.q.equalclassprobs.txt';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_print_virusobj';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Increase the number of coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_change_size_of_virusobj';
cmd{ii}.vlmax=[0 0]; %vector of new lmax values, one for each class
cmd{ii}.vpmax=[4 4]; %vector of new pmax values, one for each class
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_print_virusobj';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.fn_tilde_b='callti.out.read_by_C';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_extractdatasubset';
cmd{ii}.kmax=-1;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
cmd{ii}.use_vkminimalset_rather_than_vk=false;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%*******************Solve the linear least squares problem for a one-class spherically-symmetric homogeneous-problem.
ii=ii+1;
cmd{ii}.operator='misc_changedirectory';
cmd{ii}.dn='step0';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_sphericalsymmetry_homogeneous';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_print_virusobj';
%%%%%%%%%%%%%%%%%%%%
%ii=ii+1;
%cmd{ii}.operator='vobj_save_virusobj';
%cmd{ii}.fn_virusobj=[outputbasename '.vobj.mat'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_write_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}=[outputbasename '.eta1.clnp.txt'];
cmd{ii}.fn_clnp{2}=[outputbasename '.eta2.clnp.txt'];
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}=[outputbasename '.eta1.nu.txt'];
cmd{ii}.fn_nu{2}=[outputbasename '.eta2.nu.txt'];
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}=[outputbasename 'eta1.q.txt'];
cmd{ii}.fn_q{2}=[outputbasename 'eta2.q.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_changedirectory';
cmd{ii}.dn='..';
%*******************Done solving the linear least squares problem
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary2.txt'];
%%%%%%%%%%%%%%%%%%%%
%Define the resolution steps used in YinZhengDoerschukNatarajanJohnson JSB 2003 Table 1.
%In the sh/C code used in YinZhengDoerschukNatarajanJohnson JSB 2003, the spatial frequency vector with largest magnitude was computed by directly finding the vector in vkminimalset with maximum magnitude.  vkminimalset was computedd by vk_indexset.m which contains the function index2k.
ZhyeYinresolutionsteps=struct('lmax',[10,15,21,25,31,36,45],'pmax',[6,8,8,10,10,15,20],'Nic',[100,10,10,1,1,1,1],'kmaxpow',[3,3,3,2,1,1,1],'kmax',zeros(1,7));
%deltachi=[4.7 4.7]; %YinZhengDoerschukNatarajanJohnson JSB 2003 FHV example
NaNb=[91 91]; %YinZhengDoerschukNatarajanJohnson JSB 2003 FHV example
kbiggest=sqrt( set_kbiggest(NaNb(1),deltachi(1)).^2 + set_kbiggest(NaNb(2),deltachi(2)).^2 );
ZhyeYinresolutionsteps.kmax=kbiggest*exp((-ZhyeYinresolutionsteps.kmaxpow./2).*log(2));
clear kbiggest
%%%%%%%%%%%%%%%%%%%%
%Scale back the size of the computing for this synthetic problem.
%Do the lmax=0 case because this is a 2 class problem while the initial linear estimator is 1 class.
ZhyeYinresolutionsteps.lmax=[0 6 10];
ZhyeYinresolutionsteps.pmax=[4 5 5];
ZhyeYinresolutionsteps.Nic=[10,2,2];
ZhyeYinresolutionsteps.kmaxpow=ZhyeYinresolutionsteps.kmaxpow(1:3);
ZhyeYinresolutionsteps.kmax=ZhyeYinresolutionsteps.kmax(1:3);
%%%%%%%%%%%%%%%%%%%%
MultiplierForPixelnoisevarIC_eachstep=[10.0 4.0 1.0];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%*******************Begin computing homogeneous reconstructions at increasing resolutions.
for step=1:length(ZhyeYinresolutionsteps.lmax)
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='misc_changedirectory';
  cmd{ii}.dn=['step' num2str(step)];
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='misc_diary';
  cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='vobj_change_size_of_virusobj';
  cmd{ii}.vlmax=[ZhyeYinresolutionsteps.lmax(step) ZhyeYinresolutionsteps.lmax(step)]; %vector of new lmax values, one for each class
  cmd{ii}.vpmax=[ZhyeYinresolutionsteps.pmax(step) ZhyeYinresolutionsteps.pmax(step)]; %vector of new pmax values, one for each class
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='vobj_print_virusobj';
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='EM_read_tilde_b';
  cmd{ii}.fn_tilde_b='callti.out.read_by_C';
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='EM_extractdatasubset';
  cmd{ii}.kmax=ZhyeYinresolutionsteps.kmax(step);
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
  cmd{ii}.use_vkminimalset_rather_than_vk=false;
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='quad_read_integration_rule';
  cmd{ii}.fn_rule='/home/pd83/hetero/hetero00.71/rule_small_3_rulechopper';
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='EM_expectationmaximization';
  cmd{ii}.MC_Nrandic=ZhyeYinresolutionsteps.Nic(step);
  cmd{ii}.MC_FractionOfMeanForMinimum=0.005;
  cmd{ii}.MC_FractionOfMean=0.2;
  cmd{ii}.maxiter=200;
  cmd{ii}.maxiter4pixelnoisevarupdate=5;
  cmd{ii}.cbarftol=-1.0; %unused.
  cmd{ii}.cbarrtol=1.0e-4;
  cmd{ii}.cbarftol_dividebyNc=false;
  cmd{ii}.loglikeftol=-1.0; %unused.
  cmd{ii}.loglikertol=-1.0; %unused.
  cmd{ii}.nu_ic_FractionOfcbarForMinimum=0.005;
  cmd{ii}.nu_ic_FractionOfcbar=0.1;
  cmd{ii}.estimate_noise_var_in_homogeneous_problem=false;
  cmd{ii}.pixelnoisevar_initialcondition='from_image_statistics';
  cmd{ii}.nu_ic_always_proportional2cbar=[];
  cmd{ii}.V_TolX=NaN;
  cmd{ii}.V_MaxIter=NaN;
  cmd{ii}.fn_savehistory=[outputbasename '.history.mat'];
  cmd{ii}.verbosity=1;
  cmd{ii}.MultiplierForPixelnoisevarIC=MultiplierForPixelnoisevarIC_eachstep(step);
  cmd{ii}.MinimumClassProb=1.0e-3/Neta;
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='box_writepixelnoisevar';
  cmd{ii}.fn_pixelnoisevar=[outputbasename '.pixelnoisevar.txt'];
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='vobj_print_virusobj';
  %%%%%%%%%%%%%%%%%%%%
%  ii=ii+1;
%  cmd{ii}.operator='vobj_save_virusobj';
%  cmd{ii}.fn_virusobj=[outputbasename '.vobj.mat'];
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='vobj_write_virusobj';
  cmd{ii}.fn_clnp=cell(Neta,1);
  cmd{ii}.fn_clnp{1}=[outputbasename '.eta1.clnp.txt'];
  cmd{ii}.fn_clnp{2}=[outputbasename '.eta2.clnp.txt'];
  cmd{ii}.fn_nu=cell(Neta,1);
  cmd{ii}.fn_nu{1}=[outputbasename '.eta1.nu.txt'];
  cmd{ii}.fn_nu{2}=[outputbasename '.eta2.nu.txt'];
  cmd{ii}.fn_q=cell(Neta,1);
  cmd{ii}.fn_q{1}=[outputbasename '.eta1.q.txt'];
  cmd{ii}.fn_q{2}=[outputbasename '.eta2.q.txt'];
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='misc_diary';
  cmd{ii}.fn_diary='off';
  %%%%%%%%%%%%%%%%%%%%
  ii=ii+1;
  cmd{ii}.operator='misc_changedirectory';
  cmd{ii}.dn='..';
  %%%%%%%%%%%%%%%%%%%%
end
%*******************End computing homogeneous reconstructions at increasing resolutions.
%%%%%%%%%%%%%%%%%%%%
% ii=ii+1;
% cmd{ii}.operator='misc_save_workspace';
% cmd{ii}.fn_workspace=[outputbasename '.workspace.mat'];
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
%%%%%%%%%%%%%%%%%%%%
