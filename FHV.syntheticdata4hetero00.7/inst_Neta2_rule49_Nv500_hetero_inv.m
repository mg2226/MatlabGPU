function inst_Neta2_rule49_Nv500_hetero_inv()
%function inst_Neta2_rule49_Nv500_hetero_inv()

%Loads reciprocal space image stack.
%Because there is only one kmax in this script, it is possible to prepare y in an earlier script and not have to store both Imagestack and y.

%%%%%%%%%%%%%%%%%%%%
inputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.homo.inv';
outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.hetero.inv';
fprintf(1,'inst_Neta2_rule49_Nv500_hetero_inv: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%
NaNb=[91 91]; %image dimensions in pixels
deltachi=[4.7 4.7]; %image sampling intervals in Angstroms
Neta=2;
%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.  To do nothing, specify cmd=[];.  To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
cmd=cell(50,1); %Preallocate for 50 operators.
ii=0; %Index for loading the cmd cell array.
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_readpixelnoisevar';
cmd{ii}.fn_pixelnoisevar=['step3/' inputbasename '.pixelnoisevar.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_readImagestack';
cmd{ii}.Imagestackformat='mrc';
cmd{ii}.fn_Imagestack='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.pre4hetero.inv.Imagestack.mrc';
cmd{ii}.startSlice=1;
cmd{ii}.numSlices=500; %this is number of complex reciprocal-space images
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_setsizeof2DreciprocalspaceImagesfromImagestack';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_set2Drealspacesamplingintervals';
cmd{ii}.samplingintervals=deltachi;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_compute2Dreciprocalspaceproperties';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}=['step3/' inputbasename '.eta1.clnp.txt'];
cmd{ii}.fn_clnp{2}=['step3/' inputbasename '.eta2.clnp.txt'];
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}=['step3/' inputbasename '.eta1.nu.txt'];
cmd{ii}.fn_nu{2}=['step3/' inputbasename '.eta2.nu.txt'];
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}=['step3/' inputbasename '.eta1.q.txt'];
cmd{ii}.fn_q{2}=['step3/' inputbasename '.eta2.q.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_change_homo2hetero_in_virusobj';
cmd{ii}.homo2hetero=cell(Neta,1);
cmd{ii}.homo2hetero{1}.action=1; %make this class heterogeneous
cmd{ii}.homo2hetero{1}.FractionOfMeanForMinimum=0.0025; %only required when .action=1
cmd{ii}.homo2hetero{1}.FractionOfMean=0.02; %only required when .action=1
cmd{ii}.homo2hetero{2}.action=1; %make this class heterogeneous
cmd{ii}.homo2hetero{2}.FractionOfMeanForMinimum=0.0025; %only required when .action=1
cmd{ii}.homo2hetero{2}.FractionOfMean=0.02; %only required when .action=1
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_print_virusobj';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_changedirectory';
cmd{ii}.dn='hetero.step7';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.fn_tilde_b='callti.out.read_by_C';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_extractdatasubset';
kbiggest=sqrt( set_kbiggest(NaNb(1),deltachi(1)).^2 + set_kbiggest(NaNb(2),deltachi(2)).^2 );
cmd{ii}.kmax=kbiggest*0.5/sqrt(2);
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
cmd{ii}.use_vkminimalset_rather_than_vk=false;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='quad_read_integration_rule';
cmd{ii}.fn_rule='rule_small_3';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_expectationmaximization';
cmd{ii}.MC_Nrandic=1;
cmd{ii}.MC_FractionOfMeanForMinimum=0.005; %unused because MC_Nrandic=1 but still must be defined.
cmd{ii}.MC_FractionOfMean=0.2; %unused because MC_Nrandic=1 but still must be defined.
cmd{ii}.maxiter=200;
cmd{ii}.maxiter4pixelnoisevarupdate=5;
cmd{ii}.cbarftol=-1.0; %unused.
cmd{ii}.cbarrtol=1.0e-4;
cmd{ii}.cbarftol_dividebyNc=false;
cmd{ii}.loglikeftol=-1.0; %unused.
cmd{ii}.loglikertol=-1.0; %unused.
cmd{ii}.nu_ic_FractionOfcbarForMinimum=0.15;
cmd{ii}.nu_ic_FractionOfcbar=0.1;
cmd{ii}.estimate_noise_var_in_homogeneous_problem=false;
cmd{ii}.pixelnoisevar_initialcondition='from_pixelnoisevar';
cmd{ii}.nu_ic_always_proportional2cbar=[];
cmd{ii}.V_TolX=1e-10;
cmd{ii}.V_MaxIter=8;
cmd{ii}.fn_savehistory=[outputbasename '.history.mat'];
cmd{ii}.verbosity=1;
cmd{ii}.MultiplierForPixelnoisevarIC=1.0;
cmd{ii}.MinimumClassProb=10e-3/Neta;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_writepixelnoisevar';
cmd{ii}.fn_pixelnoisevar=[outputbasename '.pixelnoisevar.txt'];
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
%ii=ii+1;
%cmd{ii}.operator='misc_save_workspace';
%cmd{ii}.fn_workspace=[outputbasename '.workspace.mat'];
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
%%%%%%%%%%%%%%%%%%%%
