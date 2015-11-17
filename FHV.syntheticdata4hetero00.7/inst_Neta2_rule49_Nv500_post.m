function inst_Neta2_rule49_Nv500_post_inv()
%function inst_Neta2_rule49_Nv500_post_inv()

%%%%%%%%%%%%%%%%%%%%
inputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.hetero.inv';
outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.post';
fprintf(1,'inst_Neta2_rule49_Nv500_post_inv: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%
deltachi=[4.7 4.7]; %image sampling intervals in Angstroms
Neta=4;
%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.  To do nothing, specify cmd=[];.  To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
cmd=cell(20,1); %Preallocate for 20 operators.
ii=0; %Index for loading the cmd cell array.
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}=['hetero.step7/' inputbasename '.eta1.clnp.txt'];
cmd{ii}.fn_clnp{2}='/home/mitchvogel/Documents/research/imaging/Matlab/FHV.syntheticdata4hetero00.7/FHV.lmax10pmax5.clnp.txt';
cmd{ii}.fn_clnp{3}=['hetero.step7/' inputbasename '.eta2.clnp.txt'];
cmd{ii}.fn_clnp{4}='/home/mitchvogel/Documents/research/imaging/Matlab/FHV.syntheticdata4hetero00.7/FHV.lmax10pmax5.clnp.perturbed.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}=['hetero.step7/' inputbasename '.eta1.nu.txt'];
cmd{ii}.fn_nu{2}='/home/mitchvogel/Documents/research/imaging/FHV.syntheticdata4hetero00.7//FHV.lmax10pmax5.nu.txt';
cmd{ii}.fn_nu{3}=['hetero.step7/' inputbasename '.eta2.nu.txt'];
cmd{ii}.fn_nu{4}='/home/mitchvogel/Documents/research/imaging/Matlab/FHV.syntheticdata4hetero00.7/FHV.lmax10pmax5.nu.perturbed.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='FHV.Neta4.q.equalclassprobs.txt';
cmd{ii}.fn_q{2}='FHV.Neta4.q.equalclassprobs.txt';
cmd{ii}.fn_q{3}='FHV.Neta4.q.equalclassprobs.txt';
cmd{ii}.fn_q{4}='FHV.Neta4.q.equalclassprobs.txt';
%%%%%%%%%%%%%%%%%%%%
%start of the block that is repeated twice, once for Class 1 and once for Class 2
ii=ii+1;
iisave=ii;
cmd{ii}.operator='post_compute_FSC';
cmd{ii}.FSC_minmagk=0.0;
cmd{ii}.FSC_maxmagk=0.5/deltachi(1);
cmd{ii}.FSC_deltamagk=(cmd{ii}.FSC_maxmagk - cmd{ii}.FSC_minmagk)/222;
cmd{ii}.FSC_eta4classA=1; %CHANGE IS HERE
cmd{ii}.FSC_eta4classB=2; %CHANGE IS HERE
cmd{ii}.FSC_is_same_vobj=true;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_plot';
cmd{ii}.plotcutoffk=false;
cmd{ii}.plottitle='Class 1: truth versus estimate'; %CHANGE IS HERE
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_LPF_FSC';
cmd{ii}.LPFcutoff4FSC=0.2;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_cutoff';
cmd{ii}.FSCcutoff=0.5;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_plot';
cmd{ii}.plotcutoffk=true;
cmd{ii}.plottitle='Class 1: truth versus estimate'; %CHANGE IS HERE
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_write_FSC';
cmd{ii}.fn_FSC=[outputbasename '.eta4classA.' num2str(cmd{iisave}.FSC_eta4classA) '.eta4classB.' num2str(cmd{iisave}.FSC_eta4classB) '.txt'];
%end of the block that is repeated twice, once for Class 1 and once for Class 2
%%%%%%%%%%%%%%%%%%%%
%start of the block that is repeated twice, once for Class 1 and once for Class 2
ii=ii+1;
iisave=ii;
cmd{ii}.operator='post_compute_FSC';
cmd{ii}.FSC_minmagk=0.0;
cmd{ii}.FSC_maxmagk=0.5/deltachi(1);
cmd{ii}.FSC_deltamagk=(cmd{ii}.FSC_maxmagk - cmd{ii}.FSC_minmagk)/222;
cmd{ii}.FSC_eta4classA=3; %CHANGE IS HERE
cmd{ii}.FSC_eta4classB=4; %CHANGE IS HERE
cmd{ii}.FSC_is_same_vobj=true;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_plot';
cmd{ii}.plotcutoffk=false;
cmd{ii}.plottitle='Class 2: truth versus estimate'; %CHANGE IS HERE
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_LPF_FSC';
cmd{ii}.LPFcutoff4FSC=0.2;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_cutoff';
cmd{ii}.FSCcutoff=0.5;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_FSC_plot';
cmd{ii}.plotcutoffk=true;
cmd{ii}.plottitle='Class 2: truth versus estimate'; %CHANGE IS HERE
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='post_write_FSC';
cmd{ii}.fn_FSC=[outputbasename '.eta4classA.' num2str(cmd{iisave}.FSC_eta4classA) '.eta4classB.' num2str(cmd{iisave}.FSC_eta4classB) '.txt'];
%end of the block that is repeated twice, once for Class 1 and once for Class 2
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
%%%%%%%%%%%%%%%%%%%%
