function inst_Neta2_rule49_Nv500_fw()
%function inst_Neta2_rule49_Nv500_fw()

outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.fw';
fprintf(1,'inst_Neta2_rule49_Nv500_fw: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.
%To do nothing, specify cmd=[];
%To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd=cell(100,1); %Preallocate for 100 operators.
ii=0;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_setpseudorandomnumberseed';
cmd{ii}.pseudorandomnumberseed=383511;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_set2Drealspacesamplingintervals';
cmd{ii}.samplingintervals=[4.7 4.7];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_setsizeof2Drealspaceimages';
cmd{ii}.NaNb=[91 91];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_compute2Dreciprocalspaceproperties';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='quad_read_integration_rule';
cmd{ii}.fn_rule='rule_small_3_rulechopper';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
Neta=2;
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}='FHV.lmax10pmax5.clnp.txt';
cmd{ii}.fn_clnp{2}='FHV.lmax10pmax5.clnp.perturbed.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}='FHV.lmax10pmax5.nu.txt';
cmd{ii}.fn_nu{2}='FHV.lmax10pmax5.nu.perturbed.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='FHV.lmax10pmax5.q.txt';
cmd{ii}.fn_q{2}='FHV.lmax10pmax5.q.perturbed.txt';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
cmd{ii}.use_vkminimalset_rather_than_vk=true;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.fn_tilde_b='callti.out.read_by_C';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='fw_mk_synthetic_2D_realspace_images';
cmd{ii}.Nv=500;
cmd{ii}.NT=1;
cmd{ii}.SNR=5.0;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='fw_write_truevalues';
cmd{ii}.fn_truevalues=[outputbasename '.truevalues.txt'];
%%%%%%%%%%%%%%%%%%%%
%ii=ii+1;
cmd{ii}.operator='misc_write_mrc';
cmd{ii}.fn_write_mrc=[outputbasename '.imagestack.mrc'];
cmd{ii}.what2write='write_image_stack';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
