function inst_Neta2_AsymObj_rule848_Nv500_fw()
%function inst_Neta2_AsymObj_rule848_Nv500_fw()

outputbasename='AsymObj.out.lmax2pmax5.Neta2.rule8-4-8.Nv500.fw';
fprintf(1,'inst_Neta2_AsymObj_rule848_Nv500_fw: outputbasename %s\n',outputbasename);
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
%ii=ii+1;
%cmd{ii}.operator='quad_read_integration_rule';
%cmd{ii}.fn_rule='rule_small_3_rulechopper';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='quad_compute_integration_rule';
cmd{ii}.whichrule=1; %all Euler angles (except handedness uncertainty): alpha [0,2pi) uniform, beta [0,pi/2] Gauss-Legendre, gamma [0,2pi) uniform
cmd{ii}.Nabc=[8 4 8]'; %order, makes sense for alpha and gamma to have 4 times more points because the region is 4 times larger
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
Neta=2;
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}='AsymObj.ic.lmax2pmax5.clnp.eta1.txt';
cmd{ii}.fn_clnp{2}='AsymObj.ic.lmax2pmax5.clnp.eta2.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}='AsymObj.ic.lmax2pmax5.nu.eta1.txt';
cmd{ii}.fn_nu{2}='AsymObj.ic.lmax2pmax5.nu.eta2.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='AsymObj.ic.lmax2pmax5.q.eta1.txt';
cmd{ii}.fn_q{2}='AsymObj.ic.lmax2pmax5.q.eta2.txt';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
cmd{ii}.use_vkminimalset_rather_than_vk=true;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.fn_tilde_b=cell(2,1); %have 2 types of angular basis functions
cmd{ii}.fn_tilde_b{1}='_set_tilde_b_to_empty_';
cmd{ii}.fn_tilde_b{2}='callti.out.read_by_C';
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
ii=ii+1;
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
