function inst_Neta2_rule49_Nv500_pre4hetero_inv()
%function inst_Neta2_rule49_Nv500_pre4hetero_inv()

%Writes reciprocal space image stack in mrc format.

%%%%%%%%%%%%%%%%%%%%
outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.pre4hetero.inv';
fprintf(1,'inst_Neta2_rule49_Nv500_pre4hetero_inv: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%
deltachi=[4.7 4.7]; %image sampling intervals in Angstroms
R2=197.4; %outer radius R2
%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.  To do nothing, specify cmd=[];.  To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
cmd=cell(20,1); %Preallocate for 100 operators.
ii=0; %Index for loading the cmd cell array.
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary.txt'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_setpseudorandomnumberseed';
cmd{ii}.pseudorandomnumberseed=29831; %arbitrary seed
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_readimagestack';
cmd{ii}.imagestackformat='mrc';
cmd{ii}.fn_imagestack='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.fw.imagestack.mrc';
cmd{ii}.startSlice=1;
cmd{ii}.numSlices=500;
%%%%%%%%%%%%%%%%%%%%
%NwV-JSB2013 ii=ii+1;
%NwV-JSB2013 cmd{ii}.operator='box_classifyviasamplemean';
%NwV-JSB2013 cmd{ii}.classifythres=0.16; %WangMatsuiDomitrovicZhengDoerschukJohnson JSB 2013, p. 197, column 2 line -3
%NwV-JSB2013 %%%%%%%%%%%%%%%%%%%%
%NwV-JSB2013 ii=ii+1;
%NwV-JSB2013 cmd{ii}.operator='box_permute';
%NwV-JSB2013 %%%%%%%%%%%%%%%%%%%%
%NwV-JSB2013 ii=ii+1;
%NwV-JSB2013 cmd{ii}.operator='box_extractsubset';
%NwV-JSB2013 cmd{ii}.maxNv=1200; %WangMatsuiDomitrovicZhengDoerschukJohnson JSB 2013, p. 198, column 1 line 6
%NwV-JSB2013 cmd{ii}.a=1; %WangMatsuiDomitrovicZhengDoerschukJohnson JSB 2013, p. 198, column 1 line 7
%NwV-JSB2013 cmd{ii}.b=4; %WangMatsuiDomitrovicZhengDoerschukJohnson JSB 2013, p. 198, column 1 line 7
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='basic_set2Drealspacesamplingintervals';
cmd{ii}.samplingintervals=deltachi;
%%%%%%%%%%%%%%%%%%%%
%NwV-JSB2013 %lines 134-171 of /home/qw32/hetero3d/newYiliCode/newcode/cacRuns/NwV_cap_preprocess/cacfw_nwv_dV_ico.m concern the subtraction of the mean and scaling by the standard deviation.  The statistics are computed in the region > 250 Angstrom from the center of the image.  Note that Qiu Wang gives radii in terms of pixels not Angstroms.
%NwV-JSB2013 ii=ii+1;
%NwV-JSB2013 cmd{ii}.operator='box_normalize2zeroone';
%NwV-JSB2013 cmd{ii}.radius01=[250 1000]; %1000 Angstrom is outside of the image, even in the corners
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_annulusstatistics';
cmd{ii}.radius01=[R2+deltachi(1) R2+10*deltachi(1)]; %Uncertain of the correspondence with Qiu Wang's code
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='realrecip_2DFFT';
%%%%%%%%%%%%%%%%%%%%
%ii=ii+1;
%cmd{ii}.operator='misc_save_workspace';
%cmd{ii}.fn_workspace=[outputbasename '.workspace.mat'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_writeannulusstatistics';
cmd{ii}.fn_annulusstatistics=[outputbasename '.annulusstatistics.txt'];
%%%%%%%%%%%%%%%%%%%%image stack has not been modified
%ii=ii+1;
%cmd{ii}.operator='misc_write_mrc';
%cmd{ii}.fn_write_mrc=[outputbasename '.Imagestack.mrc'];
%cmd{ii}.what2write='write_image_stack';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_write_mrc';
cmd{ii}.fn_write_mrc=[outputbasename '.Imagestack.mrc'];
cmd{ii}.what2write='write_Image_stack';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
%%%%%%%%%%%%%%%%%%%%
