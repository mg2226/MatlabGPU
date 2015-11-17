function test_FFT_inst()
%function test_FFT_inst()

outputbasename='test_FFT_inst';
fprintf(1,'test_FFT_inst: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.
%To do nothing, specify cmd=[];
%To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd=cell(20,1); %Preallocate for 20 operators.
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
cmd{ii}.operator='test_FFT_pt1';
cmd{ii}.N1=20;
cmd{ii}.N2=21;
cmd{ii}.delta1=2.74;
cmd{ii}.delta2=3.24;
cmd{ii}.freqintegers=[0 0 ; 1 0 ; 0 1 ; -1 0 ; 0 -1 ; 1 1 ; -1 -1 ; 1 -1 ; -1 1];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_showimagestack';
cmd{ii}.whichimages=[1:9];
cmd{ii}.printtitle=true; %true or false for the entire set of images
cmd{ii}.fn_class='';
%%%%%%%%%%%%%%%%%%%%NOT IN pre4hetero_inv
ii=ii+1;
cmd{ii}.operator='basic_setsizeof2Drealspaceimagesfromimagestack';
%%%%%%%%%%%%%%%%%%%%NOT IN pre4hetero_inv
ii=ii+1;
cmd{ii}.operator='basic_compute2Dreciprocalspaceproperties';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='realrecip_2DFFT';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='EM_extractdatasubset';
cmd{ii}.kmax=-1;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='test_FFT_pt2';
%requires nothing but uses cmd_test_FFT_pt1 set by test_FFT_pt1
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
