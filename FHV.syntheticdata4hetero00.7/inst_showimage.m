function inst_showimage()
%function inst_showimage()

outputbasename='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.showimage';
fprintf(1,'inst_showimage: outputbasename %s\n',outputbasename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Operators are executed in the order in which they appear in cmd, using the arguments that also appear in cmd.
%To do nothing, specify cmd=[];
%To indicate a 'no-op', set that element of the cell array to an empty matrix, e.g., cmd{1}=[], which is the initialization set by the 'cell' function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cmd=cell(10,1); %Preallocate for 10 operators.
ii=0;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[outputbasename '.diary'];
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_readimagestack';
cmd{ii}.imagestackformat='mrc';
cmd{ii}.fn_imagestack='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.fw.imagestack.mrc';
cmd{ii}.startSlice=1;
cmd{ii}.numSlices=500;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_showimagestack';
cmd{ii}.whichimages=[1:2:10];
cmd{ii}.printtitle=true;
cmd{ii}.fn_class='FHV.out.lmax10pmax5.Neta2.rule49.Nv500.fw.truevalues.txt';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
