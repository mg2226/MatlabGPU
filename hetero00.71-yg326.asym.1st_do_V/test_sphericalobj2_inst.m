function test_sphericalobj2_inst()
%function test_sphericalobj2_inst()

outputbasename='test_sphericalobj2_inst';
fprintf(1,'test_sphericalobj2_inst: outputbasename %s\n',outputbasename);
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
cmd{ii}.operator='test_imagestackfromsphericalobj';
cmd{ii}.numberimages=100;
cmd{ii}.minxindex=-100;
cmd{ii}.maxxindex=100;
cmd{ii}.minyindex=-100;
cmd{ii}.maxyindex=100;
cmd{ii}.deltaxy=2.74;
cmd{ii}.rho0=1.0;
cmd{ii}.R2=200.0;
cmd{ii}.pixelnoisevariance=(50.0)^2;
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='box_showimagestack';
cmd{ii}.whichimages=[1 2 3];
cmd{ii}.printtitle=true; %true or false for the entire set of images
cmd{ii}.fn_class='';
%%%%%%%%%%%%%%%%%%%%
ii=ii+1;
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary='off';
%%%%%%%%%%%%%%%%%%%%
%Execute the operations in the cmd cell array.
hetero(cmd);
return;
