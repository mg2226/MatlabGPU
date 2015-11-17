ii=1; %set this variable so that the following statements are valid Matlab statements.
basename='base_file_name'; %set this variable so that the following statements are valid Matlab statements.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%misc operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Start the diary.
cmd{ii}.operator='misc_diary';
cmd{ii}.fn_diary=[basename '.diary'];

%Add a directory to the Matlab path.
cmd{ii}.operator='misc_addtomatlabpath';
cmd{ii}.dn='/home/pd83/Hetero';

%Set the state of the random number generator.
cmd{ii}.operator='misc_setpseudorandomnumberseed';
cmd{ii}.pseudorandomnumberseed=89732634267328;

%Save the state of the random number generator to a Matlab .mat file.
cmd{ii}.operator='misc_savepseudorandomnumberstate2file';
cmd{ii}.fn_pseudorandomnumberstate=[basename '.pseudoRVstate.mat'];

%Load the state of the random number generator from a Matlab .mat file and restore the state.
cmd{ii}.operator='misc_restorepseudorandomnumberstatefromfile';
cmd{ii}.fn_pseudorandomnumberstate=[basename '.pseudoRVstate.mat'];

%Clear a variable from the hetero workspace.  The motivation is probably to reduce storage requirements.
cmd{ii}.operator='misc_clearvariable';
cmd{ii}.variablename='imagestack';
cmd{ii}.variablename='imageindex';
cmd{ii}.variablename='Imagestack';

%Change directory.
cmd{ii}.operator='misc_changedirectory';
cmd{ii}.dn='step1';
cmd{ii}.dn='..';

%Return, ignoring any additional elements in cmd cell array.
cmd{ii}.operator='misc_return';

%Enter keyboard mode, can exit by entering the command 'return'.
cmd{ii}.operator='misc_keyboard';

%Save the workspace as a Matlab .mat file.
cmd{ii}.operator='misc_save_workspace';
cmd{ii}.fn_workspace=[basename '.workspace.mat'];

%Load the workspace from a Matlab .mat file.
cmd{ii}.operator='misc_load_workspace';
cmd{ii}.fn_workspace=[basename '.workspace.mat'];
cmd{ii}.ExistingOverReload=true;

%Write a file in MRC format.
cmd{ii}.operator='misc_write_mrc';
cmd{ii}.fn_write_mrc=[basename '.extension.mrc'];
cmd{ii}.what2write='write_image_stack';
cmd{ii}.what2write='write_Image_stack';
cmd{ii}.what2write='write_rhobar';
cmd{ii}.what2write='write_rxx';
cmd{ii}.what2write='write_scaled_rhobar';
cmd{ii}.what2write='write_scaled_rxx';

%Push a data structure.  The stack is only one deep for a particular type of data structure.
cmd{ii}.operator='misc_push';
cmd{ii}.what2push='push_image_stack'
cmd{ii}.what2push='push_virusobj'

%Pop a data structure.  The stack is only one deep for a particular type of data structure.
cmd{ii}.operator='misc_pop';
cmd{ii}.what2pop='pop_image_stack'
cmd{ii}.what2pop='pop_virusobj'

%Reduce memory fragmentation by executing matlab 'pack'.
cmd{ii}.operator='misc_pack';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basic operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set real-space image sampling intervals in the first and second coordinates.
%The units that are used here determine the units that are used throughout the software.
%Need to check if spatial frequencies are 1/r or 2\pi/r.
cmd{ii}.operator='basic_set2Drealspacesamplingintervals';
cmd{ii}.samplingintervals=[2.57 2.57];

%Set the size of a 2-D real-space image.  Sometimes this is set by reading an image.  But
%certainly in the case of computing synthetic images, this needs to be set.
cmd{ii}.operator='basic_setsizeof2Drealspaceimages';
cmd{ii}.NaNb=[101 101];

%Set the size of a 2-D real-space image by examining an image.  Requires that all images in
%the stack have the same size.
cmd{ii}.operator='basic_setsizeof2Drealspaceimagesfromimagestack';

%Set the size of a 2-D reciprocal-space image by examining an Image.  Requires that all Images in
%the stack have the same size.
cmd{ii}.operator='basic_setsizeof2DreciprocalspaceImagesfromImagestack';

%Determine a minimal subset of 2-D reciprocal space such that conjugate symmetry fills in all of reciprocal space.
cmd{ii}.operator='basic_compute2Dreciprocalspaceproperties';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real-space boxed image operators (lowercase b in box):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read the pixel noise variance from a text file.
cmd{ii}.operator='box_readpixelnoisevar';
cmd{ii}.fn_pixelnoisevar=[basename '.pixelnoisevar.txt'];

%Write the pixel noise variance to a text file.
cmd{ii}.operator='box_writepixelnoisevar';
cmd{ii}.fn_pixelnoisevar=[basename '.pixelnoisevar.txt'];

%Print the imagestack and imageindex data on the display.
cmd{ii}.operator='box_printimagesasarrays';

%Read the imagestack and imageindex data from files.
cmd{ii}.operator='box_readimagestack';
cmd{ii}.imagestackformat='fake';
cmd{ii}.imagestackformat='mrc';
cmd{ii}.imagestackformat='img';
%If mrc, give full file name.
%If img, apparently can give full file name or only base name.
%If fake, currently unused.
cmd{ii}.fn_imagestack='imagestackfilename';
%If using mrc, also give
cmd{ii}.startSlice=10;
cmd{ii}.numSlices=20;
%If using img, also give
cmd{ii}.startimage=10;
cmd{ii}.numimages=20;

%Read the Imagestack data from files.
cmd{ii}.operator='box_readImagestack';
cmd{ii}.Imagestackformat='fake';
cmd{ii}.Imagestackformat='mrc';
cmd{ii}.Imagestackformat='img';
%If mrc, give full file name.
%If img, apparently can give full file name or only base name.
%If fake, currently unused.
cmd{ii}.fn_Imagestack='Imagestackfilename';
%If using mrc, also give
cmd{ii}.startSlice=10;
cmd{ii}.numSlices=20; %this is number of complex reciprocal-space images
%If using img, also give
cmd{ii}.startImage=10;
cmd{ii}.numImages=20;

%Load the imagestack and imageindex data from a Matlab .mat file.
cmd{ii}.operator='box_loadimagestack';
cmd{ii}.fn_imagestack='imagestackfilename.mat';
cmd{ii}.only_imagestack=true; %set to true if there is no imageindex data

%Load the Imagestack and imageindex data from a Matlab .mat file.
cmd{ii}.operator='box_loadImagestack';
cmd{ii}.fn_Imagestack='Imagestackfilename.mat';

%Save the imagestack and imageindex data to a Matlab .mat file.
cmd{ii}.operator='box_saveimagestack';
cmd{ii}.fn_imagestack=[basename '.imagestack.mat'];

%Save the Imagestack and imageindex data to a Matlab .mat file.
cmd{ii}.operator='box_saveImagestack';
cmd{ii}.fn_Imagestack=[basename '.Imagestack.mat'];

%Extract the images numbered a + b*(n-1) (n=1, 2, ...) from a larger image stack up to a total of maxNv images.
cmd{ii}.operator='box_extractsubset';
cmd{ii}.maxNv=6000;
cmd{ii}.a=1;
cmd{ii}.b=4;

%Randomly permute the images in a stack.
cmd{ii}.operator='box_permute';

%Shrink each image in a stack by removing a certain number of pixels from each edge.
cmd{ii}.operator='box_shrink';
cmd{ii}.pixels2delete=[1 2 3 4]; %the order is "left", "bottom", "right", "top" when you think of the image as a matrix.

%Mask the corners (replace by the sample mean of the masked region) of each image.
cmd{ii}.operator='box_maskcorners';
cmd{ii}.radiuscorner=258.0;

%Normalize each image so that the sample mean is 0 and the sample variance is one
%in an annulus.
cmd{ii}.operator='box_normalize2zeroone';
cmd{ii}.radius01=[230.5 266.3];

%Classify images as good versus bad based on their l2 norm from the mean image.
%Extract the images that are classified as good from the stack.
cmd{ii}.operator='box_classifyviasamplemean';
cmd{ii}.classifythres=0.16;

%Compute first sample moment and second central sample moment of an annulus averaging
%over the entire real-space image stack.
cmd{ii}.operator='box_annulusstatistics';
cmd{ii}.radius01=[230.5 266.3];

%Write the annulus statistics.
cmd{ii}.operator='box_writeannulusstatistics';
cmd{ii}.fn_annulusstatistics=[basename '.annulusstatistics.txt'];

%Read the annulus statistics.
cmd{ii}.operator='box_readannulusstatistics';
cmd{ii}.fn_annulusstatistics=[basename '.annulusstatistics.txt'];

%Show selected images.
cmd{ii}.operator='box_showimagestack';
cmd{ii}.whichimages=[1:2:10];
cmd{ii}.printtitle=true; %true or false for the entire set of images
cmd{ii}.fn_class=[basename '.fw.truevalues.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real space <-> reciprocal space transformations operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transform 2-D real-space images to reciprocal space.
cmd{ii}.operator='realrecip_2DFFT';

%Transform 2-D reciprocal-space images to real space.
cmd{ii}.operator='realrecip_2DIFFT';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reciprocal-space boxed image operators (uppercase B in Box):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print the Imagestack and imageindex data on the display.
cmd{ii}.operator='Box_printImagesasarrays';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%virus object operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Print a virusobj on the display.
cmd{ii}.operator='vobj_print_virusobj';

%Read a virusobj from a set of text files -- do not read or compute the radial basis function table part of the virusobj.
cmd{ii}.operator='vobj_read_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}='V1.clnp.txt';
cmd{ii}.fn_clnp{2}='V2.clnp.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}='V1.nu.txt';
cmd{ii}.fn_nu{2}='V2.nu.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='V1.q.txt';
cmd{ii}.fn_q{2}='V2.q.txt';

%Write a virusobj to a set of text files -- do not write the radial basis function table part of the virusobj.
cmd{ii}.operator='vobj_write_virusobj';
cmd{ii}.fn_clnp=cell(Neta,1);
cmd{ii}.fn_clnp{1}='V1.clnp.txt';
cmd{ii}.fn_clnp{2}='V2.clnp.txt';
cmd{ii}.fn_nu=cell(Neta,1);
cmd{ii}.fn_nu{1}='V1.nu.txt';
cmd{ii}.fn_nu{2}='V2.nu.txt';
cmd{ii}.fn_q=cell(Neta,1);
cmd{ii}.fn_q{1}='V1.q.txt';
cmd{ii}.fn_q{2}='V2.q.txt';

%Save a virusobj to a Matlab .mat file.
cmd{ii}.operator='vobj_save_virusobj';
cmd{ii}.fn_virusobj=[basename '.vobj.mat'];

%Load a virusobj from a Matlab .mat file.
cmd{ii}.operator='vobj_load_virusobj';
cmd{ii}.fn_virusobj=[basename '.vobj.mat'];

%Change the dimension of the virusobj.
cmd{ii}.operator='vobj_change_size_of_virusobj';
cmd{ii}.vlmax=[0 0]; %vector of new lmax values, one for each class
cmd{ii}.vpmax=[4 4]; %vector of new pmax values, one for each class

%Change a homogeneous class to a heterogeneous class or visa versa.
cmd{ii}.operator='vobj_change_homo2hetero_in_virusobj';
cmd{ii}.homo2hetero=cell(Neta,1);
cmd{ii}.homo2hetero{1}.action=1; %make this class heterogeneous
cmd{ii}.homo2hetero{1}.action=0; %no change to this class
cmd{ii}.homo2hetero{1}.action=-1; %make this class homogeneous
cmd{ii}.homo2hetero{1}.FractionOfMeanForMinimum=0.005; %only required when .action=1
cmd{ii}.homo2hetero{1}.FractionOfMean=0.2; %only required when .action=1
cmd{ii}.homo2hetero{2}.action=1; %make this class heterogeneous
cmd{ii}.homo2hetero{2}.action=0; %no change to this class
cmd{ii}.homo2hetero{2}.action=-1; %make this class homogeneous
cmd{ii}.homo2hetero{2}.FractionOfMeanForMinimum=0.005; %only required when .action=1
cmd{ii}.homo2hetero{2}.FractionOfMean=0.2; %only required when .action=1

%Change the handedness of the virus.
cmd{ii}.operator='vobj_change_handedness';
cmd{ii}.changehand=[true false]; %Neta entries

%Change the sign of all c_{l,n,p} coefficients.  As of 2013, UCSF
%Chimera seems to expect that the object be positive going.
cmd{ii}.operator='vobj_change_sign';
cmd{ii}.changesign=[true false]; %Neta entries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integration rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read the integration rule from a file.
cmd{ii}.operator='quad_read_integration_rule';
cmd{ii}.fn_rule='rule_small_3';

%Compute the integration rule.
cmd{ii}.operator='quad_compute_integration_rule';
cmd{ii}.whichrule=1; %all Euler angles (except handedness uncertainty): alpha [0,2pi) uniform, beta [0,pi/2] Gauss-Legendre, gamma [0,2pi) uniform
cmd{ii}.Nabc=[32 8 32]'; %order, makes sense for alpha and gamma to have 4 times more points because the region is 4 times larger

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%forward operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute synthetic images.
cmd{ii}.operator='fw_mk_synthetic_2D_realspace_images';
cmd{ii}.Nv=500;
cmd{ii}.NT=1;
cmd{ii}.SNR=0.5;

%Compute synthetic real-space images from a 3-D model with multiple 3-D Gaussian pulses
cmd{ii}.operator='fw_mk_synthetic_GaussianPulse_2D_realspace_images';
cmd{ii}.Nv=500;
cmd{ii}.NT=1;
cmd{ii}.SNR=0.5;
%cmd{ii}.GaussianPulses{1:Neta}.q=class probability
%cmd{ii}.GaussianPulses{1:Neta}.gainminmax(1:2)=range of uniform pdf on gain for this class which is the sole source of continuous heterogeneity
%cmd{ii}.GaussianPulses{1:Neta}.means(1:Npulse,1:3)=mean vector for each pulse in the class
%cmd{ii}.GaussianPulses{1:Neta}.covars(1:Npulse,1:3,1:3)=covariance matrix for each pulse in the class
cmd{ii}.GaussianPulses=cell(1,1);
cmd{ii}.GaussianPulses{1}.q=1.0;
cmd{ii}.GaussianPulses{1}.gainminmax=[.9 ; 1.1];
cmd{ii}.GaussianPulses{1}.means=zeros(1,3);
cmd{ii}.GaussianPulses{1}.covars=zeros(1,3,3);
cmd{ii}.GaussianPulses{1}.covars=zeros(1,1,1)=1.0;
cmd{ii}.GaussianPulses{1}.covars=zeros(1,2,2)=1.5;
cmd{ii}.GaussianPulses{1}.covars=zeros(1,3,3)=2.0;

%Save the true values of the class and the projection direction and origin offset
%(described as an index into the integration rule) as a Matlab .mat file.
cmd{ii}.operator='fw_save_truevalues';
cmd{ii}.fn_truevalues=[basename '.fw_truevalues.mat'];

%Save the true values of the class and the projection direction and origin offset
%(described as an index into the integration rule) as a text file.
cmd{ii}.operator='fw_write_truevalues';
cmd{ii}.fn_truevalues=[basename '.fw_truevalues.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%expectation-maximization operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Read the coefficients $\tilde b$ for computing icosahedral harmonics from a file.
%give one file name for each type of angular basis function:
%0 = real-valued spherical harmonics
%1 = icosahedral harmonics (identity representation only)
%use '_set_tilde_b_to_empty_' to set tilde_b to the empty matrix
cmd{ii}.operator='EM_read_tilde_b';
cmd{ii}.fn_tilde_b=cell(2,1); %have 2 types of angular basis functions
cmd{ii}.fn_tilde_b{1}='_set_tilde_b_to_empty_';
cmd{ii}.fn_tilde_b{2}='callti.out.read_by_C';

%Set the range of 2-D reciprocal space that will be used and construct the y data structure.
cmd{ii}.operator='EM_extractdatasubset';
cmd{ii}.kmax=-1;

%Set the table of reciprocal-space radial basis function values in vobj.
cmd{ii}.operator='EM_set_2Dreciprocal_in_virusobj';
cmd{ii}.use_vkminimalset_rather_than_vk=false;

%Remove the table of reciprocal-space radial basis function values in vobj.
%May want to do this before saving a virusobj via EM_save_virusobj.
cmd{ii}.operator='EM_rm_2Dreciprocal_in_virusobj';

%Estimation: spherically-symmetric homogeneous model, done by least squares.
cmd{ii}.operator='EM_sphericalsymmetry_homogeneous';

%Estimation: homogeneous or heterogeneous cases, can do various symmetries or no symmetry, done by expectation maximization.
%MC_* for Monte Carlo.  Others for expectation maximization, including V_* for covariance matrix estimation.
cmd{ii}.operator='EM_expectationmaximization';
cmd{ii}.MC_Nrandic=100;
cmd{ii}.MC_FractionOfMeanForMinimum=0.005;
cmd{ii}.MC_FractionOfMean=0.2;
cmd{ii}.maxiter=200;
cmd{ii}.maxiter4pixelnoisevarupdate=5;
cmd{ii}.cbarftol=-1.0; %unused.
cmd{ii}.cbarrtol=1.0e-4;
cmd{ii}.cbarftol_dividebyNc=false;
cmd{ii}.cbarftol_dividebyNc=true;
%An extremely strict criteria is used in the Qiu Wang code with an || not an &&.
%In new code, use && not || and just turn off the test.
%cmd{ii}.loglikeftol=1.0e-10;
cmd{ii}.loglikeftol=-1.0; %unused.
cmd{ii}.loglikertol=-1.0; %unused.
cmd{ii}.nu_ic_FractionOfcbarForMinimum=0.005;
cmd{ii}.nu_ic_FractionOfcbar=0.1;
%cmd{ii}.estimate_noise_var_in_homogeneous_problem=true;
cmd{ii}.estimate_noise_var_in_homogeneous_problem=false;
cmd{ii}.pixelnoisevar_initialcondition='from_image_statistics';
cmd{ii}.pixelnoisevar_initialcondition='from_pixelnoisevar';
%nu initial condition is not reset to be proportional to the current cbar for any model.
cmd{ii}.nu_ic_always_proportional2cbar=[];
%nu initial condition is not reset to be proportional to the current cbar for any model for particular classes.
cmd{ii}.nu_ic_always_proportional2cbar=[true true]';
cmd{ii}.nu_ic_always_proportional2cbar=[false false]';
cmd{ii}.V_TolX=1e-10;
cmd{ii}.V_MaxIter=8;
cmd{ii}.fn_savehistory=[basename '.history.mat']; %use =[] to not save the history
cmd{ii}.verbosity=1; %0,1,2 have meaning
cmd{ii}.MultiplierForPixelnoisevarIC=1.0;
cmd{ii}.MinimumClassProb=10e-3/Neta;
cmd{ii}.absthres=10e-4;
cmd{ii}.homo2hetero_1st_V_MaxIter=100; %used for first nu=diag(V) optimization
cmd{ii}.homo2hetero_estpixelnoisevar=true; %used for first nu=diag(V) optimization
cmd{ii}.homo2hetero_estV=true; %used for first nu=diag(V) optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%post-processing operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Compute real space cubes, both \bar\rho and R_{\rho}(\vec x, \vec x).
cmd{ii}.operator='post_compute_real_space_cubes';
cmd{ii}.whichclass=2;
cmd{ii}.wantrhobar=true;
cmd{ii}.wantrxx=true;
cmd{ii}.mlow=[-15 -15 -15];
cmd{ii}.mhigh=[15 15 15];
cmd{ii}.deltax=4.0 .* [4.7 4.7 4.7];
cmd{ii}.EulerAngles=[0.0 0.0 0.0];

%Save the results of computing real space cubes as a Matlab .mat file.
cmd{ii}.operator='post_save_real_space_cubes';
cmd{ii}.fn_real_space_cubes=[basename '.cubes.mat'];

%Use imagesc to plot three crossections through the rhobar and rxx cubes
cmd{ii}.operator='post_crossections_real_space_cubes';

%Compute FSC between two class that are in the same virusobj cell array or in two different virusobj cell arrays.
cmd{ii}.operator='post_compute_FSC';
cmd{ii}.FSC_minmagk=0.0;
cmd{ii}.FSC_maxmagk=0.5/4.7;
cmd{ii}.FSC_deltamagk=(cmd{ii}.FSC_maxmagk - cmd{ii}.FSC_minmagk)/100;
cmd{ii}.FSC_eta4classA=1;
cmd{ii}.FSC_eta4classB=2;
cmd{ii}.FSC_is_same_vobj=true; %if 'false' not 'true' then Class A is the virus object that was pushed onto the stack.
cmd{ii}.FSC_scaling=[1 1];% scaling factors 

%Low pass filter the FSC curve.
cmd{ii}.operator='post_LPF_FSC';
cmd{ii}.LPFcutoff4FSC=0.2;

%Low pass filter the FSC curve with moving average
cmd{ii}.operator='post_movingavg_FSC';

%Determine the resolution from the FSC curve.
cmd{ii}.operator='post_FSC_cutoff';
cmd{ii}.FSCcutoff=0.5;

%Plot FSC.
cmd{ii}.operator='post_FSC_plot';
cmd{ii}.plotcutoffk=true;
cmd{ii}.plottitle='blah-blah-blah'; %if not defined or if ='' then nothing is written.
cmd{ii}.FSCtitle='FSC.eps';% filename to be saved
cmd{ii}.energylegend1='set1';
cmd{ii}.energylegend2='set2';
cmd{ii].energytitle='energy.eps';

%Save FSC results to a Matlab .mat file.
cmd{ii}.operator='post_save_FSC';
cmd{ii}.fn_FSC=[basename 'FSC.mat'];

%Write FSC results to a text file.
cmd{ii}.operator='post_write_FSC';
cmd{ii}.fn_FSC=[basename 'FSC.txt'];

%Select capsid region
cmd{ii}.operator='post_capsid4scale';
cmd{ii}.scaling_fn1='A.realspace.mat';
cmd{ii}.scaling_fn2='B.realspace.mat';
cmd{ii}.capsid_r=[195 300];
cmd{ii}.capsid_thres=[0.1 0.1];

%Compute scaling factor 
cmd{ii}.operator='post_scaling';
cmd{ii}.scaling_type=0;% 0 for rho_new=g*rho_old; 1 for rho_new=a*rho_old+b

%Compute spherical average of variance 
cmd{ii}.operator='post_compute_spherical_average';
cmd{ii}.sphavg_eta2scale=1;
cmd{ii}.sphavg_scale=1;%''means no additional input, use the scaling_fac in workspace calculated by previous operators

%Save the results of spherical average calculation 
cmd{ii}.operator='post_save_spherical_average';
cmd{ii}.fn_spherical_avg='HK97.noproset1.spherical.average.mat';%save the stdv vector and unique_radius vector for plotting

%Plot the square root of spherical average of variance
cmd{ii}.operator='post_plot_spherical_averages';
cmd{ii}.fn_avg_plot1='HK97.noproset1.spherical.average.mat';
cmd{ii}.fn_avg_plot2='HK97.withproset1.spherical.average.mat';
cmd{ii}.sphavg_legend1='no protease';
cmd{ii}.sphavg_legend2='with protease';
cmd{ii}.sphavg_innercap= 50;% inner radius for capsid
cmd{ii}.sphavg_title_full='Spherically.averaged.square.root.of.variance.full.range.eps';
cmd{ii}.sphavg_title_cap='Spherically.averaged.square.root.of.variance.capsid.eps';% zoomed-in plot for capsid region

%Plot the histogram of deviation in orientation
cmd{ii}.operator='post_hist_orientation_deviation';
cmd{ii}.fn_initial_orientation='pthetaeta1.mat';% p_theta_eta matrix at first iteration
cmd{ii}.fn_final_orientation='pthetaeta.final.mat';% p_theta_eta matrix at final iteration
cmd{ii}.hist_title='Histogram.orientation.deviation.eps';
cmd{ii}.metric_type=4;% choose from different definition of metric for two orientation

%Plot FSC generated by EMAN
cmd{ii}.fn_EMAN_FSC='e2fsc.txt';% generated by EMAN
cmd{ii}.fn_EMAN_FSC_title='EMAN.FSC.eps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create an image stack of noisy images of a spherical object with constant electron scattering intensity
cmd{ii}.operator='test_imagestackfromsphericalobj';
cmd{ii}.numberimages=100;
cmd{ii}.minxindex=-100;
cmd{ii}.maxxindex=100;
cmd{ii}.minyindex=-100;
cmd{ii}.maxyindex=100;
cmd{ii}.deltaxy=2.74;
cmd{ii}.rho0=1.0;
cmd{ii}.R2=500.0;
cmd{ii}.pixelnoisevariance=10.0;

%Test the FFT and processes for creating the reciprocal space image vectors
%First part of the test -- make real space images
cmd{ii}.operator='test_FFT_pt1';
cmd{ii}.N1=20;
cmd{ii}.N2=21;
cmd{ii}.delta1=2.74;
cmd{ii}.delta2=3.24;
cmd{ii}.freqintegers=[0 0 ; 1 0 ; 0 1 ; -1 0 ; 0 -1 ; 1 1 ; -1 -1 ; 1 -1 ; -1 1];

%Test the FFT and processes for creating the reciprocal space image vectors
%Second part of the test -- test reciprocal space image vectors
cmd{ii}.operator='test_FFT_pt2';
%requires nothing but uses cmd_test_FFT_pt1 set by test_FFT_pt1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%not yet a part of hetero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zero_one_cube_from_angularbasisfunc.m
zero_one_cube_from_angularbasisfunc.v1.m
zero_one_cube_from_angularbasisfunc_main.m
RealBasisFunctionCoeff_RealUniIrrepPhyPiover2_noStructure_1012.txt
rd_b_general.m
