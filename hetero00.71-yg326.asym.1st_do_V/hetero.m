function hetero(cmd)
%function hetero(cmd)

if isempty(cmd)
  return;
end

for ii=1:length(cmd)

  if isempty(cmd{ii})
    fprintf(1,'hetero: ii %d no-op command\n',ii);
    continue;
  end

  switch cmd{ii}.operator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%misc operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'misc_diary'
      %requires cmd{ii}.fn_diary
      fprintf(1,'hetero: ii %d misc_diary fn_diary %s\n',ii,cmd{ii}.fn_diary);
      if ~isempty(cmd{ii}.fn_diary)
        diary(cmd{ii}.fn_diary);
      end

    case 'misc_addtomatlabpath'
      %requires cmd{ii}.dn %directory name
      fprintf(1,'hetero: ii %d misc_addtomatlabpath dn %s\n',ii,cmd{ii}.dn);
      if ~isempty(cmd{ii}.dn)
        addpath(cmd{ii}.dn);
      end

    case 'misc_setpseudorandomnumberseed'
      %requires cmd{ii}.pseudorandomnumberseed
      fprintf(1,'hetero: ii %d misc_setpseudorandomnumberseed pseudorandomnumberseed %d\n', ...
              ii,cmd{ii}.pseudorandomnumberseed);
      if cmd{ii}.pseudorandomnumberseed<0
        error('hetero: cmd{ii}.pseudorandomnumberseed %d < 0\n',cmd{ii}.pseudorandomnumberseed);
      end
      rng(cmd{ii}.pseudorandomnumberseed);

    case 'misc_savepseudorandomnumberstate2file'
      %requires cmd{ii}.fn_pseudorandomnumberstate'
      fprintf(1,'hetero: ii %d misc_savepseudorandomnumberstate2file fn_pseudorandomnumberstate %s\n', ...
              ii,cmd{ii}.pseudorandomnumberseed);
      state4rng=rng;
      save(cmd{ii}.fn_pseudorandomnumberstate,'state4rng');
      clear state4rng;

    case 'misc_restorepseudorandomnumberstatefromfile'
      %requires cmd{ii}.fn_pseudorandomnumberstate'
      fprintf(1,'hetero: ii %d misc_restorepseudorandomnumberstate2file fn_pseudorandomnumberstate %s\n', ...
              ii,cmd{ii}.pseudorandomnumberseed);
      load(cmd{ii}.fn_pseudorandomnumberstate,'state4rng');
      rng(state4rng);
      clear state4rng;

    case 'misc_clearvariable'
      %requires cmd{ii}.variablename
      fprintf(1,'hetero: ii %d misc_clearvariable variablename %s\n',ii,cmd{ii}.variablename);
      clear(cmd{ii}.variablename);

    case 'misc_changedirectory'
      %requires cmd{ii}.dn directory name
      fprintf(1,'hetero: ii %d misc_changedirectory dn %s\n',ii,cmd{ii}.dn);
      ifnotexistmkdir(cmd{ii}.dn);
      cd(cmd{ii}.dn);

    case 'misc_return'
      fprintf(1,'hetero: ii %d misc_return\n',ii);
      return;

    case 'misc_keyboard'
      fprintf(1,'hetero: ii %d misc_keyboard, enter "return" to continue, "dbquit" to halt\n',ii);
      keyboard;

    case 'misc_save_workspace'
      %requires cmd{ii}.fn_workspace
      fprintf(1,'hetero: ii %d misc_save_workspace fn_workspace %s\n',ii,cmd{ii}.fn_workspace);
      state4rng=rng; %make sure that the state of the pseudorandom number generator is in the workspace
      save(cmd{ii}.fn_workspace);

    case 'misc_load_workspace'
      %requires cmd{ii}.fn_workspace
      %requires cmd{ii}.ExistingOverReload
      fprintf(1,'hetero: ii %d misc_load_workspace fn_workspace %s existingoverreload %d\n', ...
              ii,cmd{ii}.fn_workspace,cmd{ii}.existingoverreload);
      if cmd{ii}.existingoverreload
        state4rng=rng; %make sure that the state of the pseudorandom number generator is in the workspace
        save tmp_workspace;
        load(cmd{ii}.fn_workspace);
        load tmp_workspace;
        delete tmp_workspace;
        rng(state4rng);
      else
        load(cmd{ii}.fn_workspace);
        if exist('state4rng','var')
          rng(state4rng);
        else
          fprintf(1,'hetero: misc_load_workspace: state4rng does not exist\n');
        end
      end

    case 'misc_write_mrc'
      %requires cmd{ii}.fn_write_mrc
      %requires cmd{ii}.what2write
      fprintf(1,'hetero: ii %d misc_write_mrc fn_write_mrc %s what2write %s\n', ...
              ii,cmd{ii}.fn_write_mrc,cmd{ii}.what2write);

      switch cmd{ii}.what2write

        case 'write_image_stack'
          if max(deltachi)~=min(deltachi)
            error('hetero: misc_write_mrc: deltachi %g %g\n',deltachi(1),deltachi(2));
          end
          map=zeros(size(imagestack{1},1),size(imagestack{1},2),length(imagestack));
          for jj=1:length(imagestack)
            map(:,:,jj)=imagestack{jj};
          end
          WriteMRC(map,deltachi(1),cmd{ii}.fn_write_mrc);
          clear jj map;

        case 'write_Image_stack'
          if max(deltachi)~=min(deltachi)
            error('hetero: misc_write_mrc: deltachi %g %g\n',deltachi(1),deltachi(2));
          end
          map=zeros(size(Imagestack{1},1),size(Imagestack{1},2),2*length(Imagestack));
          for jj=1:2*length(Imagestack)
            if mod(jj-1,2)==0
              map(:,:,jj)=real(Imagestack{floor((jj-1)/2)+1});
            else
              map(:,:,jj)=imag(Imagestack{floor((jj-1)/2)+1});
            end
          end
          WriteMRC(map,deltachi(1),cmd{ii}.fn_write_mrc);
          clear jj map;

        case 'write_rhobar'
          if max(real_space_cubes.deltax)~=min(real_space_cubes.deltax)
            error('hetero: misc_write_mrc: deltax %g %g %g\n', ...
                  real_space_cubes.deltax(1),real_space_cubes.deltax(2),real_space_cubes.deltax(3));
          end
          WriteMRC(rhobar,real_space_cubes.deltax(1),cmd{ii}.fn_write_mrc);

        case 'write_rxx'
          if max(real_space_cubes.deltax)~=min(real_space_cubes.deltax)
            error('hetero: misc_write_mrc: deltax %g %g %g\n', ...
                  real_space_cubes.deltax(1),real_space_cubes.deltax(2),real_space_cubes.deltax(3));
          end
          WriteMRC(rxx,real_space_cubes.deltax(1),cmd{ii}.fn_write_mrc);
        case 'write_scaled_rhobar'
          if max(real_space_cubes.deltax)~=min(real_space_cubes.deltax)
            error('hetero: misc_write_mrc: deltax %g %g %g\n', ...
                  real_space_cubes.deltax(1),real_space_cubes.deltax(2),real_space_cubes.deltax(3));
          end
          WriteMRC(rhobarscaled,real_space_cubes.deltax(1),cmd{ii}.fn_write_mrc);

        case 'write_scaled_rxx'
          if max(real_space_cubes.deltax)~=min(real_space_cubes.deltax)
            error('hetero: misc_write_mrc: deltax %g %g %g\n', ...
                  real_space_cubes.deltax(1),real_space_cubes.deltax(2),real_space_cubes.deltax(3));
          end
          WriteMRC(rxxscaled,real_space_cubes.deltax(1),cmd{ii}.fn_write_mrc);

        otherwise
          error('hetero: ii %d misc_write_mrc: what2write %s\n',ii,cmd{ii}.what2write);

      end %close switch cmd{ii}.what2write

    case 'misc_push'
      %requires cmd{ii}.what2push
      fprintf(1,'hetero: ii %d misc_push what2push %s\n',ii,cmd{ii}.what2push);
      switch cmd{ii}.what2push

        case 'push_image_stack'
          imagestack2=imagestack;

        case 'push_virusobj'
          vobj2=vobj;

        otherwise
          error('hetero: ii %d misc_push: what2push %s\n',ii,cmd{ii}.what2push);

      end %close switch cmd{ii}.what2push

    case 'misc_pop'
      %requires cmd{ii}.what2pop
      fprintf(1,'hetero: ii %d misc_pop what2pop %s\n',ii,cmd{ii}.what2pop);
      switch cmd{ii}.what2pop

        case 'pop_image_stack'
          imagestack=imagestack2;
          clear imagestack2;

        case 'pop_virusobj'
          vobj=vobj2;
          clear vobj2;

        otherwise
          error('hetero: ii %d misc_pop: what2pop %s\n',ii,cmd{ii}.what2pop);

      end %close switch cmd{ii}.what2pop

    case 'misc_pack'
      fprintf(1,'hetero: ii %d misc_pack\n',ii);
      pack;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%basic operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'basic_set2Drealspacesamplingintervals'
      %requires cmd{ii}.samplingintervals(1:2)
      fprintf(1,'hetero: ii %d basic_set2Drealspacesamplingintervals deltachi(1) %g deltachi(2) %g\n', ...
              ii,cmd{ii}.samplingintervals(1),cmd{ii}.samplingintervals(2));
      deltachi=cmd{ii}.samplingintervals;

    case 'basic_setsizeof2Drealspaceimages'
      %requires cmd{ii}.NaNb(1:2)
      fprintf(1,'hetero: ii %d basic_setsizeof2Drealspaceimages Na %d Nb %d\n', ...
              ii,cmd{ii}.NaNb(1),cmd{ii}.NaNb(2));
      Na=cmd{ii}.NaNb(1);
      Nb=cmd{ii}.NaNb(2);

    case 'basic_setsizeof2Drealspaceimagesfromimagestack'
      fprintf(1,'hetero: ii %d basic_setsizeof2Drealspaceimagesfromimagestack\n',ii);
      Na=size(imagestack{1},1);
      Nb=size(imagestack{1},2);

    case 'basic_setsizeof2DreciprocalspaceImagesfromImagestack'
      fprintf(1,'hetero: ii %d basic_setsizeof2DreciprocalspaceImagesfromImagestack\n',ii);
      Na=size(Imagestack{1},1);
      Nb=size(Imagestack{1},2);

    case 'basic_compute2Dreciprocalspaceproperties'
      fprintf(1,'hetero: ii %d basic_compute2Dreciprocalspaceproperties\n',ii)
      %Determine a minimal subset of 2-D reciprocal space such that
      %conjugate symmetry fills in all of reciprocal space.  Assume
      %that all images have the same dimensions and extract the
      %dimensions from the first of the boxed 2-D real-space images.
      [Ny4minimalset,iiminimalset,vkminimalset,ixminimalset,iyminimalset]= ...
        vk_indexset(Na,Nb,deltachi(1),deltachi(2));

      %Compute magnitude of the reciprocal space vector.
      vkmag=sqrt( vkminimalset(:,1).^2 + vkminimalset(:,2).^2 );

      %Package in a structure for simplicity.
      recipobj.Ny4minimalset=Ny4minimalset;
      recipobj.iiminimalset=iiminimalset;
      recipobj.vkminimalset=vkminimalset;
      recipobj.ixminimalset=ixminimalset;
      recipobj.iyminimalset=iyminimalset;
      recipobj.vkmag=vkmag;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real-space boxed image operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'box_readpixelnoisevar'
      %requires cmd{ii}.fn_pixelnoisevar
      fprintf(1,'hetero: ii %d box_readpixelnoisevar fn_pixelnoisevar %s\n', ...
              ii,cmd{ii}.fn_pixelnoisevar);
      fid=fopen(cmd{ii}.fn_pixelnoisevar,'r');
      if fid==-1
        error('hetero: ii %d box_readpixelnoisevar: fn_pixelnoisevar %s\n', ...
              ii,cmd{ii}.fn_pixelnoisevar);
      end
      [pixelnoisevar,count]=fscanf(fid,'%g');
      if count~=1
        error('hetero: ii %d box_readpixelnoisevar count %d ~= 1\n',ii,count);
      end
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d box_readpixelnoisevar status %d ~= 0\n',ii,status);
      end
      clear fid count status;
      fprintf(1,'hetero: ii %d box_readpixelnoisevar pixelnoisevar %g\n',ii,pixelnoisevar);

    case 'box_writepixelnoisevar'
      %requires cmd{ii}.fn_pixelnoisevar
      fprintf(1,'hetero: ii %d box_writepixelnoisevar fn_pixelnoisevar %s\n',ii,cmd{ii}.fn_pixelnoisevar);
      fid=fopen(cmd{ii}.fn_pixelnoisevar,'w');
      if fid==-1
        error('hetero: ii %d box_writepixelnoisevar: fn_pixelnoisevar %s\n', ...
              ii,cmd{ii}.fn_pixelnoisevar);
      end
      fprintf(fid,'%g\n',pixelnoisevar);
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d box_writepixelnoisevar status %d ~= 0\n',ii,status);
      end
      clear fid status;

    case 'box_printimagesasarrays'
      fprintf(1,'hetero: ii %d box_printimagesasarrays: imageindex:\n',ii);
      disp(imageindex');
      for jj=1:length(imagestack)
        fprintf(1,'hetero: box_printimagesasarrays: image %d:\n',jj);
        disp(imagestack{jj});
      end
      clear jj;

    case 'box_readimagestack'
      %requires cmd{ii}.imagestackformat
      %requires cmd{ii}.fn_imagestack
      fprintf(1,'hetero: ii %d box_readimagestack: imagestackformat %s fn_imagestack %s\n', ...
              ii,cmd{ii}.imagestackformat,cmd{ii}.fn_imagestack);
      switch cmd{ii}.imagestackformat

        case 'fake'
          imagestack=readreadread(cmd{ii}.fn_imagestack);
          imageindex=[1:length(imagestack)]';

        case 'mrc'
          %requires also cmd{ii}.startSlice
          %requires also cmd{ii}.numSlices
          fprintf(1,'hetero: ii %d box_readimagestack: mrc: startSlice %d numSlices %d\n', ...
                  ii,cmd{ii}.startSlice,cmd{ii}.numSlices);
          map=ReadMRC(cmd{ii}.fn_imagestack,cmd{ii}.startSlice,cmd{ii}.numSlices);
          if size(map,3)~=cmd{ii}.numSlices
            error('hetero: ii %d box_readimagestack: mrc: size(map,3) %d\n',ii,size(map,3))
          end
          imagestack=cell(size(map,3),1);
          for jj=1:length(imagestack)
            imagestack{jj}=map(:,:,jj);
          end
          imageindex=[cmd{ii}.startSlice:cmd{ii}.startSlice+cmd{ii}.numSlices-1];
          clear map jj;

        case 'img'
          %requires also cmd{ii}.startimage
          %requires also cmd{ii}.numimages
          fprintf(1,'hetero: ii %d box_readimagestack: img: startimage %d numimages %d\n', ...
                  ii,cmd{ii}.startimage,cmd{ii}.numimages);
          [map,info]=ReadImagic(cmd{ii}.fn_imagestack,cmd{ii}.startimage,cmd{ii}.numimages);
          if size(map,3)~=cmd{ii}.numimages
            error('hetero: ii %d box_readimagestack: img: size(map,3) %d\n',ii,size(map,3))
          end
          imagestack=cell(size(map,3),1);
          for jj=1:length(imagestack)
            imagestack{jj}=map(:,:,jj);
          end
          imageindex=[cmd{ii}.startimage:cmd{ii}.startimage+cmd{ii}.numimages-1];
          clear map jj;

        otherwise
          error('hetero: ii %d unknown imagestackformat %s\n',ii,cmd{ii}.imagestackformat);

      end %end of switch cmd{ii}.imagestackformat

    case 'box_readImagestack'
      %requires cmd{ii}.Imagestackformat
      %requires cmd{ii}.fn_Imagestack
      fprintf(1,'hetero: ii %d box_readImagestack: Imagestackformat %s fn_Imagestack %s\n', ...
              ii,cmd{ii}.Imagestackformat,cmd{ii}.fn_Imagestack);
      switch cmd{ii}.Imagestackformat

        case 'fake'
          Imagestack=readreadread(cmd{ii}.fn_Imagestack);
          Imageindex=[1:length(Imagestack)]';

        case 'mrc'
          %requires also cmd{ii}.startSlice
          %requires also cmd{ii}.numSlices %Number of real images.  Need two real images for each complex reciprocal-space image.
          fprintf(1,'hetero: ii %d box_readImagestack: mrc: startSlice %d numSlices %d\n', ...
                  ii,cmd{ii}.startSlice,cmd{ii}.numSlices);
          map=ReadMRC(cmd{ii}.fn_Imagestack,cmd{ii}.startSlice,2*cmd{ii}.numSlices);
          if size(map,3)~=2*cmd{ii}.numSlices
            error('hetero: ii %d box_readImagestack: mrc: size(map,3) %d\n',ii,size(map,3))
          end
          Imagestack=cell(size(map,3)/2,1);
          for jj=1:length(Imagestack)
            Imagestack{jj}=map(:,:,2*jj-1) + sqrt(-1)*map(:,:,2*jj);
          end
          Imageindex=[cmd{ii}.startSlice:cmd{ii}.startSlice+cmd{ii}.numSlices-1];
          clear map jj;

        case 'img'
          error('hetero: ii %d box_readImagestack: img is not implemented\n');
          %requires also cmd{ii}.startImage
          %requires also cmd{ii}.numImages
          fprintf(1,'hetero: ii %d box_readImagestack: img: startImage %d numImages %d\n', ...
                  ii,cmd{ii}.startImage,cmd{ii}.numImages);
          [map,info]=ReadImagic(cmd{ii}.fn_Imagestack,cmd{ii}.startImage,cmd{ii}.numImages);
          if size(map,3)~=cmd{ii}.numImages
            error('hetero: ii %d box_readImagestack: img: size(map,3) %d\n',ii,size(map,3))
          end
          Imagestack=cell(size(map,3),1);
          for jj=1:length(Imagestack)
            Imagestack{jj}=map(:,:,jj);
          end
          Imageindex=[cmd{ii}.startImage:cmd{ii}.startImage+cmd{ii}.numImages-1];
          clear map jj;

        otherwise
          error('hetero: ii %d unknown Imagestackformat %s\n',ii,cmd{ii}.Imagestackformat);

      end %end of switch cmd{ii}.Imagestackformat

    case 'box_loadimagestack'
      %requires cmd{ii}.fn_imagestack
      %requires cmd{ii}.only_imagestack
      fprintf(1,'hetero: ii %d box_loadimagestack: fn_imagestack %s only_imagestack %d\n',ii,cmd{ii}.fn_imagestack,cmd{ii}.only_imagestack);
      if cmd{ii}.only_imagestack
        load(cmd{ii}.fn_imagestack,'imagestack');
      else
        load(cmd{ii}.fn_imagestack,'imagestack','imageindex');
      end

    case 'box_loadImagestack'
      %requires cmd{ii}.fn_Imagestack
      fprintf(1,'hetero: ii %d box_loadImagestack: fn_Imagestack %s\n',ii,cmd{ii}.fn_Imagestack);
      load(cmd{ii}.fn_Imagestack,'Imagestack','imageindex');

    case 'box_saveimagestack'
      %requires cmd{ii}.fn_imagestack
      fprintf(1,'hetero: ii %d box_saveimagestack: fn_imagestack %s\n',ii,cmd{ii}.fn_imagestack);
      save(cmd{ii}.fn_imagestack,'imagestack','imageindex');

    case 'box_saveImagestack'
      %requires cmd{ii}.fn_Imagestack
      fprintf(1,'hetero: ii %d box_saveImagestack: fn_Imagestack %s\n',ii,cmd{ii}.fn_Imagestack);
      save(cmd{ii}.fn_Imagestack,'Imagestack','imageindex');

    case 'box_extractsubset'
      %requires cmd{ii}.maxNv
      %requires cmd{ii}.a
      %requires cmd{ii}.b
      fprintf(1,'hetero: ii %d box_extractsubset: maxNv %d a %d b %d\n', ...
              ii,cmd{ii}.maxNv,cmd{ii}.a,cmd{ii}.b);
      Nv=box_extractsubset_getnewNv(cmd{ii}.maxNv,cmd{ii}.a,cmd{ii}.b,length(imagestack));
      tmpimagestack=cell(Nv,1);
      tmpimageindex=zeros(Nv,1);
      for jj=1:Nv
        kk=(cmd{ii}.a)+(cmd{ii}.b)*(jj-1);
        tmpimagestack{jj}=imagestack{kk};
        tmpimageindex(jj)=imageindex(kk);
      end
      imagestack=tmpimagestack;
      imageindex=tmpimageindex;
      clear tmpimagestack tmpimageindex Nv jj kk;

    case 'box_permute'
      fprintf(1,'hetero: ii %d box_permute\n',ii);
      p=randperm(length(imagestack));
%      fprintf(1,'hetero: box_permute: p:\n');
%      disp(p);
      tmpimagestack=cell(size(imagestack));
      tmpimageindex=zeros(length(imagestack),1);
      for jj=1:length(imagestack)
        tmpimagestack{jj}=imagestack{p(jj)};
        tmpimageindex(jj)=imageindex(p(jj));
      end
      imagestack=tmpimagestack;
      imageindex=tmpimageindex;
      clear p tmpimagestack tmpimageindex jj;

    case 'box_shrink'
      %requires cmd{ii}.pixels2delete(1:4)
      %could be done image-by-image.
      %assume all images have the same size so mask can be computed just once.
      fprintf(1,'hetero: ii %d box_shrink: pixels2delete %d %d %d %d\n',ii,cmd{ii}.pixels2delete(1), ...
              cmd{ii}.pixels2delete(2),cmd{ii}.pixels2delete(3),cmd{ii}.pixels2delete(4));
      boxmask=box_shrink_mask(size(imagestack{1}),cmd{ii}.pixels2delete);
      for jj=1:length(imagestack)
        tmp=imagestack{jj};
        imagestack{jj}=reshape(tmp(boxmask), ...
                               size(tmp,1)-cmd{ii}.pixels2delete(2)-cmd{ii}.pixels2delete(4), ...
                               size(tmp,2)-cmd{ii}.pixels2delete(1)-cmd{ii}.pixels2delete(3));
      end
      clear boxmask jj tmp;

    case 'box_maskcorners'
      %requires cmd{ii}.radiuscorner
      %could be done image-by-image.
      %assume all images have the same size so mask can be computed just once.
      fprintf(1,'hetero: ii %d box_maskcorners: radiuscorner %g\n',ii,cmd{ii}.radiuscorner);
      [corner center]=box_maskcorners_mask(size(imagestack{1}),deltachi,cmd{ii}.radiuscorner);
      for jj=1:length(imagestack)
        samplemean=mean(imagestack{jj}(center));
        imagestack{jj}(corner)=samplemean;
      end
      clear corner center jj samplemean;

    case 'box_normalize2zeroone'
      %requires cmd{ii}.radius01(1:2)
      %could be done image-by-image.
      %assume all images have the same size so mask can be computed just once.
      fprintf(1,'hetero: ii %d box_normalize2zeroone: radius01 %g %g\n', ...
              ii,cmd{ii}.radius01(1),cmd{ii}.radius01(2));
      normalizemask=box_normalize_mask(size(imagestack{1}),deltachi,cmd{ii}.radius01);
      for jj=1:length(imagestack)
        samplemean=mean(imagestack{jj}(normalizemask));
        samplevariance=cov(imagestack{jj}(normalizemask));
        imagestack{jj}=(imagestack{jj}-samplemean)./sqrt(samplevariance);
      end
      clear normalizemask jj samplemean samplevariance;

    case 'box_classifyviasamplemean'
      %requires cmd{ii}.classifythres
      %must be done collectively over all images.
      %all images must be the same size.
      fprintf(1,'hetero: ii %d box_classifyviasamplemean: classifythres %g\n',ii,cmd{ii}.classifythres);

      if cmd{ii}.classifythres<0.0 | cmd{ii}.classifythres>1.0
        error('hetero: ii %d box_classifyviasamplemean classifythres %g\n',ii,cmd{ii}.classifythres);
      end

      samplemean=imagestack{1};
      for jj=2:length(imagestack)
        samplemean=samplemean+imagestack{jj};
      end
      samplemean=samplemean./length(imagestack);

      l2norms=zeros(length(imagestack));
      for jj=1:length(imagestack)
        difference=imagestack{jj}-samplemean;
        l2norms(jj)=norm(difference(:));
      end

      clear samplemean difference;

      [values,indices]=sort(l2norms); %ascending order
      maxindex=(1.0-cmd{ii}.classifythres)*length(l2norms);
      maxindex=fix(maxindex);
      if maxindex>length(l2norms)
        error('hetero: ii %d box_classifyviasamplemean maxindex %d > length(l2norms) %d\n', ...
              ii,maxindex,length(l2norms));
      end
      if maxindex<1
        error('hetero: ii %d box_classifyviasamplemean maxindex %d < 1\n',ii,maxindex);
      end
      goodindices=indices(1:maxindex);
      goodindices=sort(goodindices); %in order to shorten image stack in place

      clear l2norms values indices maxindex;

      tmpimagestack=cell(size(goodindices));
      tmpimageindex=zeros(size(goodindices));
      for jj=1:length(goodindices)
        tmpimagestack{jj}=imagestack{goodindices(jj)};
        tmpimageindex(jj)=imageindex(goodindices(jj));
      end
      imagestack=tmpimagestack;
      imageindex=tmpimageindex;
      clear tmpimagestack tmpimageindex goodindices jj;

    case 'box_annulusstatistics'
      %requires cmd{ii}.radius01(1:2)
      %must be done collectively over all images.
      %all images must be the same size.
      fprintf(1,'hetero: ii %d box_annulusstatistics: radius01 %g %g\n', ...
              ii,cmd{ii}.radius01(1),cmd{ii}.radius01(2));
      normalizemask=box_normalize_mask(size(imagestack{1}),deltachi,cmd{ii}.radius01);
      npixels=length(find(normalizemask==true));
      fprintf(1,'hetero: ii %d box_annulusstatistics: number of pixels in the annulus %d\n', ...
              ii,npixels);
      annulussamplemean=0.0;
      for jj=1:length(imagestack)
        annulussamplemean=annulussamplemean+sum(imagestack{jj}(normalizemask));
      end
      annulussamplemean=annulussamplemean/(npixels*length(imagestack));
      annulussamplevariance=0.0;
      for jj=1:length(imagestack)
        annulussamplevariance=annulussamplevariance+sum((imagestack{jj}(normalizemask)-annulussamplemean).^2);
      end
      annulussamplevariance=annulussamplevariance/(npixels*length(imagestack)-1);
      fprintf(1,'hetero: ii %d box_annulusstatistics: annulussamplemean %g annulussamplevariance %g\n', ...
              ii,annulussamplemean,annulussamplevariance);
      clear normalizemask npixels jj;

    case 'box_writeannulusstatistics'
      %requres cmd{ii}.fn_annulusstatistics
      fprintf(1,'hetero: ii %d box_writeannulusstatistics fn_annulusstatistics %s\n', ...
              ii,cmd{ii}.fn_annulusstatistics);
      fid=fopen(cmd{ii}.fn_annulusstatistics,'w');
      if fid==-1
        error('hetero: ii %d box_writeannulusstatistics: fn_annulusstatistics %s\n', ...
              ii,cmd{ii}.fn_annulusstatistics);
      end
      fprintf(fid,'%g %g\n',annulussamplemean,annulussamplevariance);
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d box_writeannulusstatistics status %d ~= 0\n',ii,status);
      end
      clear fid status;

    case 'box_readannulusstatistics'
      %requres cmd{ii}.fn_annulusstatistics
      fprintf(1,'hetero: ii %d box_readannulusstatistics fn_annulusstatistics %s\n', ...
              ii,cmd{ii}.fn_annulusstatistics);
      fid=fopen(cmd{ii}.fn_annulusstatistics,'r');
      if fid==-1
        error('hetero: ii %d box_readannulusstatistics: fn_annulusstatistics %s\n', ...
              ii,cmd{ii}.fn_annulusstatistics);
      end
      [tmp_statistics,count]=fscanf(fid,'%g%g');
      if count~=2
        error('hetero: ii %d box_readannulusstatistics count %d ~= 2\n',ii,count);
      end
      annulussamplemean=tmp_statistics(1);
      annulussamplevariance=tmp_statistics(2);
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d box_readannulusstatistics status %d ~= 0\n',ii,status);
      end
      clear fid tmp_statistics count status;
      fprintf(1,'hetero: ii %d box_readannulusstatistics annulussamplemean %g annulussamplevariance %g\n', ...
              ii,annulussamplemean,annulussamplevariance);

    case 'box_showimagestack'
      %requires cmd{ii}.whichimages
      %requires cmd{ii}.printtitle
      %requires cmd{ii}.fn_class
      fprintf(1,'hetero: ii %d box_showimagestack: whichimages:',ii);
      fprintf(1,' %d',cmd{ii}.whichimages);
      fprintf(1,'\n');
      fprintf(1,'hetero: ii %d box_showimagestack: printtitle %d\n', ...
              ii,cmd{ii}.printtitle);
      fprintf(1,'hetero: ii %d box_showimagestack: fn_class %s\n', ...
              ii,cmd{ii}.fn_class);
      fprintf(1,'hetero: ii %d box_showimagestack: each image has a different colormap\n',ii);

      if length(cmd{ii}.fn_class)>0
        truevalues_box_showimagestack=load(cmd{ii}.fn_class);
      end

      for jj=cmd{ii}.whichimages
        figure;
        imagesc(imagestack{jj});
        colormap(gray);
        axis equal;
        axis off;
        if cmd{ii}.printtitle
          titlestr=['image ' num2str(jj)];
          if length(cmd{ii}.fn_class)>0
            titlestr=[titlestr ', class ' num2str(truevalues_box_showimagestack(jj,1))];
          end
          title(titlestr);
        end
      end
      clear jj;
      if length(cmd{ii}.fn_class)>0
        clear truevalues_box_showimagestack;
      end
      if cmd{ii}.printtitle
        clear titlestr;
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%real space <-> reciprocal space transformations operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'realrecip_2DFFT'
      fprintf(1,'hetero: ii %d basic_realrecip_2DFFT\n',ii);
      Imagestack=cell(size(imagestack));
      for ii=1:length(Imagestack)
        Imagestack{ii}=deltachi(1)*deltachi(2)*fft2(ifftshift(imagestack{ii}));
      end
      clear ii;

    case 'realrecip_2DIFFT'
      fprintf(1,'hetero: ii %d basic_realrecip_2DIFFT\n',ii);
      imagestack=cell(size(Imagestack));
      for ii=1:length(Imagestack)
        imagestack{ii}=(1/(deltachi(1)*deltachi(2)))*fftshift(ifft2(ImgAs2DComplexImage)); %maybe add 'symmetric'
      end
      clear ii;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reciprocal-space boxed image operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'Box_printImagesasarrays'
      fprintf(1,'hetero: ii %d Box_printImagesasarrays: imageindex:\n',ii);
      disp(imageindex');
      for jj=1:length(Imagestack)
        fprintf(1,'hetero: Box_printImagesasarrays: Image %d:\n',jj);
        disp(Imagestack{jj});
      end
      clear jj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%virus object operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'vobj_print_virusobj'
      fprintf(1,'hetero: ii %d vobj_print_virusobj\n',ii);
      for eta=1:length(vobj)
        fprintf(1,'vobj_print_virusobj: eta %d clnp_fn %s nu_fn %s q_fn %s\n',eta,vobj{eta}.clnp_fn,vobj{eta}.nu_fn,vobj{eta}.q_fn);
        fprintf(1,'vobj_print_virusobj: eta %d clnp.il:\n',eta);
        disp(vobj{eta}.clnp.il');
        fprintf(1,'vobj_print_virusobj: eta %d clnp.in:\n',eta);
        disp(vobj{eta}.clnp.in');
        fprintf(1,'vobj_print_virusobj: eta %d clnp.ip:\n',eta);
        disp(vobj{eta}.clnp.ip');
        fprintf(1,'vobj_print_virusobj: eta %d clnp.optflag:\n',eta);
        disp(vobj{eta}.clnp.optflag');
        fprintf(1,'vobj_print_virusobj: eta %d clnp.c:\n',eta);
        disp(vobj{eta}.clnp.c');
        fprintf(1,'vobj_print_virusobj: eta %d cbar:\n',eta);
        disp(vobj{eta}.cbar');
        fprintf(1,'vobj_print_virusobj: eta %d BasisFunctionType %d\n',eta,vobj{eta}.BasisFunctionType);
        fprintf(1,'vobj_print_virusobj: eta %d R1 %g R2 %g\n',eta,vobj{eta}.R1,vobj{eta}.R2);
        fprintf(1,'vobj_print_virusobj: eta %d nu:\n',eta);
        disp(vobj{eta}.nu');
        fprintf(1,'vobj_print_virusobj: eta %d q %g\n',eta,vobj{eta}.q);
      end

    case 'vobj_read_virusobj'
      %requires cmd{ii}.fn_clnp, cell array of Neta file names
      %requires cmd{ii}.fn_nu, cell array of Neta file names
      %requires cmd{ii}.fn_q, cell array of Neta file names
      for jj=1:length(cmd{ii}.fn_clnp)
        fprintf(1,'hetero: ii %d vobj_read_virusobj fn_clnp %s fn_nu %s fn_q %s\n', ...
                ii,cmd{ii}.fn_clnp{jj},cmd{ii}.fn_nu{jj},cmd{ii}.fn_q{jj});
      end
      vobj=virusobj_read(cmd{ii}.fn_clnp,cmd{ii}.fn_nu,cmd{ii}.fn_q);
      clear jj;

    case 'vobj_write_virusobj'
      %requires cmd{ii}.fn_clnp, cell array of Neta file names
      %requires cmd{ii}.fn_nu, cell array of Neta file names
      %requires cmd{ii}.fn_q, cell array of Neta file names
      for eta=1:length(vobj)
        fprintf(1,'hetero: ii %d vobj_write_virusobj fn_clnp %s fn_nu %s fn_q %s\n', ...
                ii,cmd{ii}.fn_clnp{eta},cmd{ii}.fn_nu{eta},cmd{ii}.fn_q{eta});
      end
      virusobj_write(cmd{ii}.fn_clnp,cmd{ii}.fn_nu,cmd{ii}.fn_q,vobj);
      clear eta;

    case 'vobj_save_virusobj'
      %requires cmd{ii}.fn_virusobj
      fprintf(1,'hetero: ii %d vobj_save_virusobj fn_virusobj %s\n',ii,cmd{ii}.fn_virusobj);
      save(cmd{ii}.fn_virusobj,'vobj');

    case 'vobj_load_virusobj'
      %requires cmd{ii}.fn_virusobj
      fprintf(1,'hetero: ii %d vobj_load_virusobj fn_virusobj %s\n',ii,cmd{ii}.fn_virusobj);
      load(cmd{ii}.fn_virusobj,'vobj');

    case 'vobj_change_size_of_virusobj'
      %requires cmd{ii}.vlmax(1:Neta)
      %requires cmd{ii}.vpmax(1:Neta)
      %Change the size of lmax and pmax in the virus model that will be used.  Does not set 2Dreciprocal.
      fprintf(1,'hetero: ii %d vobj_change_size_of_virusobj:\n',ii);
      fprintf(1,'hetero: ii %d vlmax:\n',ii);
      disp(cmd{ii}.vlmax);
      fprintf(1,'hetero: ii %d vpmax:\n',ii);
      disp(cmd{ii}.vpmax);
      vobj=virusobj_changesize(cmd{ii}.vlmax,cmd{ii}.vpmax,vobj);

    case 'vobj_change_homo2hetero_in_virusobj'
      %requires cmd{ii}.homo2hetero
      fprintf(1,'hetero: ii %d vobj_change_homo2hetero_in_virusobj:\n',ii);
      fprintf(1,'hetero: ii %d homo2hetero:\n',ii);
      for eta=1:length(cmd{ii}.homo2hetero)
        fprintf(1,'hetero: ii %d vobj_change_homo2hetero_in_virusobj eta %d action %d\n', ...
                ii,eta,cmd{ii}.homo2hetero{eta}.action);
      end
      vobj=virusobj_change_homo2hetero(vobj,cmd{ii}.homo2hetero);
      clear eta;

    case 'vobj_change_handedness'
      %requires cmd{ii}.changehand(1:Neta)
      for eta=1:length(vobj)
        if cmd{ii}.changehand(eta)
          tochange=find(mod(vobj{eta}.clnp.il,2)==1);
          vobj{eta}.clnp.c(tochange)=-vobj{eta}.clnp.c(tochange);
          vobj{eta}.cbar=vobj{eta}.clnp.c;
        end
      end
      clear eta tochange;

    case 'vobj_change_sign'
      %requires cmd{ii}.changesign(1:Neta)
      for eta=1:length(vobj)
        if cmd{ii}.changehand(eta)
          vobj{eta}.clnp.c=-vobj{eta}.clnp.c;
          vobj{eta}.cbar=vobj{eta}.clnp.c;
        end
      end
      clear eta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%integration rules:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'quad_read_integration_rule'
      %requires cmd{ii}.fn_rule file name
      fprintf(1,'hetero: ii %d quad_read_integration_rule fn_rule %s\n',ii,cmd{ii}.fn_rule);
      rule=rd_rule(cmd{ii}.fn_rule);

    case 'quad_compute_integration_rule'
      %requires cmd{ii}.whichrule
      fprintf(1,'hetero: ii %d quad_compute_integration_rule whichrule %d\n',ii,cmd{ii}.whichrule);
      switch cmd{ii}.whichrule

        case 1 %all Euler angles (except handedness uncertainty): alpha [0,2pi) uniform, beta [0,pi/2] Gauss-Legendre, gamma [0,2pi) uniform
          %requires cmd{ii}.Nabc(1:3)
          fprintf(1,'hetero: ii %d quad_compute_integration_rule Nabc %d %d %d\n',ii,cmd{ii}.Nabc(1),cmd{ii}.Nabc(2),cmd{ii}.Nabc(3));
          rule=quadrule4asym1(0,2*pi,cmd{ii}.Nabc(1),0,pi/2,cmd{ii}.Nabc(2),0,2*pi,cmd{ii}.Nabc(3));
          %test:
          fprintf(1,'hetero: ii %d quad_compute_integration_rule: sum of weights %g\n',ii,sum(rule(:,6)));

        otherwise
          error('hetero: ii %d quad_compute_integration_rule: unknown whichrule %d\n',ii,cmd{ii}.whichrule);

      end %switch cmd{ii}.whichrule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%forward operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'fw_mk_synthetic_2D_realspace_images'
      %requires cmd{ii}.Nv
      %requires cmd{ii}.NT
      %requires cmd{ii}.SNR
      fprintf(1,'hetero: ii %d fw_mk_synthetic_2D_realspace_images Nv %d NT %d SNR %g\n', ...
              ii,cmd{ii}.Nv,cmd{ii}.NT,cmd{ii}.SNR);
      [y,imagestack,truevalues,pixelnoisevar]=fw(cmd{ii}.SNR,vobj,vkminimalset,cmd{ii}.Nv,cmd{ii}.NT,Na,Nb,rule,all_tilde_b,ixminimalset,iyminimalset,deltachi); %y has no noise
      imageindex=[1:length(imagestack)]';

    case 'fw_mk_synthetic_GaussianPulse_2D_realspace_images'
      %requires cmd{ii}.Nv
      %requires cmd{ii}.NT
      %requires cmd{ii}.SNR
      %requires cmd{ii}.GaussianPulses
      fprintf(1,'hetero: ii %d fw_mk_synthetic_GaussianPulse_2D_realspace_images Nv %d NT %d SNR %g\n', ...
              ii,cmd{ii}.Nv,cmd{ii}.NT,cmd{ii}.SNR);
      fprintf(1,'hetero: ii %d fw_mk_synthetic_GaussianPulse_2D_realspace_images should print GaussianPulses\n',ii);
      [y,imagestack,truevalues,pixelnoisevar]=fw_GaussianPulses(cmd{ii}.SNR,cmd{ii}.GaussianPulses,deltachi,cmd{ii}.Nv,cmd{ii}.NT,Na,Nb,rule); %y has no noise
      imageindex=[1:length(imagestack)]';

    case 'fw_save_truevalues'
      %requires cmd{ii}.fn_truevalues
      fprintf(1,'hetero: ii %d fw_save_truevalues fn_truevalues %s\n',ii,cmd{ii}.fn_truevalues);
      save(cmd{ii}.fn_truevalues,'truevalues');

    case 'fw_write_truevalues'
      %requires cmd{ii}.fn_truevalues
      fprintf(1,'hetero: ii %d fw_write_truevalues fn_truevalues %s\n',ii,cmd{ii}.fn_truevalues);
      fid=fopen(cmd{ii}.fn_truevalues,'w');
      if fid==-1
        error('hetero: ii %d fw_write_truevalues fid %d == -1\n',ii,fid);
      end
      fprintf(fid,'%d %d\n',truevalues');
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d fw_write_truevalues status %d ~= 0\n',ii,status);
      end
      clear fid status;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%expectation-maximization operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'EM_read_tilde_b'
      %requires cmd{ii}.fn_tilde_b{.} file names, one file name for each type of angular basis function:
      %fn_tilde_b{1} = real-valued spherical harmonics
      %fn_tilde_b{2} = icosahedral harmonics (identity representation only)
      fprintf(1,'hetero: ii %d EM_read_tilde_b fn_tilde_b{1} real-valued spherical harmonics %s\n',ii,cmd{ii}.fn_tilde_b{1});
      fprintf(1,'hetero: ii %d EM_read_tilde_b fn_tilde_b{2} icosahedral harmonics (identity representation only) %s\n',ii,cmd{ii}.fn_tilde_b{2});
      all_tilde_b=cell(size(vobj));
      for eta=1:length(vobj)
        switch vobj{eta}.BasisFunctionType
          case 0 %real-valued spherical harmonics
            all_tilde_b{eta}.lmax=0;
            all_tilde_b{eta}.tilde_b=rd_b(cmd{ii}.fn_tilde_b{1},all_tilde_b{eta}.lmax);
          case 1 %icosahedral harmonics (identity representation only)
            all_tilde_b{eta}.lmax=max(vobj{eta}.clnp.il);
            all_tilde_b{eta}.tilde_b=rd_b(cmd{ii}.fn_tilde_b{2},all_tilde_b{eta}.lmax);
          otherwise
            error('hetero: EM_real_tilde_b: vobj{eta=%d}.BasisFunctionType=%d not known',eta,vobj{eta}.BasisFunctionType);
        end
      end
      clear eta;

    case 'EM_extractdatasubset'
      %requires cmd{ii}.kmax
      %Construct the reciprocal space image data structure for the range of reciprocal space $\|\vec\kappa\|<kmax$ that will be used.
      fprintf(1,'hetero: ii %d EM_extractdatasubset kmax %g\n',ii,cmd{ii}.kmax);
      [vk,y]=subset_reciprocalspace(cmd{ii}.kmax,recipobj.vkmag,recipobj.vkminimalset,Imagestack,recipobj.iiminimalset);
      fprintf(1,'hetero: ii %d EM_extractdatasubset Ny=size(vk,1)=%d\n',ii,size(vk,1));

    case 'EM_set_2Dreciprocal_in_virusobj'
      %requires cmd{ii}.use_vkminimalset_rather_than_vk
      fprintf(1,'hetero: ii %d EM_set_2Dreciprocal_in_virusobj use_vkminimalset_rather_than_vk %d\n', ...
              ii,cmd{ii}.use_vkminimalset_rather_than_vk);
      if cmd{ii}.use_vkminimalset_rather_than_vk
        vobj=virusobj_set_2Dreciprocal(vobj,vkminimalset);
      else
        vobj=virusobj_set_2Dreciprocal(vobj,vk);
      end

    case 'EM_rm_2Dreciprocal_in_virusobj'
      fprintf(1,'hetero: ii %d EM_rm_2Dreciprocal_in_virusobj\n',ii);
      for eta=1:length(vobj)
        vobj=rmfield(vobj,{'unique_lp','map_lp2unique','map_unique2lp','Htable'});
      end
      clear eta;

    case 'EM_sphericalsymmetry_homogeneous' %least squares
      fprintf(1,'hetero: ii %d EM_sphericalsymmetry_homogeneous\n',ii);
      %A spherically-symmetric object has a pure-real Fourier transform.  Therefore,
      %such a model can make only a 0 prediction of the imaginary components of the data.
      %Imaginary components of the data are removed from y and the corresponding rows of L are removed.
      meanofy=sum(y,2)/size(y,2); %compute meanofy before deleting imaginary components.
      meanofy=meanofy(1:2:end);
      Rabc=eye(3);
      for eta=1:length(vobj)
        L=setL_nord(Rabc, vobj{eta}.clnp.il, vobj{eta}.clnp.in, vk, vobj{eta}.Htable, vobj{eta}.map_unique2lp, all_tilde_b{eta}.tilde_b, vobj{eta}.Ftable, vobj{eta}.BasisFunctionType);
%        fprintf(1,'hetero: EM_sphericalsymmetry_homogeneous: L (odd rows are 0?):\n')
%        disp(L(1:min(end,10),:));
        L=L(1:2:end,:);
        vobj{eta}.cbar=pinv(L.' * L) * (L.' * meanofy);
        vobj{eta}.clnp.c=vobj{eta}.cbar;
%        disp(vobj{eta}.clnp.c);
      end
      clear meanofy Rabc eta L;

    case 'EM_expectationmaximization' %can do homogeneous or heterogeneous cases, can do various symmetries or no symmetry
      %requires cmd{ii}.MC_Nrandic
      %requires cmd{ii}.MC_FractionOfMeanForMinimum
      %requires cmd{ii}.MC_FractionOfMean
      %requires cmd{ii}.maxiter
      %requires cmd{ii}.maxiter4pixelnoisevarupdate
      %requires cmd{ii}.cbarftol
      %requires cmd{ii}.loglikeftol
      %requires cmd{ii}.nu_ic_FractionOfcbarForMinimum
      %requires cmd{ii}.nu_ic_FractionOfcbar
      %requires cmd{ii}.estimate_noise_var_in_homogeneous_problem
      %requires cmd{ii}.pixelnoisevar_initialcondition
      %requires cmd{ii}.fn_savehistory
      %requires cmd{ii}.verbosity
      %requires cmd{ii}.MultiplierForPixelnoisevarIC
      %requires cmd{ii}.MinimumClassProb
      %requires cmd{ii}.absthres
      %requires cmd{ii}.homo2hetero_1st_V_MaxIter
      %requires cmd{ii}.homo2hetero_estpixelnoisevar
      %requires cmd{ii}.homo2hetero_estV

      fprintf(1,'hetero: ii %d EM_expectationmaximization\n',ii);

      %Package the parameters related to the Monte Carlo choice of initial condition in a structure for simplicity.
      EM_MC.Nrandic=cmd{ii}.MC_Nrandic;
      EM_MC.FractionOfMeanForMinimum=cmd{ii}.MC_FractionOfMeanForMinimum;
      EM_MC.FractionOfMean=cmd{ii}.MC_FractionOfMean;

      %Package the parameters related to the Expectation-Maximization iterations in a structure for simplicity.
      EM_iter.maxiter=cmd{ii}.maxiter;
      EM_iter.maxiter4pixelnoisevarupdate=cmd{ii}.maxiter4pixelnoisevarupdate;
      EM_iter.cbarftol=cmd{ii}.cbarftol;
      EM_iter.cbarrtol=cmd{ii}.cbarrtol;
      EM_iter.cbarftol_dividebyNc=cmd{ii}.cbarftol_dividebyNc;
      EM_iter.loglikeftol=cmd{ii}.loglikeftol;
      EM_iter.loglikertol=cmd{ii}.loglikertol;
      EM_iter.nu_ic_FractionOfcbarForMinimum=cmd{ii}.nu_ic_FractionOfcbarForMinimum;
      EM_iter.nu_ic_FractionOfcbar=cmd{ii}.nu_ic_FractionOfcbar;
      EM_iter.estimate_noise_var_in_homogeneous_problem=cmd{ii}.estimate_noise_var_in_homogeneous_problem;
      EM_iter.nu_ic_always_proportional2cbar=cmd{ii}.nu_ic_always_proportional2cbar;
      EM_iter.rule=rule;
      EM_iter.Na=Na;
      EM_iter.V_TolX=cmd{ii}.V_TolX;
      EM_iter.V_MaxIter=cmd{ii}.V_MaxIter;
      EM_iter.fn_savehistory=cmd{ii}.fn_savehistory;
      EM_iter.verbosity=cmd{ii}.verbosity;
      EM_iter.MinimumClassProb=cmd{ii}.MinimumClassProb;
      EM_iter.homo2hetero_1st_V_MaxIter=cmd{ii}.homo2hetero_1st_V_MaxIter;
      EM_iter.homo2hetero_estpixelnoisevar=cmd{ii}.homo2hetero_estpixelnoisevar;
      EM_iter.homo2hetero_estV=cmd{ii}.homo2hetero_estV;

      switch cmd{ii}.pixelnoisevar_initialcondition

        case 'from_image_statistics'
          %For the following formula, please see test_noisevar.m.  The
          %fact that the reciprocal space image is complex and the
          %code treats the Re and Im parts as separate measurements
          %(independent and with equal variance) leads to the factor
          %of 0.5.  The user must have already set Na and Nb by using
          %one of 'basic_setsizeof2Drealspaceimages',
          %'basic_setsizeof2Drealspaceimagesfromimagestack', or
          %'basic_setsizeof2DreciprocalspaceImagesfromImagestack'.
          pixelnoisevar0=0.5*Na*Nb*annulussamplevariance;
          fprintf(1,'hetero: ii %d EM_expectationmaximization: pixelnoisevar0 %g annulussamplevariance %g Na %d Nb %d\n', ...
                  ii,pixelnoisevar0,annulussamplevariance,Na,Nb);

        case 'from_pixelnoisevar'
          pixelnoisevar0=pixelnoisevar;
          fprintf(1,'hetero: ii %d EM_expectationmaximization: pixelnoisevar0 %g pixelnoisevar %g\n', ...
                  ii,pixelnoisevar0,pixelnoisevar);

        otherwise
          error('hetero: ii %d EM_expectationmaximization: unknown pixelnoisevar_initialcondition %s\n', ...
                ii,cmd{ii}.pixelnoisevar_initialcondition);

      end
      fprintf(1,'hetero: ii %d EM_expectationmaximization: cmd{ii}.MultiplierForPixelnoisevarIC %g\n', ...
              ii,cmd{ii}.MultiplierForPixelnoisevarIC);
      pixelnoisevar0=pixelnoisevar0*cmd{ii}.MultiplierForPixelnoisevarIC;

      [vobj,pixelnoisevar,loglikelihood]=EM_expmax_MonteCarlo(vk,y,vobj,pixelnoisevar0,EM_MC,EM_iter,all_tilde_b,cmd{ii}.absthres);
      fprintf(1,'hetero: ii %d EM_expectationmaximization: pixelnoisevar %g\n',ii,pixelnoisevar);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%post-processing operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'post_compute_real_space_cubes'
      %requires cmd{ii}.whichclass
      %requires cmd{ii}.wantrhobar
      %requires cmd{ii}.wantrxx
      %requires cmd{ii}.mlow(1:3)
      %requires cmd{ii}.mhigh(1:3)
      %requires cmd{ii}.deltax(1:3)
      %requires cmd{ii}.EulerAngles(1:3)
      fprintf(1,'hetero: ii %d post_compute_real_space_cubes\n',ii);

      if cmd{ii}.whichclass<1 || length(vobj)<cmd{ii}.whichclass
        error('hetero: ii %d post_compute_real_space_cubes: whichclass %d\n',ii,cmd{ii}.whichclass);
      end
      [rhobar,rxx,x_rect]=plt_realspace(cmd{ii}.wantrhobar,cmd{ii}.wantrxx,vobj{cmd{ii}.whichclass}, ...
                                        all_tilde_b{cmd{ii}.whichclass}, ...
                                        cmd{ii}.mlow(1),cmd{ii}.mhigh(1),cmd{ii}.deltax(1), ...
                                        cmd{ii}.mlow(2),cmd{ii}.mhigh(2),cmd{ii}.deltax(2), ...
                                        cmd{ii}.mlow(3),cmd{ii}.mhigh(3),cmd{ii}.deltax(3), ...
                                        cmd{ii}.EulerAngles(1),cmd{ii}.EulerAngles(2),cmd{ii}.EulerAngles(3));
      %Package the parameters for simplicity.
      real_space_cubes.whichclass=cmd{ii}.whichclass;
      real_space_cubes.wantrhobar=cmd{ii}.wantrhobar;
      real_space_cubes.wantrxx=cmd{ii}.wantrxx;
      real_space_cubes.mlow=cmd{ii}.mlow;
      real_space_cubes.mhigh=cmd{ii}.mhigh;
      real_space_cubes.deltax=cmd{ii}.deltax;
      real_space_cubes.EulerAngles=cmd{ii}.EulerAngles;

    case 'post_save_real_space_cubes'
      %requires cmd{ii}.fn_real_space_cubes
      fprintf(1,'hetero: ii %d post_save_real_space_cubes fn_real_space_cubes %s\n', ...
              ii,cmd{ii}.fn_real_space_cubes);
      save(cmd{ii}.fn_real_space_cubes,'rhobar','rxx','x_rect','real_space_cubes');

    case 'post_crossections_real_space_cubes'
      %requires nothing
      fprintf(1,'hetero: ii %d post_crossections_real_space_cubes\n',ii);
      if ~isempty(rhobar)
        plt_realspace_crossections(rhobar,'rhobar');
      end
      if ~isempty(rxx)
        plt_realspace_crossections(rxx,'rxx');
      end

    case 'post_compute_FSC'
      %requires cmd{ii}.FSC_minmagk
      %requires cmd{ii}.FSC_maxmagk
      %requires cmd{ii}.FSC_deltamagk
      %requires cmd{ii}.FSC_eta4classA
      %requires cmd{ii}.FSC_eta4classB
      %requires cmd{ii}.FSC_is_same_vobj
      %requires cmd{ii}.FSC_scaling
      fprintf(1,'hetero: ii %d post_compute_FSC minmagk %g maxmagk %g deltamagk %g eta4classA %d eta4classB %d is_same_vobj %d scaling[%d,%d]\n', ...
              ii,cmd{ii}.FSC_minmagk,cmd{ii}.FSC_maxmagk,cmd{ii}.FSC_deltamagk, ...
              cmd{ii}.FSC_eta4classA,cmd{ii}.FSC_eta4classB,cmd{ii}.FSC_is_same_vobj,cmd{ii}.FSC_scaling(1),cmd{ii}.FSC_scaling(2));
      if cmd{ii}.FSC_is_same_vobj
        [magk,fsc]=get_FSC(vobj{cmd{ii}.FSC_eta4classA},vobj{cmd{ii}.FSC_eta4classB}, ...
                           cmd{ii}.FSC_minmagk,cmd{ii}.FSC_maxmagk,cmd{ii}.FSC_deltamagk,cmd{ii}.FSC_scaling);
      else
        fprintf(1,'hetero: ii %d post_compute_FSC: class A in vobj2 and class B in vobj\n',ii);
        [magk,fsc]=get_FSC(vobj2{cmd{ii}.FSC_eta4classA},vobj{cmd{ii}.FSC_eta4classB}, ...
                           cmd{ii}.FSC_minmagk,cmd{ii}.FSC_maxmagk,cmd{ii}.FSC_deltamagk,cmd{ii}.FSC_scaling);
      end
      

    case 'post_LPF_FSC'
      %requires cmd{ii}.LPFcutoff4FSC
      fprintf(1,'hetero: ii %d post_LPF_FSC LPFcutoff4FSC\n',cmd{ii}.LPFcutoff4FSC);
      tmp=fsc(:,4);
      N=length(tmp);
      M=round(N*cmd{ii}.LPFcutoff4FSC);
      Tmp=fft(tmp);
      Tmp(M:end-M+2)=0.0;
      tmp=ifft(Tmp);
      if ~isreal(tmp)
        fprintf(1,'hetero: ii %d post_LPF_FSC: isreal(tmp) is false\n',ii);
        fprintf(1,'hetero: ii %d post_LPF_FSC: DID NOT APPLY FILTER\n',ii);
      else
        fsc(:,4)=tmp;
      end
      clear tmp N M Tmp;
      
    case 'post_movingavg_FSC' 
      % moving average
      lenfsc=length(fsc(:,4));
      dfsc=zeros(1+2*(lenfsc-1));
      dmagk=dfsc;
      dfsc(lenfsc)=fsc(1,4);
      dmagk(lenfsc)=magk(1);
      for i=1:lenfsc-1
          dfsc(lenfsc+i)=fsc(i,4);
          dmagk(lenfsc+i)=magk(i);
          dfsc(lenfsc-i)=fsc(i,4);
          dmagk(lenfsc-i)=-magk(i);
      end
      for j=lenfsc:length(dfsc)
          newfsc(j-lenfsc+1)=(dfsc(j)+dfsc(j-1)+dfsc(j+1))/3;
      end
      fsc(:,4)=newfsc;
      
    case 'post_FSC_cutoff'
      %requires cmd{ii}.FSCcutoff
      fprintf(1,'hetero: ii %d post_FSC_cutoff FSCcutoff\n',cmd{ii}.FSCcutoff);
      mm=find(fsc(:,4)<cmd{ii}.FSCcutoff);
      if isempty(mm)
        fprintf(1,'hetero: ii %d post_FSC_cutoff isempty(mm) is true\n',ii);
        if exist('magkstar','var')==1 %clear incorrect value
          clear magkstar;
        end
      else
        mmstar=min(mm);
        magkstar=magk(mmstar); %keep this value for use in post_FSC_plot
        fprintf(1,'hetero: ii %d post_FSC_cutoff mmstar %d magkstar %g 1/magkstar %g\n', ...
                ii,mmstar,magkstar,1.0/magkstar);
        clear mmstar;
      end
      clear mm;

    case 'post_FSC_plot'
      %requires cmd{ii}.FSCtitle (title for the file to be printed)
      %requires cmd{ii}.energylegend1
      %requires cmd{ii}.energylegend2
      %requires cmd{ii].energytitle
      %requires cmd{ii}.plotcutoffk
      %requires cmd{ii}.plottitle (if this is not defined or if it is '' then nothing is written)
      figure;
      plot(magk,fsc(:,4),'-k');
      axis([0 0.18 0 1]);
      hold on
      plot(magk,0.5*ones(size(magk)),'--k');
      if cmd{ii}.plotcutoffk==true
          plot(magkstar,min(fsc(:,4)),'k*');
      end
      if length(cmd{ii}.plottitle)>0
          title(cmd{ii}.plottitle);
      end
      hold off      
      print('-depsc',cmd{ii}.FSCtitle);
      
      % Energy plots
      figure;
      p1=semilogy(magk,fsc(:,2),':k');
      hold on
      p2=semilogy(magk,fsc(:,3),'-k');
      hold off
      axis([0 0.18 0 10e+8]);
      l1=legend([p1,p2],cmd{ii}.energylegend1,cmd{ii}.energylegend2);
      print('-depsc',cmd{ii}.energytitle);
      
    case 'post_save_FSC'
      %requires cmd{ii}.fn_FSC
      fprintf(1,'hetero: ii %d post_save_FSC fn_FSC %s\n',ii,cmd{ii}.fn_FSC);
      save(cmd{ii}.fn_FSC,'magk','fsc');

    case 'post_write_FSC'
      %requires cmd{ii}.fn_FSC
      fprintf(1,'hetero: ii %d post_write_FSC fn_FSC %s\n',ii,cmd{ii}.fn_FSC);
      fid=fopen(cmd{ii}.fn_FSC,'w');
      if fid==-1
        error('hetero: ii %d post_write_FSC fid %d == -1\n',ii,fid);
      end
      fprintf(fid,'%g %g %g %g %g\n',[magk fsc]');
      status=fclose(fid);
      if status~=0
        error('hetero: ii %d post_write_FSC status %d ~= 0\n',ii,status);
      end
      clear fid status;

    case 'post_capsid4scale'
        %requires cmd{ii}.scaling_fn1
        %requires cmd{ii}.scaling_fn2
        %requires cmd{ii}.capsid_r
        %requires cmd{ii}.capsid_thres
        fprintf(1,'hetero: ii %d post_capsid4scale: Select capsid region based on two rhobar maps: scaling_fn1 %s and scaling_fn2 %s\n',...
            ii,cmd{ii}.scaling_fn1,cmd{ii}.scaling_fn2);
        file1=load(cmd{ii}.scaling_fn1);% require the input file to be the realspace.mat file generated by hetero system 
        file2=load(cmd{ii}.scaling_fn2);
        rhobar1=file1.rhobar;
        rhobar2=file2.rhobar;
        rxx2=file2.rxx;
        capsid_radius=cmd{ii}.capsid_r;
        capsid=select_capsid(rhobar1,rhobar2,capsid_radius,cmd{ii}.capsid_thres);

    case 'post_scaling'
        %requires cmd{ii}.scaling_type
        [scale_g,scale_a,scale_b]=scale_rhobar(rhobar1,rhobar2,capsid); %scalingfac=[g,a,b]
        fprintf(1,'hetero: ii %d post_scaling: Calculated scaling factors [g,a,b]=[%g,%g,%g]',ii,scale_g,scale_a,scale_b);
        if cmd{ii}.scaling_type==1 % 1 for a*rhobar2+b
            
            
            rhobarscaled=scale_a*rhobar2+scale_b*ones(size(rhobar2));
            rxxscaled=(scale_a^2)*rxx2;         
        elseif cmd{ii}.scaling_type==0 
            rhobarscaled=scale_g*rhobar2;
            rxxscaled=(scale_g^2)*rxx2;
        else
            error('Invalid scaling type %d',cmd{ii}.scaling_type);
        end

   case 'post_compute_spherical_average'
       %requires cmd{ii}.sphavg_eta2scale
       %requires cmd{ii}.sphavg_scale
       tempeta=cmd{ii}.sphavg_eta2scale;
       if isempty(cmd{ii}.sphavg_scale)
           scale4avg=scale_g;
       else
           scale4avg=cmd{ii}.sphavg_scale;
       end
        [vavg,unique_r]=spherical_avg_v(vobj{tempeta},all_tilde_b{tempeta}.tilde_b);
        stdv=scale4avg*sqrt(vavg);

   case 'post_save_spherical_average'
      %requires cmd{ii}.fn_spherical_avg
      fprintf(1,'hetero: ii %d post_save_spherical_average fn_spherical_avg %s\n',ii,cmd{ii}.fn_spherical_avg);
      save(cmd{ii}.fn_spherical_avg,'unique_r','stdv');

   case 'post_plot_spherical_averages'
       %requires cmd{ii}.fn_avg_plot1
       %requires cmd{ii}.fn_avg_plot2
       %requires cmd{ii}.sphavg_legend1
       %requires cmd{ii}.sphavg_legend2
       %requires cmd{ii}.sphavg_innercap
       %requires cmd{ii}.sphavg_title_full
       %requires cmd{ii}.sphavg_title_cap
       fprintf(1,'hetero: ii %d post_plot_spherical_averages fn_avg_plot1 %s fn_avg_plot2 %s\n',...
           ii,cmd{ii}.fn_avg_plot1,cmd{ii}.fn_avg_plot2);
       load(cmd{ii}.fn_avg_plot1);
       stdv1=stdv;
       unique_r1=unique_r;
       load(cmd{ii}.fn_avg_plot2);
       figure;
       p1=plot(unique_r1,stdv1,'-k');
       hold on
       p2=plot(unique_r,stdv,'--k');
       hold off
       l1=legend([p1,p2],cmd{ii}.sphavg_legend1,cmd{ii}.sphavg_legend2);
       xlim([0 capsid_radius(2)]);
       print('-depsc',cmd{ii}.sphavg_title_full);
       xlim([cmd{ii}.sphavg_innercap capsid_radius(2)]);
       print('-depsc',cmd{ii}.sphavg_title_cap);
       clear stdv1 stdv unique_r1 unique_r

    case 'post_hist_orientation_deviation'
        %requires cmd{ii}.fn_initial_orientation
        %requires cmd{ii}.fn_final_orientation
        %requires cmd{ii}.hist_title
        %requires cmd{ii}.metric_type
        fprintf(1,'hetero: ii %d post_hist_orientation_deviation fn_initial_orientation %s fn_final_orientation %s\n',...
           ii,cmd{ii}.fn_initial_orientation,cmd{ii}.fn_final_orientation);
        load(cmd{ii}.fn_initial_orientation);
        probi=p_theta_eta{1};
        for i=1:size(probi,1)
            [maxv peaki(i)]=max(probi(i,:));
        end
        load(cmd{ii}.fn_final_orientation);
        pf=p_theta_eta{1};
        for i=1:1200
            [maxv peakf(i)]=max(pf(i,:));
        end
        diffim=find(peaki-peakf);
        display(length(diffim));
        for i=1:length(diffim)
            abi=peaki(diffim(i));
            abf=peakf(diffim(i));
            euleri = rule(abi,1:3);
            eulerf = rule(abf,1:3);
            Ri=euler2R(euleri);
            Rf=euler2R(eulerf);
            metric(i)=rot_metric(Ri,Rf,cmd{ii}.metric_type);
        end
        hist(metric);
        set(get(gca,'child'),'FaceColor','k','EdgeColor','k');
        print('-depsc',cmd{ii}.hist_title);

    case 'post_plot_EMAN_FSC'
        %requires cmd{ii}.fn_EMAN_FSC
        %requires cmd{ii}.fn_EMAN_FSC_title
        fprintf(1,'hetero: ii %d post_plot_EMAN_FSC fn_EMAN_FSC %s fn_EMAN_FSC_title %s\n',...
           ii,cmd{ii}.fn_EMAN_FSC,cmd{ii}.fn_EMAN_FSC_title);
        [fid,fopenmsg]=fopen(cmd{ii}.fn_EMAN_FSC,'r');
        if fid == -1
            error(['rd_b: fid ' num2str(fid) ' fopenmsg ' fopenmsg ' opening ' ctf_fn]);
        end
        tmp=textscan(fid,'%f %f',inf);
        x=tmp{1};
        y=tmp{2};
        fclose(fid);
        
        figure;
        plot(x,y,'k');
        hold on
        plot(x,0.5*ones(size(x)),'--k');
        axis([0 0.18 0 1]);
        print('-depsc',cmd{ii}.fn_EMAN_FSC_title);
        clear tmp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test operators:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'test_imagestackfromsphericalobj'
      %requires cmd{ii}.numberimages
      %requires cmd{ii}.minxindex
      %requires cmd{ii}.maxxindex
      %requires cmd{ii}.minyindex
      %requires cmd{ii}.maxyindex
      %requires cmd{ii}.deltaxy
      %requires cmd{ii}.rho0
      %requires cmd{ii}.R2
      %requires cmd{ii}.pixelnoisevariance
      fprintf(1,'hetero: ii %d test_imagestackfromsphericalobj: numberimages %d minxindex %d maxxindex %d minyindex %d maxyindex %d deltaxy %g rho0 %g R2 %g pixelnoisevariance %g\n',ii,cmd{ii}.numberimages,cmd{ii}.minxindex,cmd{ii}.maxxindex,cmd{ii}.minyindex,cmd{ii}.maxyindex,cmd{ii}.deltaxy,cmd{ii}.rho0,cmd{ii}.R2,cmd{ii}.pixelnoisevariance);

%     could alternatively use test_sphericalobj1
      projimage=test_sphericalobj2([cmd{ii}.minxindex cmd{ii}.minyindex],[cmd{ii}.maxxindex cmd{ii}.maxyindex],cmd{ii}.deltaxy,cmd{ii}.rho0,cmd{ii}.R2);

      imagestack=cell(cmd{ii}.numberimages,1);
      for jj=1:length(imagestack)
        imagestack{jj}=projimage + sqrt(cmd{ii}.pixelnoisevariance).*randn(size(projimage));
      end

      clear projimage jj;

    case 'test_FFT_pt1'
      %requires cmd{ii}.N1
      %requires cmd{ii}.N2
      %requires cmd{ii}.delta1
      %requires cmd{ii}.delta2
      %requires cmd{ii}.freqintegers
      fprintf(1,'hetero: ii %d test_FFT_pt1: N1 %d N2 %d delta1 %g delta2 %g\n',ii,cmd{ii}.N1,cmd{ii}.N2,cmd{ii}.delta1,cmd{ii}.delta2);
      fprintf(1,'hetero: ii %d test_FFT_pt1: freqintegers:\n',ii);
      cmd{ii}.freqintegers
      imagestack=test_FFT_pt1(cmd{ii}.N1,cmd{ii}.N2,cmd{ii}.delta1,cmd{ii}.delta2,cmd{ii}.freqintegers);
      fprintf(1,'hetero: ii %d test_FFT_pt1: save cmd value\n',ii);
      cmd_test_FFT_pt1=cmd{ii}; %so test_FFT_pt2 has access to the cmd
      fprintf(1,'hetero: ii %d test_FFT_pt1: set deltachi value\n',ii);
      deltachi=[cmd{ii}.delta1 cmd{ii}.delta2];

    case 'test_FFT_pt2'
      %requires nothing but uses and then clears cmd_test_FFT_pt1 set by test_FFT_pt1
      fprintf(1,'hetero: ii %d test_FFT_pt2\n',ii);
      test_FFT_pt2(cmd_test_FFT_pt1.N1,cmd_test_FFT_pt1.N2,cmd_test_FFT_pt1.delta1,cmd_test_FFT_pt1.delta2,cmd_test_FFT_pt1.freqintegers,vk,y);
      clear cmd_test_FFT_pt1;

    otherwise
      error('hetero: ii %d unknown operator %s\n',ii,cmd{ii}.operator);

  end %close switch cmd{ii}.operator

end %close for ii=1:length(cmd)
