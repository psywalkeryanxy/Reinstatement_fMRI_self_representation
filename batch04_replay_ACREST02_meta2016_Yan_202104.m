%batch for replay analysis
%HR-self
%xinyuanyan@sanp bnu
%20190312


%If you want to use this batch created by xinyuanyan, please acknowledge in
%any format when you published your paper


clear all;
close;

spm('defaults', 'FMRI');

%% prepare the rest01
rootpath=pwd;
addpath([rootpath,filesep,'hamming_filter.m']);
datapath_rest01=[rootpath,filesep,'3_ANALYSIS_REST_02',filesep,'FunImgARCWF'];
datapath_task01=[rootpath,filesep,'2_ANALYSIS_AC_SELF'];
subfile=dir([datapath_task01,filesep,'2019*']);
maskpath=[rootpath,filesep,'rois_meta2016'];

for i=1:length(subfile)
%     %% trans 4D vol to 3D with time
%         needtransfile = dir([datapath_rest01,filesep,subfile(i).name,filesep,'Filtered_4DVolume.nii']);
%         matlabbatch{1}.spm.util.split.vol = {[datapath_rest01,filesep,subfile(i).name,filesep,needtransfile(1).name]};
%         matlabbatch{1}.spm.util.split.outdir = {[datapath_rest01,filesep,subfile(i).name]};
%         spm_jobman('serial', matlabbatch);
%         clear matlabbatch;
%     %% get ROI-based signal
%     
% %     
    %trans to *roi*.mat
    seedROI = dir([maskpath,filesep,'*.nii']);

    for numseed = 1:length(seedROI)
        %seedroi_name = seedROI(numseed).name;
        cd(maskpath);
        seedname = seedROI(numseed).name;
        o = maroi_image(struct('vol', spm_vol(seedname), 'binarize',0,...
            'func', 'img'));
        saveroi(o, [seedROI(numseed).name(1:length(seedname)-4),'_roi.mat']);
        clear o
        
    end%for numseed = 1:length(seedROI)
    
    clear numseed
%     
    seedROI = dir([maskpath,filesep,'*roi*.mat']);
    %% calculate the percent of signal within the ROI mask
    RESTFILES = spm_select('FPList',[datapath_rest01,filesep,subfile(i).name,filesep],'^Filtered_4DVolume_0*.*\.nii'); % get all beta-files
    
    fprintf('Number of rest-files (regressors): %d\n',size(RESTFILES,1));
    
    %now, xinyuan wants to get every nodes signal to do N*N
    %matrix and save the matrix for each window
    
    for seedmatrix = 1:length(seedROI)
        % cd(seedROIpath);
        thisseedroi = maroi([maskpath,filesep,seedROI(seedmatrix).name]);
        thisseedroi = spm_hold(thisseedroi,0);
        [restData, ~,~] = getdata(thisseedroi,RESTFILES,'l');
        
        for voxel=1:size(restData,2)
            tempTimeSeries=hamming_filter(restData(:,voxel));
            tempTimeSeries=((tempTimeSeries - nanmean(tempTimeSeries))./nanmean(tempTimeSeries));
            restData(:,voxel)=tempTimeSeries;
        end
        %clean voxel, mean��3SD
        for TRnumber = 2:length(restData(:,1))-1
            fileterIDX=find(abs(restData(TRnumber,:))>(mean(restData(TRnumber,:))+3*std(restData(TRnumber,:))));
            restData(TRnumber,fileterIDX)=nan;
            
        end%for TRnumber
        
        
        %% get ROI-based beta value for each voxel FROM TASK
        
        % path to beta-files
        betapath = [datapath_task01,filesep,subfile(i).name,filesep,'first_level_for_rsa_DNT256_ST'];
        
        
        % userOptions.smoothingKernel_fwhm = [5 5 5];
        
        
        trialfile = dir([betapath,filesep,'Trait*']);
        
        for whichtrial=1:length(trialfile)
            
            
            
            %define mask img for each sub
            betapath = [datapath_task01,filesep,subfile(i).name,filesep,'first_level_for_rsa_DNT256_ST'];
            
            BETAFILES = spm_select('FPList',[betapath,filesep,trialfile(whichtrial).name,filesep],'^beta_0001.*\.nii'); % get all beta-files
            
            fprintf('Number of beta-files (regressors): %d\n',size(BETAFILES,1));
            % get voxel timeseries within mask
            clear('y','vXYZ');
            
            
            [ybeta, ~, vXYZ] = getdata(thisseedroi,BETAFILES,'l');
            
            y_pool(whichtrial,:) = ybeta;
            
            
        end%for whichtrial
        
        
        %% get the replay matrix for this ROI
        %restData&y_pool
        for stimnum=1:length(trialfile)
            for whichTR=1:length(restData(:,1))
                
                corrpara=corr(restData(whichTR,:)',y_pool(stimnum,:)', 'rows','complete');
                replaymatrix(stimnum,whichTR)=corrpara;
            end
            
        end%for stimnum=1:length(trialfile)
        
        %get the 1 0 matrix
        
        allindx={};
        clean_replay = replaymatrix;
        for whichstim=1:length(replaymatrix(:,1))
            
            index=find(replaymatrix(whichstim,:)>(nanmean(replaymatrix(whichstim,:))+1.5*nanstd(replaymatrix(whichstim,:))));
            allindx{whichstim,1}=index;
            clean_replay(whichstim,index)=1;
            
        end
        clean_replay(find(clean_replay~=1))=0;
        allroi_replaynum{seedmatrix,2}=replaymatrix;
        allroi_replaynum{seedmatrix,1}=clean_replay;
        clear y_pool
    end% for seedmatrix = 1:length(seedROI)
    %plot for this subject
%     
%     for whichplot=1:length(seedROI)
%         roiname={'LAMG';'LHIP';'RAMG';'RHIP';'mPFC-3self'};
%         subplot(2,3,whichplot);
%         imagesc(allroi_replaynum{whichplot,1});
%         xlabel('TR number');
%         ylabel('Trait number');
%         title(['replay matrix in ', roiname{whichplot,1}]);
%         colormap(hot);
%     end %for whichplot=1:length(seedROI)
    %how many replay number for this subject for each trait word?

    for whichseed = 4%:length(seedROI)
        
    for nnword=1:length(allindx)
        numbermatrix=allroi_replaynum{whichseed,1};
        allword_replaynum(i,nnword)=sum(numbermatrix(nnword,:));

    end
      %allseed_allword_replaynum{whichseed,1}=allword_replaynum;
      
      filename = ['results_allword_replaynum_',seedROI(whichseed).name];%index to ROI
      save(filename,'allword_replaynum');
      
      
    end
    
end%for i

  