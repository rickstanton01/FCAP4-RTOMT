function RTOMT(sDirLogicleConvertedFCSFilesInTxtFormat, sFullPathFCAP4SortedMetaDataInExcelFormat)
% input args: sDirLogicleConvertedFCSFilesInTxtFormat = string, directory containing the 383 FCS data files from FCAP4 organizers which have been logicle converted into text format
%             sFullPathFCAP4MetaDataInExcelFormat     = string, full path and filename of FCAP4 metadata training data (MetaDataTrain.xlsx)
%                                                     = this file was provided by FCAP4 organizers except it was converted to excel and sorted 1st with respect to Status and next
%                                                       2nd with respect to Survival Time to create 191 sorted training rows followed by 192 test rows
% notes:      sFullPathFCAP4MetaDataInExcelFormat = contains the meta data for this code where
%                   columnA = Status: 1=progression to AIDS, 0=non-progresion,      columnB = Survival Time (as described in the manuscript)
%                   columnC = FCS logicle transformed text file name, stimulated    columnD = FCS file name, unstimulated
% github repository: FCAP4-RTOMT
  bBuildHistogram                     = 0; % Processing switch to save execution time once histogram is calc'd
  cMarkers                            = {'IFNg', 'TNFa', 'CD4', 'CD27',	'CD107A', 'CD154', 'CD3', 'CCR7', 'IL2', 'CD8', 'CD57', 'CD45RO', 'CD14'};
  [mTrnProg, text, alldata]           = xlsread(sFullPathFCAP4SortedMetaDataInExcelFormat, 1, 'A2:D384');      % mTrnProg=(status, survival), alldata=status, survival, stim, unstim info from organizers
  iNFilesNostTrn                      = 191;     
  iNFilesNostTst                      = 192;
  iNFilesTotal                        = 383;
  iNHistBins                          = 8192;
%%%%%%%%%%  form a histogram relative to all files and save the histogram to avoid having to create it every time you train
  if                                    (bBuildHistogram==1)
    mHistFCAP4                            = zeros(iNFilesTotal, iNHistBins);
    for                             i = 1:iNFilesTotal;
      sFN                             = [sDirLogicleConvertedFCSFilesInTxtFormat alldata{i,4}];
      mNormAll                        = dlmread([sFN(1:end-3) 'txt'], '\t', 1, 0);
      mNorm                           = mNormAll(:,4:16);  % ignore FSC_A, FSC_H and SSC columns
      vMeanNost                       = mean(mNorm); 
      iNCells                         = size(mNorm,1);    
      vMap                            = zeros(size(mNorm,1), 1);
      for                           j = 1:iNCells;  for  k=1:13;  if (mNorm(j,k)>= vMeanNost(k));  vMap(j)=(vMap(j)+ 2^(k-1));  end; end; end  % map into 2^13 boolean bins
      for                          kk = 1:iNHistBins
        vI                            = find(vMap==(kk-1));                  % NOTE!! mHistFCAP4 bin range 1:8192, where vMap bin range=0:8191 - need to subtract 1 upon interp signif bins
        mHistFCAP4(i,kk)              = ((length(vI)*10000) / iNCells ) ;    % scale by large num just to look at larger numbers for convenience
      end
      vClk                            = fix(clock);   fprintf('%d\t %d\t %d\t %d\t %s\t iNCells: %i\n', i, vClk(4), vClk(5), vClk(6), alldata{i,4}, iNCells);      
    end
    save                                ('mHistFCAP4.mat', 'mHistFCAP4');
  end
%%%%%%%%%%  train regression tree on histogram
  load                                  ('mHistFCAP4.mat');       % 'mHist');
  vbAIDS                              = zeros(iNFilesTotal, 1);
  vIAIDS                              = find(mTrnProg(:,1)==1);
  vTrnProgression                     = mTrnProg(:,2);
  vTrnTarget                          = mTrnProg(:,2);
  vTrnTarget(vIAIDS)                  = (vTrnTarget(vIAIDS)-3000)*2;
  rtreeNost5                          = RegressionTree.fit(mHistFCAP4(1:191,:),  vTrnTarget);
  vPredict5                           = round(predict(rtreeNost5,  mHistFCAP4));
  vIPPA                               = find(vPredict5<=0);    % Predicted Progression to AIDS
  vbAIDS(vIPPA)                       = 1;
  vPredict5Map                        = vPredict5;
  vPredict5Map(vIPPA)                 = round( (vPredict5Map(vIPPA)/2)+3000);
  subplot(2,1,1); plot                  (vTrnTarget(1:191),         '-k');  grid on; hold on; plot(vPredict5   (1:191), '-b');   title('black=trainTarget - AIDS scaled and mapped to negative, blue=mapped prediction'); plot([34.5 34.5], [6000 -6000], '-r'); 
  subplot(2,1,2); plot                  (vTrnProgression(1:191),    '-k');  grid on; hold on; plot(vPredict5Map(1:191), '-b');   title('black=timeToPrg,   blue=prediction'); plot([34.5 34.5], [6000 -6000], '-r'); 
  plot                                  ([34.5 34.5], [6000 -6000], '-r'); 
  mOutTrn                             = [vbAIDS vPredict5 vPredict5Map];
  xlswrite                              ('FC4_predict.xlsx', mOutTrn);
  view                                  (rtreeNost5);
  view                                  (rtreeNost5, 'mode', 'graph');
end 
