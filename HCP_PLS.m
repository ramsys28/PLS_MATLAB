function result=HCP_PLS(fmri,behav)

result = [];
num_subj_lst=473;
num_cond = 1;
method = 3;
num_perm = 10;
is_struct = 0;
num_split = 0;
num_boot = 10;
clim = 95;
bscan = 1:num_cond;
beha=csvread('/scr/namibia1/Projects/Personality/HCP_Data/behav_NEO_HCP.csv');
beha(end,:)=[]; % 474 subject out!!!!
stacked_behavdata = behav;
meancentering_type = 0;
cormode = 0;
boot_type = 'strat';

num_groups = 1;
total_rows=473;

result.method = method;
result.is_struct = is_struct;
a.Properties.Writable=true;
behavdata_reorder = [1:size(stacked_behavdata,1)]';
datamat_reorder=behavdata_reorder;

manos=fmri;
[datamatsvd, datamatsvd_unnorm, datamatcorrs_lst, stacked_smeanmat] = ...
rri_get_covcor(method, manos, stacked_behavdata, num_groups, ...
num_subj_lst, num_cond, bscan, meancentering_type, cormode, ...
[], 1, num_boot, datamat_reorder, behavdata_reorder);

result.datamatcorrs_lst = datamatcorrs_lst;

[r c] = size(datamatsvd);
  if r <= c
     [u,s,v] = svd(datamatsvd',0);
  else
     [v,s,u] = svd(datamatsvd,0);
  end
s = diag(s);
org_s = s;
org_v = v;
original_u = u * diag(s);
original_v = v * diag(s);
kk = length(bscan);
result.u = u;
result.s = s;
result.v = v;
vsc = [];
row_idx = [];
last = 0;
k=1;% num of groups
for g = 1:num_groups
 if ~iscell(num_subj_lst)
    n = num_subj_lst(g);

    %  take this advantage (having g & n) to get row_idx
    %
    tmp = 1:n*k;
    tmp = reshape(tmp, [n k]);
    tmp = tmp(:, bscan);
    row_idx = [row_idx ; tmp(:) + last];
    last = last + n*k;
 else
    n = num_subj_lst{g};

    %  get row_idx
    %
    for k1 = 1:k
       if ismember(k1, bscan)
          tmp = [1:n(k1)]';
          row_idx = [row_idx ; tmp + last];
       end

       last = last + n(k1);
    end
 end
end

if ~iscell(num_subj_lst)
 [usc, vsc, lvcorrs] = rri_get_behavscores(manos, ...
stacked_behavdata, u, v, k, num_subj_lst, cormode);
else
 [usc, vsc, lvcorrs] = ssb_rri_get_behavscores(stacked_datamat, ...
stacked_behavdata, u, v, k, num_subj_lst, cormode);
end

result.lvcorrs = lvcorrs;
result.usc = usc;
result.vsc = vsc;
result.stacked_behavdata = stacked_behavdata;
result.num_subj_lst = num_subj_lst;
result.num_conditions = k;
is_perm_splithalf=0;
if num_perm > 0
    if ~is_perm_splithalf		%%%%%% REGULAR PERMUTATION TEST %%%%%%
        sp = zeros(size(s));
        vp = zeros(size(v));
        if ismember(method, [3 5])
            %PARALLEL
             for p = 1:num_perm
                reorder(:,p) = [randperm(size(manos,1))'];
             end
        
        end
        pcntacc = fprintf('Working on %d permutations:', num_perm);
        % PARALLEL
        for p = 1:num_perm
            pcntacc = pcntacc + fprintf(' %d', p);
            datamat_reorder = [1:size(manos,1)]';
            behavdata_reorder = reorder(:,p);   
             
            if ismember(method, [3 4 5 6])
                behav_p = stacked_behavdata(behavdata_reorder,:);

                for g = 1:num_groups
                   if ~iscell(num_subj_lst)
                      n = num_subj_lst(g);
                      span = sum(num_subj_lst(1:g-1)) * k;

                      min1 = min(std(behav_p(1+span:n*k+span,:)));
                      count = 0;

                      while (min1 == 0)
                         if ismember(method, [4 6])
                            reorder(:,p) = [rri_randperm_notall(num_subj_lst, k, bscan)'];
                            behavdata_reorder = reorder(:,p);
                         elseif ismember(method, [3 5])
                            reorder(:,p) = [randperm(size(stacked_datamat,1))'];
                            behavdata_reorder = reorder(:,p);
                         end

                         behav_p = stacked_behavdata(behavdata_reorder,:);
                         min1 = min(std(behav_p(1+span:n*k+span,:)));
                         count = count + 1;

                         if count > 100
                            msg = 'Please check your behavior data, and make ';
                            msg = [msg 'sure none of the columns are all the '];
                            msg = [msg 'same for each group'];

                            disp(' '); disp(msg);
                            return;
                         end	% if count
                      end		% while min1

                     else
                      n = num_subj_lst{g};
                      span = sum([num_subj_lst{1:g-1}]);

                      min1 = min(std(behav_p(1+span:sum(n)+span,:)));
                      count = 0;

                      while (min1 == 0)
                          reorder(:,p) = [randperm(size(stacked_datamat,1))'];
                          behavdata_reorder = reorder(:,p);
                         behav_p = stacked_behavdata(behavdata_reorder,:);
                         min1 = min(std(behav_p(1+span:sum(n)+span,:)));
                         count = count + 1;

                         if count > 100
                            msg = 'Please check your behavior data, and make ';
                            msg = [msg 'sure none of the columns are all the '];
                            msg = [msg 'same for each group'];

                            disp(' '); disp(msg);
                            return;
                         end	% if count
                      end		% while min1
                   end	% if ~iscell(num_subj_lst)
                end		% for g
            end		% if ismember

     [datamatsvd, datamatsvd_unnorm] = rri_get_covcor(method, ...
    manos, stacked_behavdata, num_groups, ...
    num_subj_lst, num_cond, bscan, meancentering_type, ...
    cormode, [], 0, 0, datamat_reorder, ...
    behavdata_reorder);
        [r c] = size(datamatsvd);
        if r <= c
           [pu, sperm, pv] = svd(datamatsvd',0);
        else
           [pv, sperm, pu] = svd(datamatsvd,0);
        end

        %  rotate pv to align with the original v
        %
        rotatemat = rri_bootprocrust(v, pv);

        %  rescale the vectors
        %
        pv = pv * sperm * rotatemat;

        sperm = sqrt(sum(pv.^2));
           ptotal_s = sum(datamatsvd_unnorm(:).^2);
           per = diag(sperm).^2 / sum(diag(sperm).^2);
           sperm = sqrt(per * ptotal_s);
           pv = normalize(pv) * diag(sperm);
           sp = sp + (sperm>=org_s);
           vp = vp + (abs(pv) >= abs(org_v));
        if pcntacc > 70
           fprintf('\n');
           pcntacc = 0;
        end
     end
  end		% for num_perm


  %  Save perm_result
  %
  result.perm_result.num_perm = num_perm;
  result.perm_result.sp = sp;
  result.perm_result.sprob = sp ./ (num_perm + 1);
  result.perm_result.permsamp = reorder;

end	% if perm

if num_boot > 0
      countnewtotal=0;
      num_LowVariability_behav_boots = [];
      badbeh = [];
      incl_seq = 0;
      [reorder, new_num_boot] = rri_boot_order(num_subj_lst, k, ...
		num_boot, [1:k], incl_seq, boot_type);
      reorder_4beh = reorder;

         orig_corr = lvcorrs;

         [r1 c1] = size(orig_corr);
         distrib = zeros(r1, c1, num_boot+1);
         distrib(:, :, 1) = orig_corr;

      max_subj_per_group = 8;

    
         if (sum(num_subj_lst <= max_subj_per_group) == num_groups)
            is_boot_samples = 1;
         else
            is_boot_samples = 0;
         end

         u_sum = original_u;

      u_sq = u_sum.^2;


         %  Check min% unique values for all behavior variables
         %
         num_LowVariability_behav_boots = zeros(1, size(stacked_behavdata, 2));

         for bw = 1:size(stacked_behavdata, 2)
            for p = 1:num_boot
               vv = stacked_behavdata(reorder_4beh(:,p),bw);

               if rri_islowvariability(vv, stacked_behavdata(:,bw))
                  num_LowVariability_behav_boots(bw) = num_LowVariability_behav_boots(bw) + 1;
               end
            end
         end

         if any(num_LowVariability_behav_boots)
            disp(' ');
            disp(' ');
            disp('For at least one behavior measure, the minimum unique values of resampled behavior data does not exceed 50% of its total.');
            disp(' ');
            disp(' ');
         end
         pcntacc = fprintf('Working on %d bootstraps:', num_boot);
         

         % PARALLEL
      for p=1:num_boot

  
            pcntacc = pcntacc + fprintf(' %d', p);
         

         datamat_reorder = reorder(:,p);
            datamat_reorder_4beh = reorder_4beh(:,p);
            behavdata_reorder = reorder_4beh(:,p);

         step = 0;

         for g = 1:num_groups

            %  check reorder for badbehav
            %
            if ismember(method,[3 4 5 6])
               if ~iscell(num_subj_lst)

                   n = num_subj_lst(g);
                   span = sum(num_subj_lst(1:g-1)) * k;	% group length

                   if ~is_boot_samples

                      % the code below is mainly trying to find a proper
                      % reorder matrix

                      % init badbehav array to 0
                      % which is used to record the bad behav data caused by
                      % bad re-order. This variable is for disp only.
                      %
                      badbehav = zeros(k, size(stacked_behavdata,2));

                      % Check for upcoming NaN and re-sample if necessary.
                      % this only happened on behavior analysis, because the
                      % 'xcor' inside of 'corr_maps' contains a 'stdev', which
                      % is a divident. If it is 0, it will cause divided by 0
                      % problem.
                      % since this happend very rarely, so the speed will not
                      % be affected that much.
                      %
                      % For behavpls_boot, also need to account for multiple
                      % scans and behavs
                      %
                      for c=1:k	% traverse all conditions in this group
                          stdmat(c,:) = std(stacked_behavdata(reorder_4beh((1+ ...
				(n*(c-1))+span):(n*c+span), p), :));
                      end		% scanloop

                      % now, check to see if any are zero
                      %
                      while sum(stdmat(:)==0)>0
                          countnewtotal = countnewtotal + 1;

                          % keep track of scan & behav that force a resample
                          %
                          badbehav(find(stdmat(:)==0)) = ...
				badbehav(find(stdmat(:)==0)) + 1;

                          badbeh{g,countnewtotal} = badbehav;	% save instead of disp

                          % num_boot is just something to be picked to prevent
                          % infinite loop
                          %
                          if countnewtotal > num_boot
                              disp('Please check behavior data');
                              breakon=1;
                              break;
                          end

                          reorder_4beh(:,p) = rri_boot_order(num_subj_lst, k, ...
				1, [1:k], incl_seq, boot_type);

                          if ismember(method, [4 6])
                             reorder_4beh(:,p) = rri_boot_order(num_subj_lst, k, ...
				1, bscan, incl_seq, boot_type);
                          end

                          for c=1:k	% recalc stdmat
                              stdmat(c,:) = std(stacked_behavdata(reorder_4beh((1+ ...
                                 (n*(c-1))+span):(n*c+span), p), :));
                          end	% scanloop
                      end	% while

                      if countnewtotal>0 & exist('bootsamp_4beh','var')
                         disp(' '); disp('bootsamp_4beh changed'); disp(' ');
                      end

                      % now, we can use this proper reorder matrix to generate 
                      % datamat_reorder_4beh & behavdata_reorder, and then to
                      % calculate datamatcorrs
                      %
                      datamat_reorder_4beh = reorder_4beh(:,p);
                      behavdata_reorder = reorder_4beh(:,p);

                   end	% if ~is_boot_samples
                                       

                   

               end	% if ~iscell(num_subj_lst)
            end		% if ismember(method,[3 4 5 6])

            [datamatsvd, datamatsvd_unnorm, datamatcorrs_lst, ...
		stacked_smeanmat] = rri_get_covcor(method, ...
		manos, stacked_behavdata, num_groups, ...
		num_subj_lst, num_cond, bscan, meancentering_type, ...
		cormode, [], 1, num_boot, datamat_reorder, ...
		behavdata_reorder, datamat_reorder_4beh);

         end		% for num_groups



            %  Singular Value Decomposition
            %
            [r c] = size(datamatsvd);
            if r <= c
               [pu, sboot, pv] = svd(datamatsvd',0);
            else
               [pv, sboot, pu] = svd(datamatsvd,0);
            end

            %  rotate pv to align with the original v
            %
            rotatemat = rri_bootprocrust(v, pv);

            %  rescale the vectors
            %
            pu = pu * sboot * rotatemat;
            pv = pv * sboot * rotatemat;

               data_p = manos(datamat_reorder_4beh(row_idx),:);
               behav_p = stacked_behavdata(behavdata_reorder(row_idx),:);

               if ~iscell(num_subj_lst)
                  [brainsctmp, behavsctmp, bcorr] = ...
			rri_get_behavscores(data_p, behav_p, ...
			normalize(pu), normalize(pv), kk, num_subj_lst, cormode);
               
               end

               distrib(:, :, p+1) = bcorr;

    
               fprintf('\n');
               pcntacc = 0;


      end		% for num_boot

      fprintf('\n');

      result.boot_result.num_boot = num_boot;
      result.boot_result.clim = clim;
      result.boot_result.num_LowVariability_behav_boots = num_LowVariability_behav_boots;

      result.boot_result.boot_type = boot_type;
      

    
      %  now compare the original unstandarized saliences
      %  with the bootstrap saliences
      %
      ul=clim;
      ll=100-clim;

      % e.g. 0.05 >> 0.025 for upper & lower tails, two-tailed
      %
      climNi = 0.5*(1-(clim*0.01));
         [llcorr, ulcorr, prop, llcorr_adj, ulcorr_adj] = ...
            rri_distrib(distrib, ll, ul, num_boot, climNi, orig_corr);


  

         result.boot_result.orig_corr = orig_corr;
         result.boot_result.ulcorr = ulcorr;
         result.boot_result.llcorr = llcorr;
         result.boot_result.ulcorr_adj = ulcorr_adj;
         result.boot_result.llcorr_adj = llcorr_adj;
         result.boot_result.badbeh = badbeh;
         result.boot_result.countnewtotal = countnewtotal;
         result.boot_result.bootsamp_4beh = reorder_4beh;

         
    

      result.boot_result.prop = prop;
      result.boot_result.distrib = distrib;
      result.boot_result.bootsamp = reorder;

      % result.boot_result.u_sum2 = u_sum2;
      % result.boot_result.v_sum2 = v_sum2;



%      test_zeros_v=find(v_se<=0);
 %     if ~isempty(test_zeros_v);
  %       v_se(test_zeros_v)=1;
   %   end

      
 
         compare_u = original_u ./ u_sq;
%         compare_v = original_v ./ v_se;


      %  for zero standard errors - replace bootstrap ratios with zero
      %  since the ratio makes no sense anyway
      %

%      if ~isempty(test_zeros_v);
 %        compare_v(test_zeros_v)=0;
  %    end

      result.boot_result.compare_u = compare_u;
%      result.boot_result.compare_v = compare_v;
      result.boot_result.u_se = u_sq;
%      result.boot_result.v_se = v_se;
  

   end	% if boot
   result.other_input.meancentering_type = meancentering_type;
   result.other_input.cormode = cormode;
   result.field_descrip = [
'Result fields description:                                           '
'                                                                     '
'method:                 PLS option                                   '
'                        1. Mean-Centering Task PLS                   '
'                        2. Non-Rotated Task PLS                      '
'                        3. Regular Behavior PLS                      '
'                        4. Multiblock PLS                            '
'                        5. Non-Rotated Behavior PLS                  '
'                        6. Non-Rotated Multiblock PLS                '
'                                                                     '
'u:                      Brainlv or Salience                          '
'                                                                     '
's:                      Singular value                               '
'                                                                     '
'v:                      Designlv or Behavlv                          '
'                                                                     '
'usc:                    Brainscores or Scalpscores                   '
'                                                                     '
'vsc:                    Designscores or Behavscores                  '
'                                                                     '
'TBv:                    Store Task / Bahavior v separately           '
'                                                                     '
'TBusc:                  Store Task / Bahavior usc separately         '
'                                                                     '
'TBvsc:                  Store Task / Bahavior vsc separately         '
'                                                                     '
'datamatcorrs_lst:       Correlation of behavior data with datamat.   '
'                        Only available in behavior PLS.              '
'                                                                     '
'lvcorrs:                Correlation of behavior data with usc,       '
'                        only available in behavior PLS.              '
'                                                                     '
'perm_result:            struct containing permutation result         '
'        num_perm:       number of permutation                        '
'        sp:             permuted singular value greater than observed'
'        sprob:          sp normalized by num_perm                    '
'        permsamp:       permutation reorder sample                   '
'        Tpermsamp:      permutation reorder sample for multiblock PLS'
'        Bpermsamp:      permutation reorder sample for multiblock PLS'
'                                                                     '
'perm_splithalf:         struct containing permutation splithalf      '
'        num_outer_perm: permutation splithalf related                '
'        num_split:      permutation splithalf related                '
'        orig_ucorr:     permutation splithalf related                '
'        orig_vcorr:     permutation splithalf related                '
'        ucorr_prob:     permutation splithalf related                '
'        vcorr_prob      permutation splithalf related                '
'        ucorr_ul:       permutation splithalf related                '
'        ucorr_ll:       permutation splithalf related                '
'        vcorr_ul:       permutation splithalf related                '
'        vcorr_ll:       permutation splithalf related                '
'                                                                     '
'boot_result:            struct containing bootstrap result           '
'        num_boot:       number of bootstrap                          '
'        boot_type:      Set to ''nonstrat'' if using Natasha''s         '
'                        ''nonstrat'' bootstrap type; set to            '
'                        ''strat'' for conventional bootstrap.          '
'        nonrotated_boot: Set to 1 if using Natasha''s Non             '
'                        Rotated bootstrap; set to 0 for              '
'                        conventional bootstrap.                      '
'        bootsamp:       bootstrap reorder sample                     '
'        bootsamp_4beh:  bootstrap reorder sample for behav PLS       '
'        compare_u:      compared salience or compared brain          '
'        u_se:           standard error of salience or brainlv        '
'        clim:           confidence level between 0 and 100.          '
'        distrib:        orig_usc or orig_corr distribution           '
'        prop:           orig_usc or orig_corr probability            '
'                                                                     '
'        following boot_result only available in task PLS:            '
'                                                                     '
'        usc2:           brain scores that are obtained from the      '
'                        mean-centered datamat                        '
'        orig_usc:       same as usc, with mean-centering on subj     '
'        ulusc:          upper boundary of orig_usc                   '
'        llusc:          lower boundary of orig_usc                   '
'        ulusc_adj:      percentile of orig_usc distribution with     '
'                        upper boundary of orig_usc                   '
'        llusc_adj:      percentile of orig_usc distribution with     '
'                        lower boundary of orig_usc                   '
'                                                                     '
'        following boot_result only available in behavior PLS:        '
'                                                                     '
'        orig_corr:      same as lvcorrs                              '
'        ulcorr:         upper boundary of orig_corr                  '
'        llcorr:         lower boundary of orig_corr                  '
'        ulcorr_adj:     percentile of orig_corr distribution with    '
'                        upper boundary of orig_corr                  '
'        llcorr_adj:     percentile of orig_corr distribution with    '
'                        lower boundary of orig_corr                  '
'        num_LowVariability_behav_boots: display numbers of low       '
'                        variability resampled hehavior data in       '
'                        bootstrap test                               '
'        badbeh:         display bad behav data that is caused by     '
'                        bad re-order (with 0 standard deviation)     '
'                        which will further cause divided by 0        '
'        countnewtotal:  count the new sample that is re-ordered      '
'                        for badbeh                                   '
'                                                                     '
'is_struct:              Set to 1 if running Non-Behavior             '
'                        Structure PLS; set to 0 for other PLS.       '
'                                                                     '
'bscan:                  Subset of conditions that are selected       '
'                        for behav block, only in Multiblock PLS.     '
'                                                                     '
'num_subj_lst            Number of subject list array, containing     '
'                        the number of subjects in each group.        '
'                                                                     '
'num_cond                Number of conditions in datamat_lst.         '
'                                                                     '
'stacked_designdata:     Stacked design contrast data for all         '
'                        the groups.                                  '
'                                                                     '
'stacked_behavdata:      Stacked behavior data for all the groups.    '
'                                                                     '
'other_input:            struct containing other input data           '
'        meancentering_type: Use Natasha''s meancentering type         '
'                        if it is not 0.                              '
'        cormode:        Use Natasha''s correlation mode if it         '
'                        is not 0.                                    '
'                                                                     '
];




