function result=for_harvard(conns,behav)
num_subj_lst=max(size(behav));
k=1;
num_cond = k;
method = 3;
num_perm = 10;
is_struct = 0;
num_split = 0;
num_boot = 10;
clim = 95;
bscan = 1:num_cond;
stacked_datamat = conns;
stacked_behavdata = behav;
is_perm_splithalf = 0;
meancentering_type = 0;
cormode = 0;
boot_type = 'strat';
nonrotated_boot = 0;
total_rows = max(size(behav));
num_groups=1;
g=1;


v7 = version;
if str2num(v7(1))<7
  singleanalysis = 0;
else
  singleanalysis = 1;
end   

if str2num(v7(1:3))<7.4 & strcmp(v7(4),'.')
pc = computer;
    if singleanalysis & ( strcmp(pc,'GLNXA64') | strcmp(pc,'GLNXI64') | strcmp(pc,'PCWIN64') )
         singleanalysis = 0;
    end
end


if singleanalysis
 stacked_datamat = single(stacked_datamat);
else
 stacked_datamat = double(stacked_datamat);
end

if singleanalysis
    stacked_behavdata = single(stacked_behavdata);
else
    stacked_behavdata = double(stacked_behavdata);
end

n = num_subj_lst(g);
span = sum(num_subj_lst(1:g-1)) * k;
min1 = min(std(stacked_behavdata(1+span:n*k+span,:)));

if min1 == 0
   msg = 'Please check your behavior data, and make sure that';
   msg = [msg char(10) 'none of the columns are all the same for each group.'];
   error(msg);
end

result.method = method;
result.is_struct = is_struct;
datamat_reorder = [1:size(stacked_datamat,1)]';
behavdata_reorder = [1:size(stacked_behavdata,1)]';

[datamatsvd, datamatsvd_unnorm, datamatcorrs_lst, stacked_smeanmat] = ...
	rri_get_covcor(method, stacked_datamat, stacked_behavdata, num_groups, ...
	num_subj_lst, num_cond, bscan, meancentering_type, cormode, ...
	1, 1, num_boot, datamat_reorder, behavdata_reorder);

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
s = diag(s);
result.u = u;
result.s = s;
result.v = v;
vsc = [];
row_idx = [];
last = 0;

n = num_subj_lst(g);
tmp = 1:n*k;
tmp = reshape(tmp, [n k]);
tmp = tmp(:, bscan);
row_idx = [row_idx ; tmp(:) + last];
last = last + n*k;
[usc, vsc, lvcorrs] = rri_get_behavscores(stacked_datamat, ...
		stacked_behavdata, u, v, k, num_subj_lst, cormode);         
result.lvcorrs = lvcorrs;
result.usc = usc;
result.vsc = vsc;
result.stacked_behavdata = stacked_behavdata;
result.num_subj_lst = num_subj_lst;
result.num_conditions = k;

%PERMUTATIONS

sp = zeros(size(s));
vp = zeros(size(v));

for p = 1:num_perm
    reorder(:,p) = [randperm(size(stacked_datamat,1))'];
end

for p = 1:num_perm


datamat_reorder = [1:size(stacked_datamat,1)]';
behavdata_reorder = reorder(:,p);
behav_p = stacked_behavdata(behavdata_reorder,:);
n = num_subj_lst(g);
span = sum(num_subj_lst(1:g-1)) * k;

min1 = min(std(behav_p(1+span:n*k+span,:)));
count = 0;

while (min1 == 0)
    reorder(:,p) = [randperm(size(stacked_datamat,1))'];
    behavdata_reorder = reorder(:,p);
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

if count>0 & exist('permsamp','var')
    disp(' '); disp('permsamp changed'); disp(' ');
end

[datamatsvd, datamatsvd_unnorm] = rri_get_covcor(method, ...
stacked_datamat, stacked_behavdata, num_groups, ...
num_subj_lst, num_cond, bscan, meancentering_type, ...
cormode, 1, 0, 0, datamat_reorder, ...
behavdata_reorder);

[r,c] = size(datamatsvd);
if r <= c
   [pu, sperm, pv] = svd(datamatsvd',0);
else
   [pv, sperm, pu] = svd(datamatsvd,0);
end
rotatemat = rri_bootprocrust(v, pv);
pv = pv * sperm * rotatemat;
sperm = sqrt(sum(pv.^2));         
sp = sp + (diag(sperm)>=s);           
end

result.perm_result.num_perm = num_perm;
result.perm_result.sp = sp;
result.perm_result.sprob = sp ./ (num_perm + 1);
result.perm_result.permsamp = reorder;

% BOOTSTRAPS

countnewtotal=0;
num_LowVariability_behav_boots = [];
badbeh = [];

if method == 1 | nonrotated_boot
    incl_seq = 1;
else
    incl_seq = 0;
end
[reorder, new_num_boot] = rri_boot_order(num_subj_lst, k, ...
		num_boot, [1:k], incl_seq, boot_type);
reorder_4beh = reorder;
orig_corr = lvcorrs;

[r1 c1] = size(orig_corr);
distrib = zeros(r1, c1, num_boot+1);
distrib(:, :, 1) = orig_corr;
max_subj_per_group = 8;
is_boot_samples = 0;
u_sum = original_u;
u_sq = u_sum.^2;
num_LowVariability_behav_boots = zeros(1, size(stacked_behavdata, 2));

for bw = 1:size(stacked_behavdata, 2)
    for p = 1:num_boot
       vv = stacked_behavdata(reorder_4beh(:,p),bw);

       if rri_islowvariability(vv, stacked_behavdata(:,bw))
          num_LowVariability_behav_boots(bw) = num_LowVariability_behav_boots(bw) + 1;
       end
    end
end
single_cond_lst=1;
kk=1;
for p=1:10

    datamat_reorder = reorder(:,p);
    datamat_reorder_4beh = reorder_4beh(:,p);
    behavdata_reorder = reorder_4beh(:,p);
    step = 0;
    n = num_subj_lst(g);
    span = sum(num_subj_lst(1:g-1)) * k;	% group length
    g=1;
    n = num_subj_lst(g);
    span = sum(num_subj_lst(1:g-1)) * k;
    
    if ~is_boot_samples
        badbehav = zeros(k, size(stacked_behavdata,2));
        for c=1:k	% traverse all conditions in this group
            stdmat(c,:) = std(stacked_behavdata(reorder_4beh((1+ ...
				(n*(c-1))+span):(n*c+span), p), :));
        end		% scanloop
        while sum(stdmat(:)==0)>0
            countnewtotal = countnewtotal + 1;
            badbehav(find(stdmat(:)==0)) =badbehav(find(stdmat(:)==0)) + 1;
            badbeh{g,countnewtotal} = badbehav;	% save instead of disp
            reorder_4beh(:,p) = rri_boot_order(num_subj_lst, k, ...
				1, [1:k], incl_seq, boot_type);
            for c=1:k	% recalc stdmat
                stdmat(c,:) = std(stacked_behavdata(reorder_4beh((1+ ...
                                 (n*(c-1))+span):(n*c+span), p), :));
            end	% scanloop
        end	% while
        if countnewtotal>0 & exist('bootsamp_4beh','var')
            disp(' '); disp('bootsamp_4beh changed'); disp(' ');
        end
        datamat_reorder_4beh = reorder_4beh(:,p);
        behavdata_reorder = reorder_4beh(:,p);
    end	% if ~is_boot_samples
    
    [datamatsvd, datamatsvd_unnorm, datamatcorrs_lst, ...
		stacked_smeanmat] = rri_get_covcor(method, ...
		stacked_datamat, stacked_behavdata, num_groups, ...
		num_subj_lst, num_cond, bscan, meancentering_type, ...
		cormode, single_cond_lst, 1, num_boot, datamat_reorder, ...
		behavdata_reorder, datamat_reorder_4beh);
    
    if nonrotated_boot
        u_p = datamatsvd' * v;
        v_p =  datamatsvd * u;
        u_sum = u_sum + u_p;
        u_sq = u_sq + u_p.^2;
        data_p = stacked_datamat(datamat_reorder_4beh(row_idx),:);
        behav_p = stacked_behavdata(behavdata_reorder(row_idx),:);
        [brainsctmp, behavsctmp, bcorr] = ...
			rri_get_behavscores(data_p, behav_p, normalize(u_p), ...
			normalize(v_p), kk, num_subj_lst, cormode);
        distrib(:,:,p+1) = bcorr;
    else
        [r c] = size(datamatsvd);
        if r <= c
           [pu, sboot, pv] = svd(datamatsvd',0);
        else
           [pv, sboot, pu] = svd(datamatsvd,0);
        end
        rotatemat = rri_bootprocrust(v, pv);
        pu = pu * sboot * rotatemat;
        pv = pv * sboot * rotatemat;
        data_p = stacked_datamat(datamat_reorder_4beh(row_idx),:);
        behav_p = stacked_behavdata(behavdata_reorder(row_idx),:);
        [brainsctmp, behavsctmp, bcorr] = ...
			rri_get_behavscores(data_p, behav_p, ...
			normalize(pu), normalize(pv), kk, num_subj_lst, cormode);
        distrib(:, :, p+1) = bcorr;
    end		% if nonrotated_boot
end	

result.boot_result.num_boot = num_boot;
result.boot_result.clim = clim;
result.boot_result.num_LowVariability_behav_boots = num_LowVariability_behav_boots;

result.boot_result.boot_type = boot_type;
result.boot_result.nonrotated_boot = nonrotated_boot;
u_sum2 = (u_sum.^2) / (num_boot+1);
u_se = sqrt(abs(u_sq-u_sum2)/(num_boot));
ul=clim;
ll=100-clim;
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
test_zeros=find(u_se<=0);
if ~isempty(test_zeros);
    u_se(test_zeros)=1;
end
compare_u = original_u ./ u_se;
if ~isempty(test_zeros);
    compare_u(test_zeros)=0;
end
result.boot_result.compare_u = compare_u;
result.boot_result.u_se = u_se;
result.boot_result.zero_u_se = test_zeros;
result.other_input.meancentering_type = meancentering_type;
result.other_input.cormode = cormode;