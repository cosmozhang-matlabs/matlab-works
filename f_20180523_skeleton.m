function [  ] = f_20180523_skeleton(  )

%% read image and transform to grayscale
im = imread('f_20180523_skeleton.png');
im = rgb2gray(im);
imw = size(im, 2);
imh = size(im, 1);

%% binarize the image
fgrng = [0 165];
imbin = zeros(imh, imw);
imbin((im >= fgrng(1)) & (im <= fgrng(2))) = 1;

%% Filter out noise
imbin = medfilt2(imbin, [5 5]);

%% find the connected components
[bwl, nbw] = bwlabel(imbin);
newbwl = zeros(size(bwl));
bwarearng = [10 * 10, imw * imh];
inewbwl = 1;
for i = 1:nbw
    bwa = bwarea(bwl == i);
    if (bwa >= bwarearng(1)) && (bwa <= bwarearng(2))
        newbwl(bwl == i) = inewbwl;
        inewbwl = inewbwl + 1;
    end
end
nbw = inewbwl - 1;
bwl = newbwl;
% imagesc(bwl);

%% Withdraw skeletons

brt = get_all_branches(bwl);
brs = get_branchset(brt);
brbw = zeros(size(bwl));
for j = 1:length(brs)
    brbw(points2indexes(imh, brs{j})) = j;
end
imagesc(brbw);
colorbar;

return;

%% Morphologic skeletonize
skl = zeros(size(bwl));
for i = 1:nbw
    skli = bwmorph(bwl == i, 'skel', Inf);
    skl(skli > 0) = i;
end
imagesc(skl);

%% Analyze the skeletons
sklbw = (skl > 0);
sklnnei = zeros(imh,imw);
for r = 1:imh
    for c = 1:imw
        if sklbw(r,c) > 0
            sklnnei(r,c) = length(find((sklbw(r-1:r+1, c-1:c+1) == 1))) - 1;
        end
    end
end
% sklvtx = indexes2points(find(sklnnei == 1));
% imagesc(sklnnei);
% hold on;
% plot(sklvtx(1,:), sklvtx(2,:), 'ro');
% hold off;

%% withdraw branches
for i = 1:nbw
    sklswi = skl == i;
    vtxs = indexes2points(imh, find(sklswi));
    %nvtx = size(vtxs, 2);
    vtx0 = vtxs(:,1);
    brt = get_branch(sklswi, vtx0);
    brs = get_branchset(brt);
    brbw = zeros(size(skl));
    for j = 1:length(brs)
        brbw(points2indexes(imh, brs{j})) = j;
    end
    imagesc(brbw);
end

function [ branchtree ] = get_all_branches(bwlabels)

nbw = max(max(bwlabels));
skl = zeros(size(bwlabels));
for i = 1:nbw
    skli = bwmorph(bwlabels == i, 'skel', Inf);
    skl(skli > 0) = i;
end
branches = cell(1, nbw);
for i = 1:nbw
    skli = zeros(size(skl));
    skli(skl == i) = 1;
    branches{i} = get_branch(skli);
end

branchtree = struct;
branchtree.major = [];
branchtree.branches = branches;

function [ branchtree, leftbw ] = get_branch(skeleton_bw, start_point)

if nargin < 2
    vtxs = indexes2points(size(skeleton_bw,1), find(skeleton_bw));
    start_point = vtxs(:,1);
end

sx = start_point(1);
sy = start_point(2);

skbw = skeleton_bw > 0;
% imw = size(skbw, 2);
% imh = size(skbw, 1);
maxn = length(find(skbw > 0));

br0 = zeros(2, maxn);
ibr0 = 1;
br0(:,ibr0) = [sx;sy];

bro = cell(0,0);

while 1
    skbw(sy,sx) = 0;
    nps = indexes2points(3, find(skbw(sy-1:sy+1, sx-1:sx+1) > 0)) + [sx-2; sy-2];
    nnps = size(nps, 2);
    if nnps == 0
        break;
    elseif nnps == 1
        ibr0 = ibr0 + 1;
        br0(:,ibr0) = nps;
        sx = nps(1);
        sy = nps(2);
    else
        bro = cell(1, nnps);
        for j = 1:nnps
            skbw(nps(2,j), nps(1,j)) = 0;
        end
        for j = 1:nnps
            [brj, skbw] = get_branch(skbw, nps(:,j));
            bro{j} = brj;
        end
        break;
    end
end

br0 = br0(:,1:ibr0);

branchtree = struct;
branchtree.major = br0;
branchtree.branches = bro;

leftbw = skbw;

function [ branchset ] = get_branchset( branchtree )

nbrs = length(branchtree.branches);
bro = cell(1, nbrs);
nnbro = 0;
for j = 1:nbrs
    bro{j} = get_branchset(branchtree.branches{j});
    nnbro = nnbro + length(bro{j});
end
if isempty(branchtree.major)
    brs = cell(1, nnbro);
    ibrs = 0;
else
    brs = cell(1, nnbro + 1);
    ibrs = 1;
    brs{ibrs} = branchtree.major;
end
for j = 1:nbrs
    brj = bro{j};
    for k = 1:length(brj)
        ibrs = ibrs + 1;
        brs{ibrs} = brj{k};
    end
end

branchset = brs;

function [ pts ] = indexes2points( imh, indexes )

indexes = reshape(indexes, [1,length(indexes)]);
pts = floor((indexes - 1) / imh) + 1;
pts = [pts; indexes - (pts - 1) * imh];

function [ indexes ] = points2indexes( imh, pts )

indexes = (pts(1,:) - 1) * imh + pts(2,:);
