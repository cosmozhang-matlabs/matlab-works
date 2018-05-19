function [  ] = f_20180516_profiling(  )

% im = imread('f_20180410_profiling/door1.raw.tiff');
% im = imread('f_20180410_profiling/door.raw.tiff');
% imshow(im);

% imp_door = im(643:1627, 2474:3182, :);
% imp_door = im(1216:2201, 2020:2724, :);
name = 'door_piece.gray.raw.tiff';
% name = 'deconv_door_piece_1_0.5.tiff';
imp_door = imread(name);
% imshow(imp_door);

% lps = {'0.2', '0.25', '0.4', '0.5', '1.0'};
lps = {'0.2', '0.25', '0.4', '0.5', '0.55', '0.63', '0.71', '0.83', '1.0', '1.13'};

lines = cell(1,2);
st = struct;
st.name = 'Horizontal Profile';
% st.line = [88,614;494,492];
st.line = [27,684;509,506];
st.color = [0.8, 0.3, 0];
st.lps = lps;
lines{1} = st;
st = struct;
st.name = 'Vertical Profile';
% st.line = [351,352;229,746];
st.line = [354,363;184,830];
st.color = [0, 0.7  , 0];
st.lps = lps;
lines{2} = st;
profile_image(imp_door, lines, ['profile_' name]);
% profile_image(rgb2gray(imp_door), lines, 'f_20180516_profiling\\door_gray');

function [] = profile_image (im, lines, saveto)

if nargin < 3
    saveto = 0;
end

ginputname = '';
if saveto
    ginputname = [saveto '.mat'];
end

imw = size(im,2);
imh = size(im,1);
spwi = 8;
sphi = 6;
spimwi = sphi * imw / imh;

spr = 1;
spc = 1 + length(lines);
spgx = 0.04;
spgy = 0.07;
% spw = (1 - spgx * (spc + 1)) / spc;
sph = (1 - spgy * (spr + 1)) / spr;
spwim = (1 - spgx * (spc + 1)) * spimwi / (spimwi + spwi * (spc-1));
spw = (1 - spgx * (spc + 1)) * spwi / (spimwi + spwi * (spc-1));
figwi = spwi / spw;
fighi = sphi / sph;

subplot('Position', [spgx, spgy, spwim, sph]);
imshow(im);
hold on;
for il = 1:length(lines)
    line = lines{il}.line;
    x0 = line(1,1);
    x1 = line(1,2);
    y0 = line(2,1);
    y1 = line(2,2);
    plot([x0 x1], [y0 y1], 'Color', lines{il}.color, 'LineWidth', 2);
end
hold off;

set(gcf, 'unit', 'inches');
set(gcf, 'position', [1, 3, figwi, fighi]);

mginputi = 1;

for il = 1:length(lines)
    line = lines{il}.line;
    x0 = line(1,1);
    x1 = line(1,2);
    y0 = line(2,1);
    y1 = line(2,2);
    
    fontname = 'Times New Roman';
    fs = 18;
    fs_title = fs;

    len = sqrt((x1-x0)^2 + (y1-y0)^2);
    ilen = round(len);
    ct = (x1 - x0) / len;
    st = (y1 - y0) / len;

    posx = round((x0 + (0:ilen) * ct));
    posy = round((y0 + (0:ilen) * st));

    nchannels = size(im,3);
    prf = zeros(nchannels, length(posx));
    if nchannels == 3
        channel_colors = {'r','g','b'};
    else
        channel_colors = {'k'};
    end
    for c = 1:nchannels
        prf(c, :) = diag(im(posy, posx, c));
    end
    
    minlv = min(min(prf));
    maxlv = max(max(prf));
    tp = 10;
    plotrange = [floor(minlv / tp) * tp - tp * 3, ceil(maxlv / tp) * tp + tp * 2];
    plottickn = 6;
    plottickitv = floor((diff(plotrange) / (plottickn - 1)) / tp) * tp;
    plotticks = plotrange(1) + floor((diff(plotrange) - plottickitv * (plottickn - 1)) / 2 / tp) * tp + (0:(plottickn-1)) * plottickitv;

    spx = spgx * (il+1) + spw * (il-1) + spwim;
    spy = spgy;
    
    subplot('Position', [spx, spy, spw, sph]);
    for c = 1:nchannels
        plot(0:ilen, prf(c,:), channel_colors{c});
        hold on;
    end
    hold off;
    
    lgdw = 0.15;
    lgdh = 0.1;
    lgdx = 0.95 - lgdw;
    lgdy = 0.95 - lgdh;
    if nchannels > 1
        legend({'Red', 'Green', 'Blue'}, 'Position', [lgdx * spw + spx, lgdy * sph + spy, lgdw * spw, lgdh * sph]);
    end
    
    ylim(plotrange);
    set(gca, 'YTick', plotticks);
    xlim([0 ilen]);
    set(gca, 'XTick', 0:100:ilen);
    
    ylabel('Graylevel');
    
    set(gca, 'FontSize', fs);
    set(gca, 'FontName', fontname);
    
    title(lines{il}.name, 'FontSize', fs_title, 'FontName', fontname);
    
    datafs = fs / 72 / sphi * diff(plotrange);
    lh = 1.2;
	ggys = plotrange(1) + datafs * lh * [2, 3, 4];
    color_gray = ones(1,3) * 0.5;
    
    hold on;
    
    lps = lines{il}.lps;
    iggys = 1;
    for ilp = [1:length(lps), fliplr(1:length(lps))]
        [gx, gy, mginputi] = mginput(mginputi, ginputname);
        ggy = ggys(iggys);
        iggys = mod(iggys, length(ggys)) + 1;
        plot([gx gx], [gy, ggy + datafs * 0.6], 'Color', color_gray);
        text(gx, ggy, lps{ilp}, 'HorizontalAlignment', 'center', 'FontSize', fs, 'FontName', fontname);
    end
    [gx, ~, mginputi] = mginput(mginputi, ginputname);
    text(gx, min(ggys) - datafs * lh, 'lp/mm', 'HorizontalAlignment', 'right', 'FontSize', fs, 'FontName', fontname);
    
    hold off;
end

set(gcf,'PaperPositionMode','auto');

if saveto
    print(saveto, '-djpeg', '-r300');
    close gcf;
end

function [outx, outy, nexti] = mginput(idx, filename)

nexti = idx + 1;

notfound = 0;
ginputdata = [];
if exist(filename) == 2
    ginputdata = load(filename);
    ginputdata = ginputdata.ginputdata;
    if size(ginputdata, 2) >= idx
        outx = ginputdata(1,idx);
        outy = ginputdata(2,idx);
        if isnan(outx) || isnan(outy)
            notfound = 1;
        end
    else
        notfound = 1;
    end
else
    notfound = 1;
end

if notfound
    [outx, outy] = ginput(1);
    if isempty(ginputdata)
        ginputdata = zeros(2, idx);
    else
        maxlen = max([size(ginputdata,2), idx]);
        newginputdata = zeros(2, maxlen);
        newginputdata(:, 1:size(ginputdata,2)) = ginputdata;
        ginputdata = newginputdata;
    end
    ginputdata(:,idx) = [outx;outy];
    save(filename, 'ginputdata');
end


