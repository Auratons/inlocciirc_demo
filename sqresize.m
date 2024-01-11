function resized = sqresize(img, sqsize)
    sz = size(img);
    sz_max = max(sz(1), sz(2));
    scale = min(sz_max / sqsize, sqsize / sz_max);
    sz(1:2) = sqsize;
    resized = zeros(sz, class(img));
    rsz = imresize(img, scale);
    rsz_sz = size(rsz);
    [resized_sz_min, min_idx] = min(rsz_sz(1:2));
    offset = floor((sqsize - resized_sz_min) / 2);
    if min_idx == 1
        resized(offset + 1:offset + rsz_sz(1), :, :) = rsz;
    else
        resized(:, offset + 1:offset + rsz_sz(2), :) = rsz;
    end
end
