function res = constback_sub(img)
    sz = size(img);
    res = zeros(sz);
    for num  = 1:sz(1)
        tmp_arr = img(num,:);
        res(num,:)=tmp_arr-mean(tmp_arr(:));
    end
    res(res<0)=0;
end