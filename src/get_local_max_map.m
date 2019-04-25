function max_map = get_local_max_map(img)
    eps_0=1e-3;
    max_map = (sign(img(:,2:end-1)-img(:,1:end-2)-eps_0)+1)/2.*...
        (sign(img(:,2:end-1)-img(:,3:end)-eps_0)+1)/2;
    sub_max = zeros(size(max_map));
    sub_max(:,2:end-1)=(sign(img(:,2:end-3)-img(:,1:end-4)-eps_0)+1)/2.*...
        (sign(img(:,4:end-1)-img(:,5:end)-eps_0)+1)/2;
    max_map=max_map.*sub_max;
end