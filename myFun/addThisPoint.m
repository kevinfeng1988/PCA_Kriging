function p_plus = addThisPoint(p_orig, p0, I)
%% Description


%% Input
[ri, ci] = size(I);
[rp, cp] = size(p_orig);
[r0, c0] = size(p0);

% Check size consistency
if ri ~= rp
    error('Wrong size provided.');
elseif ci ~= 1
    error('DIM2 should be one.');
elseif cp ~= c0
    error('DIM2 of these variables should be the same.');
end


%% Main
cond = true; % While loop breaker
i = 1; % First index

while (cond)
    new_point = p_orig(I(i),:); % Get the point
    
    % Check the point does not belong to p0
    J = ~ismember(new_point, p0, 'rows');
    if J
        cond = false;
    else
        i = i + 1;
    end
end


%% Output
p_plus = new_point;


end

