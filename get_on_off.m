function [ on, off ] = get_on_off( bin_vec )
% Takes a binary vector and returns the indices of when the vector turns
% "on" (goes from 0 to 1) and "off" (goes from 1 to 0).

if bin_vec(1) == 1;
    on_status = 'on';
    on(1) = 1;
else
    on_status = 'off';
end

n = 1;
for j = 1:length(bin_vec)
     if strcmpi(on_status,'off') && bin_vec(j) == 1
         on(n) = j;
         on_status = 'on';
     elseif strcmpi(on_status,'on') && bin_vec(j) == 0
         off(n) = j;
         on_status = 'off';
         n = n + 1;
     end

end

if strcmpi(on_status,'on') && bin_vec(end) == 1
    off(n) = length(bin_vec);
end

end

