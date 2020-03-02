function [a_new,b_new] = filter_bifurcation_vals(a_valid,b_valid)
b_copy = b_valid;
for i=1:length(a_valid)-1
    if i~=1 && b_copy(i)<(b_valid(i+1)+b_valid(i-1))/2.2
        b_copy(i)=0;
        
    end

end

index = b_copy~=0;
a_new = a_valid(index);
b_new = b_valid(index);

end
