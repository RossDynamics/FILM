function monomials = getMonomials(vars,order)
%GETMONOMIALS Calculates the monomials in the chosen symbolic variables
%vars starting at order 1 and going up to the order specified. 
%The symbolic variables should be provided in a column vector.

monomials = [];
for i = 1:order
    %We get a list containing order numeric representations of 
    %each variable. Not using symbolic variables right now is much faster.
    list = reshape((1:numel(vars)).'*ones(1,i),[],1);

    %And then we run nchoosek and then sort each monomial list so that we 
    %can then filter out non-unique entries
    monomiallists = unique(sort(nchoosek(list,i),2),'rows');
    
    %We then put the symbols in
    monomialsatorder = subs(sym(monomiallists),sym(1:numel(vars)),vars.');

    %And then multiply all variables in each monomial list together to get
    %each monomial
    monomials = [monomials
                 prod(monomialsatorder,2)];
end

% monomials = sym
% %slice = zeros([size(monomials) numel(vars)]);
% for i = 1:numel(vars)
%     %slice{i} = zeros
%     monomials(monomials == i) = vars(i);
% end

end

