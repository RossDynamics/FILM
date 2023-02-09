function Lambdax = quadratlooping(x,phi)
%A looping version of quadrat which might be slower but which is compatible
%with symbolic inputs.
    Lambdax = phi{1}(:,:,end)*x;
    
    %Slicing up the identity matrix gives us the ith standard basis vector
    B = eye(numel(x));
    for i = 1:numel(x)
        coeff = x.' * squeeze(phi{2}(i,:,:,end)) * x;
        Lambdax = Lambdax + 0.5*coeff * B(:,i);
    end
end