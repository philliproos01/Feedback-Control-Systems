n = 100

%A = [1 -1 ; 1 1]
%A = [3 -2 1; 1 -3 2; -1 2 4]
%A = rand(n, n);
dominantA = diagonally_dominantMatrix(n)  
A_new = sym_pos_def(n)
B = rand(n, 1)
   
%try chol(A_new)
%    disp('Matrix is symmetric positive definite.')
%catch ME
%    disp('Matrix is not symmetric positive definite')
%end

A = dominantA
b = B
%A=[5 -2 3 0;-3 9 1 -2;2 -1 -7 1; 4 3 -5 7]
%b=[-1 2 3 0.5]'
getSize = size(b)
x=zeros(getSize)

n=size(x,1);



%lk = makeDD(dominantA)
%tf = issymmetric(lk)
%d = eig(lk)
%isposdef = all(d > 0)
normVal=Inf; 
tol=1e-5; itr=0;
while normVal>tol
    x1=x;
    
    for I=1:n
        s=0;
        for l=1:I-1
                s=s+A(I,l)*x(l);
        end
        
        for l=I+1:n
                s=s+A(I,l)*x1(l);
        end
        
        x(I)=(1/A(I,I))*(b(I)-s);
    end   
    itr=itr+1;
    normVal=norm(x1-x);
end
fprintf('Solution set:');
%fprintf('\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f',x);
x_f = x
e = eig(A)

plot(x_f)
hold on
plot(e,'o','LineWidth',3)

function A = diagonally_dominantMatrix(n)
    A = rand(n, n);
    for i = 1:n
        A(i, i) = sum(abs(A(i, :))) + 1;
    end
end

function A = sym_pos_def(n)
    A = rand(n,n); 
    A = 0.5*(A+A');
    A = A + n*eye(n);
end




