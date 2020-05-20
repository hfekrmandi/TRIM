clear all
close all
A = [1 1 0;
    1 1 1;
    0 1 1] ;
 B(:,:,1) = A;
 C(:,:,1) = (A);
 
 b(:,1)= B(:,:,1)*ones(size(A,1),1);
 c(:,1)= (1/(size(A,1))) * ones(size(A,1),size(A,1))*ones(size(A,1),1);
 for i_=2:100
    B(:,:,i_) = repmat(1./sum( B(:,:,i_-1),2),1,size(A,2)).*B(:,:,i_-1);
%     C(:,:,i_) = (1./sum(sum( tril(C(:,:,i_-1))))).*C(:,:,i_-1);
C(:,:,i_) = (1./size(C,1)).*C(:,:,i_-1);



    E = norm( (B(:,:,i_)  -  C(:,:,i_))*ones(size(C,1),1) );
    disp('C')
     disp(C(:,:,i_)*ones(size(C,1),1))
     disp('B')
      disp(B(:,:,i_)*ones(size(C,1),1))
      disp('-------------')
      E
end