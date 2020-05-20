function output_  = matrix_tests(A,test_type)
output_ = 1;
% % 
% % switch test_type
% %     case 'is_pos'
% %         if max(abs(eig(A)))<10^(-6)
% %             output_ = -1;
% % %             disp('Matrix almmost close to zero')
% %         else
% %             d = eig(A+A');
% % %             disp(['eig of the matrix is  = ', num2str(d')]);
% % %             d(abs(d)<10^(-6))=[];
% % %             disp(['discarding very close to zero elements  = ', num2str(d')]);
% %             output_ = all(d > 0);
% %         end
% %     otherwise
% %         disp('test has not beeen imnplemented or not recognized')
% %         output_ = -1;
% %         
% % end
end