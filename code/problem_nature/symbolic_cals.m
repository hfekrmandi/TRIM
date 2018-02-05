nSamples = 500
error_norm1 = [];
error_norm2 = [];
dim = 10;
for i=1:nSamples
    disp(i)
    [error_norm1(i),error_norm2(i)] = cov_int_sample(dim);
end
figure
plot(error_norm1(1:nSamples));hold on
plot(error_norm2(1:nSamples));hold on
 grid on
grid minor
percent_fmincon = numel(find(error_norm1>=0.01))/numel(error_norm1)*100
percent_closed = numel(find(error_norm2>=0.01))/numel(error_norm2)*100


legend('fmincon','closed form')
title(['norm error fmincon  = ',num2str(percent_fmincon),'%  closed form = ',num2str(percent_closed),'%']);
t = annotation('textbox');
t.FontSize = 12;
t.String = ['D1 = CI(A,CI(B,C)); ',char(10),...
'D2 = CI(CI(A,B),C);', char(10),...
'D3 = CI(CI(A,C),B);',char(10),...
'D4 = CI(B,CI(A,C));',char(10),...
'D5 = CI(CI(B,C),A);',char(10),...
'D12 = (D1-D2);',char(10),...
'D13 = (D1-D3);',char(10),...
'D14 = (D1-D4);',char(10),...
'D15 = (D1-D5);',char(10),...
'Ddif = D1 - P_{CI} ;'];


