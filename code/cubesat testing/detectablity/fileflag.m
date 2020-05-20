% FILEFLAG : Function creates and returns a 'file flag'.
%            Function returns 1 if the number appearing in the ascii file name is non-zero otherwise 0.
%            If the number in the file is non-zero it is set to zero after reading. If the file is not present
%            it is created with a zero in it.
%
%            [flag]=fileflag(name)
%
%            Input  :
%            name   : name of the ascii file including its extension
%            
%            Output :
%            flag   : 1 if non-zero number in name else 0.
%
%            GvW last updated 4-6-2001

function [flag]=fileflag(name)

% Check existence of name file. If not present create it with a zero in it
if ~exist(name,'file'); flag=0; eval(['save ' name ' flag -ascii']); end;

% Determine length variable name
ln=length(name);
for i=ln:-1:1;
  if name(i)=='.'; lne=i-1; break; end;
end

lnb=1;
for i=lne:-1:1;
  if name(i)=='\'; lnb=i+1; break; end;
end
varnam=name(lnb:lne);

% Read flag from name and if non-zero return 1 and delete the file
eval(['load ' name]); eval(['flag=' varnam ';']);
if flag; eval(['delete ' name]); flag=1;
else flag=0; end;