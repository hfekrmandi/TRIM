% PLCONV: Plot convergence results
%
%         plconv(trt,trr,epsl,endloop,tit,ylab,nx,nc);
%
%         trt    : trace array to be plotted
%         trr    : relative trace
%         epsl   : epsilon that determines convergence
%         endloop: switch that determines title
%         tit    : String for title
%         ylab   : ylabel string
%         nx     : Full-order (optional with nc)
%         nc     : Reduced-order (optional with nx)
%
%         L.G. Van Willigenburg, W.L. De Koning 01-07-96
%
  function plconv(trt,trr,epsl,endloop,tit,ylab,nx,nc,b)

  ltrt=length(trt);

%   sprm=length(epsl)>1; % Flag indicating convergence based on spectral radius
%   if sprm; epsl=epsl(1); end
  
% parameter that determines the length
% of the final part of trt which is plotted
  lplot=max(50,round(0.2*ltrt));
 
  if ltrt>1
    if ltrt>=lplot;
      xmin=ltrt-lplot+1; xmax=ltrt;
      trt=trt(xmin:xmax);
    else
      xmin=1; xmax=ltrt;
    end

    xp=0.1*(xmax-xmin)+xmin;
    ymax=max(trt); ymin=min(trt);

%     if (ymax-ymin)/(ymax+ymin)<1e-3
%       trt=trt(1)*ones(1,max(size(trt)));
%     end;

    if (ymax-ymin)/(ymax+ymin+eps)<2e-2
      ymax=ymax+0.01*abs(ymax)+1e-2;
      ymin=ymin-0.01*abs(ymin)-1e-2;
    end;

    if nargin>6 ext1=1.3; ext2=1.2; ext3=1.1;
    else ext1=1.2; ext2=1.1; end;

    plot(xmin:xmax,trt);

%    if xmin>=xmax; keyboard; end;
%    if ymin>=ymax; keyboard; end;

    axis([xmin xmax ymin ymin+ext1*(ymax-ymin)]);
    stri=['Convergence when ' num2str(trr) ' < ' num2str(epsl)];
    yp=ymin+ext2*(ymax-ymin);
    text(xp,yp,0,stri);
 
    if nargin>6;
      if length(nc)>10; nc1=nc(1:5); nc2=nc(end-4:end);
        stri=['Full-order=' num2str(nx) '    Reduced-order=' num2str(nc1) ' ... ' num2str(nc2)];
      else
        stri=['Full-order=' num2str(nx) '    Reduced-order=' num2str(nc)];
      end
      if nargin>8; stri=[stri '  b=' num2str(b)]; end;
      yp=ymin+ext3*(ymax-ymin); text(xp,yp,0,stri);
    end;

    xlabel('Iteration number'); ylabel(ylab);
    if endloop==1;
      stri=[tit ': Final convergence plot'];
    else
      stri=[tit ': Intermediate convergence plot'];
    end
    title(stri); figure(gcf); drawnow;
  end