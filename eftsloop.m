% eftsloop:
%	Run this script to generate the simulations in Section III
%
% Title:
%	Generation of Interpretable Residuals for Fault Diagnosis based on
%	Projection Techniques: Leveraging Variable Redundancy
% Authors:
%	Abel Alberto Cuadrado Vega, Ignacio Díaz Blanco, José María Enguita González,
%	Diego García Pérez and Ana González Muñiz
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

% Loops evalfaultestim
hFHR=figure;
hF1score=figure;
printaspdf=0;
pause;
alist=[1 2 3 6 12];
threslist=[5 15 30 60 100];
meth={'res','sperbc','sr3','lasso'};
lineopt={'.-b','o-r','x-g','+-c'};
set(groot,'defaultLineMarkerSize',3);
M=1000;
fmax=length(alist);
cmax=length(threslist);
mmax=length(meth);
f=1;
for a=alist
	c=1;
	for threspcent=threslist
		for m=1:mmax
			estmode=meth{m};
			evalfaultestim;

			figure(hFHR);
			subplot(fmax,cmax,c+(f-1)*cmax);
			plot(reldiffdim,fullhitrate*100,lineopt{m});
			ylim([0 100]);
			if f==fmax
				xlabel('(n-d)/n');
			end
			if c==1
				ylabel('Full-hit rate (%)');
			end
			if f==fmax && c==cmax && m==mmax
				legend('Res','SPE-RBC','SR3','LASSO','location','northeast');
			end
			title(['n=' num2str(n) ', a=' num2str(a) ', th=' num2str(threspcent) '%']);
			hold on;

			figure(hF1score);
			subplot(fmax,cmax,c+(f-1)*cmax);
			plot(reldiffdim,F1score*100,lineopt{m});
			ylim([0 100]);
			if f==fmax
				xlabel('(n-d)/n');
			end
			if c==1
				ylabel('F1-score (%)');
			end
			if f==fmax && c==cmax && m==mmax
				legend('Res','SPE-RBC','SR3','LASSO','location','northeast');
			end
			title(['n=' num2str(n) ', a=' num2str(a) ', th=' num2str(threspcent) '%']);
			hold on;

			drawnow;
		end
		c=c+1;
	end
	f=f+1;
end
if printaspdf
	figure(hFHR);
  rows=5;
  cols=5;
  for r=1:rows
      for c=1:cols
          aa=subplot(rows,cols,(r-1)*cols+c);
          pp=get(aa,'Position');
          pp(1)=pp(1)-0.04*c;
          pp(2)=pp(2)+0.04*r;
          set(aa,'Position',pp);
      end;
  end
	print('-dpdfcrop','-S1000,800','-F:6','mc-fhr.pdf')
	figure(hF1score);
  for r=1:rows
      for c=1:cols
          aa=subplot(rows,cols,(r-1)*cols+c);
          pp=get(aa,'Position');
          pp(1)=pp(1)-0.04*c;
          pp(2)=pp(2)+0.04*r;
          set(aa,'Position',pp);
      end;
  end
	print('-dpdfcrop','-S1000,800','-F:6','mc-f1.pdf')
  % print('-dpdfcrop','-S2000,1600','-F:5','mc-f1_.pdf')
end

