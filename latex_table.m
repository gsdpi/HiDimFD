function latex_table(fname,dat,head,form,form_type,align,hline,ndline)

% latex_table(dat,head,form,form_type)

%
% Author:
%	Abel Alberto Cuadrado Vega
%
% 	Grupo de Supervisión, Diagnóstico y Descubrimiento de Conocimiento en Procesos de Ingeniería
% 	(c) Universidad de Oviedo, 2023-2026
%

[m,n]=size(dat);

f=fopen(fname,'w');
if nargin>5
	fprintf(f,['\\begin{tabular}{' align '}\n']);
	fprintf(f,[hline '\n']);
else
	fprintf(f,['\\begin{tabular}{|' repmat('c|',1,n) '}\n']);
	fprintf(f,'\\hline\n');
end

fprintf(f,'%s \\\\\n',head);
if nargin>7
	fprintf(f,'\\hline\n');
else
	fprintf(f,'\\hline\n');
	fprintf(f,'\\hline\n'); 
end
if form_type=='c'
	for r=1:m
		fprintf(f,[form ' \\\\\n'],dat(r,:));
		fprintf(f,'\\hline\n');
	end
end

fprintf(f,'\\end{tabular}\n');

fclose(f);
