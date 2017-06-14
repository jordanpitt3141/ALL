LAT = pdflatex
LATFLAGS = -shell-escape

o1h:
	$(LAT) $(LATFLAGS) newo1.tex
o1u:
	$(LAT) $(LATFLAGS) newo1u.tex
o2h:
	$(LAT) $(LATFLAGS) newo2.tex
o2u:
	$(LAT) $(LATFLAGS) newo2u.tex
o3h:
	$(LAT) $(LATFLAGS) newo3.tex
o3u:
	$(LAT) $(LATFLAGS) newo3u.tex


clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5

all: o1h o1u o2h o2u o3h o3u

