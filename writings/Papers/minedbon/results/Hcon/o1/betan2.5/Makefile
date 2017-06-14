LAT = pdflatex
LATFLAGS = -shell-escape

o1:
	$(LAT) $(LATFLAGS) sto1.tex
o2:
	$(LAT) $(LATFLAGS) sto2.tex
o3:
	$(LAT) $(LATFLAGS) sto3.tex

clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5

all: o2

