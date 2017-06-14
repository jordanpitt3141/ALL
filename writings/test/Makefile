LAT = pdflatex
LATFLAGS = -shell-escape

o1:
	$(LAT) $(LATFLAGS) o1.tex
o2:
	$(LAT) $(LATFLAGS) o2.tex
o3:
	$(LAT) $(LATFLAGS) o3.tex
clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf ./*.gz ./*.md5
all: o1 o2 o3

