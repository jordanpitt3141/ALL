LAT = pdflatex
LATFLAGS = -shell-escape

p0:
	$(LAT) $(LATFLAGS) plotp0.tex

p5:
	$(LAT) $(LATFLAGS) plotp5.tex

p10:
	$(LAT) $(LATFLAGS) plotp10.tex

p15:
	$(LAT) $(LATFLAGS) plotp15.tex

p20:
	$(LAT) $(LATFLAGS) plotp20.tex


clean:
	rm -f  *~ ./*.log ./*.aux ./*.auxlock ./*.dep ./*.dpth ./*.pdf

all: p0 p5 p10 p15 p20

