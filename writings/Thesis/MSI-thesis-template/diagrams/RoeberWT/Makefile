SRC_FILES := $(wildcard *.tex)
OUT_FILES :=  $(addsuffix .pdf,$(basename $(wildcard *.tex)))

$(OUT_FILES):
	pdflatex -shell-escape -interaction=nonstopmode -enable-write18 $(addsuffix .tex,$(basename $@))

all: $(OUT_FILES)

clean:
	rm -f *.aux *.auxlock *.log *.dep *.dpth *.md5 *temp* *.pdf *.synctex.gz
