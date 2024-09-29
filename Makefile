main: rgb2greyscale.cpp
	g++ rgb2greyscale.cpp -o main.x

format:
	clang-format --style=file -i rgb2greyscale.cpp *.h

report: doc/main.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error -Werror doc/main.tex -output-directory=doc

clean:
	rm -f main.x
	cd doc && latexmk -CA && cd ..
