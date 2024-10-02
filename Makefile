CXX_FLAGS=--std=c++17 -Wall -Wextra -Wpedantic -Werror -Wshadow
DEBUG_FLAGS=-O0 -g -fsanitize=address,undefined,leak,float-divide-by-zero
OPT_FLAGS=-O3 -DNDEBUG -march=native -mtune=native

debug: main.cpp
	g++ $(DEBUG_FLAGS) $(CXX_FLAGS) main.cpp -o main.x

release: main.cpp
	g++ $(OPT_FLAGS) $(CXX_FLAGS) main.cpp -o main.x

format:
	clang-format --style=file -i *.cpp *.h

report: doc/main.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error -Werror doc/main.tex -output-directory=doc

clean:
	rm -f main.x
	cd doc && latexmk -CA && cd ..
