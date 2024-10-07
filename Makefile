CXX_FLAGS=-Wall -Wextra -Wpedantic -Wshadow -Werror
DEBUG_FLAGS=-O0 -g #-fsanitize=address,undefined,leak,float-divide-by-zero
OPT_FLAGS=-O3 -DNDEBUG -march=native -mtune=native

debug: main.cpp
	g++ --std=c++17 $(DEBUG_FLAGS) $(CXX_FLAGS) main.cpp -o main.x
	mpicc -DUSE_MPI $(CXX_FLAGS) task_8.c -o task_8.x -I$(LIS_INC) -L$(LIS_LIB) -llis -lm -fopenmp

release: main.cpp
	g++ --std=c++17 $(OPT_FLAGS) $(CXX_FLAGS) main.cpp -o main.x
	mpicc -DUSE_MPI $(OPT_FLAGS) $(CXX_FLAGS) task_8.c -o task_8.x -I$(LIS_INC) -L$(LIS_LIB) -llis -lm -fopenmp

format:
	clang-format --style=file -i *.cpp *.c *.h

report: doc/main.tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error -Werror doc/main.tex -output-directory=doc

clean:
	rm -f main.x
	rm -f task_8.x
	cd doc && latexmk -CA && cd ..
