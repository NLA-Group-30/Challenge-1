main: rgb2greyscale.cpp
	g++ rgb2greyscale.cpp -o main.x

format:
	clang-format --style=file -i rgb2greyscale.cpp *.h

clean:
	rm -f main.x
