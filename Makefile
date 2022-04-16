CXX=g++-11.2
#CXX=clang
XLOGINXX=xoleks00

OBJ=kry.c
BIN=kry

CXXFLAGS:=-Wall -Wextra -Wsuggest-override -Wnull-dereference -Wshadow -Wold-style-cast -pedantic -lgmp -std=c++20

LINK.o = $(LINK.cpp)

all: CXXFLAGS += -Ofast -march=native -flto
all: kry

debug: CXXFLAGS += -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
debug: kry

kry: $(OBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(BIN)

pack:
	zip $(XLOGINXX).zip *.cpp *.hpp  Makefile doc.pdf

dep:
	g++ *.cpp -MM >> Makefile