CXX=g++
CXXFLAGS=-std=c++11 -I. -O0
LDFLAGS+=-lutil -lboost_iostreams -lboost_system -lboost_filesystem

%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

numerics: main.o numerics.o statistics.o
	$(CXX) -o numerics main.o numerics.o statistics.o $(CXXFLAGS) $(LDFLAGS)
