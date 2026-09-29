CXX      = g++
CXXFLAGS = -std=c++17 -O3 -fno-fast-math -Wall -Wextra -Iinclude -Isrc

all: dirs bin/sovsolve.exe bin/nirnaya.exe bin/nirnaya.dll bin/checker.exe bin/unit_tests.exe bin/gen_refinery_models.exe bin/bench_runner.exe

dirs:
	@if not exist bin mkdir bin
	@if not exist data\refinery mkdir data\refinery

bin/sovsolve.exe: src/driver/cli_main.cpp src/driver/c_api.cpp
	$(CXX) $(CXXFLAGS) src/driver/cli_main.cpp src/driver/c_api.cpp -o bin/sovsolve.exe

bin/nirnaya.exe: src/driver/cli_main.cpp src/driver/c_api.cpp
	$(CXX) $(CXXFLAGS) src/driver/cli_main.cpp src/driver/c_api.cpp -o bin/nirnaya.exe

bin/nirnaya.dll: src/driver/c_api.cpp
	$(CXX) $(CXXFLAGS) -shared src/driver/c_api.cpp -o bin/nirnaya.dll

bin/checker.exe: tools/checker.cpp
	$(CXX) $(CXXFLAGS) tools/checker.cpp -o bin/checker.exe

bin/unit_tests.exe: tests/unit_tests.cpp src/driver/c_api.cpp
	$(CXX) $(CXXFLAGS) tests/unit_tests.cpp src/driver/c_api.cpp -o bin/unit_tests.exe

bin/gen_refinery_models.exe: tools/gen_refinery_models.cpp
	$(CXX) $(CXXFLAGS) tools/gen_refinery_models.cpp -o bin/gen_refinery_models.exe

bin/bench_runner.exe: tools/bench_runner.cpp
	$(CXX) $(CXXFLAGS) tools/bench_runner.cpp -o bin/bench_runner.exe

test: all
	bin/unit_tests.exe
	bin/gen_refinery_models.exe
	bin/bench_runner.exe
