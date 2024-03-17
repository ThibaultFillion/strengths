Building Strengths from source
==============================

Strengths is not a pure python package, as its main simulation engine is implemented in
C++. Engine C++ sources must be compiled as proper shared library/DLL, which will be called internally by 
the python implementation through the ctypes standard Python module. When Strenghts is installed from source, the compiled
shared library/DLL will be missing. To have strengths working properly, one would first have to compile the 
C++ engine files as a shared library (Linux/Macos) or DLL (Windows).

Currently, only the main engine needs to be compiled.
Its C++ source files are located in the "strengths/engines/strengths_engine/src" directory.
They must be compiled as a shared library in the "strengths/engines/strengths_engine" directory,
and it must be the only file in this directory. Any proper file name should be fine for this file.
Let us call it "engine.so" from now on. We thus have to build "strengths/engines/strengths_engine/src/engine.so".

Here is an example using the g++ compiler from the GGC (https://gcc.gnu.org/).
In the "strengths/engines/strengths_engine" directory, g++ can be called as

  g++ -o engine.so src/engine.cpp -shared -static-libstdc++ -static-libgcc -static -O3 -Wall -Wextra

You may need replace the "g++"" by the path where it is installed, depending on you configuration.
If everything works fine, it should produce the file "engine.so" in the "strengths/engines/strengths_engine" directory.
From this point, you should be able to use Strengths normally, as it will be able to use its main simulation engine. 
