from setuptools import setup, Extension

setup(
     ext_modules = [
         Extension(
             name = "strengths.engines.strengths_engine.engine",
             sources = ["src/strengths/engines/strengths_engine/src/engine.cpp"],
             include_dirs = ["src/strengths/engines/strengths_engine/src"],
             define_macros = [("CPYEMVER", None)],
             extra_compile_args = ["-std=c++11"],
             export_symbols = [
                "Initialize3D",
                "InitializeGraph",
                "Run",
                "IterateN",
                "Iterate",
                "GetProgress",
                "GetOutput",
                "GetState",
                "GetT",
                "GetTSample",
                "GetNSamples",
                "Sample",
                "Finalize"
                ]
             )
         ]
     )
