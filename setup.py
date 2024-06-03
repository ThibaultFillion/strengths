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
                "engineexport_initialize_grid",
                "engineexport_initialize_graph",
                "engineexport_run",
                "engineexport_iterate_n",
                "engineexport_iterate",
                "engineexport_get_progress",
                "engineexport_get_trajectory",
                "engineexport_get_state",
                "engineexport_get_time",
                "engineexport_get_tsample",
                "engineexport_get_nsamples",
                "engineexport_sample",
                "engineexport_finalize"
                ]
             )
         ]
     )
