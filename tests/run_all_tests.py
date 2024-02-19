import test_meshgrid
import test_rdnetwork
import test_rdsystem
import test_units
import test_unitarray
import test_loadrds
import test_simulate
import test_librdengine
import test_rdoutput
import test_coarsegrain
import test_kinetics

def run_all_tests() :
    test_meshgrid.run_all_tests()
    test_rdnetwork.run_all_tests()
    test_rdsystem.run_all_tests()
    test_units.run_all_tests()
    test_unitarray.run_all_tests()
    test_loadrds.run_all_tests()
    test_simulate.run_all_tests()
    test_librdengine.run_all_tests()
    test_rdoutput.run_all_tests()
    test_coarsegrain.run_all_tests()
    test_kinetics.run_all_tests()

run_all_tests()
print("All tests have passed successfully.")
