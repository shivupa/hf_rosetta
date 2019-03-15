import numpy as np
import numpy.testing as npt
import time
import subprocess
import textwrap as tw

def print_output(code,result,runtime):
    value = (float)(result.stdout.decode("utf-8"))
    npt.assert_almost_equal(-76.02679819986787,value,decimal=6)
    print("{:^79}".format("{:<15}: {:>21}".format("CODE",code)))
    print("{:^79}".format("{:<15}: {:>21}".format("RESULT / Ha", value)))
    print("{:^79}".format("{:<15}: {:>21}".format("RUNTIME / s", runtime)))


def python():
    runtime = time.time()
    result = subprocess.run("./run.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./python/")
    runtime = time.time() - runtime
    print_output("Python",result,runtime)

def cpp():
    result = subprocess.run("./build.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++/")
    runtime = time.time()
    result = subprocess.run("./run.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++/")
    runtime = time.time() - runtime
    print_output("C++",result,runtime)

def print_header():
    header = "\n______ ________                           __________        \n___  /____  __/  ___________________________  /__  /______ _\n__  __ \_  /_    __  ___/  __ \_  ___/  _ \  __/  __/  __ `/\n_  / / /  __/    _  /   / /_/ /(__  )/  __/ /_ / /_ / /_/ / \n/_/ /_//_/       /_/    \____//____/ \___/\__/ \__/ \__,_/  \n\n"
    for line in header.splitlines():
        print("{:^79}".format(line))
    print("")
    print("{:^79}".format("Shiv Upadhyay"))
    print("")
    message = "About: HF Rosetta implements the Hartree Fock method in different languages. This is primarily for pedagogical purposes rather than efficiency so don't take the runtimes too seriously."
    message_list = tw.wrap(message, 40)
    for line in message_list:
        print("{:^79}".format("{:<40}".format(line)))

if __name__ == "__main__":
    print("{:#<79}".format("")) 
    print_header()
    print("{:#<79}".format("")) 
    print("Choose a program:")
    print("\t1. Python")
    print("\t2. C++")
    print("\t3. All")
    program = (int)(input("Enter a value: "))
    print("{:#<79}".format("")) 
    if  program == 1:
        python()
    elif  program == 2:
        cpp()
    elif program == 3:
        python()
        print("{:#<79}".format("")) 
        cpp()
    print("{:#<79}".format("")) 
