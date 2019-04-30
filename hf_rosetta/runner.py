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

def cpp_eigen():
    result = subprocess.run("./build.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++_eigen/")
    runtime = time.time()
    result = subprocess.run("./run.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++_eigen/")
    runtime = time.time() - runtime
    print_output("C++ Eigen",result,runtime)

def cpp_xtensor():
    result = subprocess.run("./build.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++_xtensor/")
    runtime = time.time()
    result = subprocess.run("./run.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./c++_xtensor/")
    runtime = time.time() - runtime
    print_output("C++ xtensor",result,runtime)

def julia():
    runtime = time.time()
    result = subprocess.run("./run.sh".split(), stdout=subprocess.PIPE, shell=True, cwd="./julia/")
    runtime = time.time() - runtime
    print_output("Julia",result,runtime)

def clean():
    result = subprocess.run("./clean.sh".split(), stdout=subprocess.PIPE, shell=True)

def print_header():
    header = "\n______ ________                           __________        \n___  /____  __/  ___________________________  /__  /______ _\n__  __ \_  /_    __  ___/  __ \_  ___/  _ \  __/  __/  __ `/\n_  / / /  __/    _  /   / /_/ /(__  )/  __/ /_ / /_ / /_/ / \n/_/ /_//_/       /_/    \____//____/ \___/\__/ \__/ \__,_/  \n\n"
    for line in header.splitlines():
        print("{:^79}".format(line))
    print("")
    print("{:^79}".format("Shiv Upadhyay"))
    print("")
    message = "About: HF Rosetta implements the Hartree Fock method in different languages. This is primarily for pedagogical purposes rather than efficiency so don't take the run times too seriously."
    message_list = tw.wrap(message, 40)
    for line in message_list:
        print("{:^79}".format("{:<40}".format(line)))

def run_all():
    python()
    print("{:#<79}".format(""))
    cpp_eigen()
    print("{:#<79}".format(""))
    cpp_xtensor()
    print("{:#<79}".format(""))
    julia()

def main(languages):
    print("{:#<79}".format(""))
    print_header()
    print("{:#<79}".format(""))
    print("Choose a program:")
    for i,key in enumerate(languages):
        print("\t{}. {}".format(i+1,key))
    program = (int)(input("Enter a value: ")) - 1
    print("{:#<79}".format(""))
    try:
        program = languages[program]
    except:
        raise Exception("Selection not recognized")
    if program == "Python":
        python()
    if program == "C++ eigen":
        cpp_eigen()
    if program == "C++ xtensor":
        cpp_xtensor()
    if program == "Julia":
        julia()
    if program == "Clean":
        clean()
    if program == "All":
        run_all()
    print("{:#<79}".format(""))

if __name__ == "__main__":
    languages = ["Python", "C++ eigen", "C++ xtensor", "Julia","All","Clean"]
    main(languages)
