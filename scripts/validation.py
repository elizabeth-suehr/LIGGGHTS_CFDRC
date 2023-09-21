import sys
import os
import subprocess

directories_in_curdir = list(filter(os.path.isdir, os.listdir(os.curdir)))
directories_in_curdir.remove("__pycache__")


# Would run generation of input files for all folders in this directory, Passes the folder name as a argument incase the test.py script uses it.

for test_folder in directories_in_curdir:
    print("Generating Validation files for Test: ", test_folder)
    os.chdir("./"+test_folder)
    try:
        subprocess.run(["python3","test.py", test_folder], check=True)
    except Exception as e:
        print("Error in Folder: "+ test_folder)
        print("exception raised: ", e)
    os.chdir("../")
    print("************************************\n\n")


# Would run generation of input files for all folders in this directory, Passes the folder name as a argument incase the test.py script uses it.

# for test_folder in directories_in_curdir:
#     print("Generating Validation files for Test: ", test_folder)
#     os.chdir("./"+test_folder)
#     try:
#         subprocess.run(["python3","validate.py", test_folder], check=True)
#     except Exception as e:
#         print("Error in Folder: "+ test_folder)
#         print("exception raised: ", e)
#     os.chdir("../")
#     print("************************************\n\n")
