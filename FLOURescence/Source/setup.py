'''
Setup.py file for building mobile package using cx_Freeze package


'''

from cx_Freeze import setup, Executable
import os

os.environ['TCL_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tcl8.6"
os.environ['TK_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tk8.6"

setup(name='FLOUResence.exe',
	version='0.1',
	options = {"build_exe": {"packages":["pandas", "numpy", "scipy"]}
			   },
	executables = [Executable(script='caller.py', targetName='FLOUResence.exe', icon="icon.ico")]
	)

#"include_files": ["ConcentrationCalculator.py", "FileCombiner.py", "FormatPlateData.py", 
#			   "LinearRegression.py", "NormalizeData.py", "UserInput.py"