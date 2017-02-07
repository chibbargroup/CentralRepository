'''
Setup.py file for building mobile package using cx_Freeze package


'''

from cx_Freeze import setup, Executable
import os

os.environ['TCL_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tcl8.6"
os.environ['TK_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tk8.6"

build_options = {'packages': ['pandas', 'numpy', 'scipy', 'matplotlib', 'matplotlib.backends.backend_qt5agg'], 'excludes': ['tkinter'], 
				'include_files': [('C:\Program Files\Anaconda3\Library\plugins\platforms\qwindows.dll', 'platforms\qwindows.dll')]}

executables = [Executable(script='caller.py', targetName='FLOUResence.exe', icon="flour.ico", base='Console')]

setup(name='FLOUResence.exe',
	version='0.2',
	options = {"build_exe": build_options},
	executables = executables)
	

#"include_files": ["ConcentrationCalculator.py", "FileCombiner.py", "FormatPlateData.py", 
#			   "LinearRegression.py", "NormalizeData.py", "UserInput.py"