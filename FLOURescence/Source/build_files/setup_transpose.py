from cx_Freeze import setup, Executable
import os

os.environ['TCL_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tcl8.6"
os.environ['TK_LIBRARY'] = "C:\\Program Files\\Anaconda3\\tcl\\tk8.6"

build_options = {'includes': ['pandas', 'numpy'], 
 				'excludes': ['tkinter', 'cryptography', 'PyQt5', 'zmq', 'matplotlib', 'scipy', 'babel', 'boto', 
 				'statsmodels', 'sqlalchemy', 'sphinx', 'mpl-data', 'mpl_toolkits', 'notebook']}

options = {'build_exe': build_options}
executables = [Executable(script='Transpose.py', targetName='Transpose.exe', icon = 'transpose.ico')]

setup(name = 'Transpose', version = '0.2', options = options, executables = executables)