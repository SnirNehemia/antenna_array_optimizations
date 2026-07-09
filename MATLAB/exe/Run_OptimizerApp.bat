@echo off
REM Launcher for OptimizerApp.exe — puts the MATLAB Runtime DLLs on PATH first.
REM On a machine with full MATLAB R2026a installed, those DLLs live under
REM <matlabroot>\runtime\win64. On a machine with only the MATLAB Runtime
REM installed, change MATLAB_RUNTIME_DIR to that installation's runtime\win64.

set MATLAB_RUNTIME_DIR=D:\MATLAB\runtime\win64
set PATH=%MATLAB_RUNTIME_DIR%;%PATH%

"%~dp0OptimizerApp.exe"
