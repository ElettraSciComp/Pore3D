# PyPore3D for Windows users

This link provides installers and source code of PyPore3D for Windows uers.

PyPore3D: An open source software tool for imaging data processing and analysis of porous and multiphase media.
More Information is found in the gitlab link:
https://gitlab.elettra.eu/aboulhassan.amal/PyPore3D

# If you use this software please cite:

Aboulhassan, A.; Brun, F.; Kourousias, G.; Lanzafame, G.; Voltolini, M.; Contillo, A.; Mancini, L. PyPore3D: An Open Source Software Tool for Imaging Data Processing and Analysis of Porous and Multiphase Media. J. Imaging 2022, 8, 187. https://doi.org/10.3390/jimaging8070187

F. Brun et al., Pore3D: A software library for quantitative analysis of porous media (2010) Nuclear Instruments and Methods in Physics Research, Section A: Accelerators, Spectrometers, Detectors and Associated Equipment, 615 (3), pp. 326-332.

# PyPore3d-Installation-Script

This is the installation script for pypore3d on windows, it includes its dependecies.
The installers are developed by Al-Hassan Hesham (https://github.com/AlHassanHK)

---
**NOTE**

Make sure to save the script in a path the does not contain the following characters: [" , ' ;]
---
```
"C:\Users\alex\Workspace\Company's Stuff\PyPore3d-Installation-Script.ps1" <- will fail
"C:\Users\alex\Desktop\Workspace\Python Stuff\PyPore3d-Installation-Script.ps1" <- should work
```
### Instructions
In order to install and build pypore3d, make sure you have internet connection and launch the "PyPore3d-Installation-Script.ps1": 

- Right click on the .ps1 file
- Choose run as powershell



Launch the script, in case the powershell window closes immediately after launching, open an admin powershell window and type the following:

*In order to launch a powershell admin session:*
- *Right click on powershell icon*
- *Select run as adminstartor*



```
Get-ExecutionPolicy -List > ./execution-policy-output.txt
Set-ExecutionPolicy Unrestricted && Set-ExecutionPolicy -Scope CurrentUser Unrestricted
```

In case this prompt shows up, choose [A]Yes to All. We will revert it back to its original state after finishing the setup 
![image](https://github.com/AlHassanHK/PyPore3d-Installation-Script/assets/87674084/8b4c2cd8-aec1-414a-918b-1a5f068b0c4d)

After finishing the execution, open python and try to run the following:
```
import pypore3d.p3dFiltPy
from pypore3d.p3dFiltPy import py_p3dReadRaw8
import pypore3d.p3dSITKPy
from pypore3d.p3dSITKPy import py_p3d_Dilate
```

After finishing, go to the saved "execution-policy-output.txt" file, and type the following into a powershell admin window

```
Set-ExecutionPolicy <your policy> ; Set-ExecutionPolicy -Scope CurrentUser <your CurrentUser policy>
```
Replace <your policy> and <your CurrentUser policy> with their corresponding values from your "execution-policy-output.txt" file.


---
**NOTE: What will the installer do?**

It will Automatically downloaded requirements: 
- Python + Python launcher (Version 3.10.0) If they are not already installed. If they are ptr-installed, it will use version 3.12.
- SWIG
- SimpleITK
- Visual Studio C++ build tools

If pypore3D is compiled with version 3.10, jupyter notebook or Python scripts need to be run under the same Python version. 
If pypore3D is compiled with version 3.12, jupyter notebook or Python scripts need to be run under the same Python version. 

** One tip to make sure all is compatible is using py -version command. **

For example:
>>Py -3.12 -m notebook
this commant will run jupyter notebook under python 3.12



## Examples
<details>
<summary>View Examples</summary>
<br>

[Sand_in_brine_Casestudy.ipynb](https://gitlab.elettra.eu/aboulhassan.amal/PyPore3D/-/wikis/uploads/addbb4d349668fa94720ae5930aea673/Sand_in_brine_Casestudy.ipynb)

[Geological_Casestudy.ipynb](https://gitlab.elettra.eu/aboulhassan.amal/PyPore3D/-/wikis/uploads/35a57688003cd5f74abb8c80d5491047/Geological_Casestudy.ipynb)

[Bone_Casestudy.ipynb](https://gitlab.elettra.eu/aboulhassan.amal/PyPore3D/-/wikis/uploads/4d179ff71c90379e8af748762e81495e/Bone_Casestudy.ipynb)
</details>

